use super::aligners::{Aligner, MultiMap, MultiMapR};

use crate::adapter::trim_poly_nucleotide;
use crate::barcode::{BarcodeCorrector, OligoFrequncy, Whitelist, filter_cellular_barcodes_ordmag_advanced};
use crate::qc::{AlignQC, Metrics};
use crate::utils::rev_compl_fastq_record;
use anyhow::{bail, Result};
use bstr::BString;
use indexmap::IndexMap;
use itertools::Itertools;
use kdam::{tqdm, BarExt};
use log::{debug, info};
use noodles::{bam, fastq};
use rayon::iter::ParallelIterator;
use rayon::slice::ParallelSlice;
use seqspec::{Assay, FastqReader, Modality, Read, RegionId, SegmentInfo, SegmentInfoElem, Region, LibSpec};
use smallvec::SmallVec;
use std::collections::{HashMap, HashSet};
use std::sync::{Arc, RwLock};

pub struct FastqProcessor {
    assay: Assay,
    current_modality: Option<Modality>,
    mito_dna: HashSet<String>,
    metrics: HashMap<Modality, Metrics>,
    align_qc: HashMap<Modality, AlignQC>,
    barcode_correct_prob: f64, // if the posterior probability of a correction
    // exceeds this threshold, the barcode will be corrected.
    // cellrange uses 0.975 for ATAC and 0.9 for multiome.
    mismatch_in_barcode: usize, // The number of mismatches allowed in barcode
    expected_cells: Option<usize>, // Expected number of cells
    barcode_filtering_quantile: f64, // Quantile for barcode filtering
    barcode_bootstrap_samples: usize, // Number of bootstrap samples for filtering
}

impl FastqProcessor {
    pub fn new(assay: Assay) -> Self {
        debug!("Creating new FastqProcessor");
        Self {
            assay,
            current_modality: None,
            metrics: HashMap::new(),
            align_qc: HashMap::new(),
            mito_dna: HashSet::new(),
            barcode_correct_prob: 0.975,
            mismatch_in_barcode: 1,
            expected_cells: None,
            barcode_filtering_quantile: 0.99,
            barcode_bootstrap_samples: 100,
        }
    }

    pub fn with_barcode_correct_prob(mut self, prob: f64) -> Self {
        debug!("Setting barcode correction probability to {}", prob);
        self.barcode_correct_prob = prob;
        self
    }

    pub fn with_expected_cells(mut self, cells: usize) -> Self {
        debug!("Setting expected number of cells to {}", cells);
        self.expected_cells = Some(cells);
        self
    }
    
    pub fn with_barcode_filtering_params(mut self, quantile: f64, bootstrap_samples: usize) -> Self {
        debug!("Setting barcode filtering parameters: quantile={}, bootstrap_samples={}", 
               quantile, bootstrap_samples);
        self.barcode_filtering_quantile = quantile;
        self.barcode_bootstrap_samples = bootstrap_samples;
        self
    }

    pub fn modality(&self) -> Modality {
        self.current_modality
            .expect("modality not set, please call set_modality first")
    }

    pub fn add_mito_dna(&mut self, mito_dna: impl Into<String>) {
        self.mito_dna.insert(mito_dna.into());
    }

    pub fn with_modality(mut self, modality: Modality) -> Self {
        debug!("Setting modality to {:?}", modality);
        self.current_modality = Some(modality);
        self
    }

    pub fn get_report(&self) -> Metrics {
        debug!("Generating metrics report for modality {:?}", self.modality());
        let mut metrics = self
            .metrics
            .get(&self.modality())
            .map_or(Metrics::default(), |x| x.clone());
        if let Some(align_qc) = self.align_qc.get(&self.modality()) {
            debug!("Adding alignment QC metrics to report");
            align_qc.report(&mut metrics);
        }
        metrics
    }

    /// Align reads and return the alignments.
    /// If the fastq file is paired-end, the alignments will be returned as a tuple.
    /// Otherwise, the alignments will be returned as a single vector.
    ///
    /// # Arguments
    ///
    /// * `num_threads` - The number of threads to use for alignment.
    /// * `chunk_size` - The maximum number of bases in a chunk.
    ///
    /// # Returns
    ///
    /// An iterator of alignments. If the fastq file is paired-end, the alignments will be returned as a tuple.
    /// Otherwise, the alignments will be returned as a single vector.
    pub fn gen_barcoded_alignments<'a, A: Aligner>(
        &'a mut self,
        aligner: &'a mut A,
        num_threads: u16,
        chunk_size: usize,
    ) -> impl Iterator<Item = Vec<(Option<MultiMapR>, Option<MultiMapR>)>> + 'a {
        debug!("Starting gen_barcoded_alignments with chunk_size={}", chunk_size);
        let fq_reader = self.gen_barcoded_fastq(true).with_chunk_size(chunk_size);
        
        // Log reader state
        debug!("FastqReader created. Is paired-end: {}", fq_reader.is_paired_end());
        debug!("Number of annotators: {}", fq_reader.annotators.len());
        debug!("Number of readers: {}", fq_reader.readers.len());
        
        // Log barcode and UMI information
        let barcodes = fq_reader.get_all_barcodes();
        debug!("Barcodes found: {:?}", barcodes);
        let umis = fq_reader.get_all_umi();
        debug!("UMIs found: {:?}", umis);

        let total_reads = fq_reader.total_reads.unwrap_or(0);
        debug!("Total reads reported: {}", total_reads);

        let modality = self.modality();
        debug!("Aligning {} reads...", total_reads);
        let header = aligner.header();
        let mut qc = AlignQC::default();
        self.mito_dna.iter().for_each(|mito| {
            header
                .reference_sequences()
                .get_index_of(&BString::from(mito.as_str()))
                .map(|x| qc.mito_dna.insert(x));
        });
        self.align_qc.insert(modality, qc);

        // Track stats for debugging
        //let mut records_with_barcode = 0;
        //let mut records_without_barcode = 0;
        //let mut alignments_with_bc_tag = 0;
        //let mut alignments_without_bc_tag = 0;

        let mut progress_bar = tqdm!(total = total_reads);
        fq_reader.map(move |data| {
            let align_qc = self.align_qc.get_mut(&modality).unwrap();
            
            debug!("Processing chunk with {} records", data.len());
            // Log details of first few records in chunk
            if !data.is_empty() {
                let sample = &data[0];
                debug!("Sample record - Barcode present: {}, UMI present: {}, Read1 present: {}, Read2 present: {}", 
                    sample.barcode.is_some(),
                    sample.umi.is_some(),
                    sample.read1.is_some(),
                    sample.read2.is_some()
                );
                
                // Check if at least one read is present
                if sample.read1.is_none() && sample.read2.is_none() {
                    panic!("Neither Read1 nor Read2 is present. Please provide at least one read.");
                }
            }
            
            // Perform alignment
            debug!("Aligning {} records...", data.len());
            let results: Vec<_> = aligner.align_reads(num_threads, data);
            
            // Log alignment results
            debug!("Got {} alignment results", results.len());
            
            // Examine first few alignment results
            let limit = std::cmp::min(5, results.len());
            if limit > 0 {
                debug!("Processing first {} alignment results (of {})...", limit, results.len());
                
                for (i, (ali1, ali2)) in results.iter().take(limit).enumerate() {
                    debug!("Alignment pair #{}:", i + 1);
                    
                    match (ali1, ali2) {
                        (Some(_), Some(_)) => {
                            debug!("  Paired alignment");
                        },
                        (Some(_), None) => {
                            debug!("  Read1 only alignment");
                        },
                        (None, Some(_)) => {
                            debug!("  Read2 only alignment");
                        },
                        _ => {
                            debug!("  No alignment information available");
                        }
                    }
                }
            }
            
            // Get QC stats about valid barcode counts
            let valid_ratio = align_qc.frac_valid_barcode();
            if !results.is_empty() && valid_ratio > 0.0 {
                debug!("Current barcoded alignment statistics:");
                debug!("  Valid barcode ratio: {:.2}%", valid_ratio * 100.0);
            } else {
                debug!("No valid barcoded alignments found yet");
            }
            
            // Process results for QC metrics
            results.iter().for_each(|ali| match ali {
                (Some(ali1), Some(ali2)) => {
                    align_qc.add_pair(&header, ali1, ali2).unwrap();
                }
                (Some(ali1), None) => {
                    align_qc.add_read1(&header, ali1).unwrap();
                }
                (None, Some(ali2)) => {
                    align_qc.add_read2(&header, ali2).unwrap();
                }
                _ => {
                    debug!("No alignment found for read");
                }
            });
            
            progress_bar.update(results.len()).unwrap();
            results
        })
    }

    pub fn gen_barcoded_fastq(&mut self, correct_barcode: bool) -> AnnotatedFastqReader {
        let modality = self.modality();
        debug!("Starting gen_barcoded_fastq for modality: {:?}", modality);

        let whitelists = if correct_barcode {
            debug!("Counting barcodes...");
            debug!("Calling count_barcodes()");
            match self.count_barcodes() {
                Ok(wl) => {
                    debug!("Successfully counted barcodes. Found {} whitelists", wl.len());
                    // Log details for each whitelist
                    for (id, whitelist) in wl.iter() {
                        debug!(
                            "Whitelist '{}': exact match rate={:.2}%, unique barcodes={}",
                            id,
                            whitelist.frac_exact_match() * 100.0,
                            whitelist.num_seen_barcodes()
                        );
                    }
                    wl
                }
                Err(e) => {
                    debug!("Error counting barcodes: {}", e);
                    IndexMap::new()
                }
            }
        } else {
            debug!("Skipping barcode correction");
            IndexMap::new()
        };

        debug!("Creating BarcodeCorrector with threshold={}", self.barcode_correct_prob);
        let corrector = BarcodeCorrector::default()
            .with_max_missmatch(self.mismatch_in_barcode)
            .with_bc_confidence_threshold(self.barcode_correct_prob);

        debug!("Creating FastqAnnotators from segments");
        
        // Added detailed logging about seqspec and FASTQ files
        let segments_by_modality: Vec<_> = self.assay.get_segments_by_modality(modality).collect();
        debug!("Found {} read segments for modality {:?}", segments_by_modality.len(), modality);
        
        // Track readers creation process
        let mut readers_created = 0;
        let mut reads_with_errors = 0;
        
        let mut fq_reader: AnnotatedFastqReader = self
            .assay
            .get_segments_by_modality(modality)
            .filter(|(read, _)| {
                // Check if file can be opened and log outcome
                // Use the read_id instead of path since path() doesn't exist
                let read_id = &read.read_id;
                let open_result = read.open();
                if open_result.is_none() {
                    debug!("ERROR: Unable to open FASTQ file for read ID: {}", read_id);
                    reads_with_errors += 1;
                    false
                } else {
                    debug!("Successfully opened FASTQ file for read ID: {}", read_id);
                    true
                }
            })
            .filter_map(|(read, index)| {
                debug!(
                    "Processing read {} with {} segments, is_reverse: {}", 
                    read.read_id,
                    index.segments.len(),
                    read.is_reverse()
                );
                
                // Log each segment's type
                for seg in &index.segments {
                    debug!(
                        "Segment in read {}: type={:?}, range={:?}, id={}", 
                        read.read_id, 
                        seg.region_type, 
                        seg.range,
                        seg.region_id
                    );
                }

                let library_spec = Arc::new(RwLock::new(self.assay.library_spec.clone()));
                let annotator = match FastqAnnotator::new(read, index, &whitelists, corrector.clone(), library_spec) {
                    Some(a) => a,
                    None => {
                        debug!("ERROR: Failed to create annotator for read {}", read.read_id);
                        return None;
                    }
                };
                
                debug!(
                    "Created annotator for read {} with {} subregions",
                    annotator.id,
                    annotator.subregions.len()
                );
                
                // Check if file can be opened again specifically for reading
                let reader = match read.open() {
                    Some(r) => {
                        readers_created += 1;
                        debug!("Successfully created reader #{} for read {}", readers_created, read.read_id);
                        r
                    },
                    None => {
                        debug!("ERROR: Failed to open reader for {} when creating FastqReader", read.read_id);
                        return None;
                    }
                };
                
                Some((annotator, reader))
            })
            .collect();

        debug!("Created {} readers out of {} potential reads ({} had errors opening)",
            readers_created, segments_by_modality.len(), reads_with_errors);

        if !whitelists.is_empty() {
            fq_reader.total_reads = Some(whitelists[0].total_count);
            debug!("Set total_reads to {}", whitelists[0].total_count);
        } else {
            debug!("No whitelists available, total_reads not set");
        }

        // Add validation of reader state
        if fq_reader.readers.is_empty() {
            debug!("CRITICAL ERROR: No valid FASTQ readers were created!");
        } else {
            for (i, _reader) in fq_reader.readers.iter().enumerate() {
                debug!("Reader #{} validation", i);
                // Remove reader.path() calls since FastqReader doesn't have this method
                
                // We can't check file metadata without the path
                debug!("  Note: Cannot check file metadata (path method not available)");
            }
        }

        debug!("Completed gen_barcoded_fastq setup");
        fq_reader
    }

    fn count_barcodes(&mut self) -> Result<IndexMap<RegionId, Whitelist>> {
        let modality = self.modality();
        let mut whitelists = self.get_whitelists();
        let mut total_filtered_bcs = 0;

        // Add debug info about initial whitelists
        for (id, whitelist) in whitelists.iter() {
            debug!("Initial whitelist '{}' state:", id);
            debug!("  - has_entries: {}", !whitelist.get_barcode_counts().is_empty());
            debug!("  - num_seen_barcodes: {}", whitelist.num_seen_barcodes());
            debug!("  - total_count: {}", whitelist.total_count);
            
            // Print a sample of the whitelist entries (first 5)
            let barcode_sample: Vec<_> = whitelist.get_barcode_counts().iter()
                .take(5)
                .map(|(bc, count)| (String::from_utf8_lossy(bc), count))
                .collect();
            if !barcode_sample.is_empty() {
                debug!("  - sample entries: {:?}", barcode_sample);
            } else {
                debug!("  - sample entries: <empty>");
            }
            
            // Check if this whitelist is predefined (has entries but no counts)
            let is_empty = whitelist.get_barcode_counts().is_empty();
            let num_seen = whitelist.num_seen_barcodes();
            let has_predefined_whitelist = !is_empty;
            
            debug!("Whitelist '{}' predefined check:", id);
            debug!("  - is_empty: {}", is_empty);
            debug!("  - num_seen_barcodes: {}", num_seen);
            debug!("  - has_predefined_whitelist: {}", has_predefined_whitelist);
        }

        debug!("Counting barcodes");
        
        // Create a list of IDs that need counting
        // We need to count barcodes for ALL whitelists, regardless of whether they have predefined entries
        // This ensures we get accurate barcode counts for filtering later
        let whitelist_ids: Vec<RegionId> = whitelists.keys().cloned().collect();
        
        if !whitelist_ids.is_empty() {
            debug!("Found {} whitelists that need counting", whitelist_ids.len());
            
            // Get segments with barcode regions
            let barcode_segments: Vec<_> = self.assay
                .get_segments_by_modality(modality)
                .filter(|(_, region_index)| {
                    region_index.segments.iter().any(|x| x.region_type.is_barcode())
                })
                .collect();
            
            debug!("Processing {} segment groups with barcode regions", barcode_segments.len());
            
            // Process each read segment for barcodes
            for (read, region_index) in barcode_segments {
                debug!("Processing read {}", read.read_id);
                
                let is_reverse = read.is_reverse();
                
                // Read all records for this segment
                if let Some(mut reader) = read.open() {
                    for fastq_record_result in reader.records() {
                        let fastq_record = fastq_record_result.unwrap();
                        
                        // Process each barcode region directly
                        for info in region_index.segments.iter() {
                            if info.region_type.is_barcode() && whitelist_ids.contains(&info.region_id) {
                                let start_pos = info.range.start as usize;
                                let end_pos = info.range.end as usize;
                                
                                // Extract barcode from this region
                                let mut barcode_slice = slice_fastq_record(&fastq_record, start_pos, end_pos);
                                if is_reverse {
                                    barcode_slice = rev_compl_fastq_record(barcode_slice);
                                }
                                
                                // Find the corresponding whitelist
                                if let Some(whitelist) = whitelists.get_mut(&info.region_id) {
                                    // Count the barcode - this is essential for filtering
                                    whitelist.count_barcode(barcode_slice.sequence(), barcode_slice.quality_scores());
                                }
                            }
                        }
                    }
                }
            }
        } else {
            debug!("No whitelists found for counting");
        }

        // Add debug info after counting
        for (id, whitelist) in whitelists.iter() {
            debug!("After counting, whitelist '{}' state:", id);
            debug!("  - has_entries: {}", !whitelist.get_barcode_counts().is_empty());
            debug!("  - num_seen_barcodes: {}", whitelist.num_seen_barcodes());
            debug!("  - total_count: {}", whitelist.total_count);
            debug!("  - exact_match_rate: {:.2}%", whitelist.frac_exact_match() * 100.0);
            
            // Show counts distribution
            let counts = whitelist.get_sorted_counts();
            if !counts.is_empty() {
                debug!("  - counts distribution: [max: {}, median: {}, min: {}]", 
                    counts.first().unwrap_or(&0),
                    counts.get(counts.len() / 2).unwrap_or(&0),
                    counts.last().unwrap_or(&0)
                );
            }
            
            // Print top 5 most frequent barcodes
            let top_barcodes: Vec<_> = whitelist.get_barcode_counts().iter()
                .sorted_by(|a, b| b.1.cmp(a.1))
                .take(5)
                .map(|(bc, count)| (String::from_utf8_lossy(bc), count))
                .collect();
            if !top_barcodes.is_empty() {
                debug!("  - top barcodes: {:?}", top_barcodes);
            }
        }

        // Apply advanced filtering to each whitelist if needed
        for (id, whitelist) in whitelists.iter_mut() {
            // Check if this whitelist needs filtering (no predefined entries or has counts)
            let is_empty = whitelist.get_barcode_counts().is_empty();
            let num_seen = whitelist.num_seen_barcodes();
            let has_predefined_whitelist = !is_empty && num_seen == 0;
            
            if !has_predefined_whitelist && num_seen > 0 {
                info!("No predefined whitelist for '{}', applying order-of-magnitude filtering", id);
                
                let results = filter_cellular_barcodes_ordmag_advanced(
                    whitelist,
                    self.expected_cells,
                    None,
                    None,
                    Some(self.barcode_filtering_quantile),
                    Some(self.barcode_bootstrap_samples),
                );
                
                info!(
                    "Found {} cellular barcodes (95% CI: {}-{}) for '{}' with counts >= {}",
                    results.filtered_bcs, 
                    results.filtered_bcs_lb,
                    results.filtered_bcs_ub,
                    id,
                    results.filtered_bcs_cutoff
                );
                
                total_filtered_bcs += results.filtered_bcs;
                
                // Store filtering metrics
                self.metrics.entry(modality).or_default().insert(
                    format!("filtered_bcs_{}", id),
                    results.filtered_bcs as f64,
                );
                self.metrics.entry(modality).or_default().insert(
                    format!("filtered_bcs_cutoff_{}", id),
                    results.filtered_bcs_cutoff as f64,
                );
            } else {
                info!("Using predefined whitelist for '{}', skipping filtering", id);
            }
        }
        
        if total_filtered_bcs > 0 {
            self.metrics.entry(modality).or_default().insert(
                "total_filtered_bcs".to_string(),
                total_filtered_bcs as f64,
            );
        }

        self.metrics.entry(modality).or_default().insert(
            "frac_q30_bases_barcode".to_string(),
            whitelists.values().map(|x| x.frac_q30_bases()).sum::<f64>() / whitelists.len() as f64,
        );
        
        Ok(whitelists)
    }

    fn get_whitelists(&self) -> IndexMap<RegionId, Whitelist> {
        debug!("Getting whitelists for modality {:?}", self.modality());
        let regions = self
            .assay
            .library_spec
            .get_modality(&self.modality())
            .unwrap()
            .read()
            .unwrap();
        
        let whitelists: IndexMap<RegionId, Whitelist> = regions
            .subregions
            .iter()
            .filter_map(|r| {
                let r = r.read().unwrap();
                if r.region_type.is_barcode() {
                    let id = r.region_id.to_string();
                    let list = if let Some(onlist) = r.onlist.as_ref() {
                        debug!("Found whitelist for barcode region {}", id);
                        Whitelist::new(onlist.read().unwrap())
                    } else {
                        debug!("No whitelist found for barcode region {}, using empty list", id);
                        Whitelist::empty()
                    };
                    Some((id, list))
                } else {
                    None
                }
            })
            .collect();

        debug!("Found {} whitelists", whitelists.len());
        whitelists
    }
}

pub struct AnnotatedFastqReader {
    buffer: fastq::Record,
    total_reads: Option<usize>,
    trim_poly_a: bool,
    annotators: Vec<FastqAnnotator>,
    readers: Vec<FastqReader>,
    chunk_size: usize,
    chunk: Vec<SmallVec<[fastq::Record; 4]>>,
}

impl AnnotatedFastqReader {
    pub fn with_chunk_size(mut self, chunk_size: usize) -> Self {
        self.chunk_size = chunk_size;
        self
    }

    pub fn with_polya_trimmed(mut self) -> Self {
        self.trim_poly_a = true;
        self
    }

    pub fn get_all_barcodes(&self) -> Vec<(&str, usize)> {
        self.annotators
            .iter()
            .flat_map(|annotator| {
                annotator
                    .subregions
                    .iter()
                    .filter(|info| info.region_type.is_barcode())
                    .map(|info| (info.region_id.as_str(), info.range.len()))
            })
            .collect()
    }

    pub fn get_all_umi(&self) -> Vec<(&str, usize)> {
        self.annotators
            .iter()
            .flat_map(|annotator| {
                annotator
                    .subregions
                    .iter()
                    .filter(|info| info.region_type.is_umi())
                    .map(|info| (info.region_id.as_str(), info.range.len()))
            })
            .collect()
    }

    pub fn is_paired_end(&self) -> bool {
        let mut has_read1 = false;
        let mut has_read2 = false;
        self.annotators.iter().for_each(|x| {
            x.subregions.iter().for_each(|info| {
                if info.region_type.is_target() {
                    if x.is_reverse {
                        has_read1 = true;
                    } else {
                        has_read2 = true;
                    }
                }
            });
        });
        has_read1 && has_read2
    }

    /// Read a chunk of records from the fastq files.
    fn read_chunk(&mut self) -> usize {
        self.chunk.clear();

        let mut accumulated_length = 0;
        let mut _total_records = 0;  // Already correctly prefixed with underscore

        while accumulated_length < self.chunk_size {
            let mut max_read = 0;
            let mut min_read = usize::MAX;
            
            // Track which readers failed/succeeded for this chunk
            let mut successful_readers = 0;
            let mut empty_readers = 0;
            let mut error_readers = 0;
            
            let records: SmallVec<[_; 4]> = self
                .readers
                .iter_mut()
                .enumerate()
                .flat_map(|(_i, reader)| {  // Add underscore to i
                    let result = reader.read_record(&mut self.buffer);
                    
                    if let Err(_e) = &result {  // Add underscore to e
                        error_readers += 1;
                        return None;
                    }
                    
                    let n = result.expect("error reading fastq record");
                    
                    if n > 0 {
                        successful_readers += 1;
                        min_read = min_read.min(n);
                        max_read = max_read.max(n);
                        accumulated_length += self.buffer.sequence().len();
                        Some(self.buffer.clone())
                    } else {
                        empty_readers += 1;
                        min_read = 0;
                        None
                    }
                })
                .collect();
            
            if max_read == 0 {
                break;
            } else if min_read == 0 && successful_readers > 0 {
                // Existing commented code
            } else {
                _total_records += 1;
                
                // Check records for name consistency
                if records.len() > 0 {
                    let _first_name = records[0].name();  // Add underscore to first_name
                    
                    let names_match = records.iter().map(|r| r.name()).all_equal();
                    if !names_match {
                        debug!("WARNING: Read names don't match in this chunk!");
                    }
                    
                    assert!(
                        records.iter().map(|r| r.name()).all_equal(),
                        "read names mismatch"
                    );
                } else {
                    debug!("WARNING: Empty records vector even though max_read > 0");
                }
                
                self.chunk.push(records);
            }
        }
        
        debug!("Finished read_chunk, got {} chunk entries", self.chunk.len());
        self.chunk.len()
    }
}

impl FromIterator<(FastqAnnotator, FastqReader)> for AnnotatedFastqReader {
    fn from_iter<T: IntoIterator<Item = (FastqAnnotator, FastqReader)>>(iter: T) -> Self {
        let (annotators, readers): (Vec<_>, Vec<_>) = iter.into_iter().unzip();
        let chunk = Vec::new();
        Self {
            buffer: fastq::Record::default(),
            total_reads: None,
            annotators,
            readers,
            trim_poly_a: false,
            chunk_size: 10000000,
            chunk,
        }
    }
}

impl Iterator for AnnotatedFastqReader {
    type Item = Vec<AnnotatedFastq>;

    fn next(&mut self) -> Option<Self::Item> {
        let n = self.read_chunk();
        debug!("read_chunk returned {} entries", n);
        
        if n == 0 {
            debug!("No chunks read, returning None");
            None
        } else {
            let n = (n / 256).max(256);
            let annotators = &self.annotators;
            
            debug!("Processing chunk with {} annotators", annotators.len());
            
            if annotators.is_empty() {
                debug!("ERROR: No annotators available!");
                return None;
            }
            
            let result: Vec<AnnotatedFastq> = self
                .chunk
                .par_chunks(n)
                .flat_map_iter(|chunk| {
                    chunk.into_iter().map(move |records| {
                        // Remove potentially problematic debug statements
                        records
                            .iter()
                            .enumerate()
                            .map(|(i, record)| {
                                if i >= annotators.len() {
                                    panic!("Annotator index out of bounds");
                                }
                                annotators[i].annotate(record).unwrap()
                            })
                            .reduce(|mut this, other| {
                                this.join(other);
                                this
                            })
                            .unwrap_or_else(|| {
                                panic!("Failed to reduce annotated records");
                            })
                    })
                })
                .collect();
            
            debug!("Returning {} processed records", result.len());
            
            // Debug print the first few records
            for (i, record) in result.iter().take(5).enumerate() {
                debug!("Record #{} details:", i);
                
                // Log barcode information
                if let Some(bc) = &record.barcode {
                    let barcode_seq = std::str::from_utf8(bc.raw.sequence()).unwrap_or("<invalid UTF-8>");
                    let barcode_qual = std::str::from_utf8(bc.raw.quality_scores()).unwrap_or("<invalid UTF-8>");
                    debug!("  Barcode: {}", barcode_seq);
                    debug!("  Barcode quality: {}", barcode_qual);
                    if let Some(corrected) = &bc.corrected {
                        let corrected_str = std::str::from_utf8(corrected).unwrap_or("<invalid UTF-8>");
                        debug!("  Corrected barcode: {}", corrected_str);
                    }
                } else {
                    debug!("  No barcode present");
                }
                
                // Log UMI information
                if let Some(umi) = &record.umi {
                    let umi_seq = std::str::from_utf8(umi.sequence()).unwrap_or("<invalid UTF-8>");
                    let umi_qual = std::str::from_utf8(umi.quality_scores()).unwrap_or("<invalid UTF-8>");
                    debug!("  UMI: {}", umi_seq);
                    debug!("  UMI quality: {}", umi_qual);
                } else {
                    debug!("  No UMI present");
                }
                
                // Log read1 information
                if let Some(read1) = &record.read1 {
                    let read1_seq = std::str::from_utf8(read1.sequence()).unwrap_or("<invalid UTF-8>");
                    let _read1_qual = std::str::from_utf8(read1.quality_scores()).unwrap_or("<invalid UTF-8>");
                    let read1_name = std::str::from_utf8(read1.name()).unwrap_or("<invalid UTF-8>");
                    debug!("  Read1 name: {}", read1_name);
                    debug!("  Read1 sequence (first 30 bp): {}", &read1_seq[..read1_seq.len().min(30)]);
                    debug!("  Read1 length: {}", read1.sequence().len());
                } else {
                    debug!("  No read1 present");
                }
                
                // Log read2 information
                if let Some(read2) = &record.read2 {
                    let read2_seq = std::str::from_utf8(read2.sequence()).unwrap_or("<invalid UTF-8>");
                    let _read2_qual = std::str::from_utf8(read2.quality_scores()).unwrap_or("<invalid UTF-8>");
                    let read2_name = std::str::from_utf8(read2.name()).unwrap_or("<invalid UTF-8>");
                    debug!("  Read2 name: {}", read2_name);
                    debug!("  Read2 sequence (first 30 bp): {}", &read2_seq[..read2_seq.len().min(30)]);
                    debug!("  Read2 length: {}", read2.sequence().len());
                } else {
                    debug!("  No read2 present");
                }
                
                debug!("  Is empty: {}", record.is_empty());
                debug!("  Total length: {}", record.len());
            }
            
            Some(result)
        }
    }
}

/// A FastqAnnotator that splits the reads into subregions, e.g., barcode, UMI, and
/// return annotated reads.
#[derive(Debug)]
struct FastqAnnotator {
    whitelists: IndexMap<String, OligoFrequncy>,
    corrector: BarcodeCorrector,
    id: String,
    is_reverse: bool,
    subregions: Vec<SegmentInfoElem>,
    min_len: usize,
    max_len: usize,
    library_spec: Arc<RwLock<LibSpec>>,
}

impl FastqAnnotator {
    pub fn new(
        read: &Read,
        index: SegmentInfo,
        whitelists: &IndexMap<String, Whitelist>,
        corrector: BarcodeCorrector,
        library_spec: Arc<RwLock<LibSpec>>,
    ) -> Option<Self> {
        let subregions: Vec<_> = index
            .segments
            .into_iter()
            .filter(|x| {
                x.region_type.is_barcode() || x.region_type.is_umi() || x.region_type.is_target()
            }) // only barcode and target regions
            .collect();
        if subregions.is_empty() {
            None
        } else {
            let whitelists = subregions
                .iter()
                .flat_map(|info| {
                    let v = whitelists.get(&info.region_id)?;
                    Some((info.region_id.clone(), v.get_barcode_counts().clone()))
                })
                .collect();
            let anno = Self {
                whitelists,
                corrector,
                id: read.read_id.clone(),
                is_reverse: read.is_reverse(),
                subregions,
                min_len: read.min_len as usize,
                max_len: read.max_len as usize,
                library_spec,
            };
            Some(anno)
        }
    }

    fn annotate(&self, record: &fastq::Record) -> Result<AnnotatedFastq> {
        let n = record.sequence().len();
        if n < self.min_len || n > self.max_len {
            bail!(
                "Read length ({}) out of range: {}-{}",
                n,
                self.min_len,
                self.max_len
            );
        }
        
        let mut barcode: Option<Barcode> = None;
        let mut umi = None;
        let mut read1 = None;
        let mut read2 = None;

        // Process all regions
        for (_, info) in self.subregions.iter().enumerate() {
            let start_pos = info.range.start as usize;
            let end_pos = info.range.end as usize;
            
            // Extract the appropriate slice of the record
            let mut record_slice = slice_fastq_record(record, start_pos, end_pos);
            
            // Reverse complement if needed
            if self.is_reverse && (info.region_type.is_barcode() || info.region_type.is_umi()) {
                record_slice = rev_compl_fastq_record(record_slice);
            }
            
            // Handle based on region type
            if info.region_type.is_barcode() {
                let corrected = self.whitelists.get(&info.region_id).map_or(
                    Some(record_slice.sequence().to_vec()),
                    |counts| {
                        self.corrector
                            .correct(counts, record_slice.sequence(), record_slice.quality_scores())
                            .ok()
                            .map(|x| x.to_vec())
                    },
                );
                
                if let Some(bc) = &mut barcode {
                    bc.extend(&Barcode { raw: record_slice, corrected });
                } else {
                    barcode = Some(Barcode { raw: record_slice, corrected });
                }
            } else if info.region_type.is_umi() {
                umi = Some(record_slice);
            } else if info.region_type.is_target() {
                if read1.is_some() || read2.is_some() {
                    panic!("Both Read1 and Read2 are set");
                } else {
                    if let Some(nucl) = info.region_type.poly_nucl() {
                        if let Some(idx) = trim_poly_nucleotide(nucl, record_slice.sequence().iter().copied())
                        {
                            record_slice = slice_fastq_record(&record_slice, idx, record_slice.sequence().len());
                        }
                    }
                    // Only keep reads with length >= 8
                    if record_slice.sequence().len() >= 8 {
                        if self.is_reverse {
                            read2 = Some(record_slice);
                        } else {
                            read1 = Some(record_slice);
                        }
                    }
                }
            }
            // We simply ignore phase blocks
        }
        
        Ok(AnnotatedFastq {
            barcode,
            umi,
            read1,
            read2,
        })
    }
}

fn find_pattern(haystack: &[u8], needle: &[u8]) -> Option<usize> {
    if needle.is_empty() || haystack.is_empty() || needle.len() > haystack.len() {
        return None;
    }
    
    // Simple pattern matching algorithm
    // You may want to use a more efficient algorithm like Boyer-Moore or Knuth-Morris-Pratt
    for i in 0..=haystack.len() - needle.len() {
        let mut match_found = true;
        for j in 0..needle.len() {
            if haystack[i + j] != needle[j] {
                match_found = false;
                break;
            }
        }
        if match_found {
            return Some(i);
        }
    }
    
    None
}

#[derive(Debug)]
pub struct Barcode {
    pub raw: fastq::Record,
    pub corrected: Option<Vec<u8>>,
}

impl Barcode {
    pub fn extend(&mut self, other: &Self) {
        extend_fastq_record(&mut self.raw, &other.raw);
        if let Some(c2) = &other.corrected {
            if let Some(c1) = &mut self.corrected {
                c1.extend_from_slice(c2);
            }
        } else {
            self.corrected = None;
        }
    }
}

pub type UMI = fastq::Record;

/// An annotated fastq record with barcode, UMI, and sequence.
#[derive(Debug)]
pub struct AnnotatedFastq {
    pub barcode: Option<Barcode>,
    pub umi: Option<UMI>,
    pub read1: Option<fastq::Record>,
    pub read2: Option<fastq::Record>,
}

impl AnnotatedFastq {
    /// The total number of bases, including read1 and read2, in the record.
    pub fn len(&self) -> usize {
        self.read1.as_ref().map_or(0, |x| x.sequence().len())
            + self.read2.as_ref().map_or(0, |x| x.sequence().len())
    }
    pub fn is_empty(&self) -> bool {
        self.read1.is_none() && self.read2.is_none()
    }
}

impl AnnotatedFastq {
    pub fn join(&mut self, other: Self) {
        if let Some(bc) = &mut self.barcode {
            if let Some(x) = other.barcode.as_ref() {
                bc.extend(x)
            }
        } else {
            self.barcode = other.barcode;
        }

        if let Some(umi) = &mut self.umi {
            if let Some(x) = other.umi.as_ref() {
                extend_fastq_record(umi, x)
            }
        } else {
            self.umi = other.umi;
        }

        if self.read1.is_some() {
            if other.read1.is_some() {
                panic!("Read1 already exists");
            }
        } else {
            self.read1 = other.read1;
        }

        if self.read2.is_some() {
            if other.read2.is_some() {
                panic!("Read2 already exists");
            }
        } else {
            self.read2 = other.read2;
        }
    }
}

fn slice_fastq_record(record: &fastq::Record, start: usize, end: usize) -> fastq::Record {
    let end = end.min(record.sequence().len());
    fastq::Record::new(
        record.definition().clone(),
        &record.sequence()[start..end],
        record.quality_scores().get(start..end).unwrap(),
    )
}

pub fn extend_fastq_record(this: &mut fastq::Record, other: &fastq::Record) {
    this.sequence_mut().extend_from_slice(other.sequence());
    this.quality_scores_mut()
        .extend_from_slice(other.quality_scores());
}

pub struct NameCollatedRecords<'a, R> {
    records: bam::io::reader::Records<'a, R>,
    prev_record: Option<(BString, bam::Record)>,
    checker: HashSet<BString>,
}

impl<'a, R: std::io::Read> NameCollatedRecords<'a, R> {
    pub fn new(records: bam::io::reader::Records<'a, R>) -> Self {
        Self {
            records,
            prev_record: None,
            checker: HashSet::new(),
        }
    }

    fn check(&mut self, name: &BString) {
        assert!(
            !self.checker.contains(name),
            "bam file must be name collated or name sorted"
        );
        self.checker.insert(name.to_owned());
    }
}

impl<'a, R: std::io::Read> Iterator for NameCollatedRecords<'a, R> {
    type Item = (MultiMap<bam::Record>, MultiMap<bam::Record>);

    fn next(&mut self) -> Option<Self::Item> {
        let record = self.records.next()?.unwrap();
        let name = record.name().unwrap().to_owned();
        if let Some((prev_name, prev_record)) = self.prev_record.take() {
            if name == prev_name {
                Some((prev_record.into(), record.into()))
            } else {
                panic!(
                    "Expecting paired end reads with the same name, found {} and {}",
                    prev_name, name
                );
            }
        } else {
            self.check(&name);
            self.prev_record = Some((name, record));
            self.next()
        }
    }
}


// Helper function to get a region from a region ID and library spec
fn get_region_by_id(region_id: &str, library_spec: &Arc<RwLock<LibSpec>>) -> Option<Region> {
    let lib_spec = library_spec.read().ok()?;
    let region = lib_spec.get(region_id)?;
    let region = region.read().ok()?;
    Some(region.clone())
}

#[cfg(test)]
mod tests {
    use bwa_mem2::{AlignerOpts, BurrowsWheelerAligner, FMIndex, PairedEndStats};
    use env_logger;
    use std::sync::{Arc, RwLock};

    use super::*;

    #[test]
    fn test_seqspec_io() {
        // Initialize logging with debug level
        env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("debug"))
            .init();

        let bwa_index = "/data/Public/BWA_MEM2_index/GRCh38";
        let seqspec = "/data/kzhang/dev/PreCellar/test/seqspec.yaml";
        let spec = Assay::from_path(seqspec).unwrap();
        let mut aligner = BurrowsWheelerAligner::new(
            FMIndex::read(bwa_index).unwrap(),
            AlignerOpts::default(),
            PairedEndStats::default(),
        );
        let mut processor = FastqProcessor::new(spec).with_modality(Modality::ATAC);

        processor
            .gen_barcoded_alignments(&mut aligner, 4, 40000)
            .take(6)
            .for_each(|x| {
                println!("{:?}", x);
            });

        println!("{}", processor.get_report());
    }
    
    #[test]
    fn test_barcode_filtering_control() {
        // Create a test whitelist with no predefined barcodes
        let mut whitelist = Whitelist::empty();
        
        // Add a mix of high and low count barcodes
        for i in 0..20 {
            let count = 1000 - i * 50;
            let barcode = format!("CELL_{:03}", i).into_bytes();
            let quality = vec![b'F'; barcode.len()];
            
            for _ in 0..count {
                whitelist.count_barcode(&barcode, &quality);
            }
        }
        
        // Add background noise
        for i in 0..100 {
            let count = 10;
            let barcode = format!("BG_{:03}", i).into_bytes();
            let quality = vec![b'F'; barcode.len()];
            
            for _ in 0..count {
                whitelist.count_barcode(&barcode, &quality);
            }
        }
        
        // Check barcode count before filtering
        let initial_count = whitelist.num_seen_barcodes();
        
        // Create a copy for comparison
        let whitelist_no_filtering = whitelist.clone();
        
        // Apply the filtering
        let results = filter_cellular_barcodes_ordmag_advanced(
            &mut whitelist,
            Some(20),
            None,
            None,
            Some(0.99),
            Some(10)
        );
        
        // Verify filtering worked
        assert!(whitelist.num_seen_barcodes() < initial_count, 
            "Filtering should reduce the number of barcodes");
        assert!(whitelist.num_seen_barcodes() <= results.filtered_bcs,
            "Number of barcodes should match the filtering results");
        
        // Verify we can skip filtering by not calling the function
        assert_eq!(whitelist_no_filtering.num_seen_barcodes(), initial_count,
            "Whitelist without filtering should maintain all barcodes");
    }

    #[test]
    fn test_find_pattern() {
        // Basic pattern matching test
        let haystack = b"ACGTACGTACGT";
        let needle = b"ACGT";
        assert_eq!(find_pattern(haystack, needle), Some(0));
        assert_eq!(find_pattern(haystack, b"CGTA"), Some(1));
        assert_eq!(find_pattern(haystack, b"ACGT"), Some(0));
        assert_eq!(find_pattern(haystack, b"CGTX"), None);
        assert_eq!(find_pattern(haystack, b""), None);
        assert_eq!(find_pattern(&[], b"ACGT"), None);
    }
    
    // Note: The following phase block tests are no longer relevant since we removed phase block integration
    // They are kept here in commented form for historical reference
    /*
    #[test]
    fn test_phase_block_pattern_detection() {
        // Test removed - phase block functionality is no longer used
    }

    #[test]
    fn test_phase_block_adjustment() {
        // Test removed - phase block functionality is no longer used
    }

    #[test]
    fn test_multiple_barcodes_with_phase_block() {
        // Test removed - phase block functionality is no longer used
    }

    #[test]
    fn test_barcode_with_phase_block_adjustment() {
        // Test removed - phase block functionality is no longer used
    }
    */
}
