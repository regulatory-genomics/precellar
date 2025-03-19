use super::aligners::{Aligner, MultiMap, MultiMapR};

use crate::barcode::{
    filter_cellular_barcodes_ordmag_advanced, BarcodeCorrector, OligoFrequncy, Whitelist,
};
use crate::qc::{AlignQC, Metrics};
use crate::utils::{rev_compl, rev_compl_fastq_record};
use anyhow::{bail, Result};
use bstr::BString;
use indexmap::IndexMap;
use itertools::Itertools;
use kdam::{tqdm, BarExt};
use log::{debug, info, warn};
use noodles::{bam, fastq};
use rayon::iter::ParallelIterator;
use rayon::slice::ParallelSlice;
use seqspec::{Assay, FastqReader, Modality, Read, RegionId, SegmentInfo};
use serde_json;
use smallvec::SmallVec;
use std::collections::{HashMap, HashSet};
use std::path::PathBuf;

/// FastqProcessor manages the preprocessing of FASTQ files including barcode correction,
/// alignment, and QC metrics.
pub struct FastqProcessor {
    assay: Assay,                         // Specification of sequencing assay.
    current_modality: Option<Modality>, // Current sequencing modality being processed (e.g., RNA, ATAC).
    mito_dna: HashSet<String>, // Set of mitochondrial DNA sequence identifiers for special handling.
    metrics: HashMap<Modality, Metrics>, // Quality control metrics for each modality.
    align_qc: HashMap<Modality, AlignQC>, // Alignment QC data for each modality.
    barcode_correct_prob: f64, // if the posterior probability of a correction
    // exceeds this threshold, the barcode will be corrected.
    // cellrange uses 0.975 for ATAC and 0.9 for multiome.
    mismatch_in_barcode: usize, // The number of mismatches allowed in barcode
    expected_cells: Option<usize>, // Expected number of cells
    barcode_filtering_quantile: f64, // Quantile for barcode filtering
    barcode_bootstrap_samples: usize, // Number of bootstrap samples for filtering
    metrics_path: Option<std::path::PathBuf>, // Path to write QC metrics
}

impl FastqProcessor {
    /// Creates a new FastqProcessor with default settings.
    pub fn new(assay: Assay) -> Self {
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
            metrics_path: None,
        }
    }

    pub fn with_barcode_correct_prob(mut self, prob: f64) -> Self {
        self.barcode_correct_prob = prob;
        self
    }

    pub fn with_expected_cells(mut self, cells: usize) -> Self {
        self.expected_cells = Some(cells);
        self
    }

    pub fn with_barcode_filtering_params(
        mut self,
        quantile: f64,
        bootstrap_samples: usize,
        metrics_path: Option<&std::path::Path>,
    ) -> Self {
        self.barcode_filtering_quantile = quantile;
        self.barcode_bootstrap_samples = bootstrap_samples;
        self.metrics_path = metrics_path.map(|p| p.to_path_buf());
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
        self.current_modality = Some(modality);
        self
    }

    pub fn get_report(&self) -> Metrics {
        let mut metrics = self
            .metrics
            .get(&self.modality())
            .map_or(Metrics::default(), |x| x.clone());
        if let Some(align_qc) = self.align_qc.get(&self.modality()) {
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
        let fq_reader = self.gen_barcoded_fastq(true).with_chunk_size(chunk_size);

        let header = aligner.header();
        let mut qc = AlignQC::default();
        self.mito_dna.iter().for_each(|mito| {
            header
                .reference_sequences()
                .get_index_of(&BString::from(mito.as_str()))
                .map(|x| qc.mito_dna.insert(x));
        });
        let modality = self.modality();
        self.align_qc.insert(modality, qc);

        let total_reads = fq_reader.total_reads.unwrap_or(0);
        let mut progress_bar = tqdm!(total = total_reads);
        info!("Aligning {} reads...", total_reads);
        fq_reader.map(move |data| {
            let align_qc = self.align_qc.get_mut(&modality).unwrap();
            let results: Vec<_> = aligner.align_reads(num_threads, data);
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

        let whitelists = if correct_barcode {
            info!("Counting barcodes...");
            let whitelists = self.count_barcodes().unwrap();
            for (id, whitelist) in whitelists.iter() {
                info!(
                    "{:.2}% of sequences have an exact match in whitelist '{}'. Number of unique barcodes: {}.",
                    whitelist.frac_exact_match() * 100.0,
                    id,
                    whitelist.num_seen_barcodes(),
                );
            }
            whitelists
        } else {
            IndexMap::new()
        };

        let corrector = BarcodeCorrector::default()
            .with_max_missmatch(self.mismatch_in_barcode)
            .with_bc_confidence_threshold(self.barcode_correct_prob);

        let mut fq_reader: AnnotatedFastqReader = self
            .assay
            .get_segments_by_modality(modality)
            .filter_map(|(read, segment_info)| {
                let annotator =
                    FastqAnnotator::new(read, segment_info, &whitelists, corrector.clone())?;
                let reader = read.open()?;
                Some((annotator, reader))
            })
            .collect();

        if !whitelists.is_empty() {
            fq_reader.total_reads = Some(whitelists[0].total_count);
        }
        fq_reader
    }

    fn count_barcodes(&mut self) -> Result<IndexMap<RegionId, Whitelist>> {
        let modality = self.modality();
        let mut whitelists = self.get_whitelists();

        self.assay
            .get_segments_by_modality(modality)
            .for_each(|(read, segment_info)| {
                let is_reverse = read.is_reverse();
                if let Some(mut reader) = read.open() {
                    for fq in reader.records() {
                        let fq = fq.unwrap();
                        segment_info.split(&fq).unwrap().iter().for_each(|segment| {
                            if segment.is_barcode() {
                                let wl = whitelists.get_mut(segment.region_id()).expect(&format!(
                                    "whitelist not found for region {}",
                                    segment.region_id()
                                ));
                                if is_reverse {
                                    wl.count_barcode(
                                        &rev_compl(segment.seq),
                                        &segment.qual.iter().rev().copied().collect::<Vec<_>>(),
                                    );
                                } else {
                                    wl.count_barcode(segment.seq, segment.qual);
                                }
                            }
                        });
                    }
                }
            });

        self.metrics.entry(modality).or_default().insert(
            "frac_q30_bases_barcode".to_string(),
            whitelists.values().map(|x| x.frac_q30_bases()).sum::<f64>() / whitelists.len() as f64,
        );
        Ok(whitelists)
    }

    /*
    fn count_barcodes(&mut self) -> Result<IndexMap<RegionId, Whitelist>> {
        let mut whitelists = self.get_whitelists();
        debug!("Number of whitelists found: {}", whitelists.len());

        if whitelists.is_empty() {
            return Ok(whitelists);
        }

        let modality = self.modality();
        let mut total_filtered_bcs = 0;

        debug!(
            "Starting barcode counting process. Metrics path is: {:?}",
            self.metrics_path
        );

        self.assay
            .get_segments_by_modality(modality)
            .filter(|(_, segments)| segments.iter().any(|x| x.region_type.is_barcode()))
            .for_each(|(read, segment_info)| {
                let is_reverse = read.is_reverse();
                if let Some(mut reader) = read.open() {
                    for fq in reader.records() {
                        let fq = fq.unwrap();
                        segment_info.split(&fq).unwrap().iter().for_each(|segment| {
                            if segment.is_barcode() {
                                if let Some(wl) = whitelists.get_mut(segment.region_id()) {
                                    if is_reverse {
                                        wl.count_barcode(
                                            &rev_compl(segment.seq),
                                            &segment.qual.iter().rev().copied().collect::<Vec<_>>(),
                                        );
                                    } else {
                                        wl.count_barcode(segment.seq, segment.qual);
                                    }
                                }
                            }
                        });
                    }
                }
            });

        // Apply advanced filtering to each whitelist if needed
        for (id, whitelist) in whitelists.iter_mut() {
            // Check if this whitelist needs filtering (no predefined entries or has counts)
            let is_empty = whitelist.get_barcode_counts().is_empty();
            let has_predefined_whitelist = !is_empty;

            if !has_predefined_whitelist {
                info!(
                    "No predefined whitelist for '{}', applying order-of-magnitude filtering",
                    id
                );

                // Create metrics file path if metrics_path is provided
                let metrics_file_path = self.metrics_path.as_ref().map(|path| {
                    let file_name = format!("barcode_metrics_{}.json", id);
                    let full_path = path.join(file_name);
                    debug!(
                        "Creating metrics file path: {:?} for barcode region: {}",
                        full_path, id
                    );
                    full_path
                });

                debug!(
                    "Metrics path for barcode region {}: {:?}",
                    id, metrics_file_path
                );

                let results = filter_cellular_barcodes_ordmag_advanced(
                    whitelist,
                    self.expected_cells,
                    None,
                    None,
                    Some(self.barcode_filtering_quantile),
                    Some(self.barcode_bootstrap_samples),
                    metrics_file_path.as_deref(),
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
            } else {
                info!("Using predefined whitelist for '{}'", id);

                // Create metrics file for predefined whitelists too
                if let Some(metrics_path) = &self.metrics_path {
                    let file_name = format!("barcode_metrics_{}.json", id);
                    let full_path = metrics_path.join(file_name);
                    debug!("Creating metrics file path for predefined whitelist: {:?} for barcode region: {}", full_path, id);

                    // Create simple metrics structure
                    let mut metrics = serde_json::Map::new();
                    metrics.insert(
                        "whitelist_id".to_string(),
                        serde_json::Value::String(id.clone()),
                    );
                    metrics.insert("is_predefined".to_string(), serde_json::Value::Bool(true));
                    metrics.insert(
                        "total_barcodes".to_string(),
                        serde_json::Value::Number(serde_json::Number::from(
                            whitelist.get_barcode_counts().len(),
                        )),
                    );

                    // Add barcode-specific metrics
                    metrics.insert(
                        "frac_valid_barcode".to_string(),
                        serde_json::Value::Number(
                            serde_json::Number::from_f64(whitelist.frac_exact_match())
                                .unwrap_or(serde_json::Number::from(0)),
                        ),
                    );
                    metrics.insert(
                        "frac_q30_bases_barcode".to_string(),
                        serde_json::Value::Number(
                            serde_json::Number::from_f64(whitelist.frac_q30_bases())
                                .unwrap_or(serde_json::Number::from(0)),
                        ),
                    );

                    // Add general QC metrics from the processor's metrics if they exist
                    if let Some(modality_metrics) = self.metrics.get(&modality) {
                        // Add all existing metrics for this modality
                        for (metric_name, metric_value) in modality_metrics.iter() {
                            metrics.insert(
                                metric_name.clone(),
                                serde_json::Value::Number(
                                    serde_json::Number::from_f64(*metric_value)
                                        .unwrap_or(serde_json::Number::from(0)),
                                ),
                            );
                        }
                    }

                    // Write metrics to file
                    if let Ok(file) = std::fs::File::create(&full_path) {
                        if let Err(e) = serde_json::to_writer_pretty(file, &metrics) {
                            warn!(
                                "Failed to write metrics file for predefined whitelist '{}': {}",
                                id, e
                            );
                        } else {
                            debug!(
                                "Successfully wrote metrics for predefined whitelist '{}'",
                                id
                            );
                        }
                    } else {
                        warn!(
                            "Failed to create metrics file for predefined whitelist '{}'",
                            id
                        );
                    }
                }

                // Count predefined barcodes toward the total
                total_filtered_bcs += whitelist.get_barcode_counts().len();
            }
        }

        // After processing all barcodes, write a consolidated metrics file
        if let Some(metrics_path) = &self.metrics_path {
            // Create a consolidated metrics file with all QC metrics
            let consolidated_file_path =
                PathBuf::from(metrics_path).join("consolidated_qc_metrics.json");
            debug!(
                "Creating consolidated QC metrics file at: {:?}",
                consolidated_file_path
            );

            if let Some(modality_metrics) = self.metrics.get(&modality) {
                let mut consolidated_metrics = serde_json::Map::new();

                // Add all metrics for this modality
                for (metric_name, metric_value) in modality_metrics.iter() {
                    consolidated_metrics.insert(
                        metric_name.clone(),
                        serde_json::Value::Number(
                            serde_json::Number::from_f64(*metric_value)
                                .unwrap_or(serde_json::Number::from(0)),
                        ),
                    );
                }

                // Add total filtered barcodes count
                consolidated_metrics.insert(
                    "total_filtered_bcs".to_string(),
                    serde_json::Value::Number(serde_json::Number::from(total_filtered_bcs)),
                );

                // Add average barcode quality metrics across all whitelists
                consolidated_metrics.insert(
                    "frac_q30_bases_barcode".to_string(),
                    serde_json::Value::Number(
                        serde_json::Number::from_f64(
                            whitelists.values().map(|x| x.frac_q30_bases()).sum::<f64>()
                                / whitelists.len() as f64,
                        )
                        .unwrap_or(serde_json::Number::from(0)),
                    ),
                );

                consolidated_metrics.insert(
                    "frac_valid_barcode".to_string(),
                    serde_json::Value::Number(
                        serde_json::Number::from_f64(
                            whitelists
                                .values()
                                .map(|x| x.frac_exact_match())
                                .sum::<f64>()
                                / whitelists.len() as f64,
                        )
                        .unwrap_or(serde_json::Number::from(0)),
                    ),
                );

                // Write the consolidated metrics to file
                if let Ok(file) = std::fs::File::create(&consolidated_file_path) {
                    if let Err(e) = serde_json::to_writer_pretty(file, &consolidated_metrics) {
                        warn!("Failed to write consolidated metrics file: {}", e);
                    } else {
                        debug!(
                            "Successfully wrote consolidated QC metrics to {:?}",
                            consolidated_file_path
                        );
                    }
                } else {
                    warn!(
                        "Failed to create consolidated metrics file at {:?}",
                        consolidated_file_path
                    );
                }
            }
        }

        if total_filtered_bcs > 0 {
            self.metrics
                .entry(modality)
                .or_default()
                .insert("total_filtered_bcs".to_string(), total_filtered_bcs as f64);
        }

        self.metrics.entry(modality).or_default().insert(
            "frac_q30_bases_barcode".to_string(),
            whitelists.values().map(|x| x.frac_q30_bases()).sum::<f64>() / whitelists.len() as f64,
        );

        Ok(whitelists)
    }
    */

    fn get_whitelists(&self) -> IndexMap<RegionId, Whitelist> {
        let regions = self
            .assay
            .library_spec
            .get_modality(&self.modality())
            .unwrap()
            .read()
            .unwrap();
        regions
            .subregions
            .iter()
            .filter_map(|r| {
                let r = r.read().unwrap();
                if r.region_type.is_barcode() {
                    let id = r.region_id.to_string();
                    let list = if let Some(onlist) = r.onlist.as_ref() {
                        Whitelist::new(onlist.read().unwrap())
                    } else {
                        Whitelist::empty()
                    };
                    Some((id, list))
                } else {
                    None
                }
            })
            .collect()
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
                    .segment_info
                    .iter()
                    .filter(|info| info.region_type.is_barcode())
                    .map(|info| (info.region_id.as_str(), info.len.len()))
            })
            .collect()
    }

    pub fn get_all_umi(&self) -> Vec<(&str, usize)> {
        self.annotators
            .iter()
            .flat_map(|annotator| {
                annotator
                    .segment_info
                    .iter()
                    .filter(|info| info.region_type.is_umi())
                    .map(|info| (info.region_id.as_str(), info.len.len()))
            })
            .collect()
    }

    pub fn is_paired_end(&self) -> bool {
        let mut has_read1 = false;
        let mut has_read2 = false;
        self.annotators.iter().for_each(|x| {
            x.segment_info.iter().for_each(|info| {
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

        while accumulated_length < self.chunk_size {
            let mut max_read = 0;
            let mut min_read = usize::MAX;
            let records: SmallVec<[_; 4]> = self
                .readers
                .iter_mut()
                .flat_map(|reader| {
                    let n = reader
                        .read_record(&mut self.buffer)
                        .expect("error reading fastq record");
                    min_read = min_read.min(n);
                    max_read = max_read.max(n);
                    if n > 0 {
                        accumulated_length += self.buffer.sequence().len();
                        Some(self.buffer.clone())
                    } else {
                        None
                    }
                })
                .collect();
            if max_read == 0 {
                break;
            } else if min_read == 0 {
                panic!("Unequal number of reads in the chunk");
            } else {
                assert!(
                    records.iter().map(|r| r.name()).all_equal(),
                    "read names mismatch"
                );
                self.chunk.push(records);
            }
        }

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
        if n == 0 {
            None
        } else {
            let n = (n / 256).max(256);
            let annotators = &self.annotators;
            let result: Vec<AnnotatedFastq> = self
                .chunk
                .par_chunks(n)
                .flat_map_iter(|chunk| {
                    chunk.into_iter().map(move |records| {
                        records
                            .iter()
                            .enumerate()
                            .map(|(i, record)| annotators[i].annotate(record).unwrap())
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
    is_reverse: bool,
    segment_info: SegmentInfo,
    min_len: usize,
    max_len: usize,
}

impl FastqAnnotator {
    pub fn new(
        read: &Read,
        segment_info: SegmentInfo,
        whitelists: &IndexMap<String, Whitelist>,
        corrector: BarcodeCorrector,
    ) -> Option<Self> {
        if !segment_info.iter().any(|x| {
            x.region_type.is_barcode() || x.region_type.is_umi() || x.region_type.is_target()
        }) {
            None
        } else {
            let whitelists = segment_info
                .iter()
                .flat_map(|segment| {
                    let v = whitelists.get(&segment.region_id)?;
                    Some((segment.region_id.clone(), v.get_barcode_counts().clone()))
                })
                .collect();
            let anno = Self {
                whitelists,
                corrector,
                is_reverse: read.is_reverse(),
                segment_info,
                min_len: read.min_len as usize,
                max_len: read.max_len as usize,
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

        if let Some(segments) = self.segment_info.split(record) {
            segments.into_iter().for_each(|segment| {
                if segment.is_barcode() || segment.is_umi() {
                    let mut fq = segment.into_fq(record.definition());
                    if self.is_reverse {
                        fq = rev_compl_fastq_record(fq);
                    }

                    if segment.is_barcode() {
                        let corrected = self.whitelists.get(segment.region_id()).map_or(
                            Some(fq.sequence().to_vec()),
                            |counts| {
                                self.corrector
                                    .correct(counts, fq.sequence(), fq.quality_scores())
                                    .ok()
                                    .map(|x| x.to_vec())
                            },
                        );

                        if let Some(bc) = &mut barcode {
                            bc.extend(&Barcode { raw: fq, corrected });
                        } else {
                            barcode = Some(Barcode { raw: fq, corrected });
                        }
                    } else {
                        umi = Some(fq);
                    }
                } else if segment.contains_target() {
                    if read1.is_some() || read2.is_some() {
                        panic!("Both Read1 and Read2 are set");
                    } else {
                        let fq = segment.into_fq(record.definition());
                        // TODO: polyA and adapter trimming
                        if self.is_reverse {
                            read2 = Some(fq);
                        } else {
                            read1 = Some(fq);
                        }
                    }
                }
            });
        }

        Ok(AnnotatedFastq {
            barcode,
            umi,
            read1,
            read2,
        })
    }
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

#[cfg(test)]
mod tests {
    use super::*;

    use std::io::{BufRead, Write};

    fn show_fq(fq: &AnnotatedFastq) -> String {
        format!(
            "{}\t{}\t{}\t{}",
            fq.barcode
                .as_ref()
                .map_or("", |x| std::str::from_utf8(x.raw.sequence()).unwrap()),
            fq.umi
                .as_ref()
                .map_or("", |x| std::str::from_utf8(x.sequence()).unwrap()),
            fq.read1
                .as_ref()
                .map_or("", |x| std::str::from_utf8(x.sequence()).unwrap()),
            fq.read2
                .as_ref()
                .map_or("", |x| std::str::from_utf8(x.sequence()).unwrap())
        )
    }

    fn test_fq(input: &str, output: &str) {
        let assay = Assay::from_path(input).unwrap();
        let mut fq_proc = FastqProcessor::new(assay).with_modality(Modality::RNA);
        let file = std::fs::File::open(output).unwrap();
        let reader = std::io::BufReader::new(flate2::read::GzDecoder::new(file));
        for (fq, line) in fq_proc
            .gen_barcoded_fastq(false)
            .flatten()
            .zip(reader.lines())
        {
            assert_eq!(show_fq(&fq), line.unwrap());
        }
    }

    #[test]
    fn test_io() {
        let seqspec = "data/test4.yaml";
        let assay = Assay::from_path(seqspec).unwrap();
        let mut fq_proc = FastqProcessor::new(assay).with_modality(Modality::RNA);
        // open a file
        let file = std::fs::File::create("test.out").unwrap();
        let mut writer = std::io::BufWriter::new(file);
        for fq in fq_proc.gen_barcoded_fastq(false).flatten() {
            writeln!(writer, "{}", show_fq(&fq)).unwrap();
        }
    }

    #[test]
    fn test_fastq() {
        test_fq("data/test1.yaml", "data/test1.out.gz");
        test_fq("data/test2.yaml", "data/test2.out.gz");
        test_fq("data/test3.yaml", "data/test3.out.gz");
        test_fq("data/test4.yaml", "data/test4.out.gz");
    }
}
