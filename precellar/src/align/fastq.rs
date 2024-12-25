use super::aligners::{Aligner, MultiMap, MultiMapR};

use crate::adapter::trim_poly_nucleotide;
use crate::barcode::{BarcodeCorrector, OligoFrequncy, Whitelist};
use crate::qc::{AlignQC, Metrics};
use anyhow::{bail, Result};
use bstr::BString;
use indexmap::IndexMap;
use kdam::{tqdm, BarExt};
use log::info;
use noodles::{bam, fastq};
use rayon::iter::ParallelIterator;
use rayon::slice::ParallelSlice;
use seqspec::{Assay, FastqReader, Modality, Read, RegionId, SegmentInfo, SegmentInfoElem};
use smallvec::SmallVec;
use std::collections::{HashMap, HashSet};

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
}

impl FastqProcessor {
    pub fn new(assay: Assay) -> Self {
        Self {
            assay,
            current_modality: None,
            metrics: HashMap::new(),
            align_qc: HashMap::new(),
            mito_dna: HashSet::new(),
            barcode_correct_prob: 0.975,
            mismatch_in_barcode: 1,
        }
    }

    pub fn with_barcode_correct_prob(mut self, prob: f64) -> Self {
        self.barcode_correct_prob = prob;
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
    ) -> impl Iterator<Item = Vec<(MultiMapR, Option<MultiMapR>)>> + 'a {
        let fq_reader = self.gen_barcoded_fastq(true).with_chunk_size(chunk_size);
        let total_reads = fq_reader.total_reads.unwrap_or(0);

        let modality = self.modality();
        info!("Aligning {} reads...", total_reads);
        let header = aligner.header();
        let mut qc = AlignQC::default();
        self.mito_dna.iter().for_each(|mito| {
            header
                .reference_sequences()
                .get_index_of(&BString::from(mito.as_str()))
                .map(|x| qc.mito_dna.insert(x));
        });
        self.align_qc.insert(modality, qc);

        let mut progress_bar = tqdm!(total = total_reads);
        fq_reader.map(move |data| {
            let align_qc = self.align_qc.get_mut(&modality).unwrap();
            let results: Vec<_> = aligner.align_reads(num_threads, data);
            results
                .iter()
                .for_each(|(ali1, ali2)| align_qc.add(&header, ali1, ali2.as_ref()).unwrap());
            progress_bar.update(results.len()).unwrap();
            results
        })
    }

    pub fn gen_barcoded_fastq(&mut self, correct_barcode: bool) -> AnnotatedFastqReader {
        let modality = self.modality();

        let whitelists = if correct_barcode {
            info!("Counting barcodes...");
            // let whitelist = self.count_barcodes().unwrap();
            let whitelists = self.count_barcodes().unwrap();
            for (id, whitelist) in whitelists.iter() {
                info!(
                    "Found {} barcodes. {:.2}% of them have an exact match in whitelist {}",
                    whitelist.total_count,
                    whitelist.frac_exact_match() * 100.0,
                    id
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
            .filter_map(|(read, index)| {
                let annotator = FastqAnnotator::new(read, index, &whitelists, corrector.clone())?;
                Some((annotator, read.open().unwrap()))
            })
            .collect();
        if !whitelists.is_empty() {
            fq_reader.total_reads = Some(whitelists[0].total_count);
        }
        fq_reader
    }

    fn count_barcodes(&mut self) -> Result<IndexMap<RegionId, Whitelist>> {
        let modality = self.modality();
        let mut whitelists = self.get_whitelists()?;

        fn count(
            read: &Read,
            barcode_region_index: SegmentInfo,
            whitelist: &mut Whitelist,
        ) -> Result<()> {
            let range = &barcode_region_index
                .segments
                .iter()
                .find(|x| x.region_type.is_barcode())
                .unwrap()
                .range;
            read.open().unwrap().records().for_each(|record| {
                let mut record = record.unwrap();
                record = slice_fastq_record(&record, range.start as usize, range.end as usize);
                if read.is_reverse() {
                    record = rev_compl_fastq_record(record);
                }
                whitelist.count_barcode(record.sequence(), record.quality_scores());
            });
            Ok(())
        }

        for (i, (read, barcode_region_index)) in self
            .assay
            .get_segments_by_modality(modality)
            .filter(|(_, region_index)| {
                region_index
                    .segments
                    .iter()
                    .any(|x| x.region_type.is_barcode())
            })
            .enumerate()
        {
            count(read, barcode_region_index, &mut whitelists[i])?;
        }

        self.metrics.entry(modality).or_default().insert(
            "frac_q30_bases_barcode".to_string(),
            whitelists.values().map(|x| x.frac_q30_bases()).sum::<f64>() / whitelists.len() as f64,
        );
        Ok(whitelists)
    }

    fn get_whitelists(&self) -> Result<IndexMap<RegionId, Whitelist>> {
        let regions = self
            .assay
            .library_spec
            .get_modality(&self.modality())
            .unwrap()
            .read()
            .map_err(|_| anyhow::anyhow!("Cannot obtain lock"))?;
        regions
            .subregions
            .iter()
            .filter_map(|r| {
                let r = r.read().unwrap();
                if r.region_type.is_barcode() {
                    if let Some(onlist) = r.onlist.as_ref() {
                        let list = onlist
                            .read()
                            .map(|list| (r.region_id.to_string(), Whitelist::new(list)));
                        Some(list)
                    } else {
                        None
                    }
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

    fn read_chunk(&mut self) -> usize {
        self.chunk.clear();

        let mut accumulated_length = 0;

        while accumulated_length < self.chunk_size {
            let mut max_read = 0;
            let mut min_read = usize::MAX;
            let records = self.readers.iter_mut().flat_map(|reader| {
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
            }).collect();
            if max_read == 0 {
                break;
            } else if min_read == 0 {
                panic!("Unequal number of reads in the chunk");
            } else {
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
            let result = self
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
                            .unwrap()
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
    id: String,
    is_reverse: bool,
    subregions: Vec<SegmentInfoElem>,
    min_len: usize,
    max_len: usize,
}

impl FastqAnnotator {
    pub fn new(
        read: &Read,
        index: SegmentInfo,
        whitelists: &IndexMap<String, Whitelist>,
        corrector: BarcodeCorrector,
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
        self.subregions.iter().for_each(|info| {
            let mut fq =
                slice_fastq_record(record, info.range.start as usize, info.range.end as usize);
            if self.is_reverse && (info.region_type.is_barcode() || info.region_type.is_umi()) {
                fq = rev_compl_fastq_record(fq);
            }
            if info.region_type.is_umi() {
                umi = Some(fq);
            } else if info.region_type.is_barcode() {
                let corrected = self.whitelists.get(&info.region_id).map_or(
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
            } else if info.region_type.is_target() {
                if read1.is_some() || read2.is_some() {
                    panic!("Both Read1 and Read2 are set");
                } else {
                    if let Some(nucl) = info.region_type.poly_nucl() {
                        if let Some(idx) = trim_poly_nucleotide(nucl, fq.sequence().iter().copied())
                        {
                            fq = slice_fastq_record(&fq, idx, fq.sequence().len());
                        }
                    }
                    // Only keep reads with length >= 8
                    if fq.sequence().len() >= 8 {
                        if self.is_reverse {
                            read2 = Some(fq);
                        } else {
                            read1 = Some(fq);
                        }
                    }
                }
            }
        });
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

fn slice_fastq_record(record: &fastq::Record, start: usize, end: usize) -> fastq::Record {
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

fn rev_compl_fastq_record(mut record: fastq::Record) -> fastq::Record {
    *record.quality_scores_mut() = record.quality_scores().iter().rev().copied().collect();
    *record.sequence_mut() = record
        .sequence()
        .iter()
        .rev()
        .map(|&x| match x {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            _ => x,
        })
        .collect();
    record
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
    use bwa_mem2::{AlignerOpts, BurrowsWheelerAligner, FMIndex, PairedEndStats};

    use super::*;

    #[test]
    fn test_seqspec_io() {
        let spec = Assay::from_path("tests/data/spec.yaml").unwrap();
        let mut aligner = BurrowsWheelerAligner::new(
            FMIndex::read("tests/data/hg38").unwrap(),
            AlignerOpts::default(),
            PairedEndStats::default(),
        );
        let mut processor = FastqProcessor::new(spec).with_modality(Modality::ATAC);

        processor
            .gen_barcoded_alignments(&mut aligner, 8, 50000)
            .take(6)
            .for_each(|x| {
                println!("{:?}", x);
            });

        println!("{}", processor.get_report());
    }
}
