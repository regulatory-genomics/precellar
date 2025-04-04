use super::aligners::{Aligner, MultiMap, MultiMapR};

use crate::barcode::{BarcodeCorrector, OligoFrequncy, Whitelist};
use crate::qc::{QcAlign, QcFastq};
use crate::utils::{rev_compl, rev_compl_fastq_record};
use anyhow::Result;
use bstr::BString;
use indexmap::IndexMap;
use indicatif::{ProgressBar, ProgressDrawTarget, ProgressIterator, ProgressStyle};
use itertools::Itertools;
use log::{debug, info, warn};
use noodles::{bam, fastq};
use rayon::iter::ParallelIterator;
use rayon::slice::ParallelSlice;
use seqspec::{
    Assay, FastqReader, Modality, Read, RegionId, SegmentInfo, SequenceType, SplitError,
};
use smallvec::SmallVec;
use std::collections::{HashMap, HashSet};
use std::sync::{Arc, Mutex};

/// FastqProcessor manages the preprocessing of FASTQ files including barcode correction,
/// alignment, and QC metrics.
pub struct FastqProcessor {
    assay: Assay,                       // Specification of sequencing assay.
    current_modality: Option<Modality>, // Current sequencing modality being processed (e.g., RNA, ATAC).
    mito_dna: HashSet<String>, // Set of mitochondrial DNA sequence identifiers for special handling.
    barcode_correct_prob: f64, // if the posterior probability of a correction
    // exceeds this threshold, the barcode will be corrected.
    // cellrange uses 0.975 for ATAC and 0.9 for multiome.
    mismatch_in_barcode: usize, // The number of mismatches allowed in barcode
    qc_align: HashMap<Modality, Arc<Mutex<QcAlign>>>,
    qc_fastq: HashMap<Modality, Arc<Mutex<QcFastq>>>,
}

impl FastqProcessor {
    /// Creates a new FastqProcessor with default settings.
    pub fn new(assay: Assay) -> Self {
        Self {
            assay,
            current_modality: None,
            mito_dna: HashSet::new(),
            barcode_correct_prob: 0.975,
            mismatch_in_barcode: 1,
            qc_align: HashMap::new(),
            qc_fastq: HashMap::new(),
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

    pub fn get_fastq_qc(&self) -> &Arc<Mutex<QcFastq>> {
        self.qc_fastq
            .get(&self.modality())
            .expect("fastq qc not found")
    }

    pub fn get_align_qc(&self) -> &Arc<Mutex<QcAlign>> {
        self.qc_align
            .get(&self.modality())
            .expect("align qc not found")
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
    ) -> AlignmentResult<'a, A> {
        let fq_reader = self.gen_barcoded_fastq(true).with_chunk_size(chunk_size);

        // Initialize qc
        let mut qc = QcAlign::default();
        let header = aligner.header();
        self.mito_dna.iter().for_each(|mito| {
            header
                .reference_sequences()
                .get_index_of(&BString::from(mito.as_str()))
                .map(|x| qc.mito_dna.insert(x));
        });
        self.qc_align
            .insert(self.modality(), Arc::new(Mutex::new(qc)));

        info!("Aligning reads to reference genome...");
        AlignmentResult {
            aligner,
            fastq_reader: fq_reader,
            qc: self.get_align_qc().clone(),
            header,
            num_threads,
        }
    }

    pub fn gen_barcoded_fastq(&mut self, correct_barcode: bool) -> AnnotatedFastqReader {
        let modality = self.modality();
        // Initialize qc
        self.qc_fastq
            .insert(modality.clone(), Arc::new(Mutex::new(QcFastq::default())));

        let corrector = if correct_barcode {
            Some(
                BarcodeCorrector::default()
                    .with_max_missmatch(self.mismatch_in_barcode)
                    .with_bc_confidence_threshold(self.barcode_correct_prob),
            )
        } else {
            None
        };

        let (mut whitelists, num_reads) = self.count_barcodes().unwrap();
        for (id, whitelist) in whitelists.iter_mut() {
            if whitelist.len() > 0 {
                info!(
                    "{:.2}% of sequences have an exact match in whitelist '{}'. Number of unique barcodes: {}.",
                    whitelist.frac_exact_match() * 100.0,
                    id,
                    whitelist.num_seen_barcodes(),
                );
            } else if self
                .assay
                .library_spec
                .get(id)
                .unwrap()
                .read()
                .unwrap()
                .sequence_type
                == SequenceType::Onlist
            {
                whitelist.predict_whitelist();
            }
        }

        let readers =
            self.assay
                .get_segments_by_modality(modality)
                .filter_map(|(read, segment_info)| {
                    let annotator =
                        FastqAnnotator::new(read, segment_info, &whitelists, corrector.clone())?;
                    let reader = read.open()?;
                    Some((annotator, reader))
                });
        let mut fq_reader = AnnotatedFastqReader::new(readers, self.get_fastq_qc().clone());

        fq_reader.total_reads = Some(num_reads);

        fq_reader
    }

    fn count_barcodes(&self) -> Result<(IndexMap<RegionId, Whitelist>, usize)> {
        let spinner = ProgressBar::with_draw_target(None, ProgressDrawTarget::stderr_with_hz(1))
        .with_style(
            ProgressStyle::with_template(
                "{spinner} Processed {human_pos} reads in {elapsed} ({per_sec}) ...",
            )
            .unwrap(),
        );

        let modality = self.modality();
        let mut whitelists = self.get_whitelists();

        let mut num_reads = 0;
        self.assay
            .get_segments_by_modality(modality)
            .filter(|(_, info)| info.iter().any(|x| x.is_barcode()))
            .for_each(|(read, segment_info)| {
                let is_reverse = read.is_reverse();
                if let Some(mut reader) = read.open() {
                    info!("Counting barcodes in read {}...", read.read_id);
                    num_reads = 0;
                    for fq in reader.records().progress_with(spinner.clone())
                    {
                        num_reads += 1;
                        if let Ok(segments) = segment_info.split(&fq.unwrap()) {
                            segments.iter().for_each(|segment| {
                                if segment.is_barcode() {
                                    let wl =
                                        whitelists.get_mut(segment.region_id()).expect(&format!(
                                            "whitelist not found for region {}",
                                            segment.region_id()
                                        ));
                                    if is_reverse {
                                        wl.count_barcode(&rev_compl(segment.seq));
                                    } else {
                                        wl.count_barcode(segment.seq);
                                    }
                                }
                            })
                        }
                    }
                }
            });

        Ok((whitelists, num_reads))
    }

    fn get_whitelists(&self) -> IndexMap<RegionId, Whitelist> {
        let regions = self
            .assay
            .library_spec
            .get_modality(&self.modality())
            .expect(&format!("modality not found: {}", self.modality()))
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
                        if r.sequence_type == SequenceType::Onlist {
                            warn!("Barcode region '{}' does not have a whitelist", id);
                        }
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

pub struct AlignmentResult<'a, A> {
    aligner: &'a mut A,
    pub fastq_reader: AnnotatedFastqReader,
    qc: Arc<Mutex<QcAlign>>,
    header: noodles::sam::Header,
    num_threads: u16,
}

impl<'a, A: Aligner> Iterator for AlignmentResult<'a, A> {
    type Item = Vec<(Option<MultiMapR>, Option<MultiMapR>)>;
    fn next(&mut self) -> Option<Self::Item> {
        let data = self.fastq_reader.next()?;

        let results: Vec<_> = self.aligner.align_reads(self.num_threads, data);
        let mut qc = self.qc.lock().unwrap();
        results.iter().for_each(|ali| match ali {
            (Some(ali1), Some(ali2)) => {
                qc.add_pair(&self.header, ali1, ali2).unwrap();
            }
            (Some(ali1), None) => {
                qc.add_read1(&self.header, ali1).unwrap();
            }
            (None, Some(ali2)) => {
                qc.add_read2(&self.header, ali2).unwrap();
            }
            _ => {
                debug!("No alignment found for read");
            }
        });
        Some(results)
    }
}

pub struct AnnotatedFastqReader {
    buffer: fastq::Record,
    pub total_reads: Option<usize>,
    pub num_processed: usize,
    trim_poly_a: bool,
    annotators: Vec<FastqAnnotator>,
    readers: Vec<FastqReader>,
    chunk_size: usize,
    chunk: Vec<SmallVec<[fastq::Record; 4]>>,
    qc: Arc<Mutex<QcFastq>>,
}

impl AnnotatedFastqReader {
    fn new<T: IntoIterator<Item = (FastqAnnotator, FastqReader)>>(
        iter: T,
        qc: Arc<Mutex<QcFastq>>,
    ) -> Self {
        let (annotators, readers): (Vec<_>, Vec<_>) = iter.into_iter().unzip();
        let chunk = Vec::new();
        Self {
            buffer: fastq::Record::default(),
            total_reads: None,
            num_processed: 0,
            annotators,
            readers,
            trim_poly_a: false,
            chunk_size: 10000000,
            chunk,
            qc,
        }
    }

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

impl Iterator for AnnotatedFastqReader {
    type Item = Vec<AnnotatedFastq>;

    fn next(&mut self) -> Option<Self::Item> {
        let n = self.read_chunk();
        if n == 0 {
            None
        } else {
            self.num_processed += n;
            let n = (n / 256).max(1024);
            let annotators = &self.annotators;
            let (result, qc): (Vec<_>, Vec<_>) = self
                .chunk
                .par_chunks(n)
                .map(|chunk| process_chunk(&annotators, chunk))
                .unzip();

            self.qc.lock().unwrap().extend(qc);
            Some(result.into_iter().flatten().collect())
        }
    }
}

fn process_chunk<'a, I: IntoIterator<Item = &'a SmallVec<[fastq::Record; 4]>>>(
    annotators: &[FastqAnnotator],
    chunk: I,
) -> (Vec<AnnotatedFastq>, QcFastq) {
    let mut qc = QcFastq::default();
    let annotated = chunk
        .into_iter()
        .flat_map(|records| {
            let fq = records
                .iter()
                .enumerate()
                .flat_map(|(i, record)| {
                    let annotator = &annotators[i];
                    let id = &annotator.read_id;
                    *qc.num_reads.entry(id.clone()).or_insert(0) += 1;
                    if let Ok(anno) = annotator.annotate(record) {
                        Some(anno)
                    } else {
                        *qc.num_defect.entry(id.clone()).or_insert(0) += 1;
                        None
                    }
                })
                .reduce(|mut this, other| {
                    this.join(other);
                    this
                })?;
            qc.update(&fq);
            if fq.barcode.is_none() {
                None
            } else {
                Some(fq)
            }
        })
        .collect();
    (annotated, qc)
}

/// A FastqAnnotator that splits the reads into subregions, e.g., barcode, UMI, and
/// return annotated reads.
#[derive(Debug)]
struct FastqAnnotator {
    read_id: String,
    whitelists: IndexMap<String, OligoFrequncy>,
    corrector: Option<BarcodeCorrector>,
    is_reverse: bool,
    segment_info: SegmentInfo,
}

impl FastqAnnotator {
    pub fn new(
        read: &Read,
        segment_info: SegmentInfo,
        whitelists: &IndexMap<String, Whitelist>,
        corrector: Option<BarcodeCorrector>,
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
                read_id: read.read_id.to_string(),
                whitelists,
                corrector,
                is_reverse: read.is_reverse(),
                segment_info,
            };
            Some(anno)
        }
    }

    fn annotate(&self, record: &fastq::Record) -> Result<AnnotatedFastq, SplitError> {
        let mut barcode: Option<Barcode> = None;
        let mut umi = None;
        let mut read1 = None;
        let mut read2 = None;

        let segments = self.segment_info.split(record)?;
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
                            if let Some(corrector) = self.corrector.as_ref() {
                                corrector
                                    .correct(counts, fq.sequence(), fq.quality_scores())
                                    .ok()
                                    .map(|x| x.to_vec())
                            } else {
                                Some(fq.sequence().to_vec())
                            }
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

    use std::io::BufRead;

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
        let modality = assay.modalities[0].clone();
        let mut fq_proc = FastqProcessor::new(assay).with_modality(modality);
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
        let modality = assay.modalities[0].clone();
        let mut fq_proc = FastqProcessor::new(assay).with_modality(modality);
        // open a file
        for fq in fq_proc.gen_barcoded_fastq(false).flatten() {
            println!("{}", show_fq(&fq));
        }
    }

    #[test]
    fn test_fastq() {
        test_fq("data/test2.yaml", "data/test2.out.gz");
        test_fq("data/test3.yaml", "data/test3.out.gz");
        test_fq("data/test4.yaml", "data/test4.out.gz");
    }
}
