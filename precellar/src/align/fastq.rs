use super::aligners::{Aligner, MultiMapR};

use crate::barcode::{BarcodeAnalyzer, BarcodeCorrectOptions};
use crate::pipeline::{FastqStage, FastqStagePipeline, MiddlewareQcReport};
use crate::qc::{QcAlign, QcFastq};
use crate::utils::{rev_compl_fastq_record, PrefetchIterator};
use anyhow::Result;
use bstr::BString;
use itertools::Itertools;
use log::{debug, info};
use noodles::fastq;
use rayon::iter::ParallelIterator;
use rayon::slice::ParallelSlice;
use seqspec::{Assay, FastqReader, Modality, SegmentInfo, SplitError};
use smallvec::SmallVec;
use std::collections::HashSet;
use std::sync::Arc;

/// Configuration for standard barcode correction during annotation.
#[derive(Debug, Clone)]
pub struct BarcodeCorrectionConfig {
    pub confidence_threshold: f64,
    pub max_mismatch: usize,
}

impl Default for BarcodeCorrectionConfig {
    fn default() -> Self {
        Self {
            confidence_threshold: 0.975,
            max_mismatch: 1,
        }
    }
}

/// Immutable configuration used to construct annotated FASTQ streams.
pub struct FastqPlan {
    assays: Arc<[Assay]>,
    modality: Modality,
    barcode_config: BarcodeCorrectionConfig,
    stages: FastqStagePipeline,
}

impl FastqPlan {
    pub fn new(assays: Vec<Assay>, modality: Modality) -> Self {
        Self {
            assays: assays.into(),
            modality,
            barcode_config: BarcodeCorrectionConfig::default(),
            stages: FastqStagePipeline::default(),
        }
    }

    pub fn with_barcode_config(mut self, config: BarcodeCorrectionConfig) -> Self {
        self.barcode_config = config;
        self
    }

    pub fn with_stage<S>(mut self, stage: S) -> Self
    where
        S: FastqStage + 'static,
    {
        self.stages.push_stage(stage);
        self
    }

    pub fn build(self, correct_barcode: bool, chunk_size: usize) -> Result<FastqExecution> {
        let reader = self.build_reader(correct_barcode, chunk_size);
        let num_records = reader.num_records();
        let paired_end = reader.is_paired_end()?;
        let n_reads = reader
            .readers
            .iter()
            .map(|r| indicatif::HumanCount(r.annotation.num_reads() as u64).to_string())
            .intersperse(" + ".to_string())
            .collect();
        Ok(FastqExecution {
            source: reader,
            stages: self.stages,
            qc: QcFastq::default(),
            num_records,
            read_summary: n_reads,
            paired_end,
            finished: false,
        })
    }

    fn build_reader(&self, correct_barcode: bool, chunk_size: usize) -> MultiAnnotatedFqReader {
        let num_assays = self.assays.len();
        let readers = self
            .assays
            .iter()
            .enumerate()
            .map(|(i, assay)| {
                if num_assays > 1 {
                    info!(">>>Processing assay {}/{}<<<", i + 1, num_assays);
                }

                let mut barcode_analyzer = BarcodeAnalyzer::new(assay, self.modality);
                barcode_analyzer.summary();
                if correct_barcode {
                    barcode_analyzer.barcode_correct_options = Some(BarcodeCorrectOptions {
                        bc_confidence_threshold: self.barcode_config.confidence_threshold,
                        max_mismatch: self.barcode_config.max_mismatch,
                        ..Default::default()
                    });
                }

                let readers = assay.get_segments_by_modality(self.modality).filter_map(
                    |(read, segment_info)| {
                        let annotator = FastqAnnotator::new(&read.read_id, segment_info)?;
                        let reader = read.open()?;
                        Some((annotator, reader))
                    },
                );
                AnnotatedFastqReader::new(readers, barcode_analyzer, chunk_size)
            })
            .collect();

        MultiAnnotatedFqReader::new(readers)
    }
}

#[derive(Debug)]
pub struct FastqReport {
    pub fastq: QcFastq,
    pub middleware: Vec<MiddlewareQcReport>,
}

/// One execution of a prepared FASTQ plan.
pub struct FastqExecution {
    source: MultiAnnotatedFqReader,
    stages: FastqStagePipeline,
    qc: QcFastq,
    num_records: usize,
    read_summary: String,
    paired_end: bool,
    finished: bool,
}

impl FastqExecution {
    pub fn num_records(&self) -> usize {
        self.num_records
    }

    pub fn read_summary(&self) -> &str {
        &self.read_summary
    }

    pub fn is_paired_end(&self) -> bool {
        self.paired_end
    }

    pub fn next_batch(&mut self) -> Result<Option<Vec<AnnotatedFastq>>> {
        if self.finished {
            return Ok(None);
        }

        loop {
            let Some((batch, qc)) = self.source.next() else {
                self.stages.finish()?;
                self.finished = true;
                return Ok(None);
            };
            self.qc.extend(std::iter::once(qc));
            let batch = self.stages.process(batch)?;
            if !batch.is_empty() {
                return Ok(Some(batch));
            }
        }
    }

    pub fn finish(self) -> Result<FastqReport> {
        if !self.finished {
            anyhow::bail!("FASTQ execution has not reached end-of-input");
        }
        Ok(FastqReport {
            fastq: self.qc,
            middleware: self.stages.reports(),
        })
    }
}

/// Alignment execution configuration and stream construction.
pub struct AlignmentRunner<'a, A> {
    aligner: &'a mut A,
    num_threads: u16,
    mito_dna: HashSet<String>,
}

impl<'a, A: Aligner> AlignmentRunner<'a, A> {
    pub fn new(aligner: &'a mut A, num_threads: u16) -> Self {
        Self {
            aligner,
            num_threads,
            mito_dna: HashSet::new(),
        }
    }

    pub fn with_mito_dna<I, S>(mut self, mito_dna: I) -> Self
    where
        I: IntoIterator<Item = S>,
        S: Into<String>,
    {
        self.mito_dna.extend(mito_dna.into_iter().map(Into::into));
        self
    }

    pub fn stream(self, execution: FastqExecution) -> AlignmentResult<'a, A> {
        AlignmentResult::new(self.aligner, execution, &self.mito_dna, self.num_threads)
    }

    pub fn run<F>(self, execution: FastqExecution, mut consume: F) -> Result<RunReport>
    where
        F: FnMut(&AlignmentBatch) -> Result<()>,
    {
        let mut stream = self.stream(execution);
        while let Some(batch) = stream.next() {
            consume(&batch)?;
        }
        stream.finish()
    }
}

pub type AlignmentBatch = Vec<(Option<MultiMapR>, Option<MultiMapR>)>;
#[derive(Debug)]
pub struct RunReport {
    pub fastq: FastqReport,
    pub alignment: QcAlign,
}

/// Iterator that yields alignment results from annotated FASTQ reads with QC metrics.
pub struct AlignmentResult<'a, A> {
    aligner: &'a mut A,
    execution: FastqExecution,
    qc: QcAlign,
    header: noodles::sam::Header,
    num_threads: u16,
    num_records: usize,
    num_processed: usize,
    complete: bool,
    error: Option<anyhow::Error>,
}

impl<'a, A: Aligner> AlignmentResult<'a, A> {
    fn new(
        aligner: &'a mut A,
        execution: FastqExecution,
        mito_dna: &HashSet<String>,
        num_threads: u16,
    ) -> Self {
        let header = aligner.header();
        let num_records = execution.num_records();
        let mut qc = QcAlign::default();
        mito_dna.iter().for_each(|mito| {
            header
                .reference_sequences()
                .get_index_of(&BString::from(mito.as_str()))
                .map(|x| qc.mito_dna.insert(x));
        });

        Self {
            aligner,
            execution,
            qc,
            header,
            num_threads,
            num_records,
            num_processed: 0,
            complete: false,
            error: None,
        }
    }
}

impl<'a, A> AlignmentResult<'a, A> {
    pub fn num_records(&self) -> usize {
        self.num_records
    }

    pub fn num_processed(&self) -> usize {
        self.num_processed
    }

    pub fn finish(self) -> Result<RunReport> {
        if let Some(error) = self.error {
            return Err(error);
        }
        if !self.complete {
            anyhow::bail!("alignment stream has not reached end-of-input");
        }
        Ok(RunReport {
            fastq: self.execution.finish()?,
            alignment: self.qc,
        })
    }
}

/// Implement the Iterator trait for AlignmentResult.
/// The alignment results are yielded as a tuple of two MultiMapR.
/// If the read is unpaired, the second element is None.
impl<'a, A: Aligner> Iterator for AlignmentResult<'a, A> {
    type Item = AlignmentBatch;
    fn next(&mut self) -> Option<Self::Item> {
        let data = match self.execution.next_batch() {
            Ok(Some(data)) => data,
            Ok(None) => {
                self.num_processed = self.num_records;
                self.complete = true;
                return None;
            }
            Err(error) => {
                self.error = Some(error);
                self.complete = true;
                return None;
            }
        };
        self.num_processed += data.len();

        // Align the reads.
        let records = data.into_iter().map(AlignmentInput::from).collect();
        let results: Vec<_> = self.aligner.align_reads(self.num_threads, records);
        for alignment in &results {
            let result = match alignment {
                (Some(ali1), Some(ali2)) => self.qc.add_pair(&self.header, ali1, ali2),
                (Some(ali1), None) => self.qc.add_read1(&self.header, ali1),
                (None, Some(ali2)) => self.qc.add_read2(&self.header, ali2),
                _ => {
                    debug!("No alignment found for read");
                    Ok(())
                }
            };
            if let Err(error) = result {
                self.error = Some(error);
                self.complete = true;
                return None;
            }
        }
        Some(results)
    }
}

/// AnnotatedFastqReaders is formed by concatenating multiple AnnotatedFastqReader instances.
struct MultiAnnotatedFqReader {
    readers: Vec<AnnotatedFastqReader>,
    current: usize,
}

impl Iterator for MultiAnnotatedFqReader {
    type Item = (Vec<AnnotatedFastq>, QcFastq);

    fn next(&mut self) -> Option<Self::Item> {
        let reader = self.readers.get_mut(self.current)?;
        if let Some(chunk) = reader.next() {
            Some(chunk)
        } else {
            self.current += 1;
            self.next()
        }
    }
}

impl MultiAnnotatedFqReader {
    fn new(readers: Vec<AnnotatedFastqReader>) -> Self {
        Self {
            readers,
            current: 0,
        }
    }

    pub fn num_records(&self) -> usize {
        self.readers.iter().map(|x| x.annotation.num_reads()).sum()
    }

    pub fn is_paired_end(&self) -> Result<bool> {
        self.readers
            .iter()
            .map(|x| x.is_paired_end())
            .all_equal_value()
            .map_err(|_| anyhow::anyhow!("Not all readers are with the same paired-end status"))
    }
}

struct AnnotatedFastqReader {
    readers: PrefetchIterator<Vec<SmallVec<[fastq::Record; 4]>>>,
    annotation: AnnotationStage,
}

struct AnnotationStage {
    annotators: Vec<FastqAnnotator>,
    barcode_analyzer: BarcodeAnalyzer,
}

impl AnnotationStage {
    fn new(annotators: Vec<FastqAnnotator>, barcode_analyzer: BarcodeAnalyzer) -> Self {
        Self {
            annotators,
            barcode_analyzer,
        }
    }

    fn num_reads(&self) -> usize {
        self.barcode_analyzer.num_reads()
    }

    fn is_paired_end(&self) -> bool {
        let mut has_read1 = false;
        let mut has_read2 = false;
        self.annotators.iter().for_each(|x| {
            x.segment_info.iter().for_each(|info| {
                if info.region_type.is_target() {
                    if x.segment_info.is_reverse() {
                        has_read1 = true;
                    } else {
                        has_read2 = true;
                    }
                }
            });
        });
        has_read1 && has_read2
    }

    fn process_chunk<'a, I: IntoIterator<Item = &'a SmallVec<[fastq::Record; 4]>>>(
        &self,
        chunk: I,
    ) -> (Vec<AnnotatedFastq>, QcFastq) {
        process_chunk(&self.barcode_analyzer, &self.annotators, chunk)
    }
}

impl AnnotatedFastqReader {
    fn new<T: IntoIterator<Item = (FastqAnnotator, FastqReader)>>(
        iter: T,
        barcode_analyzer: BarcodeAnalyzer,
        chunk_size: usize,
    ) -> Self {
        let (annotators, readers): (Vec<_>, Vec<_>) = iter.into_iter().unzip();
        Self {
            annotation: AnnotationStage::new(annotators, barcode_analyzer),
            readers: PrefetchIterator::new(
                BatchedFqReader {
                    readers,
                    batch_size: chunk_size,
                },
                1,
            ),
        }
    }

    fn is_paired_end(&self) -> bool {
        self.annotation.is_paired_end()
    }
}

impl Iterator for AnnotatedFastqReader {
    type Item = (Vec<AnnotatedFastq>, QcFastq);

    fn next(&mut self) -> Option<Self::Item> {
        // Group reads of the same index from different files into a SmallVec.
        let chunk = self.readers.next()?;

        let parallel_chunk_size = (chunk.len() / 128).max(1);
        let annotation = &self.annotation;
        let processed: Vec<_> = chunk
            .par_chunks(parallel_chunk_size)
            .map(|chunk| annotation.process_chunk(chunk))
            .collect();
        let mut records = Vec::new();
        let mut qc = QcFastq::default();
        for (chunk, chunk_qc) in processed {
            records.extend(chunk);
            qc.extend(std::iter::once(chunk_qc));
        }
        Some((records, qc))
    }
}

fn process_chunk<'a, I: IntoIterator<Item = &'a SmallVec<[fastq::Record; 4]>>>(
    barcode_analyzer: &BarcodeAnalyzer,
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
                    if let Ok(anno) = annotator.annotate(record, barcode_analyzer) {
                        // annotate the read
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

/// A batched FASTQ reader that reads multiple FASTQ files in batches.
struct BatchedFqReader {
    readers: Vec<FastqReader>, // a list of open file handles
    batch_size: usize,         // target size for one chunk of data
}

impl Iterator for BatchedFqReader {
    type Item = Vec<SmallVec<[fastq::Record; 4]>>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut batch = Vec::new();
        let mut accumulated_length = 0;

        // Read records from all readers until reaching the batch size.
        // while loop for vertical iteration; readers.iter_mut() for horizontal iteration.
        while accumulated_length < self.batch_size {
            let mut max_read = 0;
            let mut min_read = usize::MAX;
            let records: SmallVec<[_; 4]> = self
                .readers
                .iter_mut() // read one record from each FASTQ file at the same position
                .flat_map(|reader| {
                    let mut buffer = fastq::Record::default();
                    let n = reader
                        .read_record(&mut buffer)
                        .expect("error reading fastq record");
                    min_read = min_read.min(n);
                    max_read = max_read.max(n);
                    if n > 0 {
                        accumulated_length += buffer.sequence().len();
                        strip_fq_suffix(&mut buffer);
                        Some(buffer)
                    } else {
                        None
                    }
                })
                .collect();
            if max_read == 0 {
                // All readers have reached EOF.
                if batch.is_empty() {
                    return None;
                } else {
                    break;
                }
            } else if min_read == 0 {
                panic!("Unequal number of reads in the chunk");
            } else {
                // Check records from all readers at the same position have the same name.
                assert!(
                    records.iter().map(|r| r.name()).all_equal(),
                    "read names mismatch"
                );
                batch.push(records);
            }
        }

        Some(batch)
    }
}

/// A FastqAnnotator that splits the reads into subregions, e.g., barcode, UMI, and
/// return annotated reads.
#[derive(Debug)]
struct FastqAnnotator {
    read_id: String,
    segment_info: SegmentInfo,
}

impl FastqAnnotator {
    pub fn new(read_id: impl Into<String>, segment_info: SegmentInfo) -> Option<Self> {
        if !segment_info.iter().any(|x| {
            x.region_type.is_barcode() || x.region_type.is_umi() || x.region_type.is_target()
        }) {
            None
        } else {
            Some(Self {
                read_id: read_id.into(),
                segment_info,
            })
        }
    }

    /// Annotate a single fastq record.
    fn annotate(
        &self,
        record: &fastq::Record,
        barcode_analyzer: &BarcodeAnalyzer,
    ) -> Result<AnnotatedFastq, SplitError> {
        let mut barcode: Option<Barcode> = None;
        let mut umi = None;
        let mut read1 = None;
        let mut read2 = None;

        let segments = self.segment_info.split(record)?;
        segments.into_iter().for_each(|segment| {
            if segment.is_barcode() || segment.is_umi() {
                let mut fq = segment.into_fq(record.definition());
                if self.segment_info.is_reverse() {
                    fq = rev_compl_fastq_record(fq);
                }

                if segment.is_barcode() {
                    let corrected = barcode_analyzer
                        .correct_barcode(segment.region_id(), fq.sequence(), fq.quality_scores())
                        .ok()
                        .map(|x| x.to_vec());
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
                    panic!("Multiple target regions found in one fastq record!");
                } else {
                    let fq = segment.into_fq(record.definition());
                    // TODO: polyA and adapter trimming
                    if self.segment_info.is_reverse() {
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

/// Metadata attached to every alignment generated from a read.
#[derive(Debug)]
pub struct ReadMetadata {
    pub barcode: Barcode,
    pub umi: Option<UMI>,
}

/// Backend-neutral input passed to an aligner.
#[derive(Debug)]
pub struct AlignmentInput {
    pub read1: Option<fastq::Record>,
    pub read2: Option<fastq::Record>,
    pub metadata: ReadMetadata,
}

/// An annotated fastq record with barcode, UMI, and sequence.
#[derive(Debug)]
pub struct AnnotatedFastq {
    pub barcode: Option<Barcode>,
    pub umi: Option<UMI>,
    pub read1: Option<fastq::Record>,
    pub read2: Option<fastq::Record>,
}

impl From<AnnotatedFastq> for AlignmentInput {
    fn from(record: AnnotatedFastq) -> Self {
        Self {
            read1: record.read1,
            read2: record.read2,
            metadata: ReadMetadata {
                barcode: record
                    .barcode
                    .expect("annotated FASTQ passed to alignment without a barcode"),
                umi: record.umi,
            },
        }
    }
}

impl AnnotatedFastq {
    /// Join another AnnotatedFastq from the same insert into self.
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

fn strip_fq_suffix(record: &mut fastq::Record) {
    let read_name = record.name();
    let n = read_name.len();
    if n > 2 {
        let suffix = &read_name[n - 2..];
        if suffix == b"/1" || suffix == b"/2" {
            record.name_mut().truncate(n - 2);
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
        let mut execution = FastqPlan::new(vec![assay], modality)
            .build(false, 500)
            .unwrap();
        let file = std::fs::File::open(output).unwrap();
        let mut lines = std::io::BufReader::new(flate2::read::GzDecoder::new(file)).lines();
        while let Some(batch) = execution.next_batch().unwrap() {
            for fq in batch {
                assert_eq!(show_fq(&fq), lines.next().unwrap().unwrap());
            }
        }
        execution.finish().unwrap();
    }

    #[test]
    fn test_io() {
        let seqspec = "data/test4.yaml";
        let assay = Assay::from_path(seqspec).unwrap();
        let modality = assay.modalities[0].clone();
        let mut execution = FastqPlan::new(vec![assay], modality)
            .build(false, 5000)
            .unwrap();
        // open a file
        while let Some(batch) = execution.next_batch().unwrap() {
            for fq in batch {
                println!("{}", show_fq(&fq));
            }
        }
        execution.finish().unwrap();
    }

    #[test]
    fn test_fastq() {
        test_fq("data/test2.yaml", "data/test2.out.gz");
        test_fq("data/test3.yaml", "data/test3.out.gz");
        test_fq("data/test4.yaml", "data/test4.out.gz");
    }
}
