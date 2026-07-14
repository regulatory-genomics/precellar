//! FASTQ annotation, middleware execution, and alignment streaming.
//!
//! This module implements the path from assay-defined FASTQ inputs to batches
//! of alignments:
//!
//! 1. [`FastqPlan`] selects a modality, barcode-correction settings, and
//!    optional [`FastqStage`] middleware.
//! 2. [`FastqPlan::build`] opens the assay inputs and creates a one-shot
//!    [`FastqExecution`].
//! 3. [`FastqExecution::next_batch`] reads synchronized FASTQ records,
//!    annotates barcode, UMI, and target segments, accumulates FASTQ QC, and
//!    applies middleware in registration order.
//! 4. [`AlignmentRunner::stream`] converts the execution into an
//!    [`AlignmentResult`] iterator that aligns each surviving batch and
//!    accumulates alignment QC.
//! 5. The execution must be drained to end-of-input and explicitly finalized
//!    with [`FastqExecution::finish`] or [`AlignmentResult::finish`] to obtain
//!    its owned report.
//!
//! Assays are processed sequentially, while annotation within each input batch
//! is parallelized with Rayon. Corresponding records from all FASTQ files in an
//! assay are synchronized by position and normalized read name. The low-level
//! reader assumes equal record counts and matching names; violations currently
//! panic rather than being returned as recoverable errors.
//!
//! Alignment iteration uses deferred error reporting because [`Iterator`]
//! cannot yield `Result` without changing downstream iterator interfaces. A
//! processing error ends the current iteration and is returned by
//! [`AlignmentResult::finish`]. Callers must therefore always call `finish`
//! after consuming the stream.

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

/// Configuration for standard barcode correction during FASTQ annotation.
///
/// This configuration is installed on each assay's [`BarcodeAnalyzer`] only
/// when barcode correction is enabled in [`FastqPlan::build`]. Annotation still
/// preserves the raw barcode when correction fails; the corrected sequence is
/// represented by [`Barcode::corrected`] and remains `None` in that case.
#[derive(Debug, Clone)]
pub struct BarcodeCorrectionConfig {
    /// Minimum posterior confidence required to accept a correction.
    pub confidence_threshold: f64,
    /// Maximum number of barcode mismatches considered by correction.
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

/// Consuming builder for a single annotated FASTQ execution.
///
/// A plan owns the configured middleware stages. Calling [`Self::build`]
/// consumes the plan and moves those stateful stages into the resulting
/// [`FastqExecution`], so a plan is intentionally one-shot rather than reusable.
/// The assays themselves are retained behind an `Arc` only to provide cheap,
/// immutable shared ownership while constructing assay-level readers.
pub struct FastqPlan {
    assays: Arc<[Assay]>,
    modality: Modality,
    barcode_config: BarcodeCorrectionConfig,
    stages: FastqStagePipeline,
}

impl FastqPlan {
    /// Creates a plan for `modality` across the supplied assays.
    ///
    /// Assays are consumed in vector order. The default barcode-correction
    /// configuration is installed, but correction is not enabled until
    /// [`Self::build`] is called with `correct_barcode = true`.
    pub fn new(assays: Vec<Assay>, modality: Modality) -> Self {
        Self {
            assays: assays.into(),
            modality,
            barcode_config: BarcodeCorrectionConfig::default(),
            stages: FastqStagePipeline::default(),
        }
    }

    /// Replaces the standard barcode-correction configuration.
    ///
    /// This has no effect when [`Self::build`] is called with barcode
    /// correction disabled.
    pub fn with_barcode_config(mut self, config: BarcodeCorrectionConfig) -> Self {
        self.barcode_config = config;
        self
    }

    /// Appends a stateful middleware stage to the annotated FASTQ pipeline.
    ///
    /// Stages run in registration order. A stage may transform, filter, or
    /// buffer records. Its [`FastqStage::finish`] method is called only after
    /// the FASTQ source reaches end-of-input.
    pub fn with_stage<S>(mut self, stage: S) -> Self
    where
        S: FastqStage + 'static,
    {
        self.stages.push_stage(stage);
        self
    }

    /// Opens the configured assay inputs and creates a one-shot execution.
    ///
    /// `chunk_size` is an approximate aggregate sequence-length target, in
    /// bases, for each synchronized input batch. It is not a record count. A
    /// batch may exceed the target by one complete synchronized group of FASTQ
    /// records.
    ///
    /// `correct_barcode` controls whether the configured correction options are
    /// installed on the assay barcode analyzers. Raw barcodes are annotated in
    /// either mode.
    ///
    /// # Requirements
    ///
    /// `chunk_size` must be greater than zero. A zero value cannot advance the
    /// current low-level batch reader.
    ///
    /// # Errors
    ///
    /// Returns an error when the assay-level readers do not agree on whether
    /// the selected modality is paired-end.
    ///
    /// FASTQ opening and low-level reading currently follow `seqspec`'s
    /// optional/panic behavior rather than returning every I/O failure here.
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

/// Completed FASTQ-side quality-control report.
///
/// The report is available only after the source has reached end-of-input and
/// every middleware stage has finalized successfully.
#[derive(Debug)]
pub struct FastqReport {
    /// Aggregated annotation and base-quality metrics from all assays.
    pub fastq: QcFastq,
    /// Stage-specific reports in middleware registration order.
    pub middleware: Vec<MiddlewareQcReport>,
}

/// Stateful, one-shot execution of a prepared [`FastqPlan`].
///
/// The execution owns the input readers, middleware stages, and FASTQ QC
/// accumulator. Repeated calls to [`Self::next_batch`] advance all of that
/// state. Empty batches produced by filtering middleware are skipped internally
/// and are never returned to the caller.
///
/// Callers must drain the execution until `next_batch` returns `Ok(None)` before
/// calling [`Self::finish`]. Reaching EOF finalizes middleware exactly once.
/// Dropping an execution does not finalize middleware or produce a report.
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
    /// Returns the estimated total number of logical records across all assays.
    ///
    /// This value is used for progress reporting and comes from the assay
    /// barcode analyzers rather than from records already consumed.
    pub fn num_records(&self) -> usize {
        self.num_records
    }

    /// Returns a human-readable per-assay read-count summary.
    ///
    /// Multiple assay counts are separated by `" + "`.
    pub fn read_summary(&self) -> &str {
        &self.read_summary
    }

    /// Returns whether every assay reader was classified as paired-end.
    pub fn is_paired_end(&self) -> bool {
        self.paired_end
    }

    /// Advances annotation and middleware processing by one non-empty batch.
    ///
    /// The method reads and annotates source records, merges their local FASTQ
    /// QC contribution, and passes the resulting batch through every configured
    /// middleware stage. If middleware removes the entire batch, processing
    /// continues until a non-empty batch or EOF is reached.
    ///
    /// At EOF, middleware is finalized before `Ok(None)` is returned. Subsequent
    /// calls also return `Ok(None)`.
    ///
    /// # Errors
    ///
    /// Returns errors raised by middleware processing or finalization. If
    /// finalization fails, the execution is not marked finished and a later
    /// call retries finalization.
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

    /// Consumes a drained execution and returns its completed FASTQ report.
    ///
    /// # Errors
    ///
    /// Returns an error if [`Self::next_batch`] has not yet observed and
    /// successfully finalized end-of-input.
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

/// Configuration used to construct the canonical alignment stream.
///
/// The runner borrows an aligner for the lifetime of the stream and owns the
/// execution settings that apply uniformly to each batch. Calling
/// [`Self::stream`] consumes the runner and moves a [`FastqExecution`] into the
/// returned [`AlignmentResult`].
pub struct AlignmentRunner<'a, A> {
    aligner: &'a mut A,
    num_threads: u16,
    mito_dna: HashSet<String>,
}

impl<'a, A: Aligner> AlignmentRunner<'a, A> {
    /// Creates a runner that uses `num_threads` for each aligner batch call.
    ///
    /// No mitochondrial references are configured initially. The accepted
    /// thread-count range is determined by the concrete [`Aligner`]
    /// implementation.
    pub fn new(aligner: &'a mut A, num_threads: u16) -> Self {
        Self {
            aligner,
            num_threads,
            mito_dna: HashSet::new(),
        }
    }

    /// Adds reference names that should contribute to mitochondrial QC.
    ///
    /// Names are deduplicated. When the stream is created, each name is resolved
    /// against the aligner's SAM header. Names absent from that header are
    /// ignored.
    pub fn with_mito_dna<I, S>(mut self, mito_dna: I) -> Self
    where
        I: IntoIterator<Item = S>,
        S: Into<String>,
    {
        self.mito_dna.extend(mito_dna.into_iter().map(Into::into));
        self
    }

    /// Creates the alignment iterator for `execution`.
    ///
    /// The stream owns the FASTQ execution and alignment QC state while
    /// borrowing the aligner. It must be consumed to EOF and finalized with
    /// [`AlignmentResult::finish`] to surface deferred errors and obtain the
    /// final [`RunReport`].
    pub fn stream(self, execution: FastqExecution) -> AlignmentResult<'a, A> {
        AlignmentResult::new(self.aligner, execution, &self.mito_dna, self.num_threads)
    }
}

/// Alignment results for one annotated FASTQ batch.
///
/// Each vector element corresponds to one [`AlignmentInput`]. The tuple stores
/// read-1 and read-2 results respectively. Either side may be `None` for
/// single-end input, missing target reads, or a backend that produced no
/// alignment. Each present value contains one primary alignment and any
/// secondary alignments reported by the backend.
pub type AlignmentBatch = Vec<(Option<MultiMapR>, Option<MultiMapR>)>;

/// Completed report for an alignment stream.
///
/// FASTQ and middleware metrics are finalized together with the accumulated
/// alignment metrics, preserving one ownership boundary for the complete run.
#[derive(Debug)]
pub struct RunReport {
    /// Annotation and middleware QC from the owned FASTQ execution.
    pub fastq: FastqReport,
    /// QC accumulated from every alignment batch yielded by the stream.
    pub alignment: QcAlign,
}

/// Stateful alignment iterator with owned QC and deferred error reporting.
///
/// Each call to [`Iterator::next`] obtains one annotated FASTQ batch, converts
/// it to backend-neutral [`AlignmentInput`] values, invokes the aligner, updates
/// alignment QC, and yields an [`AlignmentBatch`]. QC is not exposed while the
/// stream is running; it is returned by [`Self::finish`] after EOF.
///
/// FASTQ, middleware, and QC errors cannot be represented by this iterator's
/// item type. Such an error causes `next` to return `None` and is retained for
/// [`Self::finish`] to return. Consumers must stop after the first `None` and
/// always call `finish`; this type does not implement `FusedIterator`.
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
    /// Returns the estimated total logical-record count for progress reporting.
    pub fn num_records(&self) -> usize {
        self.num_records
    }

    /// Returns the number of annotated records submitted to the aligner so far.
    ///
    /// Records removed by middleware are not counted as batches are processed.
    /// On successful EOF this value is set to [`Self::num_records`] so progress
    /// displays finish at the configured total.
    pub fn num_processed(&self) -> usize {
        self.num_processed
    }

    /// Consumes a completed stream and returns its owned run report.
    ///
    /// # Errors
    ///
    /// Returns the processing error retained by the iterator, if any. Returns
    /// an error when the stream has not yet reached EOF. Successful completion
    /// also requires the underlying [`FastqExecution`] to finalize successfully.
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

/// Advances FASTQ processing, alignment, and alignment QC by one batch.
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

/// Concatenates assay-level readers without interleaving their batches.
///
/// The first reader is exhausted before the next reader is advanced. This
/// preserves assay order in both record output and QC accumulation.
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

/// Prefetches synchronized physical records and annotates each batch in parallel.
struct AnnotatedFastqReader {
    readers: PrefetchIterator<Vec<SmallVec<[fastq::Record; 4]>>>,
    annotation: AnnotationStage,
}

/// Immutable annotation configuration shared by Rayon workers.
///
/// Workers return local [`QcFastq`] values, which the owning reader merges after
/// parallel processing. The barcode analyzer is therefore read concurrently but
/// never mutated during annotation.
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
        // Synchronized record groups are already assembled by BatchedFqReader.
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

/// Annotates synchronized record groups and returns records plus local FASTQ QC.
///
/// One group contains records at the same position from every participating
/// FASTQ input. Individual physical annotations are joined into one logical
/// insert. Split failures increment the per-read defect count and drop only the
/// failed physical annotation. Logical inserts without a barcode are omitted
/// from downstream output after their available QC has been recorded.
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

/// Reads positionally synchronized records from multiple FASTQ files.
///
/// `batch_size` is an approximate aggregate sequence-length target in bases,
/// not a record count. One record is read from every input during each vertical
/// step. A complete synchronized group is always retained, so a returned batch
/// may exceed the target.
///
/// Read names are normalized by removing a trailing `/1` or `/2` and must then
/// match across every input. All readers must reach EOF at the same position.
///
/// # Panics
///
/// Panics when FASTQ reading fails, input files contain different record counts,
/// or synchronized records have different normalized names.
struct BatchedFqReader {
    readers: Vec<FastqReader>,
    batch_size: usize,
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

/// Splits one physical FASTQ record into barcode, UMI, and target regions.
///
/// Segment orientation determines whether barcode and UMI records are reverse
/// complemented and whether a target is assigned to read 1 or read 2. Multiple
/// barcode segments are concatenated. If several UMI segments are present in
/// one physical record, the final segment replaces earlier ones.
#[derive(Debug)]
struct FastqAnnotator {
    read_id: String,
    segment_info: SegmentInfo,
}

impl FastqAnnotator {
    /// Creates an annotator when the segment specification contains useful data.
    ///
    /// Specifications with no barcode, UMI, or target regions are ignored by
    /// returning `None`.
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

    /// Annotates one physical FASTQ record according to the segment definition.
    ///
    /// Barcode-correction failures are not annotation errors: the raw barcode is
    /// retained with no corrected value. Structural split failures are returned
    /// to the caller and counted as malformed reads by [`process_chunk`].
    ///
    /// # Panics
    ///
    /// Panics if one physical record produces more than one target-bearing
    /// segment. A logical read pair must instead be assembled from separate
    /// physical records through [`AnnotatedFastq::join`].
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

/// Raw and optionally corrected cell-barcode sequence.
///
/// The raw FASTQ record retains both sequence and quality scores. Corrected
/// sequences contain bases only and are present when every constituent barcode
/// segment was corrected successfully.
#[derive(Debug)]
pub struct Barcode {
    /// Concatenated raw barcode bases and their quality scores.
    pub raw: fastq::Record,
    /// Concatenated corrected bases, or `None` when correction was disabled or
    /// any constituent segment could not be corrected.
    pub corrected: Option<Vec<u8>>,
}

impl Barcode {
    /// Appends another barcode segment to this barcode.
    ///
    /// Raw sequence and quality scores are always appended. Corrected sequence
    /// is appended only when both barcodes have corrected values; otherwise the
    /// combined barcode is marked uncorrected.
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

/// UMI sequence and quality scores represented as a FASTQ record.
pub type UMI = fastq::Record;

/// Barcode and UMI metadata carried alongside target reads into an aligner.
///
/// The barcode is required at the alignment boundary. UMI metadata remains
/// optional because not every assay defines a UMI region.
#[derive(Debug)]
pub struct ReadMetadata {
    /// Required raw and optionally corrected cell barcode.
    pub barcode: Barcode,
    /// Optional UMI sequence and quality scores.
    pub umi: Option<UMI>,
}

/// Backend-neutral input passed to an aligner.
///
/// Inputs may contain read 1, read 2, or both. The metadata is independent of
/// backend-specific alignment input and is later attached to generated SAM
/// records by the aligner abstraction.
#[derive(Debug)]
pub struct AlignmentInput {
    /// Target sequence assigned to read 1, if present.
    pub read1: Option<fastq::Record>,
    /// Target sequence assigned to read 2, if present.
    pub read2: Option<fastq::Record>,
    /// Barcode and UMI metadata shared by the target reads.
    pub metadata: ReadMetadata,
}

/// One logical insert assembled from annotated physical FASTQ records.
///
/// During annotation all fields are optional because a physical record may
/// contribute only one component. After synchronized records are joined,
/// barcode-less inserts are filtered before alignment. Converting a remaining
/// value into [`AlignmentInput`] therefore requires a barcode and moves all
/// owned FASTQ records without copying them.
#[derive(Debug)]
pub struct AnnotatedFastq {
    /// Raw and optionally corrected cell barcode.
    pub barcode: Option<Barcode>,
    /// Optional UMI sequence and quality scores.
    pub umi: Option<UMI>,
    /// Optional target sequence assigned to read 1.
    pub read1: Option<fastq::Record>,
    /// Optional target sequence assigned to read 2.
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
    /// Joins another physical annotation from the same logical insert.
    ///
    /// Barcode and UMI sequence/quality fields are concatenated in input order.
    /// Missing fields are adopted from `other`. Read 1 and read 2 are moved into
    /// their corresponding empty slots.
    ///
    /// # Panics
    ///
    /// Panics when both values contain read 1 or both contain read 2. A logical
    /// insert may have at most one target record for each side.
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

/// Appends sequence and quality scores from `other` to `this`.
///
/// The definition, read name, and description of `this` are preserved. This is
/// used to concatenate segmented barcodes and UMIs while keeping their sequence
/// and quality lengths synchronized.
pub fn extend_fastq_record(this: &mut fastq::Record, other: &fastq::Record) {
    this.sequence_mut().extend_from_slice(other.sequence());
    this.quality_scores_mut()
        .extend_from_slice(other.quality_scores());
}

/// Removes a conventional `/1` or `/2` suffix before synchronization checks.
///
/// Other naming conventions are left unchanged.
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
