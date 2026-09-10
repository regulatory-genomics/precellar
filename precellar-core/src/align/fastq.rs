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
use crate::long::{BarcodeExtractor, LrBarcodeExtractionStats};
use crate::pipeline::{FastqStage, FastqStagePipeline, MiddlewareQcReport};
use crate::qc::{QcAlign, QcFastq, QcLongRead};
use crate::utils::{rev_compl_fastq_record, PrefetchIterator};
use anyhow::{Context, Result};
use bstr::BString;
use itertools::Itertools;
use log::{debug, info};
use noodles_fastq as fastq;
use rayon::iter::ParallelIterator;
use rayon::slice::ParallelSlice;
use seqspec::{Assay, AssayType, FastqReader, Modality, SegmentInfo};
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
    sequencing_type: Option<AssayType>,
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
            sequencing_type: None,
            barcode_config: BarcodeCorrectionConfig::default(),
            stages: FastqStagePipeline::default(),
        }
    }

    /// Explicitly selects the sequencing type used for FASTQ annotation.
    ///
    /// When omitted, [`Self::build`] falls back to sampling the configured
    /// FASTQ input at runtime. This override is scoped to the plan's selected
    /// modality and applies to every assay in the plan.
    pub fn with_sequencing_type(mut self, sequencing_type: AssayType) -> Self {
        self.sequencing_type = Some(sequencing_type);
        self
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
        anyhow::ensure!(chunk_size > 0, "FASTQ chunk size must be greater than zero");

        let reader = self.build_reader(correct_barcode, chunk_size)?;
        let has_long_read = reader.has_long_read();
        if has_long_read && !self.stages.is_empty() {
            anyhow::bail!("FASTQ middleware is not supported for long-read assays yet");
        }
        let num_records = reader.num_records();
        let paired_end = reader.is_paired_end()?;
        let n_reads = if has_long_read {
            "unknown".to_string()
        } else {
            reader
                .readers
                .iter()
                .map(|r| indicatif::HumanCount(r.annotation.num_reads() as u64).to_string())
                .join(" + ")
        };
        Ok(FastqExecution {
            source: reader,
            stages: self.stages,
            qc: QcFastq::default(),
            qc_long_read: has_long_read.then(QcLongRead::default),
            num_records,
            read_summary: n_reads,
            paired_end,
            finished: false,
        })
    }

    fn build_reader(
        &self,
        correct_barcode: bool,
        chunk_size: usize,
    ) -> Result<MultiAnnotatedFqReader> {
        let num_assays = self.assays.len();
        let mut readers = Vec::with_capacity(num_assays);
        let mut plan_assay_type = None;

        for (i, assay) in self.assays.iter().enumerate() {
            if num_assays > 1 {
                info!(">>>Processing assay {}/{}<<<", i + 1, num_assays);
            }

            let assay_type = match self.sequencing_type {
                Some(sequencing_type) => sequencing_type,
                None => assay.detect_assay_type(&self.modality).with_context(|| {
                    format!(
                        "failed to detect assay type for assay '{}' and modality {}",
                        assay.assay_id, self.modality
                    )
                })?,
            };
            if let Some(expected) = plan_assay_type {
                anyhow::ensure!(
                    expected == assay_type,
                    "mixing short-read and long-read assays in one FASTQ plan is not supported"
                );
            } else {
                plan_assay_type = Some(assay_type);
            }
            if self.sequencing_type.is_some() {
                info!("Using configured sequencing type: {:?}", assay_type);
            } else {
                info!("Detected assay type: {:?}", assay_type);
            }
            if assay_type == AssayType::ShortRead {
                assay.validate_short_read_strands(&self.modality)?;
            }

            let physical_fastq_reads: Vec<_> = assay
                .iter_reads(self.modality)
                .filter(|read| {
                    read.files
                        .as_ref()
                        .is_some_and(|files| files.iter().any(|file| file.filetype == "fastq"))
                })
                .map(|read| read.read_id.as_str())
                .collect();

            let inputs: Vec<(FastqAnnotator, FastqReader)> = assay
                .get_segments_by_modality(self.modality)
                .filter_map(|(read, segment_info)| {
                    FastqAnnotator::new(&read.read_id, segment_info).map(|annotator| {
                        let reader = read
                            .try_open()
                            .with_context(|| {
                                format!(
                                    "failed to open FASTQ input for assay '{}' read '{}'",
                                    assay.assay_id, read.read_id
                                )
                            })?
                            .ok_or_else(|| {
                                anyhow::anyhow!(
                                    "assay '{}' read '{}' has no FASTQ input",
                                    assay.assay_id,
                                    read.read_id
                                )
                            })?;
                        Ok((annotator, reader))
                    })
                })
                .collect::<Result<_>>()?;

            let barcode_processor = match assay_type {
                AssayType::ShortRead => {
                    let mut barcode_analyzer = BarcodeAnalyzer::new(assay, self.modality);
                    barcode_analyzer.summary();
                    if correct_barcode {
                        barcode_analyzer.barcode_correct_options = Some(BarcodeCorrectOptions {
                            bc_confidence_threshold: self.barcode_config.confidence_threshold,
                            max_mismatch: self.barcode_config.max_mismatch,
                            ..Default::default()
                        });
                    }
                    BarcodeProcessor::ShortRead(barcode_analyzer)
                }
                AssayType::LongRead => {
                    anyhow::ensure!(
                        physical_fastq_reads.len() == 1,
                        "long-read assay '{}' modality {} must have exactly one physical FASTQ Read; found {} ({})",
                        assay.assay_id,
                        self.modality,
                        physical_fastq_reads.len(),
                        physical_fastq_reads.join(", ")
                    );
                    anyhow::ensure!(
                        inputs.len() == 1,
                        "long-read assay '{}' modality {} must have exactly one full-length FASTQ Read; found {}",
                        assay.assay_id,
                        self.modality,
                        inputs.len()
                    );
                    anyhow::ensure!(
                        inputs[0].0.target_count() == 1,
                        "long-read assay '{}' modality {} must contain exactly one target segment in its full-length Read; found {}",
                        assay.assay_id,
                        self.modality,
                        inputs[0].0.target_count()
                    );

                    let whitelists = assay.get_whitelists(self.modality);
                    anyhow::ensure!(
                        !whitelists.is_empty(),
                        "long-read assay '{}' modality {} requires at least one barcode whitelist",
                        assay.assay_id,
                        self.modality
                    );
                    let barcode_extractor = BarcodeExtractor::new(
                        &assay.library_spec,
                        &self.modality,
                        whitelists,
                    )
                    .with_context(|| {
                        format!(
                            "failed to configure long-read barcode extraction for assay '{}' and modality {}",
                            assay.assay_id, self.modality
                        )
                    })?;
                    BarcodeProcessor::LongRead(barcode_extractor)
                }
            };

            readers.push(AnnotatedFastqReader::new(
                inputs,
                barcode_processor,
                chunk_size,
            ));
        }

        Ok(MultiAnnotatedFqReader::new(readers))
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
    /// Long-read barcode-extraction metrics, when this was a long-read execution.
    pub long_read: Option<QcLongRead>,
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
    qc_long_read: Option<QcLongRead>,
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

    /// Returns whether this execution processes long-read input.
    pub fn is_long_read(&self) -> bool {
        self.qc_long_read.is_some()
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
            let Some((batch, qc, qc_long_read)) = self.source.next() else {
                self.stages.finish()?;
                self.finished = true;
                return Ok(None);
            };
            self.qc.extend(std::iter::once(qc));
            if let (Some(total), Some(local)) = (&mut self.qc_long_read, qc_long_read) {
                total.extend(std::iter::once(local));
            }
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
            long_read: self.qc_long_read,
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
    header: noodles_sam::Header,
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
                if self.num_records > 0 {
                    self.num_processed = self.num_records;
                }
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

/// Assay-specific barcode annotation strategy.
#[derive(Debug)]
enum BarcodeProcessor {
    ShortRead(BarcodeAnalyzer),
    LongRead(BarcodeExtractor),
}

impl BarcodeProcessor {
    fn num_reads(&self) -> usize {
        match self {
            Self::ShortRead(analyzer) => analyzer.num_reads(),
            Self::LongRead(_) => 0,
        }
    }

    fn is_long_read(&self) -> bool {
        matches!(self, Self::LongRead(_))
    }

    fn expects_umi(&self) -> bool {
        matches!(self, Self::LongRead(extractor) if extractor.expects_umi())
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
    type Item = (Vec<AnnotatedFastq>, QcFastq, Option<QcLongRead>);

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

    pub fn has_long_read(&self) -> bool {
        self.readers.iter().any(|x| x.annotation.is_long_read())
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
    barcode_processor: BarcodeProcessor,
}

impl AnnotationStage {
    fn new(annotators: Vec<FastqAnnotator>, barcode_processor: BarcodeProcessor) -> Self {
        Self {
            annotators,
            barcode_processor,
        }
    }

    fn num_reads(&self) -> usize {
        self.barcode_processor.num_reads()
    }

    fn is_long_read(&self) -> bool {
        self.barcode_processor.is_long_read()
    }

    fn is_paired_end(&self) -> bool {
        if self.is_long_read() {
            return false;
        }
        let mut has_read1 = false;
        let mut has_read2 = false;
        self.annotators.iter().for_each(|x| {
            x.segment_info.iter().for_each(|info| {
                if info.region_type.is_target() {
                    if x.segment_info.is_reverse() {
                        has_read2 = true;
                    } else {
                        has_read1 = true;
                    }
                }
            });
        });
        has_read1 && has_read2
    }

    fn process_chunk<'a, I: IntoIterator<Item = &'a SmallVec<[fastq::Record; 4]>>>(
        &self,
        chunk: I,
    ) -> (Vec<AnnotatedFastq>, QcFastq, Option<QcLongRead>) {
        process_chunk(&self.barcode_processor, &self.annotators, chunk)
    }
}

impl AnnotatedFastqReader {
    fn new<T: IntoIterator<Item = (FastqAnnotator, FastqReader)>>(
        iter: T,
        barcode_processor: BarcodeProcessor,
        chunk_size: usize,
    ) -> Self {
        let (annotators, readers): (Vec<_>, Vec<_>) = iter.into_iter().unzip();
        Self {
            annotation: AnnotationStage::new(annotators, barcode_processor),
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
    type Item = (Vec<AnnotatedFastq>, QcFastq, Option<QcLongRead>);

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
        let mut qc_long_read = self.annotation.is_long_read().then(QcLongRead::default);
        for (chunk, chunk_qc, chunk_qc_long_read) in processed {
            records.extend(chunk);
            qc.extend(std::iter::once(chunk_qc));
            if let (Some(total), Some(local)) = (&mut qc_long_read, chunk_qc_long_read) {
                total.extend(std::iter::once(local));
            }
        }
        Some((records, qc, qc_long_read))
    }
}

/// Annotates synchronized record groups and returns records plus local FASTQ QC.
///
/// One group contains records at the same position from every participating
/// FASTQ input. Individual physical annotations are joined into one logical
/// insert. Split failures increment the per-read defect count and drop only the
/// failed physical annotation. Logical inserts without a barcode, or without an
/// expected long-read UMI, are omitted after their available QC is recorded.
fn process_chunk<'a, I: IntoIterator<Item = &'a SmallVec<[fastq::Record; 4]>>>(
    barcode_processor: &BarcodeProcessor,
    annotators: &[FastqAnnotator],
    chunk: I,
) -> (Vec<AnnotatedFastq>, QcFastq, Option<QcLongRead>) {
    let mut qc = QcFastq::default();
    let mut qc_long_read = barcode_processor.is_long_read().then(QcLongRead::default);
    let expects_umi = barcode_processor.expects_umi();
    if expects_umi {
        qc_long_read.as_mut().unwrap().enable_umi_extraction();
    }
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
                    match annotator.annotate(record, barcode_processor) {
                        Ok((anno, long_read_stats)) => {
                            if let (Some(qc_long_read), Some(stats)) =
                                (&mut qc_long_read, long_read_stats)
                            {
                                qc_long_read
                                    .record_orientation(stats.is_reverse, !stats.composite_pass);
                                qc_long_read.record_composite_alignment(stats.composite_pass);
                                qc_long_read.record_barcode_extraction(anno.barcode.is_some());
                                if expects_umi {
                                    qc_long_read.record_umi_extraction(anno.umi.is_some());
                                }
                                for (group_name, intersection_hit) in stats.consensus_results {
                                    qc_long_read.record_consensus(&group_name, intersection_hit);
                                }
                            }
                            Some(anno)
                        }
                        Err(error) => {
                            debug!("Failed to annotate read '{}': {error:#}", id);
                            *qc.num_defect.entry(id.clone()).or_insert(0) += 1;
                            None
                        }
                    }
                })
                .reduce(|mut this, other| {
                    this.join(other);
                    this
                })?;
            qc.update(&fq);
            if fq.barcode.is_none() || (expects_umi && fq.umi.is_none()) {
                None
            } else {
                Some(fq)
            }
        })
        .collect();
    (annotated, qc, qc_long_read)
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

    fn target_count(&self) -> usize {
        self.segment_info
            .iter()
            .filter(|segment| segment.region_type.is_target())
            .count()
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
        barcode_processor: &BarcodeProcessor,
    ) -> Result<(AnnotatedFastq, Option<LrBarcodeExtractionStats>)> {
        match barcode_processor {
            BarcodeProcessor::ShortRead(barcode_analyzer) => {
                Ok((self.annotate_short_read(record, barcode_analyzer)?, None))
            }
            BarcodeProcessor::LongRead(barcode_extractor) => {
                let (annotated, stats) = self.annotate_long_read(record, barcode_extractor)?;
                Ok((annotated, Some(stats)))
            }
        }
    }

    fn annotate_short_read(
        &self,
        record: &fastq::Record,
        barcode_analyzer: &BarcodeAnalyzer,
    ) -> Result<AnnotatedFastq> {
        let mut barcode: Option<Barcode> = None;
        let mut umi = None;
        let mut read1 = None;
        let mut read2 = None;

        let segments = self.segment_info.split(record)?;
        segments.into_iter().for_each(|segment| {
            if segment.is_barcode() || segment.is_umi() {
                let mut fq = segment.into_fq(record.definition());
                let should_reverse_complement = if segment.is_barcode() {
                    self.segment_info.is_reverse() ^ barcode_analyzer.should_rc(segment.region_id())
                } else {
                    self.segment_info.is_reverse()
                };
                if should_reverse_complement {
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

    fn annotate_long_read(
        &self,
        record: &fastq::Record,
        barcode_extractor: &BarcodeExtractor,
    ) -> Result<(AnnotatedFastq, LrBarcodeExtractionStats)> {
        let (barcode_result, extraction_stats) = barcode_extractor.extract_barcode(record)?;

        let barcode = barcode_result.barcode.as_ref().map(|sequence| {
            let raw = fastq::Record::new(
                record.definition().clone(),
                sequence.clone(),
                vec![b'I'; sequence.len()],
            );
            Barcode {
                raw,
                corrected: Some(sequence.clone()),
            }
        });

        let (head_trim, tail_trim) = if barcode_result.is_reverse_complemented {
            (
                barcode_result.three_prime_trim,
                barcode_result.five_prime_trim,
            )
        } else {
            (
                barcode_result.five_prime_trim,
                barcode_result.three_prime_trim,
            )
        };

        let sequence = record.sequence();
        let quality_scores = record.quality_scores();
        anyhow::ensure!(
            sequence.len() == quality_scores.len(),
            "long-read FASTQ sequence and quality lengths differ"
        );
        let total_trim = head_trim
            .checked_add(tail_trim)
            .context("long-read trim lengths overflowed")?;
        anyhow::ensure!(
            total_trim < sequence.len(),
            "long-read sequence has no target bases after trimming {} + {} bp from {} bp",
            head_trim,
            tail_trim,
            sequence.len(),
        );
        let target_end = sequence.len() - tail_trim;
        let target = fastq::Record::new(
            record.definition().clone(),
            sequence[head_trim..target_end].to_vec(),
            quality_scores[head_trim..target_end].to_vec(),
        );

        let (read1, read2) = if barcode_result.is_reverse_complemented {
            (None, Some(target))
        } else {
            (Some(target), None)
        };

        Ok((
            AnnotatedFastq {
                barcode,
                umi: barcode_result.umi,
                read1,
                read2,
            },
            extraction_stats,
        ))
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

    use crate::qc::Metric;
    use seqspec::{File, Strand, UrlType};
    use std::path::Path;

    fn write_fastq(path: &Path, name: &str, sequence: &[u8]) {
        write_fastq_records(path, &[(name, sequence)]);
    }

    fn write_fastq_records(path: &Path, records: &[(&str, &[u8])]) {
        let content = records
            .iter()
            .map(|(name, sequence)| {
                format!(
                    "@{}\n{}\n+\n{}\n",
                    name,
                    String::from_utf8_lossy(sequence),
                    "I".repeat(sequence.len())
                )
            })
            .collect::<String>();
        std::fs::write(path, content).unwrap();
    }

    fn forward_long_read() -> Vec<u8> {
        [
            b"ATGC".as_slice(),
            vec![b'A'; 24].as_slice(),
            b"ACGTTGCAACGATTCGGAATCC".as_slice(),
            vec![b'C'; 24].as_slice(),
            b"GATCTAGCGTACCTGATCGATGGCATACGTTAA".as_slice(),
            vec![b'G'; 600].as_slice(),
            b"CCGTTAAGGCTACGATTCGGCATGACTAGGTCA".as_slice(),
            vec![b'G'; 24].as_slice(),
            b"TTGACCGATGCTAGCTACGGTA".as_slice(),
            vec![b'T'; 24].as_slice(),
            b"CAGT".as_slice(),
        ]
        .concat()
    }

    fn long_read_assay(fastqs: &[&Path], directory: &Path) -> Assay {
        let whitelist = directory.join("barcodes.txt");
        std::fs::write(
            &whitelist,
            format!("{}\n{}\n", "A".repeat(24), "C".repeat(24)),
        )
        .unwrap();

        let mut assay = Assay::from_path("../seqspec_templates/scNanoATAC.yaml").unwrap();
        assay.sequence_spec.get_mut("R1").unwrap().files = Some(
            fastqs
                .iter()
                .map(|path| File::from_fastq(path, false).unwrap())
                .collect(),
        );

        let modality = assay.library_spec.get_modality(&Modality::ATAC).unwrap();
        for region in &modality.read().unwrap().subregions {
            let mut region = region.write().unwrap();
            if let Some(sequence) = match region.region_id.as_str() {
                "start_linker" => Some("ATGC"),
                "linker_5p" => Some("ACGTTGCAACGATTCGGAATCC"),
                "adapter_5p" => Some("GATCTAGCGTACCTGATCGATGGCATACGTTAA"),
                "adapter_3p" => Some("CCGTTAAGGCTACGATTCGGCATGACTAGGTCA"),
                "linker_3p" => Some("TTGACCGATGCTAGCTACGGTA"),
                "end_linker" => Some("CAGT"),
                _ => None,
            } {
                region.sequence = sequence.to_owned();
                region.min_len = sequence.len() as u32;
                region.max_len = sequence.len() as u32;
            }
            if let Some(onlist) = &mut region.onlist {
                onlist.url = whitelist.to_string_lossy().into_owned();
                onlist.filename = "barcodes.txt".to_owned();
                onlist.urltype = UrlType::Local;
            }
        }
        assay
    }

    fn forward_long_read_with_umi(umi: &[u8], unsupported_topology: bool) -> Vec<u8> {
        let spacer = unsupported_topology.then_some(b"ACGT".as_slice());
        [
            vec![b'A'; 30].as_slice(),
            b"ACTAAAGGCCATTACGGC".as_slice(),
            b"CTACACGACGCTCTTCCGATCT".as_slice(),
            b"AACCGGTTAACCGGTT".as_slice(),
            spacer.unwrap_or_default(),
            umi,
            vec![b'T'; 12].as_slice(),
            vec![b'G'; 600].as_slice(),
            b"TGTACTCTGCGTTGATACCACTGCTT".as_slice(),
        ]
        .concat()
    }

    fn long_read_umi_assay(
        fastq_path: &Path,
        directory: &Path,
        unsupported_topology: bool,
    ) -> Assay {
        use seqspec::region::Region;
        use seqspec::{RegionType, SequenceType};
        use std::sync::{Arc, RwLock};

        let whitelist = directory.join("umi-barcodes.txt");
        std::fs::write(&whitelist, "AACCGGTTAACCGGTT\n").unwrap();

        let mut assay = Assay::from_path("../seqspec_templates/10x_lr_rna_BLAZE.yaml").unwrap();
        assay.sequence_spec.get_mut("R1").unwrap().files =
            Some(vec![File::from_fastq(fastq_path, false).unwrap()]);

        let modality = assay.library_spec.get_modality(&Modality::RNA).unwrap();
        for region in &modality.read().unwrap().subregions {
            let mut region = region.write().unwrap();
            if let Some(onlist) = &mut region.onlist {
                onlist.url = whitelist.to_string_lossy().into_owned();
                onlist.filename = "umi-barcodes.txt".to_owned();
                onlist.urltype = UrlType::Local;
            }
        }
        if unsupported_topology {
            let mut modality = modality.write().unwrap();
            let umi_idx = modality
                .subregions
                .iter()
                .position(|region| region.read().unwrap().region_type.is_umi())
                .unwrap();
            modality.subregions.insert(
                umi_idx,
                Arc::new(RwLock::new(Region {
                    region_id: "umi_spacer".to_string(),
                    region_type: RegionType::Named,
                    name: "UMI spacer".to_string(),
                    sequence_type: SequenceType::Random,
                    sequence: "NNNN".to_string(),
                    min_len: 4,
                    max_len: 4,
                    onlist: None,
                    subregions: Vec::new(),
                })),
            );
        }

        assay
    }

    fn short_read_fixture(directory: &Path) -> (Assay, Vec<u8>, Vec<u8>, Vec<u8>) {
        let read1_path = directory.join("R1.fastq");
        let read2_path = directory.join("R2.fastq");
        let barcode1 = vec![b'A'; 10];
        let barcode2 = vec![b'G'; 10];
        let barcode = [barcode1.clone(), barcode2.clone()].concat();
        let umi = vec![b'T'; 8];
        let read1 = [
            barcode1.as_slice(),
            b"CAGAGC".as_slice(),
            umi.as_slice(),
            barcode2.as_slice(),
            b"TT".as_slice(),
        ]
        .concat();
        let read2 = vec![b'C'; 50];
        write_fastq(&read1_path, "read", &read1);
        write_fastq(&read2_path, "read", &read2);

        let mut assay = Assay::from_path("data/test4.yaml").unwrap();
        assay.sequence_spec.get_mut("R1").unwrap().files =
            Some(vec![File::from_fastq(&read1_path, false).unwrap()]);
        assay.sequence_spec.get_mut("R2").unwrap().files =
            Some(vec![File::from_fastq(&read2_path, false).unwrap()]);

        (assay, barcode, umi, read2)
    }

    fn assert_missing_fastq(input: &str) {
        let assay = Assay::from_path(input).unwrap();
        let modality = assay.modalities[0].clone();
        let error = FastqPlan::new(vec![assay], modality)
            .build(false, 5000)
            .err()
            .unwrap();
        assert!(format!("{error:#}").contains("Failed to open FASTQ file"));
    }

    fn show_annotated_fastq(fq: &AnnotatedFastq) -> String {
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

    fn collect_annotated_fastq(mut execution: FastqExecution) -> Result<Vec<String>> {
        let mut output = Vec::new();
        while let Some(batch) = execution.next_batch()? {
            output.extend(batch.iter().map(show_annotated_fastq));
        }
        execution.finish()?;
        Ok(output)
    }

    /// A missing FASTQ must be reported while building the execution plan.
    #[test]
    fn test_missing_fastq() {
        assert_missing_fastq("data/test2.yaml");
        assert_missing_fastq("data/test3.yaml");
        assert_missing_fastq("data/test4.yaml");
    }

    /// Short-read annotation must preserve barcode, UMI, target, and lifecycle state.
    #[test]
    fn test_short_read() {
        let directory = tempfile::tempdir().unwrap();
        let (assay, barcode, umi, read2) = short_read_fixture(directory.path());

        let mut execution = FastqPlan::new(vec![assay], Modality::RNA)
            .build(false, 5000)
            .unwrap();
        assert!(!execution.is_long_read());
        let batch = execution.next_batch().unwrap().unwrap();
        assert_eq!(batch.len(), 1);
        let annotated = &batch[0];
        assert_eq!(annotated.barcode.as_ref().unwrap().raw.sequence(), barcode);
        assert_eq!(annotated.umi.as_ref().unwrap().sequence(), umi);
        assert!(annotated.read1.is_none());
        assert_eq!(annotated.read2.as_ref().unwrap().sequence(), read2);
        assert!(execution.next_batch().unwrap().is_none());
        assert!(execution.finish().unwrap().long_read.is_none());
    }

    /// A short-read workflow must reject unstranded reads before annotation.
    #[test]
    fn test_short_read_rejects_unstranded() {
        let directory = tempfile::tempdir().unwrap();
        let (mut assay, _, _, _) = short_read_fixture(directory.path());
        assay.sequence_spec.get_mut("R1").unwrap().strand = Strand::Unstranded;

        let error = FastqPlan::new(vec![assay], Modality::RNA)
            .with_sequencing_type(AssayType::ShortRead)
            .build(false, 5000)
            .err()
            .unwrap();
        assert!(error.to_string().contains("cannot use strand=unstranded"));
    }

    /// An explicit sequencing type must take precedence over FASTQ length detection.
    #[test]
    fn test_configured_sequencing_type_overrides_detection() {
        let directory = tempfile::tempdir().unwrap();
        let (assay, _, _, _) = short_read_fixture(directory.path());

        let error = FastqPlan::new(vec![assay], Modality::RNA)
            .with_sequencing_type(AssayType::LongRead)
            .build(false, 5000)
            .err()
            .unwrap();
        assert!(
            error
                .to_string()
                .contains("must have exactly one physical FASTQ Read"),
            "{error:#}"
        );
    }

    /// Compare complete annotated FASTQ output with a deterministic golden row.
    #[test]
    fn test_fastq() {
        let directory = tempfile::tempdir().unwrap();
        let (assay, _, _, _) = short_read_fixture(directory.path());
        let execution = FastqPlan::new(vec![assay], Modality::RNA)
            .build(false, 5000)
            .unwrap();

        let actual = collect_annotated_fastq(execution).unwrap();
        let expected = vec![
            "AAAAAAAAAAGGGGGGGGGG\tTTTTTTTT\t\tCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"
                .to_owned(),
        ];
        assert_eq!(actual, expected);
    }

    /// Forward and reverse long reads must produce directional barcode QC.
    #[test]
    fn test_long_read() {
        let directory = tempfile::tempdir().unwrap();
        let forward_path = directory.path().join("forward.fastq");
        let reverse_path = directory.path().join("reverse.fastq");
        let forward = forward_long_read();
        let reverse = seqspec::utils::rev_compl(&forward);
        write_fastq(&forward_path, "forward", &forward);
        write_fastq(&reverse_path, "reverse", &reverse);

        let assay = long_read_assay(
            &[forward_path.as_path(), reverse_path.as_path()],
            directory.path(),
        );
        let mut execution = FastqPlan::new(vec![assay], Modality::ATAC)
            .build(false, 10_000)
            .unwrap();
        assert!(execution.is_long_read());
        assert!(!execution.is_paired_end());
        assert_eq!(execution.num_records(), 0);
        assert_eq!(execution.read_summary(), "unknown");

        let batch = execution.next_batch().unwrap().unwrap();
        assert_eq!(batch.len(), 2);
        let expected_barcode = [vec![b'A'; 24], vec![b'C'; 24]].concat();

        let forward = &batch[0];
        let forward_barcode = forward.barcode.as_ref().unwrap();
        assert_eq!(forward_barcode.raw.sequence(), expected_barcode);
        assert_eq!(
            forward_barcode.corrected.as_deref(),
            Some(expected_barcode.as_slice())
        );
        assert!(forward_barcode
            .raw
            .quality_scores()
            .iter()
            .all(|base| *base == b'I'));
        assert!(forward.umi.is_none());
        assert!(forward.read2.is_none());
        let forward_target = forward.read1.as_ref().unwrap();
        assert_eq!(forward_target.sequence().len(), 600);
        assert!(forward_target.sequence().iter().all(|base| *base == b'G'));
        assert_eq!(
            forward_target.sequence().len(),
            forward_target.quality_scores().len()
        );

        let reverse = &batch[1];
        assert!(reverse.umi.is_none());
        assert!(reverse.read1.is_none());
        let reverse_target = reverse.read2.as_ref().unwrap();
        assert_eq!(reverse_target.sequence().len(), 600);
        assert!(reverse_target.sequence().iter().all(|base| *base == b'C'));
        assert_eq!(
            reverse_target.sequence().len(),
            reverse_target.quality_scores().len()
        );

        assert!(execution.next_batch().unwrap().is_none());
        let report = execution.finish().unwrap();
        let long_read = report.long_read.unwrap().to_json();
        assert_eq!(long_read["total_reads"], 2);
        assert_eq!(long_read["orientation"]["forward"], 1);
        assert_eq!(long_read["orientation"]["reverse"], 1);
        assert_eq!(long_read["barcode_extraction"]["success"], 2);
        assert!(long_read.get("umi_extraction").is_none());
    }

    #[test]
    fn test_long_read_umi_failure_is_filtered_and_recorded() {
        let directory = tempfile::tempdir().unwrap();
        let fastq_path = directory.path().join("umi.fastq");
        let valid = forward_long_read_with_umi(b"TGCATGCATGCA", false);
        let invalid = forward_long_read_with_umi(b"NNNNNNNNNNNN", false);
        write_fastq_records(
            &fastq_path,
            &[("valid", valid.as_slice()), ("invalid", invalid.as_slice())],
        );
        let assay = long_read_umi_assay(&fastq_path, directory.path(), false);
        let mut execution = FastqPlan::new(vec![assay], Modality::RNA)
            .with_sequencing_type(AssayType::LongRead)
            .build(false, 10_000)
            .unwrap();

        let batch = execution.next_batch().unwrap().unwrap();
        assert_eq!(batch.len(), 1);
        assert_eq!(batch[0].read1.as_ref().unwrap().name(), b"valid");
        assert_eq!(batch[0].umi.as_ref().unwrap().sequence(), b"TGCATGCATGCA");
        assert_eq!(batch[0].read1.as_ref().unwrap().sequence().len(), 600);
        assert!(execution.next_batch().unwrap().is_none());

        let report = execution.finish().unwrap();
        assert_eq!(report.fastq.num_defect.get("R1").copied().unwrap_or(0), 0);
        let long_read = report.long_read.unwrap().to_json();
        assert_eq!(long_read["barcode_extraction"]["success"], 2);
        assert_eq!(long_read["umi_extraction"]["success"], 1);
        assert_eq!(long_read["umi_extraction"]["fail"], 1);
    }

    #[test]
    fn test_unsupported_umi_topology_keeps_reads_without_umi_qc() {
        let directory = tempfile::tempdir().unwrap();
        let fastq_path = directory.path().join("unsupported-umi.fastq");
        let read = forward_long_read_with_umi(b"TGCATGCATGCA", true);
        write_fastq(&fastq_path, "read", &read);
        let assay = long_read_umi_assay(&fastq_path, directory.path(), true);
        let mut execution = FastqPlan::new(vec![assay], Modality::RNA)
            .with_sequencing_type(AssayType::LongRead)
            .build(false, 10_000)
            .unwrap();

        let batch = execution.next_batch().unwrap().unwrap();
        assert_eq!(batch.len(), 1);
        assert!(batch[0].umi.is_none());
        assert_eq!(batch[0].read1.as_ref().unwrap().sequence().len(), 600);
        assert!(execution.next_batch().unwrap().is_none());

        let long_read = execution.finish().unwrap().long_read.unwrap().to_json();
        assert_eq!(long_read["barcode_extraction"]["success"], 1);
        assert!(long_read.get("umi_extraction").is_none());
    }

    #[test]
    fn test_long_read_filters_one_failed_record_and_continues() {
        let directory = tempfile::tempdir().unwrap();
        let fastq_path = directory.path().join("reads.fastq");
        let valid = forward_long_read();
        let failed = vec![b'N'; valid.len()];
        write_fastq_records(
            &fastq_path,
            &[("failed", failed.as_slice()), ("valid", valid.as_slice())],
        );

        let assay = long_read_assay(&[fastq_path.as_path()], directory.path());
        let mut execution = FastqPlan::new(vec![assay], Modality::ATAC)
            .build(false, 10_000)
            .unwrap();

        let batch = execution.next_batch().unwrap().unwrap();
        assert_eq!(batch.len(), 1);
        assert_eq!(batch[0].read1.as_ref().unwrap().name(), b"valid");
        assert!(execution.next_batch().unwrap().is_none());

        let report = execution.finish().unwrap().long_read.unwrap().to_json();
        assert_eq!(report["total_reads"], 2);
        assert_eq!(report["orientation"]["forward"], 1);
        assert_eq!(report["orientation"]["undetermined"], 1);
        assert_eq!(report["composite_alignment"]["pass"], 1);
        assert_eq!(report["composite_alignment"]["fail"], 1);
        assert_eq!(report["barcode_extraction"]["success"], 1);
        assert_eq!(report["barcode_extraction"]["fail"], 1);
    }

    /// End trimming that consumes the whole read must be filtered before alignment.
    #[test]
    fn test_long_read_filters_empty_target_and_continues() {
        let directory = tempfile::tempdir().unwrap();
        let fastq_path = directory.path().join("reads.fastq");
        let valid = forward_long_read();
        let mut empty_target = valid.clone();
        empty_target.drain(107..707);
        assert_eq!(empty_target.len(), 214);
        write_fastq_records(
            &fastq_path,
            &[
                ("empty-target", empty_target.as_slice()),
                ("valid", valid.as_slice()),
            ],
        );

        let assay = long_read_assay(&[fastq_path.as_path()], directory.path());
        let mut execution = FastqPlan::new(vec![assay], Modality::ATAC)
            .build(false, 10_000)
            .unwrap();

        let batch = execution.next_batch().unwrap().unwrap();
        assert_eq!(batch.len(), 1);
        assert_eq!(batch[0].read1.as_ref().unwrap().name(), b"valid");
        assert!(execution.next_batch().unwrap().is_none());

        let report = execution.finish().unwrap();
        assert_eq!(report.fastq.num_reads["R1"], 2);
        assert_eq!(report.fastq.num_defect["R1"], 1);
    }

    struct IdentityStage;

    impl FastqStage for IdentityStage {
        fn process(&mut self, batch: Vec<AnnotatedFastq>) -> Result<Vec<AnnotatedFastq>> {
            Ok(batch)
        }
    }

    struct RecordingAligner {
        layouts: Arc<std::sync::Mutex<Vec<(bool, bool)>>>,
    }

    impl Aligner for RecordingAligner {
        fn header(&self) -> noodles_sam::Header {
            noodles_sam::Header::default()
        }

        fn align_reads(
            &mut self,
            _num_threads: u16,
            records: Vec<AlignmentInput>,
        ) -> Vec<(Option<MultiMapR>, Option<MultiMapR>)> {
            let mut layouts = self.layouts.lock().unwrap();
            records
                .into_iter()
                .map(|record| {
                    layouts.push((record.read1.is_some(), record.read2.is_some()));
                    (None, None)
                })
                .collect()
        }
    }

    /// Alignment streaming must retain mixed read slots and processed counts.
    #[test]
    fn test_alignment_stream() {
        let directory = tempfile::tempdir().unwrap();
        let forward_path = directory.path().join("forward.fastq");
        let reverse_path = directory.path().join("reverse.fastq");
        let forward = forward_long_read();
        write_fastq(&forward_path, "forward", &forward);
        write_fastq(
            &reverse_path,
            "reverse",
            &seqspec::utils::rev_compl(&forward),
        );
        let assay = long_read_assay(
            &[forward_path.as_path(), reverse_path.as_path()],
            directory.path(),
        );
        let execution = FastqPlan::new(vec![assay], Modality::ATAC)
            .build(false, 10_000)
            .unwrap();
        let layouts = Arc::new(std::sync::Mutex::new(Vec::new()));
        let mut aligner = RecordingAligner {
            layouts: layouts.clone(),
        };
        let mut alignments = AlignmentRunner::new(&mut aligner, 1).stream(execution);

        assert_eq!(alignments.num_records(), 0);
        assert_eq!(alignments.next().unwrap().len(), 2);
        assert_eq!(alignments.num_processed(), 2);
        assert!(alignments.next().is_none());
        assert_eq!(alignments.num_processed(), 2);
        assert!(alignments.finish().unwrap().fastq.long_read.is_some());
        assert_eq!(*layouts.lock().unwrap(), [(true, false), (false, true)]);
    }

    /// Invalid long-read plans must fail with actionable validation errors.
    #[test]
    fn test_invalid_plan() {
        let directory = tempfile::tempdir().unwrap();
        let fastq_path = directory.path().join("long.fastq");
        write_fastq(&fastq_path, "forward", &forward_long_read());

        let assay = long_read_assay(&[fastq_path.as_path()], directory.path());
        let error = FastqPlan::new(vec![assay.clone()], Modality::ATAC)
            .with_stage(IdentityStage)
            .build(false, 10_000)
            .err()
            .unwrap();
        assert!(error
            .to_string()
            .contains("middleware is not supported for long-read"));

        let mut assay = assay;
        let mut extra_read = assay.sequence_spec.get("R1").unwrap().clone();
        extra_read.read_id = "I1".to_owned();
        assay
            .sequence_spec
            .insert(extra_read.read_id.clone(), extra_read);
        let error = FastqPlan::new(vec![assay], Modality::ATAC)
            .build(false, 10_000)
            .err()
            .unwrap();
        assert!(error
            .to_string()
            .contains("exactly one physical FASTQ Read"));

        let mut missing = long_read_assay(&[fastq_path.as_path()], directory.path());
        missing
            .sequence_spec
            .get_mut("R1")
            .unwrap()
            .files
            .as_mut()
            .unwrap()[0]
            .url = directory
            .path()
            .join("missing.fastq")
            .to_string_lossy()
            .into_owned();
        let error = FastqPlan::new(vec![missing], Modality::ATAC)
            .build(false, 10_000)
            .err()
            .unwrap();
        assert!(format!("{error:#}").contains("cannot open file"));

        let long = long_read_assay(&[fastq_path.as_path()], directory.path());
        let short_path = directory.path().join("short.fastq");
        write_fastq(&short_path, "short", &vec![b'A'; 100]);
        let short = long_read_assay(&[short_path.as_path()], directory.path());
        let error = FastqPlan::new(vec![long, short], Modality::ATAC)
            .build(false, 10_000)
            .err()
            .unwrap();
        assert!(error
            .to_string()
            .contains("mixing short-read and long-read assays"));
    }
}
