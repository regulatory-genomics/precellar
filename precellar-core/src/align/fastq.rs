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
use indexmap::IndexSet;
use itertools::Itertools;
use log::{debug, info, warn};
use noodles_fastq as fastq;
use rayon::iter::ParallelIterator;
use rayon::slice::ParallelSlice;
use regex::Regex;
use seqspec::{Assay, AssayType, FastqReader, LibSpec, Modality, Read, Region, SegmentInfo};
use smallvec::SmallVec;
use std::collections::{HashMap, HashSet};
use std::sync::{Arc, RwLock};

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

/// Runtime configuration for barcode/UMI extraction from FASTQ read names.
///
/// Use this for already-demultiplexed FASTQs whose barcode and UMI are
/// embedded in the read name (for example LR-splitpipe output) instead of the
/// sequence. The named capture groups of `pattern` map to `region_id`s of the
/// assay's library spec; the combined cell barcode is concatenated in the
/// pattern's capture-group order.
///
/// Construction validates the group-to-region mapping against every supplied
/// assay so that misconfiguration fails at the call site instead of during the
/// alignment run.
#[derive(Debug, Clone)]
pub struct ReadNameBarcodeConfig {
    pattern: Regex,
    read_id: Option<String>,
}

impl ReadNameBarcodeConfig {
    /// Compile `pattern` and validate it against every supplied assay.
    ///
    /// Non-fatal problems (currently: an omitted `read_id` for a modality
    /// with multiple FASTQ reads) are returned as warnings. The caller is
    /// responsible for logging them.
    pub fn new(
        pattern: &str,
        read_id: Option<&str>,
        assays: &[Assay],
        modality: &Modality,
    ) -> Result<(Self, Vec<String>)> {
        let pattern = Regex::new(pattern)
            .with_context(|| format!("invalid read-name barcode pattern: {pattern}"))?;
        let config = Self {
            pattern,
            read_id: read_id.map(str::to_owned),
        };
        let mut warnings = Vec::new();
        for assay in assays {
            config.validate_assay(assay, modality, &mut warnings)?;
        }
        Ok((config, warnings))
    }

    fn validate_assay(
        &self,
        assay: &Assay,
        modality: &Modality,
        warnings: &mut Vec<String>,
    ) -> Result<()> {
        let reads: Vec<&Read> = assay.sequence_spec.values().collect();
        Self::validate_parts(
            &self.pattern,
            self.read_id.as_deref(),
            &assay.library_spec,
            &reads,
            &assay.assay_id,
            modality,
            warnings,
        )
    }

    /// Validate one assay's configuration against its library spec and reads.
    ///
    /// Checks, in order:
    ///
    /// 1. Every named capture group resolves to a direct barcode/umi subregion
    ///    of the modality.
    /// 2. At least one barcode group exists and at most one umi group exists.
    /// 3. Every barcode region of the modality is covered by a group.
    /// 4. A given `read_id` exists, belongs to the modality, and carries FASTQ
    ///    files; an omitted `read_id` with multiple FASTQ reads produces a
    ///    warning.
    fn validate_parts(
        pattern: &Regex,
        read_id: Option<&str>,
        lib_spec: &LibSpec,
        reads: &[&Read],
        assay_id: &str,
        modality: &Modality,
        warnings: &mut Vec<String>,
    ) -> Result<()> {
        let modality_region = lib_spec.get_modality(modality).ok_or_else(|| {
            anyhow::anyhow!(
                "assay '{}' does not contain modality {}",
                assay_id,
                modality
            )
        })?;
        let modality_guard = modality_region.read().unwrap();
        let subregions = &modality_guard.subregions;

        let mut barcode_regions: Vec<String> = Vec::new();
        for region in subregions {
            let region = region.read().unwrap();
            if region.region_type.is_barcode() {
                barcode_regions.push(region.region_id.clone());
            }
        }

        let group_names: Vec<String> = pattern
            .capture_names()
            .flatten()
            .map(str::to_owned)
            .collect();
        let mut num_barcode_groups = 0usize;
        let mut num_umi_groups = 0usize;
        for name in &group_names {
            let region = find_field_region(subregions, name, assay_id, modality)?;
            let region = region.read().unwrap();
            anyhow::ensure!(
                region.region_type.is_barcode() || region.region_type.is_umi(),
                "read-name capture group '{}' maps to region '{}' of type {} in assay '{}'; \
                 groups must map to barcode or umi regions",
                name,
                region.region_id,
                region.region_type,
                assay_id
            );
            if region.region_type.is_barcode() {
                num_barcode_groups += 1;
            } else {
                num_umi_groups += 1;
            }
        }
        drop(modality_guard);

        anyhow::ensure!(
            num_barcode_groups >= 1,
            "read-name pattern for assay '{}' must contain at least one barcode capture group",
            assay_id
        );
        anyhow::ensure!(
            num_umi_groups <= 1,
            "read-name pattern for assay '{}' contains {} umi capture groups; at most 1 is supported",
            assay_id,
            num_umi_groups
        );

        let missing: Vec<String> = barcode_regions
            .iter()
            .filter(|id| !group_names.contains(id))
            .cloned()
            .collect();
        anyhow::ensure!(
            missing.is_empty(),
            "read-name groups do not cover all barcode regions of assay '{}' modality {}: missing {}",
            assay_id,
            modality,
            missing.join(", ")
        );

        if let Some(read_id) = read_id {
            let read = reads
                .iter()
                .find(|read| read.read_id == read_id)
                .ok_or_else(|| {
                    anyhow::anyhow!("read_id '{}' not found in assay '{}'", read_id, assay_id)
                })?;
            anyhow::ensure!(
                read.modality == *modality,
                "read '{}' of assay '{}' belongs to modality {}, not {}",
                read_id,
                assay_id,
                read.modality,
                modality
            );
            anyhow::ensure!(
                has_fastq_input(read),
                "read '{}' of assay '{}' has no fastq files",
                read_id,
                assay_id
            );
        } else {
            let num_fastq_reads = reads
                .iter()
                .filter(|read| read.modality == *modality && has_fastq_input(read))
                .count();
            if num_fastq_reads > 1 {
                warnings.push(format!(
                    "read_id omitted for assay '{}': modality {} has {} reads with fastq input; \
                     every read parses its name and the parsed barcodes are concatenated once per read",
                    assay_id, modality, num_fastq_reads
                ));
            }
        }
        Ok(())
    }
}

/// Whether the read declares at least one FASTQ input file.
fn has_fastq_input(read: &Read) -> bool {
    read.files
        .as_ref()
        .is_some_and(|files| files.iter().any(|file| file.filetype == "fastq"))
}

/// Finds the modality subregion that a read-name capture group maps to.
fn find_field_region(
    subregions: &[Arc<RwLock<Region>>],
    name: &str,
    assay_id: &str,
    modality: &Modality,
) -> Result<Arc<RwLock<Region>>> {
    subregions
        .iter()
        .find(|region| region.read().unwrap().region_id == name)
        .cloned()
        .ok_or_else(|| {
            anyhow::anyhow!(
                "read-name capture group '{}' does not match any region of assay '{}' modality {}",
                name,
                assay_id,
                modality
            )
        })
}

/// Sample size for the plan-build smoke test.
const SMOKE_TEST_SAMPLE_SIZE: usize = 1_000;

/// Structural threshold for the plan-build smoke test: the fraction of
/// sampled records that must annotate (split or extract) successfully.
const SMOKE_TEST_STRUCTURAL_THRESHOLD: f64 = 0.7;

/// Content threshold for the plan-build smoke test: the fraction of parsed
/// records whose barcode fields must hit their whitelists.
const SMOKE_TEST_CONTENT_THRESHOLD: f64 = 0.9;

/// Samples records from every modality FASTQ read that has an annotator and
/// feeds each to `check`. Best effort: reads without FASTQ input are skipped.
fn sample_reads(
    assay: &Assay,
    modality: &Modality,
    sample_size: usize,
    mut check: impl FnMut(&FastqAnnotator, &fastq::Record),
) -> Result<()> {
    for (read, segment_info) in assay.get_segments_by_modality(*modality) {
        let Some(annotator) = FastqAnnotator::new(&read.read_id, segment_info) else {
            continue;
        };
        let Some(mut reader) = read.try_open().with_context(|| {
            format!(
                "failed to open FASTQ input for assay '{}' read '{}'",
                assay.assay_id, read.read_id
            )
        })?
        else {
            continue;
        };
        let mut record = fastq::Record::default();
        for _ in 0..sample_size {
            if reader.read_record(&mut record)? == 0 {
                break;
            }
            check(&annotator, &record);
        }
    }
    Ok(())
}

/// Runs a sampling smoke test of one assay against its configured processor.
///
/// The test guards against data that does not fit the configuration: it
/// samples records from every modality FASTQ read, runs the real extraction
/// path, and returns one warning per metric below its threshold. Structural
/// metrics (annotation success) use [`SMOKE_TEST_STRUCTURAL_THRESHOLD`];
/// content metrics (parse and whitelist-hit ratios) use
/// [`SMOKE_TEST_CONTENT_THRESHOLD`]. The caller is responsible for logging
/// the returned warnings.
fn smoke_test_assay(
    assay: &Assay,
    modality: &Modality,
    barcode_processor: &BarcodeProcessor,
) -> Result<Vec<String>> {
    let mut warnings = Vec::new();
    match barcode_processor {
        BarcodeProcessor::ShortRead(_) => {
            let mut sampled = 0u64;
            let mut succeeded = 0u64;
            let mut onlist_checked = 0u64;
            let mut onlist_hit = 0u64;

            // Region onlists for exact-match content checks.
            let mut onlists: HashMap<String, IndexSet<Vec<u8>>> = HashMap::new();
            let modality_region = assay.library_spec.get_modality(modality).ok_or_else(|| {
                anyhow::anyhow!(
                    "assay '{}' does not contain modality {}",
                    assay.assay_id,
                    modality
                )
            })?;
            for region in &modality_region.read().unwrap().subregions {
                let region = region.read().unwrap();
                if let (true, Some(onlist)) = (region.region_type.is_barcode(), &region.onlist) {
                    onlists.insert(region.region_id.clone(), onlist.read()?);
                }
            }

            sample_reads(
                assay,
                modality,
                SMOKE_TEST_SAMPLE_SIZE,
                |annotator, record| {
                    sampled += 1;
                    let Ok(segments) = annotator.segment_info.split(record) else {
                        return;
                    };
                    succeeded += 1;
                    for segment in &segments {
                        if !segment.is_barcode() {
                            continue;
                        }
                        let Some(onlist) = onlists.get(segment.region_id()) else {
                            continue;
                        };
                        onlist_checked += 1;
                        let sequence = if annotator.segment_info.is_reverse() {
                            seqspec::utils::rev_compl(segment.seq)
                        } else {
                            segment.seq.to_vec()
                        };
                        if onlist.contains(&sequence) {
                            onlist_hit += 1;
                        }
                    }
                },
            )?;

            if sampled > 0 {
                let ratio = succeeded as f64 / sampled as f64;
                if ratio < SMOKE_TEST_STRUCTURAL_THRESHOLD {
                    warnings.push(format!(
                        "smoke test for assay '{}': structural ratio {:.3} is below {:.2} \
                         ({succeeded}/{sampled} records annotated); the seqspec layout may not match the data",
                        assay.assay_id, ratio, SMOKE_TEST_STRUCTURAL_THRESHOLD
                    ));
                }
            }
            if onlist_checked > 0 {
                let ratio = onlist_hit as f64 / onlist_checked as f64;
                if ratio < SMOKE_TEST_CONTENT_THRESHOLD {
                    warnings.push(format!(
                        "smoke test for assay '{}': barcode onlist hit ratio {:.3} is below {:.2} \
                         ({onlist_hit}/{onlist_checked})",
                        assay.assay_id, ratio, SMOKE_TEST_CONTENT_THRESHOLD
                    ));
                }
            }
        }
        BarcodeProcessor::LongRead(extractor) => {
            let mut sampled = 0u64;
            let mut extracted = 0u64;
            sample_reads(
                assay,
                modality,
                SMOKE_TEST_SAMPLE_SIZE,
                |_annotator, record| {
                    sampled += 1;
                    if let Ok((result, _)) = extractor.extract_barcode(record) {
                        if result.is_success() {
                            extracted += 1;
                        }
                    }
                },
            )?;

            if sampled > 0 {
                let ratio = extracted as f64 / sampled as f64;
                if ratio < SMOKE_TEST_STRUCTURAL_THRESHOLD {
                    warnings.push(format!(
                        "smoke test for assay '{}': barcode extraction ratio {:.3} is below {:.2} \
                         ({extracted}/{sampled} records)",
                        assay.assay_id, ratio, SMOKE_TEST_STRUCTURAL_THRESHOLD
                    ));
                }
            }
        }
        BarcodeProcessor::ReadName(parser) => {
            let mut sampled = 0u64;
            let mut parsed = 0u64;
            // region id -> (checked, hit)
            let mut whitelists: HashMap<String, (u64, u64)> = HashMap::new();
            sample_reads(
                assay,
                modality,
                SMOKE_TEST_SAMPLE_SIZE,
                |annotator, record| {
                    if !parser.applies_to(&annotator.read_id) {
                        return;
                    }
                    sampled += 1;
                    if let Ok((_, Some(BarcodeStats::ReadName(stats)))) =
                        annotator.annotate(record, barcode_processor)
                    {
                        if stats.parsed {
                            parsed += 1;
                        }
                        for (region_id, hit) in &stats.region_hits {
                            let entry = whitelists.entry(region_id.clone()).or_insert((0, 0));
                            entry.0 += 1;
                            if *hit {
                                entry.1 += 1;
                            }
                        }
                    }
                },
            )?;

            if sampled > 0 {
                let ratio = parsed as f64 / sampled as f64;
                if ratio < SMOKE_TEST_CONTENT_THRESHOLD {
                    warnings.push(format!(
                        "smoke test for assay '{}': read-name parse ratio {:.3} is below {:.2} \
                         ({parsed}/{sampled}); check the pattern against the read names",
                        assay.assay_id, ratio, SMOKE_TEST_CONTENT_THRESHOLD
                    ));
                }
            }
            for (region_id, (checked, hit)) in &whitelists {
                if *checked == 0 {
                    continue;
                }
                let ratio = *hit as f64 / *checked as f64;
                if ratio < SMOKE_TEST_CONTENT_THRESHOLD {
                    warnings.push(format!(
                        "smoke test for assay '{}': region '{}' whitelist hit ratio {:.3} is below {:.2} \
                         ({hit}/{checked}); the capture group may be mapped to the wrong region",
                        assay.assay_id, region_id, ratio, SMOKE_TEST_CONTENT_THRESHOLD
                    ));
                }
            }
        }
    }
    Ok(warnings)
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
    read_name_barcode: Option<ReadNameBarcodeConfig>,
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
            read_name_barcode: None,
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

    /// Installs the read-name barcode/UMI extraction configuration.
    ///
    /// When set, the plan's processor parses barcode/UMI fields from FASTQ
    /// read names instead of the sequence; see [`ReadNameBarcodeConfig`].
    /// The configuration takes precedence over the sequencing-type branch
    /// selection in [`Self::build`].
    pub fn with_read_name_barcode(mut self, config: ReadNameBarcodeConfig) -> Self {
        self.read_name_barcode = Some(config);
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

            let barcode_processor = if let Some(config) = &self.read_name_barcode {
                BarcodeProcessor::ReadName(
                    ReadNameParser::new(config, assay, &self.modality).with_context(|| {
                        format!(
                            "failed to configure read-name barcode extraction for assay '{}'",
                            assay.assay_id
                        )
                    })?,
                )
            } else {
                match assay_type {
                    AssayType::ShortRead => {
                        let mut barcode_analyzer = BarcodeAnalyzer::new(assay, self.modality);
                        barcode_analyzer.summary();
                        if correct_barcode {
                            barcode_analyzer.barcode_correct_options =
                                Some(BarcodeCorrectOptions {
                                    bc_confidence_threshold: self
                                        .barcode_config
                                        .confidence_threshold,
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
                }
            };

            for warning in smoke_test_assay(assay, &self.modality, &barcode_processor)? {
                warn!("{warning}");
            }

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
    ReadName(ReadNameParser),
}

impl BarcodeProcessor {
    fn num_reads(&self) -> usize {
        match self {
            Self::ShortRead(analyzer) => analyzer.num_reads(),
            Self::LongRead(_) | Self::ReadName(_) => 0,
        }
    }

    /// Whether this processor walks the composite-alignment extraction path.
    ///
    /// The name is historical: this flag gates behavior that only makes sense
    /// for composite-alignment extraction (paired-end handling, long-read QC,
    /// the `make_fastq` bail), not the sequencing modality of the data. A
    /// read-name processor must return `false` even for long-read data.
    fn is_long_read(&self) -> bool {
        matches!(self, Self::LongRead(_))
    }

    fn expects_umi(&self) -> bool {
        match self {
            Self::ShortRead(_) => false,
            Self::LongRead(extractor) => extractor.expects_umi(),
            Self::ReadName(parser) => parser.has_umi(),
        }
    }
}

/// One resolved named-capture group of a read-name pattern.
///
/// `region_id` is both the capture-group name and the library-spec region the
/// group maps to. Fields are stored in the pattern's capture-group order,
/// which is the order the combined cell barcode is concatenated in.
#[derive(Debug)]
struct FieldRegion {
    region_id: String,
    /// Whether the field is a UMI (false = barcode).
    is_umi: bool,
    /// Whitelist for exact-match statistics; `None` when the region declares
    /// no onlist and the field is therefore excluded from whitelist stats.
    whitelist: Option<IndexSet<Vec<u8>>>,
}

/// Extracts barcode/UMI fields from the read names of an already
/// demultiplexed assay (LR-splitpipe and similar outputs).
///
/// The named capture groups of the pattern map to barcode/umi regions of the
/// modality. The whole FASTQ sequence is the target; the pattern is matched
/// against the read name only (never the description).
#[derive(Debug)]
struct ReadNameParser {
    /// Validated user configuration.
    config: ReadNameBarcodeConfig,
    /// Resolved fields in the pattern's capture-group order.
    field_regions: Vec<FieldRegion>,
}

/// Statistics from read-name barcode extraction of one physical record.
#[derive(Debug, Default)]
pub struct ReadNameStats {
    /// Whether the pattern matched the read name and every declared group
    /// participated in the match.
    pub parsed: bool,
    /// Per-region exact whitelist hits; only regions with an onlist appear.
    pub region_hits: Vec<(String, bool)>,
}

/// Statistics produced by whichever barcode extraction strategy ran.
#[derive(Debug)]
pub enum BarcodeStats {
    Composite(LrBarcodeExtractionStats),
    ReadName(ReadNameStats),
}

impl ReadNameParser {
    /// Resolves the config's named capture groups against one assay's library
    /// spec. Field regions are stored in the pattern's capture-group order.
    fn new(config: &ReadNameBarcodeConfig, assay: &Assay, modality: &Modality) -> Result<Self> {
        let modality_region = assay.library_spec.get_modality(modality).ok_or_else(|| {
            anyhow::anyhow!(
                "assay '{}' does not contain modality {}",
                assay.assay_id,
                modality
            )
        })?;
        let modality_guard = modality_region.read().unwrap();

        let mut field_regions = Vec::new();
        for name in config.pattern.capture_names().flatten() {
            let region =
                find_field_region(&modality_guard.subregions, name, &assay.assay_id, modality)?;
            let region = region.read().unwrap();
            let whitelist = match &region.onlist {
                Some(onlist) => Some(onlist.read().with_context(|| {
                    format!(
                        "failed to read onlist for region '{}' of assay '{}'",
                        region.region_id, assay.assay_id
                    )
                })?),
                None => None,
            };
            field_regions.push(FieldRegion {
                is_umi: region.region_type.is_umi(),
                region_id: region.region_id.clone(),
                whitelist,
            });
        }
        drop(modality_guard);

        Ok(Self {
            config: config.clone(),
            field_regions,
        })
    }

    /// Whether this physical read parses its name.
    fn applies_to(&self, read_id: &str) -> bool {
        self.config
            .read_id
            .as_deref()
            .is_none_or(|id| id == read_id)
    }

    /// Whether the pattern declares a umi group.
    fn has_umi(&self) -> bool {
        self.field_regions.iter().any(|field| field.is_umi)
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
        if let Some(qc_long_read) = qc_long_read.as_mut() {
            qc_long_read.enable_umi_extraction();
        }
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
                        Ok((anno, stats)) => {
                            match stats {
                                Some(BarcodeStats::Composite(stats)) => {
                                    if let Some(qc_long_read) = qc_long_read.as_mut() {
                                        qc_long_read.record_orientation(
                                            stats.is_reverse,
                                            !stats.composite_pass,
                                        );
                                        qc_long_read
                                            .record_composite_alignment(stats.composite_pass);
                                        qc_long_read
                                            .record_barcode_extraction(anno.barcode.is_some());
                                        if expects_umi {
                                            qc_long_read.record_umi_extraction(anno.umi.is_some());
                                        }
                                        for (group_name, intersection_hit) in
                                            stats.consensus_results
                                        {
                                            qc_long_read
                                                .record_consensus(&group_name, intersection_hit);
                                        }
                                    }
                                }
                                Some(BarcodeStats::ReadName(stats)) => {
                                    qc.record_read_name_stats(id, &stats);
                                }
                                None => {}
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
    ) -> Result<(AnnotatedFastq, Option<BarcodeStats>)> {
        match barcode_processor {
            BarcodeProcessor::ShortRead(barcode_analyzer) => {
                Ok((self.annotate_short_read(record, barcode_analyzer)?, None))
            }
            BarcodeProcessor::LongRead(barcode_extractor) => {
                let (annotated, stats) = self.annotate_long_read(record, barcode_extractor)?;
                Ok((annotated, Some(BarcodeStats::Composite(stats))))
            }
            BarcodeProcessor::ReadName(parser) => {
                if parser.applies_to(&self.read_id) {
                    let (annotated, stats) = self.annotate_from_read_name(record, parser)?;
                    Ok((annotated, Some(BarcodeStats::ReadName(stats))))
                } else {
                    Ok((self.annotate_target_only(record)?, None))
                }
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

    /// The whole record as the target, assigned to read 1 or read 2 by the
    /// read's strand.
    fn annotate_target_only(&self, record: &fastq::Record) -> Result<AnnotatedFastq> {
        let (read1, read2) = if self.segment_info.is_reverse() {
            (None, Some(record.clone()))
        } else {
            (Some(record.clone()), None)
        };
        Ok(AnnotatedFastq {
            barcode: None,
            umi: None,
            read1,
            read2,
        })
    }

    /// Annotates a record whose barcode/UMI live in the read name.
    ///
    /// The whole sequence is the target, assigned to read 1 or read 2 by the
    /// read's strand. On a pattern miss — or when an optional group drops out
    /// of the match — the barcode and UMI are `None` and the record is dropped
    /// by the usual barcode-less filter; the target is still produced so QC
    /// sees the full record, mirroring the long-read failure path.
    ///
    /// Barcode and UMI quality scores are synthetic (`I`) because the read
    /// name carries no base qualities, and the corrected barcode is the raw
    /// one: the header barcode is the demultiplexer's canonical output.
    fn annotate_from_read_name(
        &self,
        record: &fastq::Record,
        parser: &ReadNameParser,
    ) -> Result<(AnnotatedFastq, ReadNameStats)> {
        let miss_stats = ReadNameStats {
            parsed: false,
            region_hits: Vec::new(),
        };
        // Read names are `&BStr`; a non-UTF-8 name cannot match the pattern.
        let Ok(name) = std::str::from_utf8(record.name()) else {
            return Ok((self.annotate_target_only(record)?, miss_stats));
        };
        let Some(captures) = parser.config.pattern.captures(name) else {
            return Ok((self.annotate_target_only(record)?, miss_stats));
        };

        let mut barcode_seq: Vec<u8> = Vec::new();
        let mut umi: Option<fastq::Record> = None;
        let mut region_hits = Vec::new();
        for field in &parser.field_regions {
            let Some(sequence) = captures.name(&field.region_id) else {
                return Ok((self.annotate_target_only(record)?, miss_stats));
            };
            let sequence = sequence.as_str();
            if let Some(whitelist) = &field.whitelist {
                let hit = whitelist.contains(sequence.as_bytes());
                region_hits.push((field.region_id.clone(), hit));
            }
            if field.is_umi {
                umi = Some(fastq::Record::new(
                    record.definition().clone(),
                    sequence.as_bytes().to_vec(),
                    vec![b'I'; sequence.len()],
                ));
            } else {
                barcode_seq.extend_from_slice(sequence.as_bytes());
            }
        }

        let barcode = Barcode {
            raw: fastq::Record::new(
                record.definition().clone(),
                barcode_seq.clone(),
                vec![b'I'; barcode_seq.len()],
            ),
            corrected: Some(barcode_seq),
        };
        let mut annotated = self.annotate_target_only(record)?;
        annotated.barcode = Some(barcode);
        annotated.umi = umi;

        Ok((
            annotated,
            ReadNameStats {
                parsed: true,
                region_hits,
            },
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

    /// Tests for [`ReadNameBarcodeConfig`] validation (read-name barcode
    /// feature, configuration step).
    mod read_name_config {
        use anyhow::Result;
        use regex::Regex;
        use seqspec::{File, LibSpec, Modality, Read, Region, RegionType, SequenceType, UrlType};
        use std::sync::{Arc, RwLock};

        use crate::align::fastq::ReadNameBarcodeConfig;

        fn region(
            id: &str,
            region_type: RegionType,
            min_len: u32,
            max_len: u32,
        ) -> Arc<RwLock<Region>> {
            let sequence_type = match region_type {
                RegionType::Barcode | RegionType::Umi => SequenceType::Onlist,
                _ => SequenceType::Random,
            };
            Arc::new(RwLock::new(Region {
                region_id: id.to_string(),
                region_type,
                name: id.to_string(),
                sequence_type,
                sequence: String::new(),
                min_len,
                max_len,
                onlist: None,
                subregions: vec![],
            }))
        }

        /// Library layout: bc1(8) -> bc2(8) -> umi(10) -> umi2(10) -> cdna.
        fn lib_spec() -> LibSpec {
            LibSpec::new(vec![Region {
                region_id: "rna".to_string(),
                region_type: RegionType::Modality(Modality::RNA),
                name: "RNA".to_string(),
                sequence_type: SequenceType::Joined,
                sequence: String::new(),
                min_len: 0,
                max_len: 0,
                onlist: None,
                subregions: vec![
                    region("bc1", RegionType::Barcode, 8, 8),
                    region("bc2", RegionType::Barcode, 8, 8),
                    region("umi", RegionType::Umi, 10, 10),
                    region("umi2", RegionType::Umi, 10, 10),
                    region("cdna", RegionType::Cdna, 100, 1000),
                ],
            }])
            .unwrap()
        }

        fn fastq_file(name: &str) -> File {
            File {
                file_id: name.to_string(),
                filename: name.to_string(),
                filetype: "fastq".to_string(),
                filesize: 0,
                url: name.to_string(),
                urltype: UrlType::Local,
                md5: String::new(),
            }
        }

        fn read(id: &str, modality: Modality, with_fastq: bool) -> Read {
            Read {
                read_id: id.to_string(),
                modality,
                files: with_fastq.then(|| vec![fastq_file(&format!("{id}.fastq.gz"))]),
                ..Default::default()
            }
        }

        fn validate(pattern: &str, read_id: Option<&str>, reads: &[Read]) -> Result<Vec<String>> {
            let pattern = Regex::new(pattern).unwrap();
            let lib_spec = lib_spec();
            let reads: Vec<&Read> = reads.iter().collect();
            let mut warnings = Vec::new();
            ReadNameBarcodeConfig::validate_parts(
                &pattern,
                read_id,
                &lib_spec,
                &reads,
                "test_assay",
                &Modality::RNA,
                &mut warnings,
            )?;
            Ok(warnings)
        }

        #[test]
        fn valid_config_passes_silently() {
            let warnings = validate(
                ":(?P<bc1>[ACGTN]{8})_(?P<bc2>[ACGTN]{8})_(?P<umi>[ACGTN]{10})$",
                None,
                &[read("R1", Modality::RNA, true)],
            )
            .unwrap();
            assert!(warnings.is_empty());
        }

        /// Group names must resolve to barcode/umi regions of the modality.
        #[test]
        fn invalid_group_mapping_fails() {
            let unknown = validate(
                "(?P<bcx>[ACGTN]{8})",
                None,
                &[read("R1", Modality::RNA, true)],
            )
            .unwrap_err();
            assert!(unknown.to_string().contains("does not match any region"));

            let non_barcode = validate(
                "(?P<cdna>[ACGTN]{8})",
                None,
                &[read("R1", Modality::RNA, true)],
            )
            .unwrap_err();
            assert!(non_barcode
                .to_string()
                .contains("must map to barcode or umi regions"));
        }

        /// All barcode regions must be covered, and at most one umi group is
        /// allowed.
        #[test]
        fn structural_rule_violations_fail() {
            let uncovered = validate(
                "(?P<bc1>[ACGTN]{8})",
                None,
                &[read("R1", Modality::RNA, true)],
            )
            .unwrap_err();
            let message = uncovered.to_string();
            assert!(message.contains("do not cover all barcode regions"));
            assert!(message.contains("bc2"));

            let two_umis = validate(
                "(?P<bc1>[ACGTN]{8})_(?P<bc2>[ACGTN]{8})_(?P<umi>[ACGTN]{10})_(?P<umi2>[ACGTN]{10})",
                None,
                &[read("R1", Modality::RNA, true)],
            )
            .unwrap_err();
            assert!(two_umis.to_string().contains("2 umi capture groups"));
        }

        /// A given read_id must exist, belong to the modality, and carry fastq
        /// files. An omitted read_id with multiple fastq reads warns.
        #[test]
        fn read_id_rules() {
            let not_found = validate(
                ":(?P<bc1>[ACGTN]{8})_(?P<bc2>[ACGTN]{8})$",
                Some("R9"),
                &[read("R1", Modality::RNA, true)],
            )
            .unwrap_err();
            assert!(not_found.to_string().contains("read_id 'R9' not found"));

            let other_modality = validate(
                ":(?P<bc1>[ACGTN]{8})_(?P<bc2>[ACGTN]{8})$",
                Some("R_atac"),
                &[
                    read("R1", Modality::RNA, true),
                    read("R_atac", Modality::ATAC, true),
                ],
            )
            .unwrap_err();
            assert!(other_modality.to_string().contains("belongs to modality"));

            let no_fastq = validate(
                ":(?P<bc1>[ACGTN]{8})_(?P<bc2>[ACGTN]{8})$",
                Some("R_nofq"),
                &[
                    read("R1", Modality::RNA, true),
                    read("R_nofq", Modality::RNA, false),
                ],
            )
            .unwrap_err();
            assert!(no_fastq.to_string().contains("has no fastq files"));

            let warnings = validate(
                ":(?P<bc1>[ACGTN]{8})_(?P<bc2>[ACGTN]{8})$",
                None,
                &[
                    read("R1", Modality::RNA, true),
                    read("R2", Modality::RNA, true),
                ],
            )
            .unwrap();
            assert_eq!(warnings.len(), 1);
            assert!(warnings[0].contains("read_id omitted"));
        }
    }

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

    /// Attaches `records` (read name, target sequence) to the read-name assay
    /// (`data/test_readname.yaml`) as read R1.
    fn read_name_assay(directory: &Path, records: &[(&str, &[u8])]) -> Assay {
        let path = directory.join("R1.fastq");
        write_fastq_records(&path, records);
        let mut assay = Assay::from_path("data/test_readname.yaml").unwrap();
        assay.sequence_spec.get_mut("R1").unwrap().files =
            Some(vec![File::from_fastq(&path, false).unwrap()]);
        assay
    }

    /// Read-name extraction must concatenate barcodes in the pattern's
    /// capture-group order, and must pass the whole sequence through as the
    /// target. PacBio, ONT, and short prefixes are all handled by the same
    /// pattern.
    #[test]
    fn test_read_name_extraction() {
        let directory = tempfile::tempdir().unwrap();
        let bc1 = b"AAAAAAAA".as_slice();
        let bc2 = b"CCCCCCCC".as_slice();
        let umi = b"ACGTACGTAC".as_slice();
        let target: Vec<u8> = vec![b'C'; 600];
        // The pattern's capture order (umi, bc2, bc1) differs from the
        // library-spec order (bc1, bc2, umi).
        let records = [
            (
                "m54284U_210909_065506/102302155/ccs:ACGTACGTAC_CCCCCCCC_AAAAAAAA",
                &target[..],
            ),
            (
                "be0ebe90-9d3a-4cc6-8a5e-ce17ba8b4094:ACGTACGTAC_CCCCCCCC_AAAAAAAA",
                &target[..],
            ),
            ("read1:ACGTACGTAC_CCCCCCCC_AAAAAAAA", &target[..]),
        ];
        let assay = read_name_assay(directory.path(), &records);
        let (config, warnings) = ReadNameBarcodeConfig::new(
            ":(?P<umi>[ACGTN]{10})_(?P<bc2>[ACGTN]{8})_(?P<bc1>[ACGTN]{8})$",
            Some("R1"),
            std::slice::from_ref(&assay),
            &Modality::RNA,
        )
        .unwrap();
        assert!(warnings.is_empty());

        let mut execution = FastqPlan::new(vec![assay], Modality::RNA)
            .with_read_name_barcode(config)
            .build(false, 5000)
            .unwrap();
        assert!(!execution.is_long_read());
        let batch = execution.next_batch().unwrap().unwrap();
        assert_eq!(batch.len(), 3);
        // Pattern capture order: umi first, then bc2, then bc1. The barcode
        // therefore concatenates bc2 before bc1.
        let expected_barcode: Vec<u8> = [bc2, bc1].concat();
        for annotated in &batch {
            let barcode = annotated.barcode.as_ref().unwrap();
            assert_eq!(barcode.raw.sequence(), &expected_barcode);
            assert_eq!(
                barcode.corrected.as_deref(),
                Some(expected_barcode.as_slice())
            );
            assert_eq!(barcode.raw.quality_scores(), vec![b'I'; 16]);
            assert_eq!(annotated.umi.as_ref().unwrap().sequence(), umi);
            assert_eq!(
                annotated.umi.as_ref().unwrap().quality_scores(),
                vec![b'I'; 10]
            );
            assert!(annotated.read2.is_none());
            assert_eq!(annotated.read1.as_ref().unwrap().sequence(), &target);
        }
        assert!(execution.next_batch().unwrap().is_none());
        execution.finish().unwrap();
    }

    /// Parse failures are filtered by the usual barcode-less gate, and the QC
    /// counters separate name parsing from whitelist membership.
    #[test]
    fn test_read_name_qc_and_filtering() {
        let directory = tempfile::tempdir().unwrap();
        // data/bc1_readname.txt contains AAAAAAAA and GGGGGGGG but not CCCCCCCC.
        let target: Vec<u8> = vec![b'C'; 600];
        let records = [
            ("read1:AAAAAAAA_GGGGGGGG_ACGTACGTAC", &target[..]),
            ("read2:CCCCCCCC_GGGGGGGG_ACGTACGTAC", &target[..]),
            ("read3_no_payload", &target[..]),
        ];
        let assay = read_name_assay(directory.path(), &records);
        let (config, _) = ReadNameBarcodeConfig::new(
            ":(?P<bc1>[ACGTN]{8})_(?P<bc2>[ACGTN]{8})_(?P<umi>[ACGTN]{10})$",
            Some("R1"),
            std::slice::from_ref(&assay),
            &Modality::RNA,
        )
        .unwrap();
        let mut execution = FastqPlan::new(vec![assay], Modality::RNA)
            .with_read_name_barcode(config)
            .build(false, 5000)
            .unwrap();
        let batch = execution.next_batch().unwrap().unwrap();
        // The unparsed record is dropped; the two parsed records survive.
        assert_eq!(batch.len(), 2);
        assert!(execution.next_batch().unwrap().is_none());
        let report = execution.finish().unwrap();
        assert_eq!(report.fastq.num_reads["R1"], 3);
        assert_eq!(report.fastq.num_readname_parsed["R1"], 2);
        assert_eq!(report.fastq.num_readname_in_whitelist["bc1"], 1);
        // bc2 declares no onlist and is absent from the whitelist stats.
        assert!(!report.fastq.num_readname_in_whitelist.contains_key("bc2"));
        let json = report.fastq.to_json();
        assert_eq!(
            json["frac_readname_parsed"]["R1"],
            serde_json::json!(2.0 / 3.0)
        );
        assert_eq!(
            json["frac_readname_in_whitelist"]["bc1"],
            serde_json::json!(0.5)
        );
        assert!(json["frac_readname_in_whitelist"].get("bc2").is_none());
    }

    /// A pattern that matches nothing must still surface a 0.0 parse ratio —
    /// the misconfiguration signal must not vanish from the QC report.
    #[test]
    fn test_read_name_unparsed_ratio_is_visible() {
        let directory = tempfile::tempdir().unwrap();
        let target: Vec<u8> = vec![b'C'; 600];
        let records = [("read1:AAAAAAAA_GGGGGGGG_ACGTACGTAC", &target[..])];
        let assay = read_name_assay(directory.path(), &records);
        // The trailing literal `_ZZ` never occurs, so no name matches.
        let (config, _) = ReadNameBarcodeConfig::new(
            ":(?P<bc1>[ACGTN]{8})_(?P<bc2>[ACGTN]{8})_(?P<umi>[ACGTN]{10})_ZZ$",
            Some("R1"),
            std::slice::from_ref(&assay),
            &Modality::RNA,
        )
        .unwrap();
        let mut execution = FastqPlan::new(vec![assay], Modality::RNA)
            .with_read_name_barcode(config)
            .build(false, 5000)
            .unwrap();
        assert!(execution.next_batch().unwrap().is_none());
        let json = execution.finish().unwrap().fastq.to_json();
        assert_eq!(json["frac_readname_parsed"]["R1"], serde_json::json!(0.0));
        // No record was parsed, so whitelist ratios cannot exist.
        assert!(json.get("frac_readname_in_whitelist").is_none());
    }

    /// Section 4: the plan-build smoke test must warn when short-read data
    /// does not match the seqspec layout.
    #[test]
    fn test_smoke_test_short_read_warns_on_wrong_layout() {
        let directory = tempfile::tempdir().unwrap();
        let (mut assay, _, _, _) = short_read_fixture(directory.path());
        let garbage: Vec<u8> = vec![b'A'; 36];
        let path = directory.path().join("garbage.fastq");
        write_fastq_records(
            &path,
            &[
                ("read1", &garbage),
                ("read2", &garbage),
                ("read3", &garbage),
            ],
        );
        assay.sequence_spec.get_mut("R1").unwrap().files =
            Some(vec![File::from_fastq(&path, false).unwrap()]);

        let analyzer = BarcodeAnalyzer::new(&assay, Modality::RNA);
        let warnings = smoke_test_assay(
            &assay,
            &Modality::RNA,
            &BarcodeProcessor::ShortRead(analyzer),
        )
        .unwrap();
        assert!(
            warnings
                .iter()
                .any(|warning| warning.contains("structural ratio")),
            "{warnings:?}"
        );
    }

    /// Section 4: the plan-build smoke test must warn when long-read data
    /// carries none of the fixed anchors the extractor aligns against.
    #[test]
    fn test_smoke_test_long_read_warns_on_wrong_data() {
        let directory = tempfile::tempdir().unwrap();
        let garbage: Vec<u8> = vec![b'A'; 600];
        let path = directory.path().join("garbage.fastq");
        write_fastq(&path, "read1", &garbage);
        let assay = long_read_assay(&[path.as_path()], directory.path());

        let extractor = BarcodeExtractor::new(
            &assay.library_spec,
            &Modality::ATAC,
            assay.get_whitelists(Modality::ATAC),
        )
        .unwrap();
        let warnings = smoke_test_assay(
            &assay,
            &Modality::ATAC,
            &BarcodeProcessor::LongRead(extractor),
        )
        .unwrap();
        assert!(
            warnings
                .iter()
                .any(|warning| warning.contains("extraction ratio")),
            "{warnings:?}"
        );
    }

    /// Section 4: the plan-build smoke test must warn when the pattern never
    /// matches (parse ratio) and when a capture group is mapped to the wrong
    /// region (whitelist ratio collapses while parsing stays intact).
    #[test]
    fn test_smoke_test_read_name_warns_on_mismatch() {
        let directory = tempfile::tempdir().unwrap();
        let target: Vec<u8> = vec![b'C'; 600];
        // bc1's field holds CCCCCCCC, which is not in bc1's whitelist.
        let records = [("read1:CCCCCCCC_GGGGGGGG_ACGTACGTAC", &target[..])];
        let assay = read_name_assay(directory.path(), &records);

        // Mapping error: parsing stays intact while bc1's whitelist ratio
        // collapses to zero.
        let (config, _) = ReadNameBarcodeConfig::new(
            ":(?P<bc1>[ACGTN]{8})_(?P<bc2>[ACGTN]{8})_(?P<umi>[ACGTN]{10})$",
            Some("R1"),
            std::slice::from_ref(&assay),
            &Modality::RNA,
        )
        .unwrap();
        let parser = ReadNameParser::new(&config, &assay, &Modality::RNA).unwrap();
        let warnings =
            smoke_test_assay(&assay, &Modality::RNA, &BarcodeProcessor::ReadName(parser)).unwrap();
        assert!(
            !warnings
                .iter()
                .any(|warning| warning.contains("parse ratio")),
            "{warnings:?}"
        );
        assert!(
            warnings
                .iter()
                .any(|warning| warning.contains("bc1") && warning.contains("whitelist hit ratio")),
            "{warnings:?}"
        );

        // A pattern that matches nothing produces the parse-ratio warning.
        let (config, _) = ReadNameBarcodeConfig::new(
            ":(?P<bc1>[ACGTN]{8})_(?P<bc2>[ACGTN]{8})_(?P<umi>[ACGTN]{10})_ZZ$",
            Some("R1"),
            std::slice::from_ref(&assay),
            &Modality::RNA,
        )
        .unwrap();
        let parser = ReadNameParser::new(&config, &assay, &Modality::RNA).unwrap();
        let warnings =
            smoke_test_assay(&assay, &Modality::RNA, &BarcodeProcessor::ReadName(parser)).unwrap();
        assert!(
            warnings
                .iter()
                .any(|warning| warning.contains("parse ratio")),
            "{warnings:?}"
        );
    }

    /// Paired-end read-name assays: with `read_id`, exactly one read parses
    /// and the other contributes read 2. Omitting `read_id` makes every read
    /// parse, doubling the barcode and UMI while the targets stay single.
    #[test]
    fn test_read_name_paired_end() {
        let directory = tempfile::tempdir().unwrap();
        let target1: Vec<u8> = vec![b'C'; 600];
        let target2: Vec<u8> = vec![b'G'; 600];
        let r1_path = directory.path().join("R1.fastq");
        let r2_path = directory.path().join("R2.fastq");
        write_fastq(&r1_path, "read1:AAAAAAAA_GGGGGGGG_ACGTACGTAC", &target1);
        write_fastq(&r2_path, "read1:AAAAAAAA_GGGGGGGG_ACGTACGTAC", &target2);
        let make_assay = |attach_r2: bool| -> Assay {
            let mut assay = Assay::from_path("data/test_readname.yaml").unwrap();
            assay.sequence_spec.get_mut("R1").unwrap().files =
                Some(vec![File::from_fastq(&r1_path, false).unwrap()]);
            if attach_r2 {
                let mut r2 = assay.sequence_spec.get("R1").unwrap().clone();
                r2.read_id = "R2".to_string();
                r2.strand = Strand::Neg;
                // A neg-strand read walks the library layout in reverse, so
                // its primer is the 3'-side sequencing primer.
                r2.primer_id = "p7".to_string();
                r2.files = Some(vec![File::from_fastq(&r2_path, false).unwrap()]);
                assay.sequence_spec.insert("R2".to_string(), r2);
            }
            assay
        };

        // Default: only R1 parses its name.
        let assay = make_assay(true);
        let (config, warnings) = ReadNameBarcodeConfig::new(
            ":(?P<bc1>[ACGTN]{8})_(?P<bc2>[ACGTN]{8})_(?P<umi>[ACGTN]{10})$",
            Some("R1"),
            std::slice::from_ref(&assay),
            &Modality::RNA,
        )
        .unwrap();
        assert!(warnings.is_empty());
        let mut execution = FastqPlan::new(vec![assay], Modality::RNA)
            .with_read_name_barcode(config)
            .build(false, 5000)
            .unwrap();
        assert!(execution.is_paired_end());
        let batch = execution.next_batch().unwrap().unwrap();
        assert_eq!(batch.len(), 1);
        let annotated = &batch[0];
        assert_eq!(
            annotated.barcode.as_ref().unwrap().raw.sequence(),
            b"AAAAAAAAGGGGGGGG"
        );
        assert_eq!(annotated.umi.as_ref().unwrap().sequence(), b"ACGTACGTAC");
        assert_eq!(annotated.read1.as_ref().unwrap().sequence(), &target1);
        assert_eq!(annotated.read2.as_ref().unwrap().sequence(), &target2);
        assert!(execution.next_batch().unwrap().is_none());
        execution.finish().unwrap();

        // Omitted read_id: both reads parse and the barcode is concatenated
        // once per read; a warning is emitted at configuration time.
        let assay = make_assay(true);
        let (config, warnings) = ReadNameBarcodeConfig::new(
            ":(?P<bc1>[ACGTN]{8})_(?P<bc2>[ACGTN]{8})_(?P<umi>[ACGTN]{10})$",
            None,
            std::slice::from_ref(&assay),
            &Modality::RNA,
        )
        .unwrap();
        assert_eq!(warnings.len(), 1);
        let mut execution = FastqPlan::new(vec![assay], Modality::RNA)
            .with_read_name_barcode(config)
            .build(false, 5000)
            .unwrap();
        let batch = execution.next_batch().unwrap().unwrap();
        assert_eq!(batch.len(), 1);
        let annotated = &batch[0];
        assert_eq!(
            annotated.barcode.as_ref().unwrap().raw.sequence(),
            b"AAAAAAAAGGGGGGGGAAAAAAAAGGGGGGGG"
        );
        assert_eq!(
            annotated.umi.as_ref().unwrap().sequence(),
            b"ACGTACGTACACGTACGTAC"
        );
        assert_eq!(annotated.read1.as_ref().unwrap().sequence(), &target1);
        assert_eq!(annotated.read2.as_ref().unwrap().sequence(), &target2);
        assert!(execution.next_batch().unwrap().is_none());
        execution.finish().unwrap();
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
