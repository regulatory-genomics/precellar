use crate::aligners::AlignerRef;
use crate::pyseqspec::extract_assays;
use crate::sinks::{
    AlignmentContext, AlignmentSink, BamSink, FragmentsSink, GeneQuantificationSink,
};
use anyhow::{bail, Result};
use indicatif::{ProgressBar, ProgressFinish, ProgressStyle};
use log::info;
use precellar::align::{Aligner, AlignmentResult, AlignmentRunner};
use precellar::align::{BarcodeCorrectionConfig, FastqPlan, MultiMapR};
use precellar::qc::Metric;
use pyo3::prelude::*;
use pyo3::types::PyDict;
use pyo3::BoundObject;
use seqspec::Modality;
use serde_json::{Map, Value};
use std::{path::PathBuf, str::FromStr};

/// Create a genome index from a fasta file.
///
/// Parameters
/// ----------
///
/// fasta: Path
///    File path to the fasta file.
/// genome_prefix: Path
///   File path to the genome index.
#[pyfunction]
pub fn make_bwa_mem2_index(fasta: PathBuf, genome_prefix: PathBuf) -> Result<()> {
    bwa_mem2::FMIndex::new(fasta, genome_prefix)?;
    Ok(())
}

/// Create a minibwa index from a FASTA file.
///
/// Parameters
/// ----------
///
/// fasta: Path
///    File path to the FASTA file.
/// index_prefix: Path
///    File path prefix for the minibwa index files.
/// num_threads: int
///    The number of threads to use when building the index.
/// methylation: bool
///    Whether to build a methylation-aware index.
#[pyfunction]
#[pyo3(
    signature = (fasta, index_prefix, *, num_threads=8, methylation=false),
    text_signature = "(fasta, index_prefix, *, num_threads=8, methylation=False)",
)]
pub fn make_minibwa_index(
    fasta: PathBuf,
    index_prefix: PathBuf,
    num_threads: i32,
    methylation: bool,
) -> Result<()> {
    minibwa::build_index(fasta, index_prefix, num_threads, methylation)?;
    Ok(())
}

/// Create a minimap2 index from a FASTA file.
///
/// `preset` selects the minimap2 read or alignment model. See the Python API
/// documentation for the list of supported presets.
#[pyfunction]
#[pyo3(
    signature = (fasta, output_index, *, preset="map-ont"),
    text_signature = "(fasta, output_index, *, preset='map-ont')",
)]
pub fn make_minimap2_index(fasta: PathBuf, output_index: PathBuf, preset: &str) -> Result<()> {
    let preset = match preset.to_lowercase().as_str() {
        "map-ont" => minimap2::Preset::MapOnt,
        "map-pb" => minimap2::Preset::MapPb,
        "map-hifi" => minimap2::Preset::MapHifi,
        "lr:hq" => minimap2::Preset::LrHq,
        "splice" => minimap2::Preset::Splice,
        "splice:hq" => minimap2::Preset::SpliceHq,
        "splice:sr" => minimap2::Preset::SpliceSr,
        "asm5" => minimap2::Preset::Asm5,
        "asm10" => minimap2::Preset::Asm10,
        "asm20" => minimap2::Preset::Asm20,
        "short" => minimap2::Preset::Short,
        "sr" => minimap2::Preset::Sr,
        "ava-pb" => minimap2::Preset::AvaPb,
        "ava-ont" => minimap2::Preset::AvaOnt,
        _ => bail!(
            "Invalid preset '{}'. Valid presets: map-ont, map-pb, map-hifi, lr:hq, splice, splice:hq, splice:sr, asm5, asm10, asm20, short, sr, ava-pb, ava-ont",
            preset,
        ),
    };

    info!(
        "Creating minimap2 index for fasta: {:?} with preset: {:?}",
        fasta, preset
    );
    minimap2::Aligner::builder()
        .preset(preset)
        .with_index(
            fasta.to_str().unwrap(),
            Some(output_index.to_str().unwrap()),
        )
        .map_err(|error| anyhow::anyhow!("Failed to create minimap2 index: {}", error))?;
    Ok(())
}

/// A reusable configuration for reading and annotating assay FASTQs.
///
/// Use this builder to configure the input assay, modality, barcode correction,
/// and optional middleware before attaching an aligner with
/// [`FastqPipeline::align_with`]. The constructor does not open FASTQ files.
/// Processing starts only when a terminal method on the returned
/// [`AlignmentJob`] is called.
///
/// Anti-Patterns
/// -------------
/// - Do not pass an aligner to the constructor. Attach it with `align_with`.
/// - Do not put fragment-only options such as `shift_left` here. Pass them to
///   `precellar.sinks.FragmentsSink`.
/// - Do not expect this object to execute a workflow by itself. Select one
///   terminal output after calling `align_with`.
#[pyclass]
pub struct FastqPipeline {
    assays: Vec<seqspec::Assay>,
    modality: Modality,
    barcode_config: BarcodeCorrectionConfig,
    middleware: Option<Py<PyAny>>,
}

#[pymethods]
impl FastqPipeline {
    /// Create a lazy FASTQ pipeline from an assay object, assay path, or list.
    ///
    /// Parameters
    /// ----------
    /// assay : precellar.Assay | pathlib.Path | str | list
    ///     A seqspec assay object, a path to a seqspec file, or a list of
    ///     assays. Assays are processed in the supplied order.
    /// modality : str, optional
    ///     Modality to annotate, such as `"rna"` or `"atac"`. If omitted,
    ///     use the modality declared by the first assay.
    ///
    /// Returns
    /// -------
    /// FastqPipeline
    ///     A lazy configuration object that can be further configured.
    ///
    /// Raises
    /// ------
    /// ValueError
    ///     If the assay list is empty or the modality is invalid.
    ///
    /// Examples
    /// --------
    /// The following complete script creates a pipeline and selects the
    /// modality explicitly before attaching an aligner:
    ///
    ///     import precellar
    ///
    ///     assay = precellar.Assay("assay.yaml")
    ///     assay.update_read("R1", fastq="R1.fastq.gz")
    ///     assay.update_read("R2", fastq="R2.fastq.gz")
    ///     pipeline = precellar.FastqPipeline(assay, modality="atac")
    ///     job = pipeline.align_with(
    ///         precellar.aligners.BwaMem2("/references/GRCh38"),
    ///         num_threads=4,
    ///     )
    ///     report = job.run_sink(precellar.sinks.BamSink("aligned.bam"))
    #[new]
    #[pyo3(signature = (assay, *, modality=None))]
    pub fn new(assay: Bound<'_, PyAny>, modality: Option<&str>) -> Result<Self> {
        let assays = extract_assays(assay)?;
        if assays.is_empty() {
            bail!("assay must contain at least one assay");
        }
        let modality = match modality {
            Some(modality) => Modality::from_str(modality)?,
            None => {
                let modality = assays[0].modality()?;
                log::info!("Using default modality: {}", modality);
                modality
            }
        };
        Ok(Self {
            assays,
            modality,
            barcode_config: BarcodeCorrectionConfig {
                confidence_threshold: 0.9,
                max_mismatch: 1,
            },
            middleware: None,
        })
    }

    /// Configure standard barcode correction for FASTQ annotation.
    ///
    /// Parameters
    /// ----------
    /// confidence_threshold : float, keyword-only, default=0.9
    ///     Minimum confidence required to accept a barcode correction.
    /// max_mismatch : int, keyword-only, default=1
    ///     Maximum number of barcode mismatches considered during correction.
    ///
    /// Returns
    /// -------
    /// FastqPipeline
    ///     The same pipeline object, allowing method chaining.
    ///
    /// Raises
    /// ------
    /// ValueError
    ///     If `confidence_threshold` is not finite or is outside `[0, 1]`.
    ///
    /// Anti-Patterns
    /// -------------
    /// - Do not use this method to configure TFseq floating barcodes. Use
    ///   `with_middleware` for floating-barcode extraction.
    ///
    /// Examples
    /// --------
    ///     import precellar
    ///
    ///     assay = precellar.Assay("assay.yaml")
    ///     assay.update_read("R1", fastq="R1.fastq.gz")
    ///     assay.update_read("R2", fastq="R2.fastq.gz")
    ///     job = (
    ///         precellar.FastqPipeline(assay, modality="rna")
    ///         .with_barcode_correction(
    ///             confidence_threshold=0.95,
    ///             max_mismatch=1,
    ///         )
    ///         .align_with(precellar.aligners.Star("/references/STAR"))
    ///     )
    ///     report = job.run_sink(
    ///         precellar.sinks.GeneQuantificationSink(
    ///             "matrix.h5ad", strandedness="reverse"
    ///         )
    ///     )
    #[pyo3(signature = (*, confidence_threshold=0.9, max_mismatch=1))]
    pub fn with_barcode_correction<'py>(
        mut slf: PyRefMut<'py, Self>,
        confidence_threshold: f64,
        max_mismatch: usize,
    ) -> Result<PyRefMut<'py, Self>> {
        if !confidence_threshold.is_finite() || !(0.0..=1.0).contains(&confidence_threshold) {
            bail!("confidence_threshold must be finite and between zero and one");
        }
        slf.barcode_config = BarcodeCorrectionConfig {
            confidence_threshold,
            max_mismatch,
        };
        Ok(slf)
    }

    /// Add Rust-backed TFseq floating-barcode middleware.
    ///
    /// Parameters
    /// ----------
    /// middleware : precellar.middleware.tfseq.FloatingBarcodeStage
    ///     Middleware configured with floating-barcode templates and a barcode
    ///     table.
    ///
    /// Returns
    /// -------
    /// FastqPipeline
    ///     The same pipeline object, allowing method chaining.
    ///
    /// Raises
    /// ------
    /// TypeError
    ///     If `middleware` is not a `FloatingBarcodeStage`.
    ///
    /// Anti-Patterns
    /// -------------
    /// - Do not pass a generic Python callable. Middleware must be the Rust
    ///   `FloatingBarcodeStage` object.
    /// - Do not call a terminal output method on the middleware object. Attach
    ///   it here, then terminate the returned alignment job.
    ///
    /// Examples
    /// --------
    ///     import precellar
    ///
    ///     assay = precellar.Assay("assay.yaml")
    ///     assay.update_read("R1", fastq="R1.fastq.gz")
    ///     assay.update_read("R2", fastq="R2.fastq.gz")
    ///     floating = precellar.middleware.tfseq.FloatingBarcodeStage(
    ///         output="floating.tsv.gz",
    ///         barcode_table="barcodes.tsv",
    ///         flanks=["ACGT", "TGCA", "GACT"],
    ///     )
    ///     job = (
    ///         precellar.FastqPipeline(assay, modality="rna")
    ///         .with_middleware(floating)
    ///         .align_with(precellar.aligners.Star("/references/STAR"))
    ///     )
    ///     report = job.run_sink(
    ///         precellar.sinks.GeneQuantificationSink(
    ///             "matrix.h5ad", strandedness="reverse"
    ///         )
    ///     )
    pub fn with_middleware<'py>(
        mut slf: PyRefMut<'py, Self>,
        middleware: Bound<'py, PyAny>,
    ) -> Result<PyRefMut<'py, Self>> {
        middleware
            .extract::<PyRef<'_, crate::middleware::FloatingBarcodeStage>>()
            .map_err(|_| {
                anyhow::anyhow!(
                    "middleware must be a precellar.middleware.tfseq.FloatingBarcodeStage"
                )
            })?;
        slf.middleware = Some(middleware.unbind());
        Ok(slf)
    }

    /// Attach an aligner and create a lazy, one-shot alignment job.
    ///
    /// Parameters
    /// ----------
    /// aligner : precellar.aligners.Star | BwaMem2 | MiniBwa | Minimap2
    ///     Configured aligner object.
    /// num_threads : int, keyword-only, default=8
    ///     Number of threads used for FASTQ processing and alignment.
    /// chunk_size : int, keyword-only, default=10000000
    ///     Approximate number of bases per synchronized FASTQ batch per thread.
    /// mito_dna : list[str], keyword-only, default=["chrM", "M"]
    ///     Reference names treated as mitochondrial DNA for alignment QC.
    /// temp_dir : pathlib.Path | str, optional
    ///     Temporary directory used by fragment and quantification consumers.
    ///
    /// Returns
    /// -------
    /// AlignmentJob
    ///     A job that must be consumed by exactly one terminal output method.
    ///
    /// Raises
    /// ------
    /// ValueError
    ///     If `num_threads` or `chunk_size` is zero.
    /// TypeError
    ///     If the object is not a supported aligner. The check is performed
    ///     when the terminal operation begins because the job is lazy.
    ///
    /// Anti-Patterns
    /// -------------
    /// - Do not call `run_sink` twice on the same job. A job is one-shot.
    /// - Do not pass `shift_left`, `shift_right`, or `strandedness` here. Those
    ///   options belong to their terminal consumers.
    ///
    /// Examples
    /// --------
    ///     import precellar
    ///
    ///     assay = precellar.Assay("assay.yaml")
    ///     assay.update_read("R1", fastq="R1.fastq.gz")
    ///     assay.update_read("R2", fastq="R2.fastq.gz")
    ///     pipeline = precellar.FastqPipeline(assay, modality="atac")
    ///     job = pipeline.align_with(
    ///         precellar.aligners.BwaMem2("/references/GRCh38"),
    ///         num_threads=8,
    ///         chunk_size=10_000_000,
    ///     )
    ///     report = job.run_sink(
    ///         precellar.sinks.FragmentsSink("fragments.tsv.gz")
    ///     )
    #[pyo3(signature = (
        aligner, *, num_threads=8, chunk_size=10000000,
        mito_dna=vec!["chrM".to_owned(), "M".to_owned()], temp_dir=None
    ))]
    pub fn align_with(
        slf: PyRef<'_, Self>,
        py: Python<'_>,
        aligner: Bound<'_, PyAny>,
        num_threads: u16,
        chunk_size: usize,
        mito_dna: Vec<String>,
        temp_dir: Option<PathBuf>,
    ) -> Result<AlignmentJob> {
        if num_threads == 0 {
            bail!("num_threads must be greater than zero");
        }
        if chunk_size == 0 {
            bail!("chunk_size must be greater than zero");
        }
        Ok(AlignmentJob {
            state: Some(AlignmentJobState {
                assays: slf.assays.clone(),
                modality: slf.modality,
                barcode_config: slf.barcode_config.clone(),
                middleware: slf
                    .middleware
                    .as_ref()
                    .map(|middleware| middleware.clone_ref(py)),
                aligner: aligner.unbind(),
                num_threads,
                chunk_size,
                mito_dna,
                temp_dir,
            }),
        })
    }
}

/// A lazy, one-shot alignment workflow.
///
/// Call [`AlignmentJob::run_sink`] with one sink from `precellar.sinks` to
/// execute the workflow. The active core stream is created only after a sink
/// is selected.
///
/// Anti-Patterns
/// -------------
/// - Do not reuse a job after a terminal method returns. Create another job
///   from the `FastqPipeline` when a second output is required.
/// - Do not expect `run_sink` to return alignment records. It returns a QC
///   dictionary after the selected sink consumes the output.
#[pyclass]
pub struct AlignmentJob {
    state: Option<AlignmentJobState>,
}

struct AlignmentJobState {
    assays: Vec<seqspec::Assay>,
    modality: Modality,
    barcode_config: BarcodeCorrectionConfig,
    middleware: Option<Py<PyAny>>,
    aligner: Py<PyAny>,
    num_threads: u16,
    chunk_size: usize,
    mito_dna: Vec<String>,
    temp_dir: Option<PathBuf>,
}

#[pymethods]
impl AlignmentJob {
    /// Run one configured output sink and return the complete QC report.
    ///
    /// Parameters
    /// ----------
    /// sink : precellar.sinks.BamSink | precellar.sinks.FragmentsSink | precellar.sinks.GeneQuantificationSink
    ///     Configured output sink.
    ///
    /// Returns
    /// -------
    /// dict
    ///     QC sections for FASTQ annotation, alignment, and the selected sink.
    ///
    /// Raises
    /// ------
    /// RuntimeError
    ///     If the job has already been consumed.
    /// TypeError
    ///     If `sink` is not a supported sink type.
    ///
    /// Anti-Patterns
    /// -------------
    /// - Do not pass an output path directly. Configure it on the sink.
    /// - Do not call `run_sink` twice on the same job.
    ///
    /// Examples
    /// --------
    ///     import precellar
    ///
    ///     assay = precellar.Assay("assay.yaml")
    ///     assay.update_read("R1", fastq="R1.fastq.gz")
    ///     assay.update_read("R2", fastq="R2.fastq.gz")
    ///     job = precellar.FastqPipeline(assay, modality="rna").align_with(
    ///         precellar.aligners.Star("/references/STAR")
    ///     )
    ///     report = job.run_sink(
    ///         precellar.sinks.BamSink("aligned.bam")
    ///     )
    pub fn run_sink<'py>(&mut self, py: Python<'py>, sink: Bound<'py, PyAny>) -> Result<Py<PyAny>> {
        if let Ok(sink) = sink.extract::<PyRef<'_, BamSink>>() {
            return self.take_state()?.run_sink(py, &*sink);
        }
        if let Ok(sink) = sink.extract::<PyRef<'_, FragmentsSink>>() {
            return self.take_state()?.run_sink(py, &*sink);
        }
        if let Ok(sink) = sink.extract::<PyRef<'_, GeneQuantificationSink>>() {
            return self.take_state()?.run_sink(py, &*sink);
        }
        Err(pyo3::exceptions::PyTypeError::new_err(
            "sink must be precellar.sinks.BamSink, precellar.sinks.FragmentsSink, or precellar.sinks.GeneQuantificationSink",
        )
        .into())
    }
}

impl AlignmentJob {
    fn take_state(&mut self) -> Result<AlignmentJobState> {
        self.state
            .take()
            .ok_or_else(|| anyhow::anyhow!("this alignment job has already been consumed"))
    }
}

impl AlignmentJobState {
    fn run_sink<'py, S: AlignmentSink>(self, py: Python<'py>, sink: &S) -> Result<Py<PyAny>> {
        let strandedness = sink.strandedness(&self.assays, self.modality)?;
        let mut aligner = AlignerRef::try_from(self.aligner.bind(py).clone())?;
        let header = aligner.header();
        let transcript_annotator = if sink.needs_transcriptome() {
            aligner.transcript_annotator(strandedness)
        } else {
            None
        };

        let mut plan =
            FastqPlan::new(self.assays, self.modality).with_barcode_config(self.barcode_config);
        if let Some(middleware) = self.middleware {
            let middleware = middleware
                .bind(py)
                .extract::<PyRef<'_, crate::middleware::FloatingBarcodeStage>>()
                .map_err(|_| {
                    anyhow::anyhow!(
                        "middleware must be a precellar.middleware.tfseq.FloatingBarcodeStage"
                    )
                })?;
            plan = middleware.configure_plan(plan, self.num_threads as usize)?;
        }

        let execution = plan.build(true, self.num_threads as usize * self.chunk_size)?;
        log::info!(
            "Aligning {} reads to reference genome ...",
            execution.read_summary()
        );
        let alignments = AlignmentRunner::new(&mut aligner, self.num_threads)
            .with_mito_dna(self.mito_dna.iter().cloned())
            .stream(execution);
        let mut alignments = AlignProgressBar::new(alignments);
        let mut context = AlignmentContext {
            num_threads: self.num_threads,
            mito_dna: self.mito_dna,
            temp_dir: self.temp_dir,
            transcript_annotator,
        };
        let output_qc = sink.consume(py, &header, &mut alignments, &mut context)?;
        let report = alignments.finish()?;
        let report = assemble_report(
            report,
            output_qc.map(|report| (report.name, report.metrics)),
        );
        Ok(value_into_pyobject(report, py).unbind())
    }
}

pub(crate) struct AlignProgressBar<'a, A> {
    pb: ProgressBar,
    alignments: AlignmentResult<'a, A>,
}

impl<'a, A> AlignProgressBar<'a, A> {
    fn new(alignments: AlignmentResult<'a, A>) -> Self {
        let pb = ProgressBar::new(alignments.num_records() as u64);
        let style = ProgressStyle::with_template(
            "{percent}%|{wide_bar:.cyan/blue}| {human_pos:>}/{human_len:} [{elapsed}<{eta}, {per_sec}]",
        )
        .unwrap();
        pb.set_style(style);
        Self {
            pb: pb.with_finish(ProgressFinish::Abandon),
            alignments,
        }
    }

    fn finish(self) -> Result<precellar::align::RunReport> {
        self.alignments.finish()
    }
}

impl<A: Aligner> Iterator for AlignProgressBar<'_, A> {
    type Item = Vec<(Option<MultiMapR>, Option<MultiMapR>)>;

    fn next(&mut self) -> Option<Self::Item> {
        let item = self.alignments.next();
        self.pb.set_position(self.alignments.num_processed() as u64);
        item
    }
}

fn assemble_report(
    report: precellar::align::RunReport,
    output_qc: Option<(String, Value)>,
) -> Value {
    let mut qc_metrics = Map::new();
    if let Some((name, metrics)) = output_qc {
        qc_metrics.insert(name, metrics);
    }
    for middleware in report.fastq.middleware {
        qc_metrics.insert(middleware.name, middleware.metrics);
    }
    qc_metrics.insert("fastq".to_owned(), report.fastq.fastq.to_json());
    qc_metrics.insert("alignment".to_owned(), report.alignment.to_json());
    Value::Object(qc_metrics)
}

fn value_into_pyobject<'py>(val: Value, py: Python<'py>) -> Bound<'py, PyAny> {
    match val {
        Value::Null => py.None().into_bound(py),
        Value::Bool(value) => pyo3::types::PyBool::new(py, value)
            .into_pyobject(py)
            .unwrap()
            .into_bound()
            .into_any(),
        Value::Number(number) => {
            if let Some(value) = number.as_i64() {
                value.into_pyobject(py).unwrap().into_any()
            } else if let Some(value) = number.as_u64() {
                value.into_pyobject(py).unwrap().into_any()
            } else if let Some(value) = number.as_f64() {
                value.into_pyobject(py).unwrap().into_any()
            } else {
                panic!("invalid number type")
            }
        }
        Value::String(value) => value.into_pyobject(py).unwrap().into_any(),
        Value::Array(values) => {
            let list = pyo3::types::PyList::empty(py);
            for value in values {
                list.append(value_into_pyobject(value, py)).unwrap();
            }
            list.into_any()
        }
        Value::Object(values) => {
            let dict = PyDict::new(py);
            for (key, value) in values {
                dict.set_item(key, value_into_pyobject(value, py)).unwrap();
            }
            dict.into_any()
        }
    }
}
