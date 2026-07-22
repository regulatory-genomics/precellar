use crate::align::AlignProgressBar;
use crate::aligners::AlignerRef;
use anyhow::{bail, Result};
use noodles_sam::{self as sam, alignment::io::Write};
use precellar::{
    align::MultiMapR,
    fragment::{IntoFragOpts, IntoFragments},
    qc::QcFragment,
    transcriptome::Quantifier,
};
use pyo3::prelude::*;
use pyo3::types::PyDict;
use seqspec::{
    utils::{create_file, Compression},
    ChemistryStrandedness, Modality,
};
use std::{path::PathBuf, str::FromStr};

pub(crate) struct AlignmentContext {
    pub(crate) num_threads: u16,
    pub(crate) mito_dna: Vec<String>,
    pub(crate) temp_dir: Option<PathBuf>,
    pub(crate) transcript_annotator: Option<precellar::transcriptome::TxAligner>,
}

pub(crate) struct OutputReport {
    pub(crate) name: String,
    pub(crate) metrics: serde_json::Value,
}

pub(crate) trait AlignmentSink {
    fn needs_transcriptome(&self) -> bool {
        false
    }

    fn strandedness(
        &self,
        _assays: &[seqspec::Assay],
        _modality: Modality,
    ) -> Result<Option<ChemistryStrandedness>> {
        Ok(None)
    }

    fn consume<'py, 'stream>(
        &self,
        py: Python<'py>,
        header: &sam::Header,
        alignments: &mut AlignProgressBar<'stream, AlignerRef<'py>>,
        context: &mut AlignmentContext,
    ) -> Result<Option<OutputReport>>;
}

fn parse_strandedness(
    assays: &[seqspec::Assay],
    modality: Modality,
    value: Option<&str>,
) -> Result<Option<ChemistryStrandedness>> {
    if let Some(value) = value {
        return match value.to_lowercase().as_str() {
            "forward" => Ok(Some(ChemistryStrandedness::Forward)),
            "reverse" => Ok(Some(ChemistryStrandedness::Reverse)),
            "unstranded" => Ok(Some(ChemistryStrandedness::Unstranded)),
            "auto" => Ok(None),
            _ => bail!("strandedness must be 'unstranded', 'forward', 'reverse' or 'auto'"),
        };
    }
    if let Some(strandedness) = assays[0].chemistry_strandedness {
        Ok(Some(strandedness))
    } else if modality == Modality::RNA {
        bail!("strandedness must be provided for RNA gene quantification")
    } else {
        Ok(None)
    }
}

fn fragment_file_header(compute_snv: bool, shift_left: i64, shift_right: i64) -> String {
    let header = if compute_snv {
        "# chromosome\tstart\tend\tbarcode\tcount\tstrand\talignment1_start\talignment1_snv\talignment2_start\talignment2_snv"
    } else {
        "# chromosome\tstart\tend\tbarcode\tcount\tstrand"
    };
    [
        &format!(
            "# This file contains unique fragments generated using precellar-v{}",
            env!("CARGO_PKG_VERSION")
        ),
        "#",
        "# Parameters",
        &format!("# shift_left = {}", shift_left),
        &format!("# shift_right = {}", shift_right),
        "#",
        "# Each line represents a unique fragment with the following fields:",
        header,
    ]
    .join("\n")
}

fn write_alignments<'a>(
    py: Python<'a>,
    output: PathBuf,
    header: &'a sam::Header,
    alignments: impl Iterator<Item = Vec<(Option<MultiMapR>, Option<MultiMapR>)>> + 'a,
) -> Result<()> {
    let mut writer = noodles_bam::io::writer::Builder::default().build_from_path(output)?;
    writer.write_header(header)?;
    for data in alignments {
        py.check_signals()?;
        for (read1, read2) in data {
            if let Some(read1) = read1 {
                for alignment in read1.iter() {
                    writer.write_alignment_record(header, alignment)?;
                }
            }
            if let Some(read2) = read2 {
                for alignment in read2.iter() {
                    writer.write_alignment_record(header, alignment)?;
                }
            }
        }
    }
    Ok(())
}

#[pyclass(module = "precellar.sinks", name = "NullSink")]
pub(crate) struct NullSink;

#[pymethods]
impl NullSink {
    /// Create a sink that discards all alignment output.
    #[new]
    fn new() -> Self {
        Self
    }
}

impl AlignmentSink for NullSink {
    fn consume<'py, 'stream>(
        &self,
        py: Python<'py>,
        _header: &sam::Header,
        alignments: &mut AlignProgressBar<'stream, AlignerRef<'py>>,
        _context: &mut AlignmentContext,
    ) -> Result<Option<OutputReport>> {
        for _ in alignments {
            py.check_signals()?;
        }
        Ok(None)
    }
}

#[pyclass(module = "precellar.sinks", name = "BamSink")]
pub(crate) struct BamSink {
    output: PathBuf,
}

#[pymethods]
impl BamSink {
    /// Create a sink that writes SAM alignments to a BAM file.
    ///
    /// Parameters
    /// ----------
    /// output : pathlib.Path | str
    ///     Destination BAM path.
    #[new]
    fn new(output: PathBuf) -> Self {
        Self { output }
    }
}

impl AlignmentSink for BamSink {
    fn consume<'py, 'stream>(
        &self,
        py: Python<'py>,
        header: &sam::Header,
        alignments: &mut AlignProgressBar<'stream, AlignerRef<'py>>,
        _context: &mut AlignmentContext,
    ) -> Result<Option<OutputReport>> {
        write_alignments(py, self.output.clone(), header, alignments)?;
        Ok(None)
    }
}

#[pyclass(module = "precellar.sinks", name = "FragmentsSink")]
pub(crate) struct FragmentsSink {
    output: PathBuf,
    shift_left: i64,
    shift_right: i64,
    compute_snv: bool,
    compression: Option<String>,
    compression_level: Option<u32>,
}

#[pymethods]
impl FragmentsSink {
    /// Create a sink that writes deduplicated fragments.
    ///
    /// Parameters
    /// ----------
    /// output : pathlib.Path | str
    ///     Destination fragment path.
    /// shift_left : int, keyword-only, default=4
    ///     Shift applied to the left end of paired-end fragments.
    /// shift_right : int, keyword-only, default=-5
    ///     Shift applied to the right end of paired-end fragments.
    /// compute_snv : bool, keyword-only, default=False
    ///     Include alignment SNV fields in the output.
    /// compression : str, keyword-only, optional
    ///     Explicit compression algorithm. Infer it from the output suffix if
    ///     omitted.
    /// compression_level : int, keyword-only, optional
    ///     Compression level passed to the writer.
    #[new]
    #[pyo3(signature = (
        output, *, shift_left=4, shift_right=-5, compute_snv=false,
        compression=None, compression_level=None
    ))]
    fn new(
        output: PathBuf,
        shift_left: i64,
        shift_right: i64,
        compute_snv: bool,
        compression: Option<String>,
        compression_level: Option<u32>,
    ) -> Self {
        Self {
            output,
            shift_left,
            shift_right,
            compute_snv,
            compression,
            compression_level,
        }
    }
}

impl AlignmentSink for FragmentsSink {
    fn consume<'py, 'stream>(
        &self,
        py: Python<'py>,
        header: &sam::Header,
        alignments: &mut AlignProgressBar<'stream, AlignerRef<'py>>,
        context: &mut AlignmentContext,
    ) -> Result<Option<OutputReport>> {
        let compression = self
            .compression
            .as_deref()
            .map(|value| Compression::from_str(value).map_err(|error| anyhow::anyhow!(error)))
            .transpose()?;
        let compression = compression.or_else(|| (&self.output).try_into().ok());
        let mut writer = create_file(
            self.output.clone(),
            compression,
            self.compression_level,
            context.num_threads as u32,
        )?;
        writeln!(
            writer,
            "{}",
            fragment_file_header(self.compute_snv, self.shift_left, self.shift_right)
        )?;
        let opts = IntoFragOpts {
            shift_left: self.shift_left,
            shift_right: self.shift_right,
            temp_dir: context.temp_dir.clone(),
            compute_snv: self.compute_snv,
            ..Default::default()
        };
        let mut fragment_qc = QcFragment::default();
        for mito in &context.mito_dna {
            fragment_qc.add_mito_dna(mito);
        }
        for fragments in &(&mut *alignments).into_fragments(header, opts) {
            py.check_signals()?;
            for fragment in fragments {
                fragment_qc.update(&fragment);
                writeln!(writer, "{}", fragment)?;
            }
        }
        Ok(Some(OutputReport {
            name: "fragment".to_owned(),
            metrics: fragment_qc.into(),
        }))
    }
}

#[pyclass(module = "precellar.sinks", name = "GeneQuantificationSink")]
pub(crate) struct GeneQuantificationSink {
    output: PathBuf,
    strandedness: Option<String>,
}

#[pymethods]
impl GeneQuantificationSink {
    /// Create a sink that writes gene counts to an H5AD file.
    ///
    /// Parameters
    /// ----------
    /// output : pathlib.Path | str
    ///     Destination H5AD path.
    /// strandedness : {"forward", "reverse", "unstranded", "auto"}, optional
    ///     Strand model. If omitted, use the assay declaration.
    #[new]
    #[pyo3(signature = (output, *, strandedness=None))]
    fn new(output: PathBuf, strandedness: Option<String>) -> Self {
        Self {
            output,
            strandedness,
        }
    }
}

impl AlignmentSink for GeneQuantificationSink {
    fn needs_transcriptome(&self) -> bool {
        true
    }

    fn strandedness(
        &self,
        assays: &[seqspec::Assay],
        modality: Modality,
    ) -> Result<Option<ChemistryStrandedness>> {
        parse_strandedness(assays, modality, self.strandedness.as_deref())
    }

    fn consume<'py, 'stream>(
        &self,
        _py: Python<'py>,
        header: &sam::Header,
        alignments: &mut AlignProgressBar<'stream, AlignerRef<'py>>,
        context: &mut AlignmentContext,
    ) -> Result<Option<OutputReport>> {
        let annotator = context
            .transcript_annotator
            .take()
            .ok_or_else(|| anyhow::anyhow!("gene quantification requires a STAR aligner"))?;
        let mut quantifier = Quantifier::new(annotator)?;
        quantifier.num_threads = context.num_threads as usize;
        quantifier.temp_dir = context.temp_dir.clone();
        let quant_qc = quantifier.quantify(header, alignments, self.output.clone())?;
        Ok(Some(OutputReport {
            name: "gene_quantification".to_owned(),
            metrics: quant_qc.into(),
        }))
    }
}

pub(crate) fn register_sinks(parent: &Bound<'_, PyModule>) -> PyResult<()> {
    let sinks = PyModule::new(parent.py(), "sinks")?;
    sinks.add_class::<NullSink>()?;
    sinks.add_class::<BamSink>()?;
    sinks.add_class::<FragmentsSink>()?;
    sinks.add_class::<GeneQuantificationSink>()?;
    parent.add_submodule(&sinks)?;

    let sys = PyModule::import(parent.py(), "sys")?;
    let modules = sys.getattr("modules")?.cast_into::<PyDict>()?;
    modules.set_item("precellar.sinks", &sinks)?;
    Ok(())
}
