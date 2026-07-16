use anyhow::{bail, Result};
use precellar::{
    align::FastqPlan, barcode::BarcodeCorrectOptions, middleware::FloatingBarcodeStage,
    utils::insertion_extractor::InsertionExtractor,
};
use pyo3::prelude::*;
use pyo3::types::PyDict;
use seqspec::utils::{create_file, open_file, Compression};
use std::{
    io::{BufRead, BufReader},
    path::PathBuf,
};

/// Rust-backed middleware for extracting floating barcodes before alignment.
#[pyclass(module = "precellar.middleware")]
pub struct FloatingBarcodeExtracter {
    output: PathBuf,
    valid_barcodes: Vec<Vec<Vec<u8>>>,
    flanks: Vec<String>,
    expected_lens: Vec<usize>,
    use_read1: bool,
    kmer_size: usize,
    max_fixed_edits: usize,
    max_barcode_mismatch: usize,
    confidence_threshold: f64,
    max_expected_errors: f64,
}

#[pymethods]
impl FloatingBarcodeExtracter {
    #[new]
    #[pyo3(
        signature = (
            output, valid_barcodes, flanks, expected_lens, *,
            use_read1=true, kmer_size=15, max_fixed_edits=1,
            max_barcode_mismatch=1, confidence_threshold=0.9,
            max_expected_errors=None,
        )
    )]
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        output: PathBuf,
        valid_barcodes: Bound<'_, PyAny>,
        flanks: Vec<String>,
        expected_lens: Vec<usize>,
        use_read1: bool,
        kmer_size: usize,
        max_fixed_edits: usize,
        max_barcode_mismatch: usize,
        confidence_threshold: f64,
        max_expected_errors: Option<f64>,
    ) -> Result<Self> {
        if kmer_size == 0 {
            bail!("kmer_size must be greater than zero");
        }
        if flanks.len() < 2 {
            bail!("flanks must contain at least two fixed sequences");
        }
        if expected_lens.len() != flanks.len() - 1 {
            bail!(
                "expected_lens must contain one length per insertion gap (expected {}, got {})",
                flanks.len() - 1,
                expected_lens.len()
            );
        }
        if !(0.0..=1.0).contains(&confidence_threshold) {
            bail!("confidence_threshold must be between zero and one");
        }
        if let Some(value) = max_expected_errors {
            if value.is_sign_negative() {
                bail!("max_expected_errors must not be negative");
            }
        }
        if max_barcode_mismatch > 2 {
            bail!("max_barcode_mismatch must not exceed 2");
        }

        let valid_barcodes = valid_barcodes
            .extract::<Vec<Bound<'_, PyAny>>>()?
            .into_iter()
            .map(read_valid_barcodes)
            .collect::<Result<Vec<_>>>()?;
        if valid_barcodes.len() != expected_lens.len() {
            bail!(
                "valid_barcodes must contain one whitelist per insertion gap (expected {}, got {})",
                expected_lens.len(),
                valid_barcodes.len()
            );
        }
        for (index, barcodes) in valid_barcodes.iter().enumerate() {
            if let Some(length) = barcodes.first().map(Vec::len) {
                if barcodes.iter().any(|barcode| barcode.len() != length) {
                    bail!("all valid barcodes within a whitelist must have the same length");
                }
                if length != expected_lens[index] {
                    bail!(
                        "whitelist {index} contains {length}-base barcodes but expected_lens[{index}] is {}",
                        expected_lens[index]
                    );
                }
            }
        }

        Ok(Self {
            output,
            valid_barcodes,
            flanks,
            expected_lens,
            use_read1,
            kmer_size,
            max_fixed_edits,
            max_barcode_mismatch,
            confidence_threshold,
            max_expected_errors: max_expected_errors.unwrap_or(f64::MAX),
        })
    }
}

fn read_valid_barcodes(value: Bound<'_, PyAny>) -> Result<Vec<Vec<u8>>> {
    let barcodes = if let Ok(barcodes) = value.extract::<Vec<String>>() {
        barcodes
    } else {
        let path = value.extract::<PathBuf>().map_err(|_| {
            anyhow::anyhow!("valid_barcodes must be a sequence of strings or a file path")
        })?;
        let reader = BufReader::new(open_file(path)?);
        reader
            .lines()
            .map(|line| line.map(|line| line.trim().to_owned()))
            .filter_map(|line| match line {
                Ok(line) if !line.is_empty() => Some(Ok(line)),
                Ok(_) => None,
                Err(error) => Some(Err(error.into())),
            })
            .collect::<Result<Vec<_>>>()?
    };

    if barcodes.is_empty() {
        bail!("valid_barcodes must contain at least one barcode");
    }

    barcodes
        .into_iter()
        .map(|barcode| {
            let barcode = barcode.to_ascii_uppercase();
            if barcode.is_empty()
                || !barcode
                    .bytes()
                    .all(|base| matches!(base, b'A' | b'C' | b'G' | b'T'))
            {
                bail!("valid barcodes must be non-empty DNA sequences");
            }
            Ok(barcode.into_bytes())
        })
        .collect::<Result<_>>()
}

impl FloatingBarcodeExtracter {
    pub(crate) fn configure_plan(&self, plan: FastqPlan) -> Result<FastqPlan> {
        let compression = Compression::try_from(&self.output).ok();
        let output = create_file(&self.output, compression, None, 1)?;
        let extractor = InsertionExtractor::new(
            self.flanks.clone(),
            self.expected_lens.clone(),
            self.kmer_size,
            self.max_fixed_edits,
        )?;
        let valid_barcodes = self.valid_barcodes.clone();
        let correction_options = BarcodeCorrectOptions {
            max_expected_errors: self.max_expected_errors,
            bc_confidence_threshold: self.confidence_threshold,
            max_mismatch: self.max_barcode_mismatch,
        };
        let use_read1 = self.use_read1;

        let stage = FloatingBarcodeStage::new(
            use_read1,
            extractor,
            valid_barcodes,
            correction_options,
            output,
        )?;
        Ok(plan.with_stage(stage))
    }
}

pub(crate) fn register_middleware(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let middleware = PyModule::new(parent_module.py(), "middleware")?;
    middleware.add_class::<FloatingBarcodeExtracter>()?;
    parent_module.add_submodule(&middleware)?;

    let sys = PyModule::import(parent_module.py(), "sys")?;
    let modules = sys.getattr("modules")?.cast_into::<PyDict>()?;
    modules.set_item("precellar.middleware", middleware)
}
