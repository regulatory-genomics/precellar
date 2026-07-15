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
    valid_barcodes: Vec<Vec<u8>>,
    flank_5: String,
    flank_3: String,
    expected_len: usize,
    use_read1: bool,
    kmer_size: usize,
    mismatch_tolerance: usize,
    max_barcode_mismatch: usize,
    confidence_threshold: f64,
    max_expected_errors: f64,
}

#[pymethods]
impl FloatingBarcodeExtracter {
    #[new]
    #[pyo3(
        signature = (
            output, valid_barcodes, flank_5, flank_3, expected_len, *,
            use_read1=true, kmer_size=15, mismatch_tolerance=1,
            max_barcode_mismatch=1, confidence_threshold=0.9,
            max_expected_errors=None,
        )
    )]
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        output: PathBuf,
        valid_barcodes: Bound<'_, PyAny>,
        flank_5: String,
        flank_3: String,
        expected_len: usize,
        use_read1: bool,
        kmer_size: usize,
        mismatch_tolerance: usize,
        max_barcode_mismatch: usize,
        confidence_threshold: f64,
        max_expected_errors: Option<f64>,
    ) -> Result<Self> {
        if kmer_size == 0 {
            bail!("kmer_size must be greater than zero");
        }
        if !(0.0..=1.0).contains(&confidence_threshold) {
            bail!("confidence_threshold must be between zero and one");
        }
        if let Some(value) = max_expected_errors {
            if value.is_sign_negative() {
                bail!("max_expected_errors must not be negative");
            }
        }

        let valid_barcodes = read_valid_barcodes(valid_barcodes)?;
        if let Some(length) = valid_barcodes.first().map(Vec::len) {
            if valid_barcodes.iter().any(|barcode| barcode.len() != length) {
                bail!("all valid barcodes must have the same length");
            }
        }

        Ok(Self {
            output,
            valid_barcodes,
            flank_5,
            flank_3,
            expected_len,
            use_read1,
            kmer_size,
            mismatch_tolerance,
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
            &self.flank_5,
            &self.flank_3,
            self.kmer_size,
            self.expected_len,
            self.mismatch_tolerance,
        );
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
