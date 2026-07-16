use anyhow::{bail, Result};
use precellar::{
    align::FastqPlan,
    barcode::BarcodeCorrectOptions,
    middleware::tfseq::{FloatingBarcodeEntry, FloatingBarcodeStage, FloatingBarcodeTable},
    utils::insertion_extractor::InsertionExtractor,
};
use pyo3::prelude::*;
use seqspec::utils::{create_file, open_file, Compression};
use std::{
    io::{BufRead, BufReader},
    path::PathBuf,
};

/// Rust-backed middleware for extracting floating barcodes before alignment.
#[pyclass(module = "precellar.middleware.tfseq")]
pub struct FloatingBarcodeFinder {
    output: PathBuf,
    barcode_table: FloatingBarcodeTable,
    flanks: Vec<String>,
    expected_lens: Vec<usize>,
    use_read1: bool,
    kmer_size: usize,
    max_fixed_edit_rate: f64,
    max_barcode_mismatch: usize,
    confidence_threshold: f64,
    max_expected_errors: f64,
}

#[pymethods]
impl FloatingBarcodeFinder {
    #[new]
    #[pyo3(
        signature = (
            output, barcode_table, flanks, expected_lens, *,
            use_read1=true, kmer_size=15, max_fixed_edit_rate=0.1,
            max_barcode_mismatch=1, confidence_threshold=0.9,
            max_expected_errors=None,
        )
    )]
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        output: PathBuf,
        barcode_table: PathBuf,
        flanks: Vec<String>,
        expected_lens: Vec<usize>,
        use_read1: bool,
        kmer_size: usize,
        max_fixed_edit_rate: f64,
        max_barcode_mismatch: usize,
        confidence_threshold: f64,
        max_expected_errors: Option<f64>,
    ) -> Result<Self> {
        if kmer_size == 0 {
            bail!("kmer_size must be greater than zero");
        }
        if !max_fixed_edit_rate.is_finite() || !(0.0..=1.0).contains(&max_fixed_edit_rate) {
            bail!("max_fixed_edit_rate must be finite and between zero and one");
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

        if expected_lens.len() != 2 {
            bail!("barcode_table requires exactly two insertion gaps");
        }
        let barcode_table =
            FloatingBarcodeTable::new(read_barcode_table(barcode_table)?, &expected_lens)?;

        Ok(Self {
            output,
            barcode_table,
            flanks,
            expected_lens,
            use_read1,
            kmer_size,
            max_fixed_edit_rate,
            max_barcode_mismatch,
            confidence_threshold,
            max_expected_errors: max_expected_errors.unwrap_or(f64::MAX),
        })
    }
}

fn read_barcode_table(path: PathBuf) -> Result<Vec<FloatingBarcodeEntry>> {
    let reader = BufReader::new(open_file(path)?);
    let mut entries = Vec::new();
    for (line_index, line) in reader.lines().enumerate() {
        let line = line?;
        if line.trim().is_empty() {
            continue;
        }
        let fields: Vec<_> = line.split('\t').map(str::trim).collect();
        if fields.len() != 4 {
            bail!(
                "barcode table line {} must contain four tab-separated columns",
                line_index + 1
            );
        }
        entries.push(FloatingBarcodeEntry::new(
            fields[0].as_bytes(),
            fields[1].as_bytes(),
            fields[2].as_bytes(),
            fields[3].as_bytes(),
        ));
    }
    if entries.is_empty() {
        bail!("barcode_table must contain at least one row");
    }
    Ok(entries)
}

impl FloatingBarcodeFinder {
    pub(crate) fn configure_plan(&self, plan: FastqPlan, num_threads: usize) -> Result<FastqPlan> {
        let compression = Compression::try_from(&self.output).ok();
        let output = create_file(&self.output, compression, None, 1)?;
        let extractor = InsertionExtractor::new(
            self.flanks.clone(),
            self.expected_lens.clone(),
            self.kmer_size,
            self.max_fixed_edit_rate,
        )?;
        let barcode_table = self.barcode_table.clone();
        let correction_options = BarcodeCorrectOptions {
            max_expected_errors: self.max_expected_errors,
            bc_confidence_threshold: self.confidence_threshold,
            max_mismatch: self.max_barcode_mismatch,
        };
        let use_read1 = self.use_read1;

        let stage = FloatingBarcodeStage::new(
            use_read1,
            extractor,
            barcode_table,
            correction_options,
            output,
        )?;
        Ok(plan.with_stage(stage.with_num_threads(num_threads)?))
    }
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<FloatingBarcodeFinder>()
}
