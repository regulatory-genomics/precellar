use anyhow::{bail, Result};
use precellar::{
    align::FastqPlan,
    barcode::BarcodeCorrectOptions,
    middleware::tfseq::{
        FloatingBarcodeEntry, FloatingBarcodeStage as CoreFloatingBarcodeStage,
        FloatingBarcodeTable,
    },
    utils::insertion_extractor::InsertionExtractor,
};
use pyo3::prelude::*;
use seqspec::utils::{create_file, open_file, Compression};
use std::{
    collections::HashMap,
    io::{BufRead, BufReader},
    path::PathBuf,
};

/// Rust-backed middleware for extracting floating barcodes before alignment.
///
/// Valid insertion-length profiles are inferred from barcode-table rows.
#[pyclass(module = "precellar.middleware.tfseq", name = "FloatingBarcodeStage")]
pub struct FloatingBarcodeStage {
    output: PathBuf,
    barcode_table: FloatingBarcodeTable,
    flanks: Vec<String>,
    use_read1: bool,
    kmer_size: usize,
    max_fixed_edit_rate: f64,
    max_barcode_mismatch: usize,
    confidence_threshold: f64,
    max_expected_errors: f64,
}

#[pymethods]
impl FloatingBarcodeStage {
    #[new]
    #[pyo3(
        signature = (
            output, barcode_table, flanks, *,
            use_read1=false, kmer_size=15, max_fixed_edit_rate=0.1,
            max_barcode_mismatch=1, confidence_threshold=0.9,
            max_expected_errors=None,
        )
    )]
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        output: PathBuf,
        barcode_table: PathBuf,
        flanks: Vec<String>,
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
        if flanks.len() != 3 {
            bail!("flanks must contain exactly three fixed sequences");
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

        let barcode_table = FloatingBarcodeTable::new(read_barcode_table(barcode_table)?)?;

        Ok(Self {
            output,
            barcode_table,
            flanks,
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
    let mut name_lines = HashMap::new();
    for (line_index, line) in reader.lines().enumerate() {
        let line = line?;
        if line.trim().is_empty() {
            continue;
        }
        let fields: Vec<_> = line.split('\t').map(str::trim).collect();
        if fields.len() != 3 {
            bail!(
                "barcode table line {} must contain three tab-separated columns",
                line_index + 1
            );
        }
        if let Some(first_line) = name_lines.insert(fields[0].to_owned(), line_index + 1) {
            bail!(
                "barcode table line {} has duplicate name '{}' first seen on line {}",
                line_index + 1,
                fields[0],
                first_line
            );
        }
        entries.push(FloatingBarcodeEntry::new(
            fields[0].as_bytes(),
            fields[1].as_bytes(),
            fields[2].as_bytes(),
        ));
    }
    if entries.is_empty() {
        bail!("barcode_table must contain at least one row");
    }
    Ok(entries)
}

impl FloatingBarcodeStage {
    pub(crate) fn configure_plan(&self, plan: FastqPlan, num_threads: usize) -> Result<FastqPlan> {
        let compression = Compression::try_from(&self.output).ok();
        let output = create_file(&self.output, compression, None, 1)?;
        let extractor = InsertionExtractor::new(
            self.flanks.clone(),
            self.barcode_table.length_profiles().to_vec(),
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

        let stage = CoreFloatingBarcodeStage::new(
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
    module.add_class::<FloatingBarcodeStage>()
}

#[cfg(test)]
mod tests {
    use super::read_barcode_table;
    use std::io::Write;

    #[test]
    fn reads_three_column_barcode_table() {
        let mut file = tempfile::NamedTempFile::new().unwrap();
        writeln!(file, "TF1\tTT\tAAA").unwrap();

        let entries = read_barcode_table(file.path().to_path_buf()).unwrap();
        assert_eq!(entries.len(), 1);
        assert_eq!(entries[0].name, b"TF1");
        assert_eq!(entries[0].barcode_1, b"TT");
        assert_eq!(entries[0].barcode_2, b"AAA");
    }

    #[test]
    fn rejects_legacy_four_column_barcode_table() {
        let mut file = tempfile::NamedTempFile::new().unwrap();
        writeln!(file, "TF1\tSITE1\tTT\tAAA").unwrap();

        let error = read_barcode_table(file.path().to_path_buf()).unwrap_err();
        assert!(error
            .to_string()
            .contains("must contain three tab-separated columns"));
    }

    #[test]
    fn rejects_duplicate_trimmed_barcode_table_names() {
        let mut file = tempfile::NamedTempFile::new().unwrap();
        writeln!(file, "TF1\tTT\tAAA").unwrap();
        writeln!(file, "  TF1  \tGG\tCCC").unwrap();

        let error = read_barcode_table(file.path().to_path_buf()).unwrap_err();
        assert!(error
            .to_string()
            .contains("line 2 has duplicate name 'TF1'"));
        assert!(error.to_string().contains("first seen on line 1"));
    }
}
