use std::{io::{BufReader, BufWriter}, path::PathBuf, str::FromStr};
use precellar::utils::strip_barcode_from_read_name;
use seqspec::utils::{open_file_for_read, open_file_for_write, Compression};
use noodles::fastq::{self, Reader, io::Writer};
use anyhow::Result;
use pyo3::prelude::*;
use regex::Regex;

/// Remove barcode from the read names of fastq records.
/// 
/// The old practice of storing barcodes in read names is not recommended. This
/// function extracts barcodes from read names using regular expressions and
/// writes them to a separate fastq file.
/// 
/// Parameters
/// ----------
/// in_fq: Path
///     File path to the input fastq file.
/// out_fq: Path
///     File path to the output fastq file.
/// regex: str
///     Extract barcodes from read names of BAM records using regular expressions.
///     Reguler expressions should contain exactly one capturing group 
///     (Parentheses group the regex between them) that matches
///     the barcodes. For example, `barcode_regex="(..:..:..:..):\\w+$"`
///     extracts `bd:69:Y6:10` from
///     `A01535:24:HW2MMDSX2:2:1359:8513:3458:bd:69:Y6:10:TGATAGGTTG`.
///     You can test your regex on this website: https://regex101.com/.
/// out_barcode: Path | None
///     File path to the output fastq file containing the extracted barcodes.
///     If None, the barcodes will not be written to a file.
/// left_add: int
///     Additional bases to strip from the left side of the barcode.
/// right_add: int
///     Additional bases to strip from the right side of the barcode.
/// compression: str | None
///     Compression algorithm to use. If None, the compression algorithm will be inferred from the file extension.
/// compression_level: int | None
///     Compression level to use.
/// num_threads: int
///     The number of threads to use.
#[pyfunction]
#[pyo3(
    signature = (in_fq, out_fq,
        *, regex, out_barcode=None, left_add=0, right_add=0,
        compression=None, compression_level=None, num_threads=8,
    ),
    text_signature = "(in_fq, out_fq,
        *, regex, out_barcode, left_add=0, right_add=0,
        compression=None, compression_level=None, num_threads=8)",
)]
fn strip_barcode_from_fastq(
    py: Python<'_>,
    in_fq: PathBuf,
    out_fq: PathBuf,
    regex: &str,
    out_barcode: Option<PathBuf>,
    left_add: usize,
    right_add: usize,
    compression: Option<&str>,
    compression_level: Option<u32>,
    num_threads: u32,
) -> Result<()> {
    let mut reader = Reader::new(BufReader::new(open_file_for_read(in_fq)));

    let mut fq_writer = {
        let compression = compression.map(|x| Compression::from_str(x).unwrap())
            .or((&out_fq).try_into().ok());
        Writer::new(BufWriter::new(open_file_for_write(out_fq, compression, compression_level, num_threads)?))
    };
    let mut barcode_writer = out_barcode.as_ref().map(|output| {
        let compression = compression.map(|x| Compression::from_str(x).unwrap())
            .or(output.try_into().ok());
        Writer::new(BufWriter::new(open_file_for_write(output, compression, compression_level, num_threads).unwrap()))
    });

    let re = Regex::new(regex)?;
    reader.records().enumerate().try_for_each(|(i, record)| {
        if i % 1000000 == 0 {
            py.check_signals().unwrap();
        }
        let record = record?;
        let (record, barcode) = strip_barcode_from_read_name(record, &re, left_add, right_add)?;

        fq_writer.write_record(&record)?;
        if let Some(barcode_writer) = &mut barcode_writer {
            let qual_score = vec![b'~'; barcode.len()];
            let fq = fastq::Record::new(
                fastq::record::Definition::new(format!("BARCODE.{}", i), ""),
                barcode.into_bytes(),
                qual_score,
            );
            barcode_writer.write_record(&fq)?;
        }

        anyhow::Ok(())
    })?;

    Ok(())
}

#[pymodule]
pub(crate) fn register_submodule(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let utils = PyModule::new_bound(parent_module.py(), "utils")?;

    utils.add_function(wrap_pyfunction!(strip_barcode_from_fastq, &utils)?)?;

    parent_module.add_submodule(&utils)
}
