use std::{io::{BufReader, BufWriter}, path::PathBuf, str::FromStr};
use precellar::{io::{open_file_for_read, open_file_for_write, Compression}, utils::strip_barcode_from_read_name};
use noodles::fastq::{self, Reader, io::Writer};
use anyhow::Result;
use pyo3::prelude::*;
use regex::Regex;

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
