mod align;
mod aligners;
mod examples;
mod middleware;
mod pyseqspec;
mod sinks;
mod utils;

use anyhow::Result;
use noodles_fastq as fastq;
use pyo3::prelude::*;
use std::io::Write;
use std::{io::BufWriter, path::PathBuf, str::FromStr};

use ::precellar::align::{extend_fastq_record, Barcode, BarcodeCorrectionConfig, FastqPlan};
use pyseqspec::extract_assays;
use seqspec::{
    utils::{create_file, Compression},
    Modality,
};

#[cfg(not(target_env = "msvc"))]
use tikv_jemallocator::Jemalloc;

#[cfg(not(target_env = "msvc"))]
#[global_allocator]
static GLOBAL: Jemalloc = Jemalloc;

/// Generate consolidated fastq files from the sequencing specification.
/// The barcodes and UMIs are concatenated to the read 1 sequence.
/// Fixed sequences and linkers are removed.
#[pyfunction]
#[pyo3(
    signature = (assay, *, modality, out_dir, correct_barcode=false),
    text_signature = "(assay, *, modality, out_dir, corect_barcode=False)",
)]
fn make_fastq(
    py: Python<'_>,
    assay: Bound<'_, PyAny>,
    modality: &str,
    out_dir: PathBuf,
    correct_barcode: bool,
) -> Result<()> {
    let modality = Modality::from_str(modality).unwrap();
    let assay = extract_assays(assay)?;

    let mut execution = FastqPlan::new(assay, modality)
        .with_barcode_config(BarcodeCorrectionConfig::default())
        .build(correct_barcode, 1000000)?;
    let paired_end = execution.is_paired_end();

    std::fs::create_dir_all(&out_dir)?;
    let read1_fq = out_dir.join("R1.fq.zst");
    let read1_writer = create_file(read1_fq, Some(Compression::Zstd), None, 8)?;
    let mut read1_writer = fastq::io::Writer::new(BufWriter::new(read1_writer));
    let mut read2_writer = if paired_end {
        let read2_fq = out_dir.join("R2.fq.zst");
        let read2_writer = create_file(read2_fq, Some(Compression::Zstd), None, 8)?;
        let read2_writer = fastq::io::Writer::new(BufWriter::new(read2_writer));
        Some(read2_writer)
    } else {
        None
    };

    let mut i = 0;
    while let Some(record_batch) = execution.next_batch()? {
        for record in record_batch {
            if i % 1000000 == 0 {
                py.check_signals().unwrap();
            }
            let Barcode { mut raw, corrected } = record.barcode.unwrap();
            if !correct_barcode || corrected.is_some() {
                if let Some(corrected) = corrected {
                    *raw.sequence_mut() = corrected;
                }
                if let Some(umi) = record.umi {
                    extend_fastq_record(&mut raw, &umi);
                }
                extend_fastq_record(&mut raw, &record.read1.unwrap());

                read1_writer.write_record(&raw)?;
                if let Some(writer) = &mut read2_writer {
                    writer.write_record(&record.read2.unwrap())?;
                }
            }
            i += 1;
        }
    }
    execution.finish()?;

    Ok(())
}

/// A Python module implemented in Rust.
#[pymodule]
fn precellar(m: &Bound<'_, PyModule>) -> PyResult<()> {
    env_logger::builder()
        .format(|buf, record| {
            let timestamp = buf.timestamp();
            let style = buf.default_level_style(record.level());
            writeln!(
                buf,
                "[{timestamp} {style}{}{style:#}] {}",
                record.level(),
                record.args()
            )
        })
        .filter_level(log::LevelFilter::Info)
        .try_init()
        .unwrap();

    m.add("__version__", env!("CARGO_PKG_VERSION"))?;

    m.add_class::<pyseqspec::Assay>()?;

    m.add_function(wrap_pyfunction!(align::make_bwa_mem2_index, m)?)?;
    m.add_function(wrap_pyfunction!(align::make_minibwa_index, m)?)?;
    m.add_function(wrap_pyfunction!(align::make_minimap2_index, m)?)?;
    m.add_class::<align::FastqPipeline>()?;
    m.add_class::<align::AlignmentJob>()?;
    m.add_function(wrap_pyfunction!(make_fastq, m)?)?;

    utils::register_utils(m)?;
    middleware::register_middleware(m)?;
    sinks::register_sinks(m)?;
    aligners::register_aligners(m)?;
    examples::register_examples(m)?;

    Ok(())
}
