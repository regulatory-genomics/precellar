mod align;
mod aligners;
mod examples;
mod pyseqspec;
mod utils;

use anyhow::Result;
use noodles::fastq;
use pyo3::prelude::*;
use std::io::Write;
use std::{io::BufWriter, path::PathBuf, str::FromStr};

use ::precellar::align::{extend_fastq_record, Barcode, FastqProcessor};
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

/*
#[pyfunction]
#[pyo3(
    signature = (
        input,
        output,
        *,
        mito_dna=vec!["chrM".to_owned(), "M".to_owned()],
        chunk_size=50000000,
        compression=None,
        compression_level=None,
        temp_dir=None,
        num_threads=8,
    ),
    text_signature = "(input, output, *,
        mito_dna=['chrM', 'M'],
        chunk_size=50000000,
        compression=None, compression_level=None,
        temp_dir=None, num_threads=8)",
)]
fn make_fragment(
    py: Python<'_>,
    input: PathBuf,
    output: PathBuf,
    mito_dna: Vec<String>,
    chunk_size: usize,
    compression: Option<&str>,
    compression_level: Option<u32>,
    temp_dir: Option<PathBuf>,
    num_threads: u32,
) -> Result<HashMap<String, f64>> {
    let file = std::fs::File::open(input)?;
    let decoder = bgzf::MultithreadedReader::with_worker_count(
        (num_threads as usize).try_into().unwrap(),
        file,
    );
    let mut reader = bam::io::Reader::from(decoder);
    let header = reader.read_header()?;

    let mut fragment_qc = FragmentQC::default();
    let mut align_qc = AlignQC::default();
    mito_dna.into_iter().for_each(|mito| {
        fragment_qc.add_mito_dna(&mito);
        header
            .reference_sequences()
            .get_index_of(&bstr::BString::from(mito))
            .map(|x| align_qc.add_mito_dna(x));
    });

    let chunks = NameCollatedRecords::new(reader.records())
        .map(|x| {
            align_qc.add(&header, &x.0, Some(&x.1)).unwrap();
            x
        })
        .chunks(chunk_size);
    let alignments = chunks
        .into_iter()
        .map(|chunk| Either::Right(chunk.collect_vec()));

    let compression = compression
        .map(|x| Compression::from_str(x).unwrap())
        .or((&output).try_into().ok());
    let mut writer = create_file(output, compression, compression_level, num_threads)?;

    let mut fragment_generator = FragmentGenerator::default();
    if let Some(dir) = temp_dir {
        fragment_generator.set_temp_dir(dir)
    };

    fragment_generator
        .gen_unique_fragments(&header, alignments)
        .into_iter()
        .for_each(|fragments| {
            py.check_signals().unwrap();
            fragments.into_iter().for_each(|frag| {
                fragment_qc.update(&frag);
                writeln!(writer.as_mut(), "{}", frag).unwrap();
            })
        });

    let mut report = Metrics::default();
    align_qc.report(&mut report);
    fragment_qc.report(&mut report);
    Ok(report.into())
}
    */

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

    let fq_reader = FastqProcessor::new(assay)
        .with_modality(modality)
        .gen_barcoded_fastq(correct_barcode, 1000000);

    std::fs::create_dir_all(&out_dir)?;
    let read1_fq = out_dir.join("R1.fq.zst");
    let read1_writer = create_file(read1_fq, Some(Compression::Zstd), None, 8)?;
    let mut read1_writer = fastq::io::Writer::new(BufWriter::new(read1_writer));
    let mut read2_writer = if fq_reader.is_paired_end()? {
        let read2_fq = out_dir.join("R2.fq.zst");
        let read2_writer = create_file(read2_fq, Some(Compression::Zstd), None, 8)?;
        let read2_writer = fastq::io::Writer::new(BufWriter::new(read2_writer));
        Some(read2_writer)
    } else {
        None
    };

    for (i, record) in fq_reader.flatten().enumerate() {
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
    }

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

    m.add_function(wrap_pyfunction!(align::make_bwa_index, m)?)?;
    m.add_function(wrap_pyfunction!(align::make_minimap2_index, m)?)?;
    m.add_function(wrap_pyfunction!(align::align, m)?)?;
    //m.add_function(wrap_pyfunction!(make_fragment, m)?)?;
    m.add_function(wrap_pyfunction!(make_fastq, m)?)?;

    utils::register_utils(m)?;
    aligners::register_aligners(m)?;
    examples::register_examples(m)?;

    Ok(())
}
