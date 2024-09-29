mod utils;

use std::{collections::HashMap, path::PathBuf, str::FromStr};
use bwa_mem2::{AlignerOpts, BurrowsWheelerAligner, FMIndex, PairedEndStats};
use either::Either;
use pyo3::prelude::*;
use anyhow::Result;
use noodles::{bgzf, sam::alignment::io::Write};
use noodles::bam;
use itertools::Itertools;

use ::precellar::{
    align::{Alinger, FastqProcessor, NameCollatedRecords},
    fragment::FragmentGenerator,
    io::{open_file_for_write, Compression},
    qc::{FragmentQC, Metrics, AlignQC}, seqspec::{Assay, Modality},
};

#[cfg(not(target_env = "msvc"))]
use tikv_jemallocator::Jemalloc;

#[cfg(not(target_env = "msvc"))]
#[global_allocator]
static GLOBAL: Jemalloc = Jemalloc;

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
fn make_genome_index(
    fasta: PathBuf,
    genome_prefix: PathBuf,
) -> Result<()> {
    FMIndex::new(fasta, genome_prefix).unwrap();
    Ok(())
}

/// Align fastq reads to the reference genome and generate unique fragments.
///
/// Parameters
/// ----------
///
/// seqspec: Path
///     File path to the sequencing specification, see https://github.com/pachterlab/seqspec.
/// genom_index: Path
///     File path to the genome index. The genome index can be created by the `make_genome_index` function.
/// modality: str
///     The modality of the sequencing data, e.g., "rna" or "atac".
/// output_bam: Path | None
///     File path to the output bam file. If None, the bam file will not be generated.
/// output_fragment: Path | None
///     File path to the output fragment file. If None, the fragment file will not be generated.
/// mito_dna: list[str]
///     List of mitochondrial DNA names.
/// shift_left: int
///     The number of bases to shift the left end of the fragment.
/// shift_right: int
///     The number of bases to shift the right end of the fragment.
/// compression: str | None
///     The compression algorithm to use for the output fragment file.
///     If None, the compression algorithm will be inferred from the file extension.
/// compression_level: int | None
///     The compression level to use for the output fragment file.
/// temp_dir: Path | None
///     The temporary directory to use.
/// num_threads: int
///     The number of threads to use.
/// 
/// Returns
/// -------
/// dict
///    A dictionary containing the QC metrics of the alignment and fragment generation.
#[pyfunction]
#[pyo3(
    signature = (
        seqspec, genome_index, *,
        modality, output_bam=None, output_fragment=None,
        mito_dna=vec!["chrM".to_owned(), "M".to_owned()],
        shift_left=4, shift_right=-5,
        compression=None, compression_level=None,
        temp_dir=None, num_threads=8,
    ),
    text_signature = "(seqspec, genome_index, *,
        modality, output_bam=None, output_fragment=None,
        mito_dna=['chrM', 'M'],
        shift_left=4, shift_right=-5,
        compression=None, compression_level=None,
        temp_dir=None, num_threads=8)",
)]
fn align(
    py: Python<'_>,
    seqspec: PathBuf,
    genome_index: PathBuf,
    modality: &str,
    output_bam: Option<PathBuf>,
    output_fragment: Option<PathBuf>,
    mito_dna: Vec<String>,
    shift_left: i64,
    shift_right: i64,
    compression: Option<&str>,
    compression_level: Option<u32>,
    temp_dir: Option<PathBuf>,
    num_threads: u32,
) -> Result<HashMap<String, f64>> {
    let modality = Modality::from_str(modality).unwrap();

    assert!(output_bam.is_some() || output_fragment.is_some(), "either output_bam or output_fragment must be provided");

    let spec = Assay::from_path(&seqspec).unwrap();
    let aligner = BurrowsWheelerAligner::new(
        FMIndex::read(genome_index).unwrap(),
        AlignerOpts::default().with_n_threads(num_threads as usize),
        PairedEndStats::default()
    );
    let header = aligner.header();
    let mut processor = FastqProcessor::new(spec, aligner).with_modality(modality)
        .with_barcode_correct_prob(0.9)
        .with_base_dir(seqspec.parent().unwrap());
    let mut fragment_qc = FragmentQC::default();
    mito_dna.into_iter().for_each(|x| {
        processor.add_mito_dna(&x);
        fragment_qc.add_mito_dna(x);
    });

    {
        let mut bam_writer = output_bam.map(|output| {
            let mut writer = noodles::bam::io::writer::Builder::default().build_from_path(output)?;
            writer.write_header(&header)?;
            anyhow::Ok(writer)
        }).transpose()?;
        let alignments = processor.gen_barcoded_alignments().map(|data| {
            py.check_signals().unwrap();
            if let Some(writer) = &mut bam_writer {
                match data.as_ref() {
                    Either::Left(chunk) => chunk.iter().for_each(|x| writer.write_alignment_record(&header, x).unwrap()),
                    Either::Right(chunk) => chunk.iter().for_each(|(a, b)| {
                        writer.write_alignment_record(&header, a).unwrap();
                        writer.write_alignment_record(&header, b).unwrap();
                    }),
                };
            }
            data
        });
        
        let fragment_writer = output_fragment.as_ref().map(|output| {
            let compression = compression.map(|x| Compression::from_str(x).unwrap())
                .or(output.try_into().ok());
            open_file_for_write(output, compression, compression_level, num_threads)
        }).transpose()?;
        if let Some(mut writer) = fragment_writer {
            let mut fragment_generator = FragmentGenerator::default();
            if let Some(dir) = temp_dir {
                fragment_generator.set_temp_dir(dir);
                fragment_generator.set_shift_left(shift_left);
                fragment_generator.set_shift_right(shift_right);
            }
            fragment_generator.gen_unique_fragments(&header, alignments).into_iter().for_each(|fragments| {
                py.check_signals().unwrap();
                fragments.into_iter().for_each(|frag| {
                    fragment_qc.update(&frag);
                    writeln!(writer.as_mut(), "{}", frag).unwrap();
                })
            });
        }
    }

    let mut report = processor.get_report();
    if output_fragment.is_some() {
        fragment_qc.report(&mut report);
    }
    Ok(report.into())
}

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
    text_signature = "(seqspec, genome_index, *, modality, output_bam=None, output_fragment=None, num_threads=8)",
)]
fn make_fragment(
    py: Python<'_>,
    input: PathBuf,
    output : PathBuf,
    mito_dna: Vec<String>,
    chunk_size: usize,
    compression: Option<&str>,
    compression_level: Option<u32>,
    temp_dir: Option<PathBuf>,
    num_threads: u32,
) -> Result<HashMap<String, f64>> {
    let file = std::fs::File::open(input)?;
    let decoder = bgzf::MultithreadedReader::with_worker_count((num_threads as usize).try_into().unwrap(), file);
    let mut reader = bam::io::Reader::from(decoder);
    let header = reader.read_header()?;

    let mut fragment_qc = FragmentQC::default();
    let mut align_qc = AlignQC::default();
    mito_dna.into_iter().for_each(|mito| {
        fragment_qc.add_mito_dna(&mito);
        header.reference_sequences().get_index_of(&bstr::BString::from(mito)).map(|x| align_qc.add_mito_dna(x));
    });

    let chunks = NameCollatedRecords::new(reader.records()).map(|x| {
        align_qc.update(&x.0, &header);
        align_qc.update(&x.1, &header);
        x
    }).chunks(chunk_size);
    let alignments = chunks.into_iter().map(|chunk| Either::Right(chunk.collect_vec()));

    let compression = compression.map(|x| Compression::from_str(x).unwrap())
        .or((&output).try_into().ok());
    let mut writer = open_file_for_write(output, compression, compression_level, num_threads)?;

    let mut fragment_generator = FragmentGenerator::default();
    if let Some(dir) = temp_dir {
        fragment_generator.set_temp_dir(dir)
    };

    fragment_generator.gen_unique_fragments(&header, alignments).into_iter().for_each(|fragments| {
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


/// A Python module implemented in Rust.
#[pymodule]
fn precellar(m: &Bound<'_, PyModule>) -> PyResult<()> {
    env_logger::builder()
        .filter_level(log::LevelFilter::Info)
        .try_init().unwrap();

    m.add("__version__", env!("CARGO_PKG_VERSION"))?;

    m.add_function(wrap_pyfunction!(make_genome_index, m)?)?;
    m.add_function(wrap_pyfunction!(align, m)?)?;
    m.add_function(wrap_pyfunction!(make_fragment, m)?)?;

    utils::register_submodule(m)?;
    Ok(())
}
