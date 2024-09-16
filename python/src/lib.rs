use std::{collections::HashMap, path::PathBuf, str::FromStr};

use bwa::{AlignerOpts, BurrowsWheelerAligner, FMIndex, PairedEndStats};
use either::Either;
use log::info;
use pyo3::prelude::*;
use anyhow::Result;
use noodles::{bgzf, sam::alignment::io::Write};
use noodles::bam;
use itertools::Itertools;

use ::precellar::{
    align::{Alinger, FastqProcessor, NameCollatedRecords},
    fragment::FragmentGenerator,
    io::{open_file_for_write, Compression},
    qc::{FragmentQC, Metrics, AlignQC}, seqspec::SeqSpec,
};

#[cfg(not(target_env = "msvc"))]
use tikv_jemallocator::Jemalloc;

#[cfg(not(target_env = "msvc"))]
#[global_allocator]
static GLOBAL: Jemalloc = Jemalloc;

#[pyfunction]
#[pyo3(
    signature = (
        seqspec,
        genome_index,
        *,
        modality,
        output_bam=None,
        output_fragment=None,
        mito_dna=vec!["chrM".to_owned(), "M".to_owned()],
        compression=None,
        compression_level=None,
        temp_dir=None,
        n_jobs=8,
    ),
    text_signature = "(seqspec, genome_index, *, modality, output_bam=None, output_fragment=None, n_jobs=8)",
)]
fn align(
    py: Python<'_>,
    seqspec: PathBuf,
    genome_index: PathBuf,
    modality: &str,
    output_bam: Option<PathBuf>,
    output_fragment: Option<PathBuf>,
    mito_dna: Vec<String>,
    compression: Option<&str>,
    compression_level: Option<u32>,
    temp_dir: Option<PathBuf>,
    n_jobs: usize,
) -> Result<HashMap<String, f64>> {
    assert!(output_bam.is_some() || output_fragment.is_some(), "either output_bam or output_fragment must be provided");

    let spec = SeqSpec::from_path(seqspec).unwrap();
    let aligner = BurrowsWheelerAligner::new(
        FMIndex::read(genome_index).unwrap(),
        AlignerOpts::default().set_n_threads(n_jobs),
        PairedEndStats::default()
    );
    let header = aligner.header();
    let mut processor = FastqProcessor::new(spec, aligner).set_modality(modality);
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
            open_file_for_write(output, compression, compression_level)
        }).transpose()?;
        if let Some(mut writer) = fragment_writer {
            info!("Writing fragments to {}", output_fragment.as_ref().unwrap().to_string_lossy());
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
        compression=None,
        compression_level=None,
        temp_dir=None,
        n_jobs=8,
    ),
    text_signature = "(seqspec, genome_index, *, modality, output_bam=None, output_fragment=None, n_jobs=8)",
)]
fn make_fragment(
    py: Python<'_>,
    input: PathBuf,
    output : PathBuf,
    mito_dna: Vec<String>,
    compression: Option<&str>,
    compression_level: Option<u32>,
    temp_dir: Option<PathBuf>,
    n_jobs: usize,
) -> Result<HashMap<String, f64>> {
    let file = std::fs::File::open(input)?;
    let decoder = bgzf::MultithreadedReader::with_worker_count(n_jobs.try_into().unwrap(), file);
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
    }).chunks(5000000);
    let alignments = chunks.into_iter().map(|chunk| Either::Right(chunk.collect_vec()));

    let compression = compression.map(|x| Compression::from_str(x).unwrap())
        .or((&output).try_into().ok());
    let mut writer = open_file_for_write(output, compression, compression_level)?;

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

    m.add_function(wrap_pyfunction!(align, m)?)?;
    m.add_function(wrap_pyfunction!(make_fragment, m)?)?;
    Ok(())
}
