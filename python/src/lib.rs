use std::{collections::HashMap, path::PathBuf, str::FromStr};

use bwa::{AlignerOpts, BurrowsWheelerAligner, FMIndex, PairedEndStats};
use either::Either;
use pyo3::prelude::*;
use anyhow::Result;
use noodles::sam::alignment::io::Write;

use ::precellar::{align::{Alinger, FastqProcessor}, fragment::FragmentGenerator, io::{open_file_for_write, Compression}, qc::FragmentQC, seqspec::SeqSpec};

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
            if let Some(writer) = &mut bam_writer {
                py.check_signals().unwrap();
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
            let compression = compression.map(|x| Compression::from_str(x).unwrap());
            open_file_for_write(output, compression, compression_level)
        }).transpose()?;
        if let Some(mut writer) = fragment_writer {
            let mut fragment_generator = FragmentGenerator::default();
            if let Some(dir) = temp_dir {
                fragment_generator.set_temp_dir(dir)
            };
            fragment_generator.gen_unique_fragments(&header, alignments).into_iter().for_each(|fragments|
                fragments.into_iter().for_each(|frag| {
                    fragment_qc.update(&frag);
                    writeln!(writer.as_mut(), "{}", frag).unwrap();
                })
            );
        }
    }

    let mut report = processor.get_report();
    if output_fragment.is_some() {
        fragment_qc.report(&mut report);
    }
    Ok(report.into())
}

/// A Python module implemented in Rust.
#[pymodule]
fn precellar(m: &Bound<'_, PyModule>) -> PyResult<()> {
    //pyo3_log::init();

    m.add_function(wrap_pyfunction!(align, m)?)?;
    Ok(())
}
