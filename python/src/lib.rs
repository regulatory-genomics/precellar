use std::{collections::HashMap, ops::Deref, path::PathBuf};

use bwa::{AlignerOpts, BurrowsWheelerAligner, FMIndex, PairedEndStats};
use either::Either;
use pyo3::prelude::*;
use anyhow::Result;
use noodles::sam::alignment::io::Write;

use ::precellar::{align::{Alinger, FastqProcessor}, seqspec::SeqSpec};

#[cfg(not(target_env = "msvc"))]
use tikv_jemallocator::Jemalloc;

#[cfg(not(target_env = "msvc"))]
#[global_allocator]
static GLOBAL: Jemalloc = Jemalloc;


#[pyfunction]
#[pyo3(
    signature = (seqspec, genome_index, *, modality, output, n_jobs=8),
    text_signature = "(seqspec, genome_index, *, modality, output, n_jobs=8)",
)]
fn align(
    py: Python<'_>,
    seqspec: PathBuf,
    genome_index: PathBuf,
    modality: &str,
    output: PathBuf,
    n_jobs: usize,
) -> Result<HashMap<String, f64>> {
    let spec = SeqSpec::from_path(seqspec).unwrap();
    let aligner = BurrowsWheelerAligner::new(
        FMIndex::read(genome_index).unwrap(),
        AlignerOpts::default().set_n_threads(n_jobs),
        PairedEndStats::default()
    );
    let header = aligner.header();
    let mut processor = FastqProcessor::new(spec, aligner).set_modality(modality);

    let mut writer = noodles::bam::io::writer::Builder::default().build_from_path(output)?;
    writer.write_header(&header)?;
    match processor.gen_barcoded_alignments() {
        Either::Left(mut data) => data.try_for_each(|chunk| {
            py.check_signals()?;
            chunk.iter().try_for_each(|x| writer.write_alignment_record(&header, x))?;
            anyhow::Ok(())
        }),
        Either::Right(mut data) => data.try_for_each(|chunk| {
            py.check_signals()?;
            chunk.iter().try_for_each(|(a, b)| {
                writer.write_alignment_record(&header, a)?;
                writer.write_alignment_record(&header, b)?;
                anyhow::Ok(())
            })?;
            Ok(())
        }),
    }?;

    Ok(processor.get_report().unwrap().deref().clone())
}

/// A Python module implemented in Rust.
#[pymodule]
fn precellar(m: &Bound<'_, PyModule>) -> PyResult<()> {
    //pyo3_log::init();

    m.add_function(wrap_pyfunction!(align, m)?)?;
    Ok(())
}
