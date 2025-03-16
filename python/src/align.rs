use crate::aligners::AlignerRef;
use crate::pyseqspec::Assay;

use anyhow::{bail, Result};
use noodles::sam::{self, alignment::io::Write};
use pyo3::prelude::*;
use std::ops::DerefMut;
use std::{collections::HashMap, path::PathBuf, str::FromStr, time::Instant};
use log::{info, warn, debug};

use precellar::{
    align::{FastqProcessor, MultiMapR},
    fragment::FragmentGenerator,
    qc::FragmentQC,
    transcript::Quantifier,
};
use seqspec::{
    utils::{create_file, Compression},
    Modality,
};

pub enum OutputType {
    Alignment,
    GeneQuantification,
    Fragment,
}

impl TryFrom<&str> for OutputType {
    type Error = anyhow::Error;

    fn try_from(value: &str) -> Result<Self> {
        match value {
            "alignment" => Ok(OutputType::Alignment),
            "gene_quantification" => Ok(OutputType::GeneQuantification),
            "fragment" => Ok(OutputType::Fragment),
            x => bail!("invalid output type: {}", x),
        }
    }
}

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
pub fn make_bwa_index(fasta: PathBuf, genome_prefix: PathBuf) -> Result<()> {
    bwa_mem2::FMIndex::new(fasta, genome_prefix)?;
    Ok(())
}

/// Align fastq reads to the reference genome and generate unique fragments.
///
/// Parameters
/// ----------
///
/// assay: Assay | Path
///     A Assay object or file path to the yaml sequencing specification file, see
///     https://github.com/pachterlab/seqspec.
/// aligner: STAR | BWAMEM2
///     The aligner to use for the alignment. Available aligners can be found at
///     `precellar.aligners` submodule.
/// output: Path
///     File path to the output file. The type of the output file is determined by the `output_type` parameter (see below).
/// modality: str | None
///     The modality of the sequencing data, e.g., "rna" or "atac".
/// output_type: Literal["alignment", "fragment", "gene_quantification"]
///     The type of the output file. If "alignment", the output will be a BAM file containing the alignments.
///     If "fragment", the output will be a fragment file containing the unique fragments.
///     If "gene_quantification", the output will be a h5ad file containing the gene quantification.
/// mito_dna: list[str]
///     List of mitochondrial DNA names.
/// shift_left: int
///     The number of bases to shift the left end of the fragment.
///     Available only when `output_type='fragment'`.
/// shift_right: int
///     The number of bases to shift the right end of the fragment.
///     Available only when `output_type='fragment'`.
/// compression: str | None
///     The compression algorithm to use for the output fragment file.
///     If None, the compression algorithm will be inferred from the file extension.
/// compression_level: int | None
///     The compression level to use for the output fragment file.
/// temp_dir: Path | None
///     The temporary directory to use.
/// num_threads: int
///     The number of threads to use.
/// chunk_size: int
///     This parameter is used to control the number of bases processed in each chunk per thread.
///     The total number of bases in each chunk is determined by: chunk_size * num_threads.
/// expected_cells: int | None
///     The expected number of cells in the sample. If None, the number of cells will be inferred.
/// barcode_filtering_quantile: float
///     The quantile to use for barcode filtering. Default is 0.99.
/// barcode_bootstrap_samples: int
///     The number of bootstrap samples to use for barcode filtering. Default is 100.
///
/// Returns
/// -------
/// dict
///    A dictionary containing the QC metrics of the alignment and fragment generation.
/// 
/// See Also
/// --------
/// aligners.BWAMEM2
/// aligners.STAR
#[pyfunction]
#[pyo3(
    signature = (
        assay, aligner, *,
        output, modality=None, output_type="alignment",
        mito_dna=vec!["chrM".to_owned(), "M".to_owned()],
        shift_left=4, shift_right=-5,
        compression=None, compression_level=None,
        temp_dir=None, num_threads=8, chunk_size=10000000,
        expected_cells=None, 
        barcode_filtering_quantile=0.99, barcode_bootstrap_samples=100,
    ),
    text_signature = "(assay, aligner, *,
        output, modality=None, output_type='alignment',
        mito_dna=['chrM', 'M'],
        shift_left=4, shift_right=-5,
        compression=None, compression_level=None,
        temp_dir=None, num_threads=8, chunk_size=10000000,
        expected_cells=None, 
        barcode_filtering_quantile=0.99, barcode_bootstrap_samples=100)",
)]
pub fn align(
    py: Python<'_>,
    assay: Bound<'_, PyAny>,
    aligner: Bound<'_, PyAny>,
    output: PathBuf,
    modality: Option<&str>,
    output_type: &str,
    mito_dna: Vec<String>,
    shift_left: i64,
    shift_right: i64,
    compression: Option<&str>,
    compression_level: Option<u32>,
    temp_dir: Option<PathBuf>,
    num_threads: u16,
    chunk_size: usize,
    expected_cells: Option<u32>,
    barcode_filtering_quantile: f64,
    barcode_bootstrap_samples: u32,
) -> Result<HashMap<String, f64>> {
    let start_time = Instant::now();
    debug!("Starting alignment process");
    
    let assay = match assay.extract::<PathBuf>() {
        Ok(p) => {
            debug!("Loading assay from path: {:?}", p);
            seqspec::Assay::from_path(&p).unwrap()
        },
        _ => {
            debug!("Using provided Assay object");
            assay.extract::<PyRef<Assay>>()?.0.clone()
        },
    };
    
    let modality = modality.map(Modality::from_str).unwrap_or(assay.modality())?;
    debug!("Using modality: {:?}", modality);
    
    let mut aligner = AlignerRef::try_from(aligner)?;
    debug!("Initialized aligner: {}", match &aligner {
        AlignerRef::STAR(_) => "STAR",
        AlignerRef::BWA(_) => "BWA-MEM2",
    });
    
    let header = aligner.header();
    let transcript_annotator = aligner.transcript_annotator();

    info!("Initializing FastqProcessor with {} threads and chunk size {}", num_threads, chunk_size);
    let mut processor = FastqProcessor::new(assay)
        .with_modality(modality)
        .with_barcode_correct_prob(0.9);
    
    // Add expected cells if provided
    if let Some(cells) = expected_cells {
        processor = processor.with_expected_cells(cells as usize);
    }
    
    // Configure barcode filtering parameters
    processor = processor
        .with_barcode_filtering_params(barcode_filtering_quantile, barcode_bootstrap_samples as usize);
    
    if !mito_dna.is_empty() {
        debug!("Adding mitochondrial DNA references: {:?}", mito_dna);
        mito_dna.iter().for_each(|x| {
            processor.add_mito_dna(x);
        });
    }

    info!("Generating alignments");
    let alignments: Box<dyn Iterator<Item = _>> = match &mut aligner {
        AlignerRef::STAR(ref mut aligner) => {
            info!("Using STAR aligner");
            Box::new(processor.gen_barcoded_alignments(
                aligner.deref_mut().deref_mut(),
                num_threads,
                num_threads as usize * chunk_size,
            ))
        },
        AlignerRef::BWA(ref mut aligner) => {
            info!("Using BWA-MEM2 aligner");
            Box::new(processor.gen_barcoded_alignments(
                aligner.deref_mut().deref_mut(),
                num_threads,
                num_threads as usize * chunk_size,
            ))
        },
    };

    debug!("Processing output type: {}", output_type);
    let result = match OutputType::try_from(output_type)? {
        OutputType::Alignment => {
            debug!("Writing alignments to BAM file: {:?}", output);
            write_alignments(py, output, &header, alignments)?;
            Ok(processor.get_report().into())
        }
        OutputType::Fragment => {
            debug!("Generating fragments");
            let compression = compression
                .map(|x| Compression::from_str(x).unwrap())
                .or((&output).try_into().ok());
            debug!("Using compression: {:?} with level: {:?}", compression, compression_level);
            
            let mut writer = create_file(output.clone(), compression, compression_level, num_threads as u32)?;
            info!("Created output file: {:?}", output);

            let mut fragment_generator = FragmentGenerator::default();
            if let Some(dir) = temp_dir.clone() {
                debug!("Using temporary directory: {:?}", dir);
                fragment_generator.set_temp_dir(dir);
                fragment_generator.set_shift_left(shift_left);
                fragment_generator.set_shift_right(shift_right);
            }
            
            let frag_qc = write_fragments(
                py,
                &mut writer,
                &header,
                &mito_dna,
                fragment_generator,
                alignments,
            )?;
            let mut qc = processor.get_report();
            frag_qc.report(&mut qc);
            Ok(qc.into())
        }
        OutputType::GeneQuantification => {
            info!("Starting gene quantification");
            let mut quantifier = Quantifier::new(transcript_annotator.unwrap());
            mito_dna.iter().for_each(|x| quantifier.add_mito_dna(x));
            let quant_qc = quantifier.quantify(&header, alignments, output.clone())?;
            info!("Completed gene quantification, writing to: {:?}", output);
            let mut qc = processor.get_report();
            quant_qc.report(&mut qc);
            Ok(qc.into())
        }
    };

    let elapsed = start_time.elapsed();
    info!("Alignment process completed in {:.2?}", elapsed);
    result
}

#[inline]
fn write_alignments<'a>(
    py: Python<'a>,
    output: PathBuf,
    header: &'a sam::Header,
    alignments: impl Iterator<Item = Vec<(Option<MultiMapR>, Option<MultiMapR>)>> + 'a,
) -> Result<()> {
    let mut writer = noodles::bam::io::writer::Builder::default().build_from_path(output)?;
    writer.write_header(&header)?;

    alignments.for_each(move |data| {
        py.check_signals().unwrap();
        data.iter().for_each(|(a, b)| {
            a.as_ref().map(|x| {
                x.iter()
                    .for_each(|x| writer.write_alignment_record(&header, x).unwrap())
            });
            b.as_ref().map(|x| {
                x.iter()
                    .for_each(|x| writer.write_alignment_record(&header, x).unwrap())
            });
        });
    });

    Ok(())
}

#[inline]
fn write_fragments<'a>(
    py: Python<'a>,
    writer: &mut Box<dyn std::io::Write + Send>,
    header: &'a sam::Header,
    mito_dna: &Vec<String>,
    fragment_generator: FragmentGenerator,
    alignments: impl Iterator<Item = Vec<(Option<MultiMapR>, Option<MultiMapR>)>> + 'a,
) -> Result<FragmentQC> {
    let mut fragment_qc = FragmentQC::default();
    mito_dna.iter().for_each(|x| {
        fragment_qc.add_mito_dna(x);
    });
    fragment_generator
        .gen_unique_fragments(&header, alignments)
        .into_iter()
        .for_each(|fragments| {
            py.check_signals().unwrap();
            fragments.into_iter().for_each(|frag| {
                fragment_qc.update(&frag);
                writeln!(writer, "{}", frag).unwrap();
            })
        });
    Ok(fragment_qc)
}
