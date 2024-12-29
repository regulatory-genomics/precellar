use crate::pyseqspec::Assay;

use anyhow::{bail, Result};
use noodles::sam::{self, alignment::io::Write};
use pyo3::prelude::*;
use std::path::Path;
use std::{collections::HashMap, path::PathBuf, str::FromStr};

use precellar::align::{Aligner, BurrowsWheelerAligner, StarAligner};
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

pub enum AlignerType {
    STAR(StarAligner),
    BWA(BurrowsWheelerAligner),
}

impl AlignerType {
    pub fn from_name<P: AsRef<Path>>(name: &str, path: P) -> Self {
        match name.to_lowercase().as_str() {
            "star" => AlignerType::STAR(StarAligner::from_path(path)),
            "bwa" => AlignerType::BWA(BurrowsWheelerAligner::from_path(path)),
            _ => unimplemented!(),
        }
    }

    pub fn from_modality<P: AsRef<Path>>(modality: Modality, path: P) -> Self {
        match modality {
            Modality::RNA => AlignerType::STAR(StarAligner::from_path(path)),
            Modality::ATAC => AlignerType::BWA(BurrowsWheelerAligner::from_path(path)),
            _ => unimplemented!(),
        }
    }

    pub fn header(&self) -> sam::Header {
        match self {
            AlignerType::STAR(aligner) => aligner.header(),
            AlignerType::BWA(aligner) => aligner.header(),
        }
    }
}

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
/// genom_index: Path
///     File path to the genome index. The genome index can be created by the `make_genome_index` function.
/// modality: str
///     The modality of the sequencing data, e.g., "rna" or "atac".
/// output: Path
///     File path to the output file. The type of the output file is determined by the `output_type` parameter (see below).
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
/// aligner: str | None
///     The aligner to use for the alignment. If None, the aligner will be inferred from the modality.
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
///
/// Returns
/// -------
/// dict
///    A dictionary containing the QC metrics of the alignment and fragment generation.
#[pyfunction]
#[pyo3(
    signature = (
        assay, genome_index, *,
        modality, output, output_type="alignment",
        mito_dna=vec!["chrM".to_owned(), "M".to_owned()],
        shift_left=4, shift_right=-5, aligner=None,
        compression=None, compression_level=None,
        temp_dir=None, num_threads=8, chunk_size=10000000,
    ),
    text_signature = "(assay, genome_index, *,
        modality, output, output_type='alignment',
        mito_dna=['chrM', 'M'],
        shift_left=4, shift_right=-5, aligner=None,
        compression=None, compression_level=None,
        temp_dir=None, num_threads=8, chunk_size=10000000)",
)]
pub fn align(
    py: Python<'_>,
    assay: Bound<'_, PyAny>,
    genome_index: PathBuf,
    modality: &str,
    output: PathBuf,
    output_type: &str,
    mito_dna: Vec<String>,
    shift_left: i64,
    shift_right: i64,
    aligner: Option<&str>,
    compression: Option<&str>,
    compression_level: Option<u32>,
    temp_dir: Option<PathBuf>,
    num_threads: u16,
    chunk_size: usize,
) -> Result<HashMap<String, f64>> {
    let modality = Modality::from_str(modality).unwrap();
    let assay = match assay.extract::<PathBuf>() {
        Ok(p) => seqspec::Assay::from_path(&p).unwrap(),
        _ => assay.extract::<PyRef<Assay>>()?.0.clone(),
    };

    let mut aligner = if let Some(name) = aligner {
        AlignerType::from_name(name, &genome_index)
    } else {
        AlignerType::from_modality(modality, &genome_index)
    };
    let header = aligner.header();
    let mut processor = FastqProcessor::new(assay)
        .with_modality(modality)
        .with_barcode_correct_prob(0.9);
    mito_dna.iter().for_each(|x| {
        processor.add_mito_dna(x);
    });

    let mut transcript_annotator = None;
    let alignments: Box<dyn Iterator<Item = _>> = match aligner {
        AlignerType::STAR(ref mut aligner) => {
            let transcriptome = ::precellar::align::read_transcriptome_star(&genome_index)?;
            transcript_annotator = Some(::precellar::transcript::AlignmentAnnotator::new(
                transcriptome,
            ));
            Box::new(processor.gen_barcoded_alignments(
                aligner,
                num_threads,
                num_threads as usize * chunk_size,
            ))
        }
        AlignerType::BWA(ref mut aligner) => Box::new(processor.gen_barcoded_alignments(
            aligner,
            num_threads,
            num_threads as usize * chunk_size,
        )),
    };

    match OutputType::try_from(output_type)? {
        OutputType::Alignment => {
            write_alignments(py, output, &header, alignments)?;
            Ok(processor.get_report().into())
        }
        OutputType::Fragment => {
            let compression = compression
                .map(|x| Compression::from_str(x).unwrap())
                .or((&output).try_into().ok());
            let mut writer =
                create_file(output, compression, compression_level, num_threads as u32)?;

            let mut fragment_generator = FragmentGenerator::default();
            if let Some(dir) = temp_dir {
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
            let mut quantifier = Quantifier::new(transcript_annotator.unwrap());
            mito_dna.iter().for_each(|x| quantifier.add_mito_dna(x));
            quantifier.quantify(&header, alignments, output)?;
            Ok(processor.get_report().into())
        }
    }
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
