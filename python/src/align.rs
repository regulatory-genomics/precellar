use crate::aligners::AlignerRef;
use crate::pyseqspec::extract_assays;

use anyhow::{bail, Result};
use core::panic;
use indicatif::{ProgressBar, ProgressFinish, ProgressStyle};
use log::info;
use noodles::sam::{self, alignment::io::Write};
use precellar::align::{Aligner, AlignmentResult};
use precellar::qc::Metric;
use pyo3::types::PyDict;
use pyo3::{prelude::*, BoundObject};
use seqspec::ChemistryStrandedness;
use serde_json::Value;
use std::{path::PathBuf, str::FromStr};

use precellar::{
    align::{FastqProcessor, MultiMapR},
    fragment::{IntoFragOpts, IntoFragments},
    qc::QcFragment,
    transcriptome::Quantifier,
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

/// Create a minimap2 index from a FASTA file.
///
/// This function creates a `.mmi` index file that can be used with the MINIMAP2 aligner.
/// The index is created with k-mer and window sizes optimized for the selected preset.
///
/// Parameters
/// ----------
/// fasta: Path
///    File path to the FASTA file containing reference sequences.
/// output_index: Path
///   File path for the output minimap2 index (.mmi file).
/// preset: str
///    Optional preset to optimize index for specific read types:
///    - Long Reads DNA Mapping:
///      - 'map-ont': Oxford Nanopore reads (default)
///      - 'map-pb': PacBio CLR reads
///      - 'map-hifi': PacBio HiFi reads
///      - 'lr:hq': Long reads, high quality
///    - Spliced / RNA-seq Alignment:
///      - 'splice': RNA-seq long reads
///      - 'splice:hq': High-quality RNA-seq long reads
///      - 'splice:sr': Short-read RNA-seq
///    - Long Assembly to Reference Mapping:
///      - 'asm5', 'asm10', 'asm20': Assembly alignment (5%, 10%, 20% divergence)
///    - Short Reads Mapping:
///      - 'short': Short single-end reads
///      - 'sr': Short paired-end reads
///    - All-vs-All Overlap Mapping:
///      - 'ava-pb': PacBio all-vs-all overlap
///      - 'ava-ont': ONT all-vs-all overlap
///
/// Examples
/// --------
/// >>> from precellar import make_minimap2_index
/// >>> make_minimap2_index("genome.fa", "genome.mmi", preset="map-ont")
/// >>> # For RNA-seq
/// >>> make_minimap2_index("transcriptome.fa", "transcriptome.mmi", preset="splice")
///
/// See Also
/// --------
/// aligners.MINIMAP2 : The aligner class that uses these indices
/// make_bwa_index : Create BWA-MEM2 index for short reads
#[pyfunction]
#[pyo3(
    signature = (fasta, output_index, *, preset="map-ont"),
    text_signature = "(fasta, output_index, *, preset='map-ont')",
)]
pub fn make_minimap2_index(fasta: PathBuf, output_index: PathBuf, preset: &str) -> Result<()> {
    let preset = match preset.to_lowercase().as_str() {
        // Long Reads DNA Mapping
        "map-ont" => minimap2::Preset::MapOnt,
        "map-pb" => minimap2::Preset::MapPb,
        "map-hifi" => minimap2::Preset::MapHifi,
        "lr:hq" => minimap2::Preset::LrHq,
        // Spliced / RNA-seq Alignment
        "splice" => minimap2::Preset::Splice,
        "splice:hq" => minimap2::Preset::SpliceHq,
        "splice:sr" => minimap2::Preset::SpliceSr,
        // Long Assembly to Reference Mapping
        "asm5" => minimap2::Preset::Asm5,
        "asm10" => minimap2::Preset::Asm10,
        "asm20" => minimap2::Preset::Asm20,
        // Short Reads Mapping
        "short" => minimap2::Preset::Short,
        "sr" => minimap2::Preset::Sr,
        // All-vs-All overlap Mapping
        "ava-pb" => minimap2::Preset::AvaPb,
        "ava-ont" => minimap2::Preset::AvaOnt,
        _ => bail!(
            "Invalid preset '{}'. Valid presets: map-ont, map-pb, map-hifi, lr:hq, splice, splice:hq, splice:sr, asm5, asm10, asm20, short, sr, ava-pb, ava-ont",
            preset,
        ),
    };

    info!(
        "Creating minimap2 index for fasta: {:?} with preset: {:?}",
        fasta, preset
    );
    // Build index from FASTA and save to output
    // with_index(input, Some(output)) reads FASTA from input and saves .mmi to output
    minimap2::Aligner::builder()
        .preset(preset)
        .with_index(
            fasta.to_str().unwrap(),
            Some(output_index.to_str().unwrap()),
        )
        .map_err(|e| anyhow::anyhow!("Failed to create minimap2 index: {}", e))?;

    Ok(())
}

/// Align fastq reads to the reference genome and generate unique fragments.
///
/// Parameters
/// ----------
///
/// assay: Assay | Path | list[Assay | Path]
///     An Assay object or file path to the yaml sequencing specification file, see
///     https://github.com/pachterlab/seqspec. The assay can also be a list of
///     Assay objects or file paths. In this case, the results will be
///     concatenated into a single output file.
/// aligner: STAR | BWAMEM2 | MINIMAP2
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
///     For example, in ATAC-seq, this is usually set to 4 to account for the Tn5 transposase insertion bias.
///     Available only when `output_type='fragment'`.
/// shift_right: int
///     The number of bases to shift the right end of the fragment.
///     For example, in ATAC-seq, this is usually set to -5 to account for the Tn5 transposase insertion bias.
///     Available only when `output_type='fragment'`.
/// compute_snv: bool
///     Whether to compute single nucleotide variants (SNVs) from the alignments.
///     If True, the SNVs will be computed and added to the fragment file.
/// strandedness: Literal['forward', 'reverse', 'unstranded', 'auto'] | None
///     The strand specificity of the assay. Can be "forward", "reverse", "unstranded", "auto".
///     "forward" means that the read 1 is expected to be aligned to the same strand as the original RNA molecule.
///     For example, in 10x Genomics 3' scRNA-seq, the read 1 is aligned to the reverse strand of the original RNA molecule,
///     so the strandedness should be set to "reverse".
///     If "auto", the strand specificity will be inferred from the data. Note that automatic
///     inference may not always be accurate, especially for in the samples where the
///     antisense transcription is prevalent, e.g., brain samples.
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
///
/// See Also
/// --------
/// aligners.BWAMEM2
/// aligners.STAR
/// aligners.MINIMAP2
#[pyfunction]
#[pyo3(
    signature = (
        assay, aligner, *,
        output, modality=None, output_type="alignment",
        mito_dna=vec!["chrM".to_owned(), "M".to_owned()],
        shift_left=4, shift_right=-5, compute_snv=false,
        strandedness=None,
        compression=None, compression_level=None,
        temp_dir=None, num_threads=8, chunk_size=10000000,
    ),
    text_signature = "(assay, aligner, *,
        output, modality=None, output_type='alignment',
        mito_dna=['chrM', 'M'],
        shift_left=4, shift_right=-5, compute_snv=False,
        strandedness=None,
        compression=None, compression_level=None,
        temp_dir=None, num_threads=8, chunk_size=10000000)"
)]
pub fn align<'py>(
    py: Python<'py>,
    assay: Bound<'py, PyAny>,
    aligner: Bound<'py, PyAny>,
    output: PathBuf,
    modality: Option<&str>,
    output_type: &str,
    mito_dna: Vec<String>,
    shift_left: i64,
    shift_right: i64,
    compute_snv: bool,
    strandedness: Option<&str>,
    compression: Option<&str>,
    compression_level: Option<u32>,
    temp_dir: Option<PathBuf>,
    num_threads: u16,
    chunk_size: usize,
) -> Result<Bound<'py, PyAny>> {
    let assay = extract_assays(assay)?;

    let modality = match modality {
        Some(m) => Modality::from_str(m),
        None => {
            let m = assay[0].modality()?;
            info!("Using default modality: {}", m);
            Ok(m)
        }
    }?;

    let mut aligner = AlignerRef::try_from(aligner)?;
    let header = aligner.header();
    let strandedness = if let Some(s) = strandedness {
        match s.to_lowercase().as_str() {
            "forward" => Some(ChemistryStrandedness::Forward),
            "reverse" => Some(ChemistryStrandedness::Reverse),
            "unstranded" => Some(ChemistryStrandedness::Unstranded),
            "auto" => None,
            _ => panic!("strandedness must be 'unstranded', 'forward', 'reverse' or 'auto'"),
        }
    } else {
        if assay[0].chemistry_strandedness.is_some() {
            assay[0].chemistry_strandedness
        } else if modality == Modality::RNA {
            panic!("strandedness must be provided if not specified in the assay. Possible values are 'unstranded', 'forward', 'reverse' or 'auto'")
        } else {
            None
        }
    };

    let transcript_annotator = aligner.transcript_annotator(strandedness);
    let mut processor = FastqProcessor::new(assay)
        .with_modality(modality)
        .with_barcode_correct_prob(0.9);
    if !mito_dna.is_empty() {
        mito_dna.iter().for_each(|x| {
            processor.add_mito_dna(x);
        });
    }

    let mut qc_metrics = serde_json::Map::new();

    let alignments = processor.gen_barcoded_alignments(
        &mut aligner,
        num_threads,
        num_threads as usize * chunk_size,
    );
    let alignments = AlignProgressBar::new(alignments);

    match OutputType::try_from(output_type)? {
        OutputType::Alignment => {
            write_alignments(py, output, &header, alignments)?;
        }
        OutputType::Fragment => {
            let compression = compression
                .map(|x| Compression::from_str(x).unwrap())
                .or((&output).try_into().ok());

            let mut writer = create_file(
                output.clone(),
                compression,
                compression_level,
                num_threads as u32,
            )?;
            writeln!(
                writer,
                "{}",
                fragment_file_header(compute_snv, shift_left, shift_right)
            )?;

            let opts = IntoFragOpts {
                shift_left,
                shift_right,
                temp_dir: temp_dir.clone(),
                compute_snv,
                ..Default::default()
            };

            let mut frag_qc = QcFragment::default();
            mito_dna.iter().for_each(|x| {
                frag_qc.add_mito_dna(x);
            });
            alignments
                .into_fragments(&header, opts)
                .into_iter()
                .for_each(|fragments| {
                    py.check_signals().unwrap();
                    fragments.into_iter().for_each(|frag| {
                        frag_qc.update(&frag);
                        writeln!(writer, "{}", frag).unwrap();
                    })
                });
            qc_metrics.insert("fragment".to_owned(), frag_qc.into());
        }
        OutputType::GeneQuantification => {
            let mut quantifier = Quantifier::new(transcript_annotator.unwrap())?;
            quantifier.num_threads = num_threads as usize;
            quantifier.temp_dir = temp_dir;
            let quant_qc = quantifier.quantify(&header, alignments, output.clone())?;
            qc_metrics.insert("gene_quantification".to_owned(), quant_qc.into());
        }
    };

    qc_metrics.insert(
        "fastq".to_owned(),
        processor.get_fastq_qc().lock().unwrap().to_json(),
    );
    qc_metrics.insert(
        "alignment".to_owned(),
        processor.get_align_qc().lock().unwrap().to_json(),
    );

    Ok(value_into_pyobject(qc_metrics.into(), py))
}

fn fragment_file_header(compute_snv: bool, shift_left: i64, shift_right: i64) -> String {
    let header = if compute_snv {
        "# chromosome\tstart\tend\tbarcode\tcount\tstrand\talignment1_start\talignment1_snv\talignment2_start\talignment2_snv"
    } else {
        "# chromosome\tstart\tend\tbarcode\tcount\tstrand"
    };
    [
        &format!(
            "# This file contains unique fragments generated using precellar-v{}",
            env!("CARGO_PKG_VERSION")
        ),
        "#",
        "# Parameters",
        &format!("# shift_left = {}", shift_left),
        &format!("# shift_right = {}", shift_right),
        "#",
        "# Each line represents a unique fragment with the following fields:",
        header,
    ]
    .join("\n")
}

struct AlignProgressBar<'a, A> {
    pb: ProgressBar,
    alignments: AlignmentResult<'a, A>,
}

impl<'a, A> AlignProgressBar<'a, A> {
    fn new(alignments: AlignmentResult<'a, A>) -> AlignProgressBar<'a, A> {
        let pb = ProgressBar::new(alignments.num_records() as u64);
        let sty = ProgressStyle::with_template(
            "{percent}%|{wide_bar:.cyan/blue}| {human_pos:>}/{human_len:} [{elapsed}<{eta}, {per_sec}]",
        )
        .unwrap();
        pb.set_style(sty);
        AlignProgressBar {
            pb: pb.with_finish(ProgressFinish::Abandon),
            alignments,
        }
    }
}

impl<A: Aligner> Iterator for AlignProgressBar<'_, A> {
    type Item = Vec<(Option<MultiMapR>, Option<MultiMapR>)>;

    fn next(&mut self) -> Option<Self::Item> {
        let item = self.alignments.next();
        self.pb.set_position(self.alignments.num_processed() as u64);
        item
    }
}

fn value_into_pyobject<'py>(val: Value, py: Python<'py>) -> Bound<'py, PyAny> {
    match val {
        Value::Null => py.None().into_bound(py),
        Value::Bool(b) => pyo3::types::PyBool::new(py, b)
            .into_pyobject(py)
            .unwrap()
            .into_bound()
            .into_any(),
        Value::Number(num) => {
            if let Some(n) = num.as_i64() {
                n.into_pyobject(py).unwrap().into_any()
            } else if let Some(n) = num.as_u64() {
                n.into_pyobject(py).unwrap().into_any()
            } else if let Some(n) = num.as_f64() {
                pyo3::types::PyFloat::new(py, n)
                    .into_pyobject(py)
                    .unwrap()
                    .into_any()
            } else {
                panic!("invalid number type")
            }
        }
        Value::String(s) => pyo3::types::PyString::new(py, &s)
            .into_pyobject(py)
            .unwrap()
            .into_any(),
        Value::Array(vec) => {
            let list = pyo3::types::PyList::empty(py);
            for v in vec {
                list.append(value_into_pyobject(v, py)).unwrap();
            }
            list.into_pyobject(py).unwrap().into_any()
        }
        Value::Object(map) => {
            let dict = PyDict::new(py);
            for (key, value) in map {
                dict.set_item(key, value_into_pyobject(value, py)).unwrap();
            }
            dict.into_pyobject(py).unwrap().into_any()
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
