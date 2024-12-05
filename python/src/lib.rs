mod pyseqspec;
mod utils;

use anyhow::Result;
use either::Either;
use itertools::Itertools;
use log::info;
use noodles::{bam, bgzf, fastq, sam::{self, alignment::io::Write}};
use pyo3::prelude::*;
use std::{collections::HashMap, io::BufWriter, path::{Path, PathBuf}, str::FromStr};

use ::precellar::{
    align::{
        extend_fastq_record, Aligner, AlignerBuilder, Barcode, BurrowsWheelerAligner, DummyAligner, FastqProcessor, MultiMapR, NameCollatedRecords, StarAligner
    },
    fragment::FragmentGenerator,
    qc::{AlignQC, FragmentQC, Metrics},
};
use pyseqspec::Assay;
use seqspec::{
    utils::{create_file, Compression},
    Modality,
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
fn make_genome_index(fasta: PathBuf, genome_prefix: PathBuf) -> Result<()> {
    //FMIndex::new(fasta, genome_prefix).unwrap();
    todo!()
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
///     This parameter is used to control the number of bases processed in each chunk.
///     The actual value is determined by: chunk_size * num_threads.
///
/// Returns
/// -------
/// dict
///    A dictionary containing the QC metrics of the alignment and fragment generation.
#[pyfunction]
#[pyo3(
    signature = (
        assay, genome_index, *,
        modality, output_bam=None, output_fragment=None,
        mito_dna=vec!["chrM".to_owned(), "M".to_owned()],
        shift_left=4, shift_right=-5, aligner=None,
        compression=None, compression_level=None,
        temp_dir=None, num_threads=8, chunk_size=10000000,
    ),
    text_signature = "(assay, genome_index, *,
        modality, output_bam=None, output_fragment=None,
        mito_dna=['chrM', 'M'],
        shift_left=4, shift_right=-5, aligner=None,
        compression=None, compression_level=None,
        temp_dir=None, num_threads=8, chunk_size=10000000)",
)]
fn align(
    py: Python<'_>,
    assay: Bound<'_, PyAny>,
    genome_index: PathBuf,
    modality: &str,
    output_bam: Option<PathBuf>,
    output_fragment: Option<PathBuf>,
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
    assert!(
        output_bam.is_some() || output_fragment.is_some(),
        "one of the following parameters must be provided: output_bam, output_fragment"
    );

    let modality = Modality::from_str(modality).unwrap();
    let assay = match assay.extract::<PathBuf>() {
        Ok(p) => seqspec::Assay::from_path(&p).unwrap(),
        _ => assay.extract::<PyRef<Assay>>()?.0.clone(),
    };

    let aligner = if let Some(name) = aligner {
        get_aligner_by_name(name, genome_index)
    } else {
        get_aligner_by_modality(modality, genome_index)
    };
    let header = aligner.header();
    let mut processor = FastqProcessor::new(assay, aligner)
        .with_modality(modality)
        .with_barcode_correct_prob(0.9);
    let mut fragment_qc = FragmentQC::default();
    mito_dna.into_iter().for_each(|x| {
        processor.add_mito_dna(&x);
        fragment_qc.add_mito_dna(x);
    });

    {
        let mut bam_writer = output_bam
            .map(|output| {
                let mut writer =
                    noodles::bam::io::writer::Builder::default().build_from_path(output)?;
                writer.write_header(&header)?;
                anyhow::Ok(writer)
            })
            .transpose()?;
        let alignments = processor
            .gen_barcoded_alignments(num_threads, chunk_size);
        let alignments = write_alignments(py, &mut bam_writer, &header, alignments);

        let fragment_writer = output_fragment
            .as_ref()
            .map(|output| {
                let compression = compression
                    .map(|x| Compression::from_str(x).unwrap())
                    .or(output.try_into().ok());
                create_file(output, compression, compression_level, num_threads as u32)
            })
            .transpose()?;
        if let Some(mut writer) = fragment_writer {
            let mut fragment_generator = FragmentGenerator::default();
            if let Some(dir) = temp_dir {
                fragment_generator.set_temp_dir(dir);
                fragment_generator.set_shift_left(shift_left);
                fragment_generator.set_shift_right(shift_right);
            }
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
        } else {
            // alignments is a lazy iterator, so we need to consume it if no other
            // output is generated.
            alignments.for_each(drop);
        }
    }

    let mut report = processor.get_report();
    if output_fragment.is_some() {
        fragment_qc.report(&mut report);
    }
    Ok(report.into())
}

fn get_aligner_by_modality<P: AsRef<Path>>(modality: Modality, path: P) -> Box<dyn Aligner> {
    match modality {
        Modality::RNA => Box::new(StarAligner::from_path(path)),
        Modality::ATAC => Box::new(BurrowsWheelerAligner::from_path(path)),
        _ => unimplemented!(),
    }
}

fn get_aligner_by_name<P: AsRef<Path>>(name: &str, path: P) -> Box<dyn Aligner> {
    match name.to_lowercase().as_str() {
        "star" => Box::new(StarAligner::from_path(path)),
        "bwa" => Box::new(BurrowsWheelerAligner::from_path(path)),
        _ => unimplemented!(),
    }
}

fn write_alignments<'a>(
    py: Python<'a>,
    bam_writer: &'a mut Option<bam::io::Writer<impl std::io::Write>>,
    header: &'a sam::Header,
    alignments: impl Iterator<Item = Either<Vec<MultiMapR>, Vec<(MultiMapR, MultiMapR)>>> + 'a,
) -> impl Iterator<Item = Either<Vec<MultiMapR>, Vec<(MultiMapR, MultiMapR)>>> + 'a
{
    alignments.map(move |data| {
        py.check_signals().unwrap();
        if let Some(writer) = bam_writer.as_mut() {
            match data.as_ref() {
                Either::Left(chunk) => chunk
                    .iter()
                    .for_each(|a| a.iter().for_each(|x|
                        writer.write_alignment_record(&header, x).unwrap()
                    )),
                Either::Right(chunk) => chunk.iter().for_each(|(a, b)| {
                    a.iter().for_each(|x| writer.write_alignment_record(&header, x).unwrap());
                    b.iter().for_each(|x| writer.write_alignment_record(&header, x).unwrap());
                }),
            };
        }
        data
    })
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
    let spec = match assay.extract::<PathBuf>() {
        Ok(p) => seqspec::Assay::from_path(&p).unwrap(),
        _ => assay.extract::<PyRef<Assay>>()?.0.clone(),
    };

    let aligner = DummyAligner;
    let mut processor = FastqProcessor::new(spec, Box::new(aligner)).with_modality(modality);
    let fq_reader = processor.gen_barcoded_fastq(correct_barcode);

    info!(
        "Adding these to the start of Read 1: {}",
        fq_reader
            .get_all_barcodes()
            .into_iter()
            .chain(fq_reader.get_all_umi())
            .map(|(x, n)| format!("{} ({})", x, n))
            .join(" + "),
    );

    std::fs::create_dir_all(&out_dir)?;
    let read1_fq = out_dir.join("R1.fq.zst");
    let read1_writer = create_file(read1_fq, Some(Compression::Zstd), None, 8)?;
    let mut read1_writer = fastq::Writer::new(BufWriter::new(read1_writer));
    let mut read2_writer = if fq_reader.is_paired_end() {
        let read2_fq = out_dir.join("R2.fq.zst");
        let read2_writer = create_file(read2_fq, Some(Compression::Zstd), None, 8)?;
        let read2_writer = fastq::Writer::new(BufWriter::new(read2_writer));
        Some(read2_writer)
    } else {
        None
    };

    for (i, record) in fq_reader.enumerate() {
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
        .filter_level(log::LevelFilter::Info)
        .try_init()
        .unwrap();

    m.add("__version__", env!("CARGO_PKG_VERSION"))?;

    m.add_class::<pyseqspec::Assay>().unwrap();

    m.add_function(wrap_pyfunction!(make_genome_index, m)?)?;
    m.add_function(wrap_pyfunction!(align, m)?)?;
    m.add_function(wrap_pyfunction!(make_fragment, m)?)?;
    m.add_function(wrap_pyfunction!(make_fastq, m)?)?;

    utils::register_submodule(m)?;
    Ok(())
}
