mod align;
mod aligners;
mod examples;
mod middleware;
mod pyseqspec;
mod sinks;
mod utils;

use anyhow::{bail, Result};
use noodles_fastq as fastq;
use pyo3::prelude::*;
use std::io::Write;
use std::{io::BufWriter, path::PathBuf, str::FromStr};

use ::precellar::align::{extend_fastq_record, Barcode, BarcodeCorrectionConfig, FastqPlan};
use pyseqspec::extract_assays;
use seqspec::{
    utils::{create_file, Compression},
    Assay, Modality,
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
    let modality = Modality::from_str(modality)?;
    let assay = extract_assays(assay)?;

    make_fastq_from_assays(assay, modality, out_dir, correct_barcode, || {
        py.check_signals()?;
        Ok(())
    })
}

fn make_fastq_from_assays<F>(
    assay: Vec<Assay>,
    modality: Modality,
    out_dir: PathBuf,
    correct_barcode: bool,
    mut check_signals: F,
) -> Result<()>
where
    F: FnMut() -> Result<()>,
{
    let mut execution = FastqPlan::new(assay, modality)
        .with_barcode_config(BarcodeCorrectionConfig::default())
        .build(correct_barcode, 1000000)?;
    if execution.is_long_read() {
        bail!(
            "make_fastq does not support long-read assays yet; use FastqPipeline.align_with(...)"
        );
    }
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
                check_signals()?;
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

#[cfg(test)]
mod tests {
    use super::*;
    use seqspec::File as SeqFile;
    use std::io::BufReader;
    use std::path::Path;

    fn write_fastq(path: &Path, name: &str, sequence: &[u8]) {
        std::fs::write(
            path,
            format!(
                "@{}\n{}\n+\n{}\n",
                name,
                String::from_utf8_lossy(sequence),
                "I".repeat(sequence.len())
            ),
        )
        .unwrap();
    }

    fn short_read_assay(fastq_path: &Path) -> Assay {
        let yaml = r#"
!Assay
seqspec_version: 0.3.0
assay_id: make-fastq-short
name: make-fastq-short
doi: none
date: none
description: synthetic short-read make_fastq test
modalities: [rna]
lib_struct: synthetic
library_protocol: synthetic
library_kit: synthetic
sequence_protocol: synthetic
sequence_kit: synthetic
sequence_spec:
- !Read
  read_id: R1
  name: Read 1
  modality: rna
  primer_id: primer
  min_len: 12
  max_len: 12
  strand: pos
library_spec:
- !Region
  region_id: rna
  region_type: rna
  name: RNA
  sequence_type: joined
  sequence: NNNNNNNNNNNNNNNN
  min_len: 16
  max_len: 16
  onlist: null
  regions:
  - !Region
    region_id: primer
    region_type: truseq_read1
    name: primer
    sequence_type: fixed
    sequence: AAAA
    min_len: 4
    max_len: 4
    onlist: null
    regions: null
  - !Region
    region_id: barcode
    region_type: barcode
    name: barcode
    sequence_type: random
    sequence: NNNN
    min_len: 4
    max_len: 4
    onlist: null
    regions: null
  - !Region
    region_id: umi
    region_type: umi
    name: UMI
    sequence_type: random
    sequence: NN
    min_len: 2
    max_len: 2
    onlist: null
    regions: null
  - !Region
    region_id: target
    region_type: cdna
    name: target
    sequence_type: random
    sequence: NNNNNN
    min_len: 6
    max_len: 6
    onlist: null
    regions: null
"#;
        let mut assay: Assay = serde_yaml::from_str(yaml).unwrap();
        assay.sequence_spec.get_mut("R1").unwrap().files =
            Some(vec![SeqFile::from_fastq(fastq_path, false).unwrap()]);
        assay
    }

    fn long_read_assay(fastq_path: &Path, directory: &Path) -> Assay {
        let whitelist = directory.join("barcodes.txt");
        std::fs::write(&whitelist, "AACC\n").unwrap();
        let yaml = format!(
            r#"
!Assay
seqspec_version: 0.3.0
assay_id: make-fastq-long
name: make-fastq-long
doi: none
date: none
description: synthetic long-read make_fastq test
modalities: [atac]
lib_struct: synthetic
library_protocol: synthetic
library_kit: synthetic
sequence_protocol: synthetic
sequence_kit: synthetic
sequence_spec:
- !Read
  read_id: R1
  name: Read 1
  modality: atac
  primer_id: primer
  min_len: 636
  max_len: 636
  strand: unstranded
library_spec:
- !Region
  region_id: atac
  region_type: atac
  name: ATAC
  sequence_type: joined
  sequence: synthetic
  min_len: 640
  max_len: 640
  onlist: null
  regions:
  - !Region
    region_id: primer
    region_type: truseq_read1
    name: primer
    sequence_type: fixed
    sequence: AAAA
    min_len: 4
    max_len: 4
    onlist: null
    regions: null
  - !Region
    region_id: anchor1
    region_type: linker
    name: anchor1
    sequence_type: fixed
    sequence: ACGTTGCAACGATTCG
    min_len: 16
    max_len: 16
    onlist: null
    regions: null
  - !Region
    region_id: barcode
    region_type: barcode
    name: barcode
    sequence_type: onlist
    sequence: NNNN
    min_len: 4
    max_len: 4
    onlist: !Onlist
      file_id: barcodes.txt
      filename: barcodes.txt
      filetype: txt
      filesize: 0
      url: "{}"
      urltype: local
      md5: ""
      location: local
    regions: null
  - !Region
    region_id: anchor2
    region_type: linker
    name: anchor2
    sequence_type: fixed
    sequence: GATCTAGCGTACCTGA
    min_len: 16
    max_len: 16
    onlist: null
    regions: null
  - !Region
    region_id: target
    region_type: gdna
    name: target
    sequence_type: random
    sequence: N
    min_len: 600
    max_len: 600
    onlist: null
    regions: null
"#,
            whitelist.display()
        );
        let mut assay: Assay = serde_yaml::from_str(&yaml).unwrap();
        assay.sequence_spec.get_mut("R1").unwrap().files =
            Some(vec![SeqFile::from_fastq(fastq_path, false).unwrap()]);
        assay
    }

    #[test]
    fn make_fastq_writes_short_read_output() {
        let directory = tempfile::tempdir().unwrap();
        let fastq_path = directory.path().join("short.fastq");
        let sequence = b"ACGTTTGGGGGG";
        write_fastq(&fastq_path, "short", sequence);
        let assay = short_read_assay(&fastq_path);
        let output = directory.path().join("output");

        make_fastq_from_assays(vec![assay], Modality::RNA, output.clone(), false, || Ok(()))
            .unwrap();

        let reader = seqspec::utils::open_file(output.join("R1.fq.zst")).unwrap();
        let mut reader = fastq::io::Reader::new(BufReader::new(reader));
        let mut record = fastq::Record::default();
        assert_ne!(reader.read_record(&mut record).unwrap(), 0);
        assert_eq!(record.name(), b"short");
        assert_eq!(record.sequence(), sequence);
        assert_eq!(record.quality_scores(), b"IIIIIIIIIIII");
        assert_eq!(reader.read_record(&mut record).unwrap(), 0);
        assert!(!output.join("R2.fq.zst").exists());
    }

    #[test]
    fn make_fastq_rejects_long_read_before_writing_output() {
        let directory = tempfile::tempdir().unwrap();
        let fastq_path = directory.path().join("long.fastq");
        let sequence = [
            b"ACGTTGCAACGATTCG".as_slice(),
            b"AACC".as_slice(),
            b"GATCTAGCGTACCTGA".as_slice(),
            vec![b'G'; 600].as_slice(),
        ]
        .concat();
        write_fastq(&fastq_path, "long", &sequence);
        let assay = long_read_assay(&fastq_path, directory.path());
        let output = directory.path().join("output");

        let error = make_fastq_from_assays(
            vec![assay],
            Modality::ATAC,
            output.clone(),
            false,
            || Ok(()),
        )
        .unwrap_err();

        assert!(error
            .to_string()
            .contains("make_fastq does not support long-read assays"));
        assert!(!output.exists());
    }
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
