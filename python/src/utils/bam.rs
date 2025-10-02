use anyhow::Result;
use itertools::Itertools;
use noodles::bam;
use noodles::fastq::{self, io::Writer};
use noodles::sam::alignment::record::data::field::{Tag, Value};
use precellar::utils::rev_compl_fastq_record;
use pyo3::prelude::*;
use rayon::slice::ParallelSliceMut;
use seqspec::utils::{create_file, Compression};
use std::{io::BufWriter, path::PathBuf};

/// Convert a BAM file to a FASTQ file.
///
/// This function reads a BAM file and convert the alignments back to a FASTQ file.
/// The function will ignore secondary/supplementary alignments.
///
/// Parameters
/// ----------
/// input: str
///    File path or url to the input BAM file.
/// output_dir: Path
///    Directory containing the output FASTQ files.
/// compression: Literal['gzip', 'zst'] | None
///    Compression algorithm to use. If None, the compression algorithm will be inferred from the file extension.
/// compression_level: int | None
///    Compression level to use.
#[pyfunction]
#[pyo3(
    signature = (input, output_dir, *, compression_level=None),
    text_signature = "(input, output_dir, *, compression_level=None)",
)]
pub fn bam_to_fastq(
    py: Python<'_>,
    input: &str,
    output_dir: PathBuf,
    compression_level: Option<u32>,
) -> Result<()> {
    fn bam_to_fq(bam: &bam::Record) -> fastq::Record {
        let name = bam.name().unwrap();
        let seq: Vec<u8> = bam.sequence().iter().collect();
        let qual: Vec<u8> = bam
            .quality_scores()
            .iter()
            .map(|x| x + 33)
            .collect();
        let fq = fastq::Record::new(fastq::record::Definition::new(name, ""), seq, qual);
        if bam.flags().is_reverse_complemented() {
            rev_compl_fastq_record(fq)
        } else {
            fq
        }
    }

    fn get_data(bam: &bam::Record, tag: &Tag) -> Option<Vec<u8>> {
        if let Some(val) = bam.data().get(tag) {
            match val.unwrap() {
                Value::String(s) => Some(s.to_vec()),
                _ => None,
            }
        } else {
            None
        }
    }

    fn extract_barcode(bam: &bam::Record) -> fastq::Record {
        let (sequence, qual) = if let Some(bc) = get_data(bam, &Tag::CELL_BARCODE_SEQUENCE) {
            let qual = get_data(bam, &Tag::CELL_BARCODE_QUALITY_SCORES).unwrap();
            (bc, qual)
        } else {
            let bc = get_data(bam, &Tag::CELL_BARCODE_ID).unwrap();
            let qual = vec![b'~'; bc.len()];
            (bc, qual)
        };

        let name = bam.name().unwrap();
        fastq::Record::new(fastq::record::Definition::new(name, ""), sequence, qual)
    }

    let read1_file = output_dir.join("R1.fq.zst");
    let read2_file = output_dir.join("R2.fq.zst");
    let index_file = output_dir.join("I1.fq.zst");

    let mut is_paired = false;
    if !output_dir.exists() {
        std::fs::create_dir_all(&output_dir)?;
    }
    let mut index_writer = Writer::new(BufWriter::new(create_file(
        index_file,
        Some(Compression::Zstd),
        compression_level,
        8,
    )?));
    let mut writer1 = Writer::new(BufWriter::new(create_file(
        read1_file,
        Some(Compression::Zstd),
        compression_level,
        8,
    )?));
    let mut writer2: Option<Writer<BufWriter<Box<dyn std::io::Write + Send>>>> = None;

    let mut bams = Vec::new();
    let reader = std::fs::File::open(input)?;
    let mut reader = bam::io::Reader::new(reader);
    reader.read_header()?;

    let mut records = reader.records();
    let mut n = 0u32;
    while let Some(record) = records.next() {
        let record = record?;
        if n % 100000 == 0 {
            py.check_signals()?;
        }

        let flag = record.flags();
        if flag.is_unmapped() || !(flag.is_secondary() || flag.is_supplementary()) {
            bams.push(record);
        }

        n += 1;
    }

    bams.par_sort_by(|a, b| a.name().cmp(&b.name()));
    bams.into_iter()
        .chunk_by(|x| x.name().unwrap().to_string())
        .into_iter()
        .for_each(|(_, group)| {
            let mut read1 = None;
            let mut read2 = None;
            group.into_iter().for_each(|r| {
                let flag = r.flags();
                if flag.is_first_segment() {
                    read1 = Some(r);
                } else if flag.is_last_segment() {
                    read2 = Some(r);
                }
            });
            if let Some(r1) = read1 {
                writer1.write_record(&bam_to_fq(&r1)).unwrap();
                index_writer.write_record(&extract_barcode(&r1)).unwrap();
            }
            if let Some(r2) = read2 {
                is_paired = true;

                if let Some(writer2) = &mut writer2 {
                    writer2.write_record(&bam_to_fq(&r2)).unwrap();
                } else {
                    writer2 = Some(Writer::new(BufWriter::new(
                        create_file(
                            read2_file.clone(),
                            Some(Compression::Zstd),
                            compression_level,
                            8,
                        )
                        .unwrap(),
                    )));
                    writer2
                        .as_mut()
                        .unwrap()
                        .write_record(&bam_to_fq(&r2))
                        .unwrap();
                }
            } else if is_paired {
                panic!("The input BAM file contains unpaired reads.");
            }
        });

    Ok(())
}