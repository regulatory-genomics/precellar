use anyhow::{bail, Result};
use futures::StreamExt;
use itertools::Itertools;
use noodles::fastq::{self, io::Writer};
use precellar::utils::strip_fastq;
use pyo3::{prelude::*, types::PyDict};
use regex::Regex;
use seqspec::utils::{create_file, is_url, open_file, open_file_async, Compression};
use std::{
    io::{BufReader, BufWriter},
    path::PathBuf,
    str::FromStr,
};
use tokio::runtime::Runtime;

/// Remove barcode from the read names of fastq records.
///
/// The old practice of storing barcodes in read names is not recommended. This
/// function extracts barcodes from read names using regular expressions and
/// writes them to a separate fastq file.
///
/// Parameters
/// ----------
/// in_fq: str
///     File path or URL to the input fastq file.
/// out_fq: Path
///     File path to the output fastq file.
/// regex: str
///     Extract barcodes from read names or descriptions from the input fastq file using regular expressions.
///     We use [capturing groups](https://en.wikipedia.org/wiki/Regular_expression#POSIX_basic_and_extended)
///     to extract and remove the barcodes. For more details, please
///     read the examples below. You can test your regex on this website: https://regex101.com/.
/// out_barcode: dict[int | str, Path] | Path | None
///     A dictionary that maps the index or name of the capturing group to the output fastq file path.
///     If a single file path is provided, the first capturing group will be written to that file.
///     If None, none of the barcodes will be written to a file.
/// from_description: bool
///    If True, the barcode will be extracted from the description instead of the name.
/// compression: Literal['gzip', 'zst'] | None
///     Compression algorithm to use. If None, the compression algorithm will be inferred from the file extension.
/// compression_level: int | None
///     Compression level to use.
/// input_compression: Liter['gzip', 'zst'] | None
///     Compression algorithm of the input fastq file. This has to be specified
///     if the input fastq file is compressed and fetched from a URL.
/// num_threads: int
///     The number of threads to use.
///
/// Examples
/// --------
/// Suppose the read names in the fastq file contain barcodes and UMIs in the following format:
/// `A01535:24:HW2MMDSX2:2:1359:8513:3458:bd:69:Y6:10:TGATAGGTTG`,
/// where `bd:69:Y6:10` corresponds to the barcode and `TGATAGGTTG` corresponds to the UMI.
/// We can extract the barcode using the following regular expression:
/// `"(..:..:..:..)(:)\\w+$"`. Here we use two capturing groups. The first group captures the barcode,
/// and the second group captures the colon that separates the barcode from the UMI.
/// Texts matched by the capturing groups are removed from the read names.
/// So the resultant read name will be `A01535:24:HW2MMDSX2:2:1359:8513:3458:TGATAGGTTG`.
/// Here we use `out_barcode="test_out.fq.zst"` to write the first capturing group to a separate fastq file.
///
/// >>> from precellar.utils import strip_barcode_from_fastq
/// >>> strip_barcode_from_fastq(
/// ...     "test.fq.gz",
/// ...     "test_out.fq.zst",
/// ...     regex="(..:..:..:..)(:)\\w+$",
/// ...     out_barcode="barcode.fq.zst",
/// ... )
///
/// It is possible to extract both the barcode and UMI at the same time, using the following regular expression:
/// `"(..:..:..:..)(:)([ACGTN]+)$"`. Here we use three capturing groups. The first group captures the barcode,
/// the second group captures the colon that separates the barcode from the UMI, and the third group captures the UMI.
/// To write the barcode and UMI to separate files, we can use a dictionary:
/// `out_barcode={0: "barcode.fq.zst", 2: "umi.fq.zst"}`. The integer keys correspond
/// to the index of the capturing groups.
///
/// >>> strip_barcode_from_fastq(
/// ...     "test.fq.gz",
/// ...     "test_out.fq.zst",
/// ...     regex="(..:..:..:..)(:)([ACGTN]+)$",
/// ...     out_barcode={0: "barcode.fq.zst", 2: "umi.fq.zst"},
/// ... )
///
/// We can use names to refer to the capturing groups instead of indices:
/// `"(?<barcode>..:..:..:..)(:)(?<umi>[ACGTN]+)$"`.
///
/// >>> strip_barcode_from_fastq(
/// ...     "test.fq.gz",
/// ...     "test_out.fq.zst",
/// ...     regex="(?<barcode>..:..:..:..)(:)(?<umi>[ACGTN]+)$",
/// ...     out_barcode={'barcode': "barcode.fq.zst", 'umi': "umi.fq.zst"},
/// ... )
///
/// If the barcode is stored in the description instead of the name, we can set `from_description=True`.
///
#[pyfunction]
#[pyo3(
    signature = (in_fq, out_fq,
        *, regex, out_barcode=None, from_description=false,
        compression=None, compression_level=None, input_compression=None, num_threads=16
    ),
    text_signature = "(in_fq, out_fq,
        *, regex, out_barcode, from_description=False,
        compression=None, compression_level=None, input_compression=None, num_threads=16)",
)]
pub fn strip_barcode_from_fastq(
    py: Python<'_>,
    in_fq: &str,
    out_fq: PathBuf,
    regex: &str,
    out_barcode: Option<Bound<'_, PyAny>>,
    from_description: bool,
    compression: Option<&str>,
    compression_level: Option<u32>,
    input_compression: Option<&str>,
    num_threads: u32,
) -> Result<()> {
    if is_url(in_fq) && input_compression.is_none() {
        log::warn!("The input source is URL and input_compression is None");
    }

    let re = Regex::new(regex)?;

    let mut fq_writer = {
        let compression = compression
            .map(|x| Compression::from_str(x).unwrap())
            .or((&out_fq).try_into().ok());
        Writer::new(BufWriter::new(create_file(
            out_fq,
            compression,
            compression_level,
            num_threads,
        )?))
    };

    let (barcode_keys, mut barcode_writers) = if let Some(output) = out_barcode {
        let mut keys = Vec::new();
        let mut files = Vec::new();
        match output.extract::<PathBuf>() {
            Ok(p) => {
                keys.push(0.into());
                files.push(p);
            }
            _ => {
                let dict = output.cast::<PyDict>().unwrap();
                dict.iter().for_each(|(k, v)| {
                    if let Ok(k) = k.extract::<usize>() {
                        keys.push(k.into());
                    } else {
                        keys.push(k.extract::<String>().unwrap().into());
                    }

                    files.push(v.extract::<PathBuf>().unwrap());
                });
            }
        };
        let writers = files
            .into_iter()
            .map(|file| {
                let compression = compression
                    .map(|x| Compression::from_str(x).unwrap())
                    .or((&file).try_into().ok());
                Writer::new(BufWriter::new(
                    create_file(file, compression, compression_level, num_threads).unwrap(),
                ))
            })
            .collect();
        (Some(keys), writers)
    } else {
        (None, Vec::new())
    };

    let rt = Runtime::new()?;
    rt.block_on(async {
        let reader = open_file_async(
            in_fq,
            input_compression.map(|x| Compression::from_str(x).unwrap()),
        )
        .await?;
        let mut reader = fastq::r#async::io::Reader::new(tokio::io::BufReader::new(reader));

        let mut i = 0u32;
        reader
            .records()
            .for_each(|record| {
                if i % 1000000 == 0 {
                    py.check_signals().unwrap();
                }
                i += 1;

                let record = record.unwrap();
                let (record, barcodes) = strip_fastq(
                    record,
                    &re,
                    barcode_keys.as_ref().map(|x| x.as_slice()),
                    from_description,
                )
                .unwrap();

                fq_writer.write_record(&record).unwrap();
                if let Some(barcodes) = barcodes {
                    barcodes.into_iter().enumerate().for_each(|(i, barcode)| {
                        let qual_score = vec![b'~'; barcode.len()];
                        let fq = fastq::Record::new(
                            fastq::record::Definition::new(record.name(), ""),
                            barcode.into_bytes(),
                            qual_score,
                        );
                        barcode_writers
                            .get_mut(i)
                            .unwrap()
                            .write_record(&fq)
                            .unwrap();
                    })
                }

                futures::future::ready(())
            })
            .await;

        Ok(())
    })
}

/// Remove barcode from the read names of fastq records.
///
/// The old practice of storing barcodes in read names is not recommended. This
/// function extracts barcodes from read names and
/// writes them to a separate fastq file.
///
/// Parameters
/// ----------
/// in_fq: str
///     File path or URL to the input fastq file.
/// out_fq: Path
///     File path to the output fastq file.
/// rename: Callable[str, str] | None
///     A function that takes in the original read name and returns the new read name.   
/// get_barcode: Callable[str, str] | None
///     A function that takes in the original read name and returns the barcode sequence.
/// barcode_fq: Path | None
///     File path to the output barcode fastq file.
/// compression: Literal['gzip', 'zst'] | None
///     Compression algorithm to use. If None, the compression algorithm will be inferred from the file extension.
/// compression_level: int | None
///     Compression level to use.
/// num_threads: int
///     The number of threads to use.
///
#[pyfunction]
#[pyo3(
    signature = (in_fq, *, out_fq=None, rename=None, get_barcode=None, barcode_fq=None,
        compression=None, compression_level=None, num_threads=16
    ),
    text_signature = "(in_fq, *, out_fq=None, rename=None, get_barcode=None, barcode_fq=None,
        compression=None, compression_level=None, num_threads=16)",
)]
pub fn extract_barcode_from_name(
    py: Python<'_>,
    in_fq: PathBuf,
    out_fq: Option<PathBuf>,
    rename: Option<Bound<'_, PyAny>>,
    get_barcode: Option<Bound<'_, PyAny>>,
    barcode_fq: Option<PathBuf>,
    compression: Option<&str>,
    compression_level: Option<u32>,
    num_threads: u32,
) -> Result<()> {
    let mut fq_writer = out_fq.map(|fq| {
        let compression = compression
            .map(|x| Compression::from_str(x).unwrap())
            .or((&fq).try_into().ok());
        Writer::new(BufWriter::new(create_file(
            fq,
            compression,
            compression_level,
            num_threads,
        ).unwrap()))
    });
    let mut barcode_writer = barcode_fq.map(|fq| {
        let compression = compression
            .map(|x| Compression::from_str(x).unwrap())
            .or((&fq).try_into().ok());
        Writer::new(BufWriter::new(create_file(
            fq,
            compression,
            compression_level,
            num_threads,
        ).unwrap()))
    });

    let reader = open_file(in_fq)?;
    let mut reader = fastq::io::Reader::new(BufReader::new(reader));
    let mut i = 0u32;
    reader.records().for_each(|record| {
        if i % 1000000 == 0 {
            py.check_signals().unwrap();
        }
        i += 1;

        let mut record = record.unwrap();
        let text = record.name().to_string() + " " + &record.definition().description().to_string();

        let new_name: String = if let Some(ref rename) = rename {
            rename
                .call1((&text,))
                .unwrap()
                .extract()
                .unwrap()
        } else {
            text.clone()
        };

        if let Some(ref mut fq_writer) = fq_writer {
            *record.definition_mut() = fastq::record::Definition::new(new_name.clone(), "");
            fq_writer.write_record(&record).unwrap();
        }

        if let Some(ref mut barcode_writer) = barcode_writer {
            let barcode: String = get_barcode.as_ref().unwrap()
                .call1((text,))
                .unwrap()
                .extract()
                .unwrap();
            let qual_score = vec![b'~'; barcode.len()];
            let fq = fastq::Record::new(
                fastq::record::Definition::new(new_name, ""),
                barcode.into_bytes(),
                qual_score,
            );
            barcode_writer.write_record(&fq).unwrap();
        }
    });

    Ok(())
}

/// Create multiplexed FASTQ files by merging multiple FASTQ files into single output files,
/// while generating barcodes for each read based on their file of origin.
///
/// This function is useful for methods such as SMART-seq.
///
/// Parameters
/// ----------
/// input_read1: list[str]
///     List of paths to the first read (R1) FASTQ files.
/// output_read1: str
///     Path to the output file for merged R1 reads.
/// output_barcode: str
///     Path to the output barcode file.
/// input_read2: list[str] | None
///     Optional list of paths to the second read (R2) FASTQ files.
/// output_read2: str | None
///     Optional path to the output file for merged R2 reads.
/// compression: Literal['gzip', 'zst'] | None
///     Compression algorithm to use. If None, the compression algorithm will be inferred from the file extension.
/// compression_level: int | None
///     Compression level to use.
/// num_threads: int
///     Number of threads to use for compression.
#[pyfunction]
#[pyo3(
    signature = (
        input_read1,
        output_read1,
        output_barcode,
        *,
        input_read2=None,
        output_read2=None,
        compression=None,
        compression_level=None,
        num_threads=16,
    ),
    text_signature = "(input_read1, output_read1, output_barcode, *, input_read2=None, output_read2=None, compression=None, compression_level=None, num_threads=16)",
)]
pub fn multiplex_fastq(
    py: Python<'_>,
    input_read1: Vec<String>,
    output_read1: PathBuf,
    output_barcode: PathBuf,
    input_read2: Option<Vec<String>>,
    output_read2: Option<PathBuf>,
    compression: Option<&str>,
    compression_level: Option<u32>,
    num_threads: u32,
) -> Result<Vec<String>> {
    let n_files = input_read1.len();
    if let Some(ref input2) = input_read2 {
        if input2.len() != n_files {
            bail!("The number of input read2 files must match the number of input read1 files");
        }
    }

    let compression = compression
        .map(|x| Compression::from_str(x).unwrap())
        .or((&output_read1).try_into().ok());

    let mut read1_writer = Writer::new(BufWriter::new(create_file(
        output_read1,
        compression,
        compression_level,
        num_threads,
    )?));
    let mut read2_writer = output_read2.as_ref().map(|output2| {
        Writer::new(BufWriter::new(
            create_file(output2, compression, compression_level, num_threads).unwrap(),
        ))
    });
    let mut barcode_writer = Writer::new(BufWriter::new(create_file(
        output_barcode,
        compression,
        compression_level,
        num_threads,
    )?));

    let barcodes = gen_unique_barcodes(n_files);
    Runtime::new()?.block_on(async {
        for (idx, (input1_file, barcode)) in
            input_read1.into_iter().zip_eq(barcodes.iter()).enumerate()
        {
            let mut count = 0;

            let reader1 = open_file_async(&input1_file, None).await?;
            let mut reader1 = fastq::r#async::io::Reader::new(tokio::io::BufReader::new(reader1));

            if let Some(ref input2) = input_read2 {
                let reader2 = open_file_async(&input2[idx], None).await?;
                let mut reader2 =
                    fastq::r#async::io::Reader::new(tokio::io::BufReader::new(reader2));
                reader1
                    .records()
                    .zip(reader2.records())
                    .for_each(|(record1, record2)| {
                        if count % 100000 == 0 {
                            py.check_signals().unwrap();
                        }

                        let record1 = record1.unwrap();
                        let record2 = record2.unwrap();
                        let barcode_record = fastq::Record::new(
                            record1.definition().clone(),
                            barcode.clone(),
                            vec![b'~'; barcode.len()],
                        );

                        assert_eq!(
                            record1.name(),
                            record2.name(),
                            "Read names do not match: {} != {}",
                            record1.name(),
                            record2.name()
                        );
                        read1_writer.write_record(&record1).unwrap();
                        read2_writer
                            .as_mut()
                            .unwrap()
                            .write_record(&record2)
                            .unwrap();
                        barcode_writer.write_record(&barcode_record).unwrap();

                        count += 1;
                        futures::future::ready(())
                    })
                    .await;
            } else {
                reader1
                    .records()
                    .for_each(|record| {
                        if count % 100000 == 0 {
                            py.check_signals().unwrap();
                        }

                        let record = record.unwrap();
                        let barcode_record = fastq::Record::new(
                            record.definition().clone(),
                            barcode.clone(),
                            vec![b'~'; barcode.len()],
                        );
                        read1_writer.write_record(&record).unwrap();
                        barcode_writer.write_record(&barcode_record).unwrap();

                        count += 1;
                        futures::future::ready(())
                    })
                    .await;
            }
        }

        anyhow::Ok(())
    })?;

    Ok(barcodes
        .into_iter()
        .map(|b| String::from_utf8_lossy(&b).to_string())
        .collect())
}

/// Generate n unique barcodes
fn gen_unique_barcodes(n: usize) -> Vec<Vec<u8>> {
    let bases = [b'A', b'C', b'G', b'T'];
    // Decide the length of the barcode needed
    let length = (n as f64).log(4.0).ceil() as usize;
    let mut barcodes = Vec::new();

    for mut i in 0..n {
        let mut kmer = vec![b'A'; length];
        let mut j = 1;
        while j <= length && i != 0 {
            kmer[length - j] = bases[i % 4];

            i = i / 4;
            j += 1;
        }
        barcodes.push(kmer);
    }
    barcodes
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gen_unique_barcodes() {
        let barcodes = gen_unique_barcodes(10);
        let expected = vec![
            b"AA".to_vec(),
            b"AC".to_vec(),
            b"AG".to_vec(),
            b"AT".to_vec(),
            b"CA".to_vec(),
            b"CC".to_vec(),
            b"CG".to_vec(),
            b"CT".to_vec(),
            b"GA".to_vec(),
            b"GC".to_vec(),
        ];
        assert_eq!(barcodes, expected);

        assert_eq!(gen_unique_barcodes(977).len(), 977);
    }
}
