use anyhow::Result;
use futures::StreamExt;
use itertools::Itertools;
use noodles::bam;
use noodles::fastq::{self, io::Writer};
use noodles::sam::alignment::record::data::field::{Tag, Value};
use noodles::sam::alignment::record::QualityScores;
use precellar::utils::{rev_compl_fastq_record, strip_fastq};
use pyo3::{prelude::*, types::{PyDict, PyList}};
use rayon::slice::ParallelSliceMut;
use regex::Regex;
use seqspec::utils::{create_file, is_url, open_file_async, Compression};
use std::{io::BufWriter, path::PathBuf, str::FromStr};
use tokio::runtime::Runtime;
use std::iter::repeat;
use rand::Rng;
use glob::glob;

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
        compression=None, compression_level=None, input_compression=None, num_threads=8
    ),
    text_signature = "(in_fq, out_fq,
        *, regex, out_barcode, from_description=False,
        compression=None, compression_level=None, input_compression=None, num_threads=8)",
)]
fn strip_barcode_from_fastq(
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
                let dict = output.downcast::<PyDict>().unwrap();
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
        let mut reader = fastq::AsyncReader::new(tokio::io::BufReader::new(reader));

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
fn bam_to_fastq(
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
            .map(|x| x.unwrap() + 33)
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

/// Convert input path to a list of file paths, supporting glob patterns
fn expand_path_pattern(pattern: &String) -> PyResult<Vec<String>> {
    let paths: Vec<String> = glob(pattern)
        .map_err(|e| {
            pyo3::exceptions::PyValueError::new_err(format!(
                "Invalid glob pattern '{}': {}", 
                pattern, e
            ))
        })?
        .filter_map(Result::ok)
        .map(|path| path.to_string_lossy().into_owned())
        .collect();

    if paths.is_empty() {
        Err(pyo3::exceptions::PyValueError::new_err(format!(
            "No files found matching pattern: {}", 
            pattern
        )))
    } else {
        Ok(paths)
    }
}

/// Merge multiple FASTQ files into a single output file and generate a barcode file.
///
/// Parameters
/// ----------
/// input1_files: list[str]
///     List of paths to the first read (R1) FASTQ files.
/// input2_files: list[str] | None
///     Optional list of paths to the second read (R2) FASTQ files.
/// output1_file: str
///     Path to the output file for merged R1 reads.
/// output2_file: str | None
///     Optional path to the output file for merged R2 reads.
/// barcode_file: str
///     Path to the output barcode file.
/// compression: Literal['gzip', 'zst'] | None
///     Compression algorithm to use. If None, the compression algorithm will be inferred from the file extension.
/// compression_level: int | None
///     Compression level to use.
/// num_threads: int
///     Number of threads to use for compression.
/// verbose: bool
///     Controls the amount of output information.
#[pyfunction]
#[pyo3(
    signature = (
        input1_files,
        output1_file,
        barcode_file,
        input2_files=None,
        output2_file=None,
        *,
        compression=None,
        compression_level=None,
        num_threads=4,
        verbose=false
    ),
    text_signature = "(input1_files, output1_file, barcode_file, input2_files=None, output2_file=None, *, compression=None, compression_level=None, num_threads=4, verbose=False)",
)]
fn merge_fastq_files(
    py: Python<'_>,
    input1_files: PyObject,
    output1_file: String,
    barcode_file: String,
    input2_files: Option<PyObject>,
    output2_file: Option<String>,
    compression: Option<&str>,
    compression_level: Option<u32>,
    num_threads: u32,
    verbose: bool,
) -> PyResult<()> {
    // Convert input1_files to Vec<String> with glob expansion
    let input1_vec = if input1_files.downcast_bound::<PyList>(py).is_ok() {
        let patterns: Vec<String> = input1_files.extract(py)?;
        patterns
            .into_iter()
            .try_fold(Vec::new(), |mut acc, pattern| -> PyResult<_> {
                acc.extend(expand_path_pattern(&pattern)?);
                Ok(acc)
            })?
    } else {
        let pattern: String = input1_files.extract(py)?;
        expand_path_pattern(&pattern)?
    };

    // Sort the files to ensure consistent ordering
    let mut input1_vec = input1_vec;
    input1_vec.sort();

    // Convert input2_files to Vec<String> if provided, with glob expansion
    let input2_vec = if let Some(input2) = input2_files {
        if input2.downcast_bound::<PyList>(py).is_ok() {
            let patterns: Vec<String> = input2.extract(py)?;
            let mut files = patterns
                .into_iter()
                .try_fold(Vec::new(), |mut acc, pattern| -> PyResult<_> {
                    acc.extend(expand_path_pattern(&pattern)?);
                    Ok(acc)
                })?;
            files.sort();
            Some(files)
        } else {
            let pattern: String = input2.extract(py)?;
            let mut files = expand_path_pattern(&pattern)?;
            files.sort();
            Some(files)
        }
    } else {
        None
    };


    // Print the files that will be processed only if verbose
    if verbose {
        println!("Files to process:");
        println!("Read 1 files:");
        for file in &input1_vec {
            println!("  {}", file);
        }
        if let Some(ref files) = input2_vec {
            println!("Read 2 files:");
            for file in files {
                println!("  {}", file);
            }
        }
    }

    // Determine compression type
    let compression = compression
        .map(|x| Compression::from_str(x).unwrap())
        .or_else(|| detect_compression_from_path(&output1_file));

    // Create runtime for async operations
    let rt = Runtime::new()?;
    
    rt.block_on(async {
        // 1. Create output files
        let mut writer1 = create_file(
            &output1_file,
            compression,
            compression_level,
            num_threads,
        )?;

        let mut writer2 = if let Some(output2) = &output2_file {
            Some(create_file(
                output2,
                compression,
                compression_level,
                num_threads,
            )?)
        } else {
            None
        };

        // 2. Create barcode file
        let mut barcode_writer = create_file(
            &barcode_file,
            compression,
            compression_level,
            num_threads,
        )?;

        // 3. Process each pair of input files
        for (idx, input1_file) in input1_vec.iter().enumerate() {
            py.check_signals()?; // Check for Python interrupt signals

            if verbose {
                println!("Processing file pair {}:", idx + 1);
                println!("  Read 1: {}", input1_file);
            }
            
            // Generate a unique barcode for this input file pair
            let barcode = generate_random_barcode();
            
            // Process first input file
            let compression1 = detect_compression_from_path(input1_file);
            let file1 = open_file_async(input1_file, compression1).await?;
            let buf_reader1 = tokio::io::BufReader::new(file1);
            let mut reader1 = fastq::AsyncReader::new(buf_reader1);
            let mut records1 = reader1.records();

            // Process second input file if provided
            let mut reader2_opt = if let Some(input2_files) = &input2_vec {
                let input2_file = &input2_files[idx];
                if verbose {
                    println!("  Read 2: {}", input2_file);
                }
                
                let compression2 = detect_compression_from_path(input2_file);
                let file2 = open_file_async(input2_file, compression2).await?;
                let buf_reader2 = tokio::io::BufReader::new(file2);
                let reader2 = fastq::AsyncReader::new(buf_reader2);
                Some(reader2)
            } else {
                None
            };

            // Process records
            let mut count = 0;
            while let Some(result1) = records1.next().await {
                if count % 100000 == 0 {
                    py.check_signals()?; // Check for Python interrupt signals periodically
                }

                match result1 {
                    Ok(record1) => {
                        // Write record1 to output file
                        writeln!(writer1, "@{}", record1.name())?;
                        writeln!(writer1, "{}", String::from_utf8_lossy(record1.sequence()))?;
                        writeln!(writer1, "+")?;
                        writeln!(writer1, "{}", String::from_utf8_lossy(record1.quality_scores()))?;

                        // Process record2 if available
                        if let Some(ref mut reader2) = reader2_opt {
                            if let Some(result2) = reader2.records().next().await {
                                match result2 {
                                    Ok(record2) => {
                                        if let Some(ref mut writer2) = writer2 {
                                            writeln!(writer2, "@{}", record2.name())?;
                                            writeln!(writer2, "{}", String::from_utf8_lossy(record2.sequence()))?;
                                            writeln!(writer2, "+")?;
                                            writeln!(writer2, "{}", String::from_utf8_lossy(record2.quality_scores()))?;
                                        }
                                    }
                                    Err(e) => {
                                        return Err(pyo3::exceptions::PyIOError::new_err(
                                            format!("Error reading record from read 2 file: {}", e)
                                        ));
                                    }
                                }
                            } else {
                                return Err(pyo3::exceptions::PyValueError::new_err(
                                    "Read 2 file has fewer records than read 1 file"
                                ));
                            }
                        }
                        
                        // Write barcode
                        writeln!(barcode_writer, "@{}", record1.name())?;
                        writeln!(barcode_writer, "{}", barcode)?;
                        writeln!(barcode_writer, "+")?;
                        writeln!(barcode_writer, "{}", repeat("I").take(barcode.len()).collect::<String>())?;
                        
                        count += 1;
                        if verbose && count % 1_000_000 == 0 {
                            println!("  Processed {} record pairs", count);
                        }
                    }
                    Err(e) => {
                        return Err(pyo3::exceptions::PyIOError::new_err(
                            format!("Error reading record from read 1 file: {}", e)
                        ));
                    }
                }
            }

            if verbose {
                println!("Completed {} records from file pair {}", count, idx + 1);
            }
        }

        Ok(())
    })?;

    Ok(())
}

// Helper function to detect compression from file extension
fn detect_compression_from_path(path: &str) -> Option<Compression> {
    if path.ends_with(".gz") {
        Some(Compression::Gzip)
    } else if path.ends_with(".zst") {
        Some(Compression::Zstd)
    } else {
        None
    }
}

fn generate_random_barcode() -> String {
    let bases = ['A', 'T', 'C', 'G'];
    let mut rng = rand::thread_rng();
    let length = 8; // Set default length to 8 bp
    (0..length)
        .map(|_| bases[rng.gen_range(0..4)])
        .collect()
}

#[pymodule]
pub(crate) fn register_utils(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let utils = PyModule::new(parent_module.py(), "utils")?;

    utils.add_function(wrap_pyfunction!(strip_barcode_from_fastq, &utils)?)?;
    utils.add_function(wrap_pyfunction!(bam_to_fastq, &utils)?)?;
    utils.add_function(wrap_pyfunction!(merge_fastq_files, &utils)?)?;

    parent_module.add_submodule(&utils)
}
