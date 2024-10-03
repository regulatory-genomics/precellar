use std::{fs::File, io::{BufWriter, Write}, path::{Path, PathBuf}, str::FromStr};
use anyhow::{Context, Result, anyhow};

/// Open a file, possibly compressed. Supports gzip and zstd.
pub fn open_file_for_read<P: AsRef<Path>>(file: P) -> Box<dyn std::io::Read> {
    match detect_compression(file.as_ref()) {
        Some(Compression::Gzip) => Box::new(flate2::read::MultiGzDecoder::new(File::open(file.as_ref()).unwrap())),
        Some(Compression::Zstd) => {
            let r = zstd::stream::read::Decoder::new(File::open(file.as_ref()).unwrap()).unwrap();
            Box::new(r)
        },
        None => Box::new(File::open(file.as_ref()).unwrap()),
    }
}

/// Determine the file compression type. Supports gzip and zstd.
fn detect_compression<P: AsRef<Path>>(file: P) -> Option<Compression> {
    if flate2::read::MultiGzDecoder::new(File::open(file.as_ref()).unwrap()).header().is_some() {
        Some(Compression::Gzip)
    } else if let Some(ext) = file.as_ref().extension() {
        if ext == "zst" {
            Some(Compression::Zstd)
        } else {
            None
        }
    } else {
        None
    }
}

#[derive(Debug, Clone, Copy)]
pub enum Compression {
    Gzip,
    Zstd,
}

impl FromStr for Compression {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "gzip" => Ok(Compression::Gzip),
            "zstd" | "zstandard" => Ok(Compression::Zstd),
            _ => Err(format!("unsupported compression: {}", s)),
        }
    }
}

impl TryFrom<&PathBuf> for Compression {
    type Error = anyhow::Error;

    fn try_from(path: &PathBuf) -> Result<Self> {
        let ext = path.extension().unwrap_or(std::ffi::OsStr::new(""));
        if ext == "gz" {
            Ok(Compression::Gzip)
        } else if ext == "zst" {
            Ok(Compression::Zstd)
        } else {
            Err(anyhow!("unsupported compression: {:?}", path))
        }
    }
}

pub fn open_file_for_write<P: AsRef<Path>>(
    filename: P,
    compression: Option<Compression>,
    compression_level: Option<u32>,
    num_threads: u32,
) -> Result<Box<dyn Write + Send>> {
    let buffer = BufWriter::new(
        File::create(&filename).with_context(|| format!("cannot create file: {}", filename.as_ref().display()))?
    );
    let writer: Box<dyn Write + Send> = match compression {
        None => Box::new(buffer),
        Some(Compression::Gzip) => Box::new(flate2::write::GzEncoder::new(buffer, flate2::Compression::new(compression_level.unwrap_or(6)))),
        Some(Compression::Zstd) => {
            let mut zstd = zstd::stream::Encoder::new(buffer, compression_level.unwrap_or(9) as i32)?;
            zstd.multithread(num_threads)?;
            Box::new(zstd.auto_finish())
        },
    };
    Ok(writer)
}