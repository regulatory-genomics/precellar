use std::{fs::File, io::{BufWriter, Write}, path::Path, str::FromStr};
use anyhow::{Context, Result};

/// Open a file, possibly compressed. Supports gzip and zstd.
pub fn open_file_for_read<P: AsRef<Path>>(file: P) -> Box<dyn std::io::Read> {
    if is_gzipped(file.as_ref()) {
        Box::new(flate2::read::MultiGzDecoder::new(File::open(file.as_ref()).unwrap()))
    } else {
        Box::new(File::open(file.as_ref()).unwrap())
    }
}

/// Determine the file compression type. Supports gzip and zstd.
fn is_gzipped<P: AsRef<Path>>(file: P) -> bool {
    flate2::read::MultiGzDecoder::new(File::open(file.as_ref()).unwrap()).header().is_some()
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

pub fn open_file_for_write<P: AsRef<Path>>(
    filename: P,
    compression: Option<Compression>,
    compression_level: Option<u32>,
) -> Result<Box<dyn Write + Send>> {
    let buffer = BufWriter::new(
        File::create(&filename).with_context(|| format!("cannot create file: {}", filename.as_ref().display()))?
    );
    let writer: Box<dyn Write + Send> = match compression {
        None => Box::new(buffer),
        Some(Compression::Gzip) => Box::new(flate2::write::GzEncoder::new(buffer, flate2::Compression::new(compression_level.unwrap_or(6)))),
        Some(Compression::Zstd) => {
            let mut zstd = zstd::stream::Encoder::new(buffer, compression_level.unwrap_or(3) as i32)?;
            zstd.multithread(8)?;
            Box::new(zstd.auto_finish())
        },
    };
    Ok(writer)
}