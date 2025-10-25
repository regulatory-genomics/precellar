use anyhow::{anyhow, bail, Context, Result};
use async_compression::tokio::bufread::GzipDecoder;
use async_compression::tokio::bufread::ZstdDecoder;
use futures::StreamExt;
use md5::{Md5, Digest};
use std::io::Read;
use std::{
    fs::File,
    io::{BufWriter, Write},
    path::{Path, PathBuf},
    str::FromStr,
};
use tokio::io::AsyncRead;
use tokio_util::io::StreamReader;
use url::Url;

pub fn create_file<P: AsRef<Path>>(
    filename: P,
    compression: Option<Compression>,
    compression_level: Option<u32>,
    num_threads: u32,
) -> Result<Box<dyn Write + Send>> {
    let buffer = BufWriter::new(
        File::create(&filename)
            .with_context(|| format!("cannot create file: {}", filename.as_ref().display()))?,
    );
    let writer: Box<dyn Write + Send> = match compression {
        None => Box::new(buffer),
        Some(Compression::Gzip) => Box::new(flate2::write::GzEncoder::new(
            buffer,
            flate2::Compression::new(compression_level.unwrap_or(6)),
        )),
        Some(Compression::Zstd) => {
            let mut zstd =
                zstd::stream::Encoder::new(buffer, compression_level.unwrap_or(9) as i32)?;
            zstd.multithread(num_threads)?;
            Box::new(zstd.auto_finish())
        }
    };
    Ok(writer)
}

/// Open a file, possibly compressed. Supports gzip and zstd.
pub fn open_file<P: AsRef<Path>>(file: P) -> Result<Box<dyn std::io::Read + Send + Sync>> {
    let reader: Box<dyn std::io::Read + Send + Sync> = match detect_compression(file.as_ref()) {
        Some(Compression::Gzip) => Box::new(flate2::read::MultiGzDecoder::new(File::open(
            file.as_ref(),
        )?)),
        Some(Compression::Zstd) => {
            let r = zstd::stream::read::Decoder::new(File::open(file.as_ref())?)?;
            Box::new(r)
        }
        None => Box::new(File::open(file.as_ref())?),
    };
    anyhow::Ok(reader).with_context(|| format!("cannot open file: {:?}", file.as_ref()))
}

/// Fetch content from a URL or a file, possibly compressed. Supports gzip and zstd.
pub async fn open_file_async(
    src: &str,
    compression: Option<Compression>,
) -> Result<Box<dyn AsyncRead + Send + Unpin>> {
    if is_url(src) {
        return Ok(Box::new(open_url_async(Url::parse(src).unwrap()).await?));
    }

    let reader = tokio::io::BufReader::new(tokio::fs::File::open(src).await?);

    let compression = compression.or_else(|| detect_compression(src));
    let reader: Box<dyn AsyncRead + Send + Unpin> = match compression {
        Some(Compression::Gzip) => Box::new(GzipDecoder::new(reader)),
        Some(Compression::Zstd) => Box::new(ZstdDecoder::new(reader)),
        None => Box::new(reader),
    };
    Ok(reader)
}

pub fn is_url(src: &str) -> bool {
    if let Ok(url) = Url::parse(src) {
        if !url.cannot_be_a_base() {
            return true;
        }
    }
    return false;
}

async fn open_url_async(url: Url) -> Result<impl AsyncRead> {
    // Create a client with automatic redirect following
    let client = reqwest::Client::builder()
        .redirect(reqwest::redirect::Policy::limited(10)) // Follow up to 10 redirects
        .build()?;

    let response = client.get(url).send().await?;
    if !response.status().is_success() {
        bail!("Request failed with status: {}", response.status());
    }

    // Convert the response into a byte stream
    let byte_stream = response
        .bytes_stream()
        .map(|result| result.map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e)));
    Ok(StreamReader::new(byte_stream))
}

/// Determine the file compression type. Supports gzip and zstd.
fn detect_compression<P: AsRef<Path>>(path: P) -> Option<Compression> {
    let file = File::open(path.as_ref())
        .with_context(|| format!("cannot open file: {:?}", path.as_ref()))
        .unwrap();
    if flate2::read::MultiGzDecoder::new(file).header().is_some() {
        Some(Compression::Gzip)
    } else if let Some(ext) = path.as_ref().extension() {
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

impl TryFrom<PathBuf> for Compression {
    type Error = anyhow::Error;

    fn try_from(path: PathBuf) -> Result<Self> {
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

/// Return reverse complement of a DNA sequence.
pub fn rev_compl(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&x| match x {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            _ => x,
        })
        .collect()
}

pub fn normalize_path<P1: AsRef<Path>, P2: AsRef<Path>>(work_dir: P1, path: P2) -> Result<PathBuf> {
    let path = path.as_ref();
    let result = if path.is_absolute() {
        path.to_path_buf()
    } else {
        work_dir.as_ref().join(path)
    };
    result
        .canonicalize()
        .with_context(|| format!("Failed to normalize path: {:?}", &result))
}

pub fn unnormalize_path<P1: AsRef<Path>, P2: AsRef<Path>>(
    work_dir: P1,
    path: P2,
) -> Result<PathBuf> {
    let work_dir = work_dir
        .as_ref()
        .canonicalize()
        .with_context(|| format!("failed to make absolute path: {:?}", work_dir.as_ref()))?;
    let path = path
        .as_ref()
        .canonicalize()
        .with_context(|| format!("failed to make absolute path: {:?}", path.as_ref()))?;
    Ok(relative_path(work_dir.as_path(), path.as_path()))
}

fn relative_path(from: &Path, to: &Path) -> PathBuf {
    let mut from_iter = from.components();
    let mut to_iter = to.components();
    let mut result = PathBuf::new();

    let mut fr;
    let mut to;
    loop {
        fr = from_iter.next();
        to = to_iter.next();
        if !(fr == to && fr.is_some()) {
            break;
        }
    }

    if fr.is_some() {
        // Now, calculate the path to "go up" from `from` to the common ancestor
        result.push("..");
        for _ in from_iter {
            result.push("..");
        }
    }

    // Append the remaining part of `to`
    if let Some(to) = to {
        result.push(to);
        for t in to_iter {
            result.push(t);
        }
    }

    result
}

pub fn md5sum<P: AsRef<Path>>(path: P) -> Result<String> {
    // Open the file and create a buffered reader
    let file = File::open(path)?;
    let mut reader = std::io::BufReader::new(file);

    let mut hasher = Md5::new();

    // Buffer for reading chunks
    let mut buffer = [0u8; 64 * 65535];

    // Read the file in a loop, updating the hasher
    loop {
        let bytes_read = reader.read(&mut buffer)?;
        if bytes_read == 0 {
            break;
        }
        hasher.update(&buffer[..bytes_read]);
    }

    let digest = hasher.finalize();
    Ok(base16ct::lower::encode_string(&digest))
}

pub fn hamming_distance(seq1: &[u8], seq2: &[u8]) -> Result<usize> {
    if seq1.len() != seq2.len() {
        bail!(
            "sequences must be of the same length: {} != {}",
            seq1.len(),
            seq2.len()
        );
    }
    Ok(seq1
        .iter()
        .zip(seq2.iter())
        .filter(|(a, b)| a != b)
        .count())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_relative_path() {
        let from = Path::new("/a/b/c");
        let to = Path::new("/a/b/d/e");
        assert_eq!(relative_path(from, to), PathBuf::from("../d/e"));

        let from = Path::new("/a/b/c");
        let to = Path::new("/a/b/c/d");
        assert_eq!(relative_path(from, to), PathBuf::from("d"));

        let from = Path::new("/a/b/c");
        let to = Path::new("/a/b/c/d/e");
        assert_eq!(relative_path(from, to), PathBuf::from("d/e"));

        let from = Path::new("/a/b/c");
        let to = Path::new("/a/b");
        assert_eq!(relative_path(from, to), PathBuf::from(".."));

        let from = Path::new("/a/b/c");
        let to = Path::new("/a");
        assert_eq!(relative_path(from, to), PathBuf::from("../.."));

        let from = Path::new("/a/b/c");
        let to = Path::new("/a/b/c");
        assert_eq!(relative_path(from, to), PathBuf::from(""));
    }
}
