mod segment;
use crate::region::Region;
use crate::Modality;
pub use segment::{Segment, SegmentInfo, SegmentInfoElem, SplitError};

use anyhow::Result;
use cached_path::Cache;
use indexmap::IndexMap;
use noodles::fastq;
use serde::{Deserialize, Serialize, Serializer};
use std::io::{BufRead, BufReader};
use std::ops::{Deref, DerefMut};
use std::path::Path;

/// Specification of a sequencing library.
#[derive(Debug, Clone, PartialEq)]
pub struct SeqSpec(IndexMap<String, Read>);

impl Deref for SeqSpec {
    type Target = IndexMap<String, Read>;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl DerefMut for SeqSpec {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

impl Serialize for SeqSpec {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        self.0.values().collect::<Vec<_>>().serialize(serializer)
    }
}

impl<'de> Deserialize<'de> for SeqSpec {
    fn deserialize<D: serde::Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        let reads = Vec::<Read>::deserialize(deserializer)?;
        let mut sequences = IndexMap::new();
        for read in reads {
            let read_id = read.read_id.clone();
            if sequences.insert(read_id.clone(), read).is_some() {
                return Err(serde::de::Error::custom(format!(
                    "Duplicate read id: {}",
                    &read_id
                )));
            }
        }
        Ok(Self(sequences))
    }
}

#[derive(Deserialize, Serialize, Debug, Clone, PartialEq)]
pub struct Read {
    pub read_id: String,
    pub name: Option<String>,
    pub modality: Modality,
    pub primer_id: String,
    pub min_len: u32,
    pub max_len: u32,
    pub strand: Strand,
    pub files: Option<Vec<File>>,
}

impl Default for Read {
    fn default() -> Self {
        Self {
            read_id: "".to_string(),
            name: None,
            modality: Modality::DNA,
            primer_id: "".to_string(),
            min_len: 0,
            max_len: 0,
            strand: Strand::Pos,
            files: None,
        }
    }
}

pub struct FastqReader {
    reader: fastq::Reader<Box<dyn BufRead>>,
    min_len: u32,
    max_len: u32,
}

impl FastqReader {
    pub fn new(reader: Box<dyn BufRead>, min_len: u32, max_len: u32) -> Self {
        Self {
            reader: fastq::Reader::new(reader),
            min_len,
            max_len,
        }
    }

    pub fn read_record(&mut self, record: &mut fastq::Record) -> Result<usize> {
        let n = self.reader.read_record(record)?;
        if self.max_len > 0 {
            record.quality_scores_mut().truncate(self.max_len as usize);
            record.sequence_mut().truncate(self.max_len as usize);
        }
        Ok(n)
    }

    pub fn records(&mut self) -> FastqRecords {
        FastqRecords { inner: self }
    }
}

pub struct FastqRecords<'a> {
    inner: &'a mut FastqReader,
}

impl<'a> Iterator for FastqRecords<'a> {
    type Item = Result<fastq::Record>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut buf = fastq::Record::default();

        match self.inner.read_record(&mut buf) {
            Ok(0) => None,
            Ok(_) => Some(Ok(buf)),
            Err(e) => Some(Err(e)),
        }
    }
}

impl Read {
    /// Open the fastq files for reading, and return a fastq reader.
    /// If the read has multiple fastq files, they will be concatenated.
    /// If the read has no fastq files, return None.
    pub fn open(&self) -> Option<FastqReader> {
        let files = self
            .files
            .clone()
            .unwrap_or_default()
            .into_iter()
            .filter(|file| file.filetype == "fastq")
            .collect::<Vec<_>>();
        if files.is_empty() {
            return None;
        }
        let reader =
            multi_reader::MultiReader::new(files.into_iter().map(move |file| file.open().unwrap()));
        Some(FastqReader::new(
            Box::new(BufReader::new(reader)),
            self.min_len,
            self.max_len,
        ))
    }

    /// Get the actual length of the read by reading the first record from the fastq file.
    pub fn actual_len(&self) -> Result<usize> {
        let mut reader = self.open().expect("No fastq files found.").reader;
        let mut record = fastq::Record::default();
        reader.read_record(&mut record)?;
        Ok(record.sequence().len())
    }

    /// Check if the read is reverse.
    pub fn is_reverse(&self) -> bool {
        match self.strand {
            Strand::Neg => true,
            Strand::Pos => false,
        }
    }

    pub(crate) fn get_segments<'a>(&'a self, region: &'a Region) -> Option<SegmentInfo> {
        if !region.sequence_type.is_joined() {
            return None;
        }

        let subregions = region.subregions.iter();
        let subregions: Box<dyn Iterator<Item = _>> = if self.is_reverse() {
            Box::new(subregions.rev())
        } else {
            Box::new(subregions)
        };

        let segment_info: SegmentInfo = subregions
            .skip_while(|region| {
                let region = region.read().unwrap();
                let found =
                    region.region_type.is_sequencing_primer() && region.region_id == self.primer_id;
                !found
            })
            .skip(1)
            .map(|x| x.read().unwrap().deref().clone())
            .collect();
        if segment_info.len() > 0 {
            Some(segment_info)
        } else {
            None
        }
    }
}

#[derive(Deserialize, Serialize, Debug, Copy, Clone, PartialEq)]
#[serde(rename_all = "lowercase")]
pub enum Strand {
    Pos,
    Neg,
}

#[derive(Deserialize, Serialize, Debug, Clone, PartialEq)]
pub struct File {
    pub file_id: String,
    pub filename: String,
    pub filetype: String,
    pub filesize: u64,
    pub url: String,
    pub urltype: UrlType,
    pub md5: String,
}

impl File {
    pub fn from_fastq<P: AsRef<Path>>(path: P, compute_md5: bool) -> Result<Self> {
        let file = std::fs::File::open(&path)?;
        let filename = path.as_ref().file_name().unwrap().to_str().unwrap();
        let md5 = if compute_md5 {
            crate::utils::md5sum(&path)?
        } else {
            "0".to_string()
        };
        Ok(File {
            file_id: filename.to_string(),
            filename: filename.to_string(),
            filetype: "fastq".to_string(),
            filesize: file.metadata()?.len(),
            url: path.as_ref().to_str().unwrap().to_string(),
            urltype: UrlType::Local,
            md5,
        })
    }

    pub(crate) fn normalize_path<P: AsRef<Path>>(&mut self, work_dir: P) -> Result<()> {
        if self.urltype == UrlType::Local {
            self.url = crate::utils::normalize_path(work_dir, &mut self.url)?
                .to_str()
                .unwrap()
                .to_owned();
        }
        Ok(())
    }

    pub(crate) fn unnormalize_path<P: AsRef<Path>>(&mut self, work_dir: P) -> Result<()> {
        if self.urltype == UrlType::Local {
            self.url = crate::utils::unnormalize_path(work_dir, &mut self.url)?
                .to_str()
                .unwrap()
                .to_owned();
        }
        Ok(())
    }

    /// Open the file for reading.
    /// If the file is remote, it will be downloaded to the cache directory.
    /// If the file is local, it will be opened directly.
    /// The base_dir is used to resolve relative paths.
    pub fn open(&self) -> Result<Box<dyn std::io::Read>> {
        match self.urltype {
            UrlType::Local => Ok(Box::new(crate::utils::open_file(&self.url)?)),
            _ => {
                let mut cache = Cache::new().unwrap();
                cache.dir = home::home_dir().unwrap().join(".cache/seqspec");
                let file = cache.cached_path(&self.url).unwrap();
                Ok(Box::new(crate::utils::open_file(file)?))
            }
        }
    }
}

#[derive(Deserialize, Serialize, Debug, Copy, Clone, PartialEq)]
#[serde(rename_all = "lowercase")]
pub enum UrlType {
    Local,
    Ftp,
    Http,
    Https,
}
