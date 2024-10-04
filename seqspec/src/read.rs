use crate::{Modality, RegionType};
use crate::region::{Region, SequenceType};

use cached_path::Cache;
use noodles::fastq;
use indexmap::IndexMap;
use serde::{Deserialize, Serialize, Serializer};
use std::ops::{Deref, DerefMut};
use std::sync::Arc;
use std::{io::{BufRead, BufReader}, ops::Range, path::PathBuf};
use anyhow::Result;
use std::path::Path;

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

impl Serialize for SeqSpec{
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
                return Err(serde::de::Error::custom(format!("Duplicate read id: {}", &read_id)));
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

impl Read {
    pub fn open<P: AsRef<Path>>(&self, base_dir: P) -> Option<fastq::Reader<impl BufRead>> {
        let files = self.files.clone().unwrap_or(Vec::new())
            .into_iter().filter(|file| file.filetype == "fastq").collect::<Vec<_>>();
        if files.is_empty() {
            return None;
        }
        let base_dir = base_dir.as_ref().to_path_buf();
        let reader = multi_reader::MultiReader::new(
            files.into_iter().map(move |file| file.open(&base_dir))
        );
        Some(fastq::Reader::new(BufReader::new(reader)))
    }

    pub fn actual_len<P: AsRef<Path>>(&self, base_dir: P) -> Result<usize> {
        let mut reader = self.open(base_dir).unwrap();
        let mut record = fastq::Record::default();
        reader.read_record(&mut record)?;
        Ok(record.sequence().len())
    }

    pub(crate) fn get_index<'a>(&'a self, region: &'a Region) -> Option<RegionIndex> {
        if region.sequence_type != SequenceType::Joined {
            return None;
        }

        let mut found_primer = false;

        let result = if self.is_reverse() {
            self.get_read_span(
                region.subregions.iter().rev()
                    .skip_while(|region| {
                        let found = region.region_type.is_sequencing_primer() && region.region_id == self.primer_id;
                        if found {
                            found_primer = true;
                        }
                        !found
                    }).skip(1)
            )
        } else {
            self.get_read_span(
                region.subregions.iter()
                    .skip_while(|region| {
                        let found = region.region_type.is_sequencing_primer() && region.region_id == self.primer_id;
                        if found {
                            found_primer = true;
                        }
                        !found
                    }).skip(1)
            )
        };
        
        if found_primer {
            Some(result)
        } else {
            None
        }
    }

    /// Get the regions of the read.
    fn get_read_span<'a, I>(&self, mut regions: I) -> RegionIndex
    where
        I: Iterator<Item = &'a Arc<Region>>,
    {
        let mut index = Vec::new();
        let read_len = self.max_len;
        let mut cur_pos = 0;
        let mut readlen_info = ReadSpan::Covered;
        while let Some(region) = regions.next() {
            let region_id = region.region_id.clone();
            let region_type = region.region_type;
            if region.is_fixed_length() {  // Fixed-length region
                let end = cur_pos + region.min_len;
                if end >= read_len {
                    index.push((region_id, region_type, cur_pos..read_len));
                    if end > read_len {
                        readlen_info = ReadSpan::NotEnough;
                    }
                    break;
                } else {
                    index.push((region_id, region_type, cur_pos..end));
                    cur_pos = end;
                }
            } else if cur_pos + region.min_len >= read_len {  // Variable-length region and read is shorter
                index.push((region_id, region_type, cur_pos..read_len));
                readlen_info = ReadSpan::Span((read_len - cur_pos) as usize);
                break;
            } else if cur_pos + region.max_len < read_len {  // Variable-length region and read is longer than max length
                index.push((region_id, region_type, cur_pos..cur_pos + region.max_len));
                if let Some(next_region) = regions.next() {
                    readlen_info = ReadSpan::ReadThrough(next_region.region_id.clone());
                }
                break;
            } else {  // Variable-length region and read is within the length range
                index.push((region_id, region_type, cur_pos..read_len));
                if let Some(next_region) = regions.next() {
                    readlen_info = ReadSpan::MayReadThrough(next_region.region_id.clone());
                }
                break;
            }
        }
        RegionIndex { index, readlen_info }
    }

    pub fn is_reverse(&self) -> bool {
        match self.strand {
            Strand::Neg => true,
            Strand::Pos => false,
        }
    }
}

#[derive(Debug, Clone)]
pub struct RegionIndex {
    pub index: Vec<(String, RegionType, Range<u32>)>,
    pub readlen_info: ReadSpan,
}

#[derive(Debug, Clone)]
pub enum ReadSpan {
    Covered,  // The read is fully contained within the target region
    Span(usize),  // The read spans the target region
    NotEnough,  // Read is too short to reach the target region
    ReadThrough(String),  // Read is longer than the target region
    MayReadThrough(String),  // Read may be longer than the target region
}

#[derive(Deserialize, Serialize, Debug, Clone, PartialEq)]
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
    pub fn from_fastq<P: AsRef<Path>>(path: P) -> Result<Self> {
        let file = std::fs::File::open(&path)?;
        let filename = path.as_ref().file_name().unwrap().to_str().unwrap();
        Ok(File {
            file_id: filename.to_string(),
            filename: filename.to_string(),
            filetype: "fastq".to_string(),
            filesize: file.metadata()?.len(),
            url: path.as_ref().to_str().unwrap().to_string(),
            urltype: UrlType::Local,
            md5: "0".to_string(),
        })
    }

    /// Open the file for reading.
    /// If the file is remote, it will be downloaded to the cache directory.
    /// If the file is local, it will be opened directly.
    /// The base_dir is used to resolve relative paths.
    pub fn open<P: AsRef<Path>>(&self, base_dir: P) -> Box<dyn std::io::Read> {
        match self.urltype {
            UrlType::Local => {
                let mut path = PathBuf::from(&self.url);
                path = if path.is_absolute() {
                    path
                } else {
                    base_dir.as_ref().join(path)
                };
                Box::new(crate::utils::open_file_for_read(path))
            }
            _ => {
                let cache = Cache::new().unwrap();
                let file = cache.cached_path(&self.url).unwrap();
                Box::new(crate::utils::open_file_for_read(file))
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