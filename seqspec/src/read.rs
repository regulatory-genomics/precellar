use crate::Modality;
use crate::region::{Region, SequenceType};

use cached_path::Cache;
use log::{debug, warn};
use noodles::fastq;
use indexmap::IndexMap;
use serde::{Deserialize, Serialize, Serializer};
use std::ops::{Deref, DerefMut};
use std::{io::{BufRead, BufReader}, ops::Range, path::PathBuf};
use anyhow::Result;
use std::path::Path;

#[derive(Debug, Clone, PartialEq)]
pub struct Sequences(IndexMap<String, Read>);

impl Deref for Sequences {
    type Target = IndexMap<String, Read>;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl DerefMut for Sequences {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

impl Serialize for Sequences {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        self.0.values().collect::<Vec<_>>().serialize(serializer)
    }
}

impl<'de> Deserialize<'de> for Sequences {
    fn deserialize<D: serde::Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        let reads = Vec::<Read>::deserialize(deserializer)?;
        let mut sequences = IndexMap::new();
        for read in reads {
            let read_id = read.read_id.clone();
            if sequences.insert(read_id.clone(), read).is_some() {
                return Err(serde::de::Error::custom(format!("Duplicate read id: {}", &read_id)));
            }
        }
        Ok(Sequences(sequences))
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

    pub(crate) fn get_index<'a>(&'a self, region: &'a Region) -> Option<Vec<(&'a Region, Range<u32>)>> {
        if region.sequence_type != SequenceType::Joined {
            return None;
        }

        let mut found_primer = false;

        let result = if self.is_reverse() {
            self.get_read_span(
                region.regions.as_ref().unwrap().iter().rev()
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
                region.regions.as_ref().unwrap().iter()
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
    fn get_read_span<'a, I: Iterator<Item = &'a Region>>(&self, regions: I) -> Vec<(&'a Region, Range<u32>)> {
        let mut result = Vec::new();
        let read_len = self.max_len;
        let mut cur_pos = 0;
        for region in regions {
            if region.min_len == region.max_len {
                let end = (cur_pos + region.min_len).min(read_len);
                result.push((region, cur_pos..end));
                if end == read_len {
                    break;
                }
                cur_pos = end;
            } else if cur_pos + region.min_len >= read_len {
                result.push((region, cur_pos..read_len));
                break;
            } else if cur_pos + region.max_len < read_len {
                warn!("Read ({}) length exceeds maximum length of the variable-length region (insertion), \
                    truncating the reads to the maximum length of the region. \
                    If this is not the desired behavior, please adjust the region lengths.", self.read_id);
                result.push((region, cur_pos..cur_pos + region.max_len));
                break;
            } else {
                debug!("Reads ({}) may contain additional bases downstream of the variable-length region, e.g., adapter sequences.", self.read_id);
                result.push((region, cur_pos..read_len));
                break;
            }
        }
        result
    }

    pub fn is_reverse(&self) -> bool {
        match self.strand {
            Strand::Neg => true,
            Strand::Pos => false,
        }
    }
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