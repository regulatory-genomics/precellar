use anyhow::{bail, Context, Result};
use std::{collections::HashMap, fs::File, path::{Path, PathBuf}};
use yaml_rust::{Yaml, YamlLoader};
use noodles::fastq;

pub type Modality = String;

#[derive(Debug)]
pub struct SeqSpec {
    pub version: String,
    pub id: String,
    pub name: String,
    pub doi: String,
    pub date: String,
    pub description: String,
    pub lib_struct: String,
    pub library_protocol: Option<String>,
    pub library_kit: Option<String>,
    pub sequence_protocol: String,
    pub sequence_kit: Option<String>,
    pub library_spec: HashMap<Modality, Region>,
    pub sequence_spec: HashMap<Modality, Read>,
}

impl TryFrom<Yaml> for SeqSpec {
    type Error = anyhow::Error;

    fn try_from(yaml: Yaml) -> Result<Self> {
        let version = yaml["seqspec_version"].as_str().unwrap().to_string();
        let id = yaml["assay_id"].as_str().unwrap().to_string();
        let name = yaml["name"].as_str().unwrap().to_string();
        let doi = yaml["doi"].as_str().unwrap().to_string();
        let date = yaml["date"].as_str().unwrap().to_string();
        let description = yaml["description"].as_str().unwrap().to_string();
        let lib_struct = yaml["lib_struct"].as_str().unwrap().to_string();
        let library_protocol = yaml["library_protocol"].as_str().map(|x| x.to_string());
        let library_kit = yaml["library_kit"].as_str().map(|x| x.to_string());
        let sequence_protocol = yaml["sequence_protocol"].as_str().unwrap().to_string();
        let sequence_kit = yaml["sequence_kit"].as_str().map(|x| x.to_string());
        let library_spec = yaml["library_spec"].clone().into_iter().map(|region|
            (region["region_type"].as_str().unwrap().to_string(), Region::try_from(region.clone()).unwrap())
        ).collect();
        let sequence_spec = yaml["sequence_spec"].clone().into_iter().map(|read|
            (read["modality"].as_str().unwrap().to_string(), Read::try_from(read.clone()).unwrap())
        ).collect();
        Ok(Self {
            version, id, name, doi, date, description, lib_struct, library_protocol,
            library_kit, sequence_protocol, sequence_kit, library_spec, sequence_spec,
        })
    }
}

impl SeqSpec {
    pub fn from_path<P: AsRef<Path>>(path: P) -> Result<Self> {
        let text = std::fs::read_to_string(&path)?;
        YamlLoader::load_from_str(&text)?
            .into_iter().next().with_context(|| format!("Empty YAML file: {:?}", path.as_ref()))?
            .try_into()
    }

    pub fn modality(&self, name: &str) -> Option<&Region> {
        self.library_spec.get(name)
    }
}

#[derive(Debug)]
pub struct Read {
    pub id: String,
    pub name: String,
    pub modality: Modality,
    pub primer_id: String,  // the primer_id should maps to the correct region id off
                        // of which the sequencing read is generated in the library_spec.
    pub length: Length,
    pub strand: bool,
}

impl TryFrom<Yaml> for Read {
    type Error = anyhow::Error;

    fn try_from(yaml: Yaml) -> Result<Self> {
        let id = yaml["read_id"].as_str().unwrap().to_string();
        let name = yaml["name"].as_str().unwrap().to_string();
        let modality = yaml["modality"].as_str().unwrap().to_string();
        let primer_id = yaml["primer_id"].as_str().unwrap().to_string();
        let min_length = yaml["min_len"].as_i64().unwrap() as usize;
        let max_length = yaml["max_len"].as_i64().unwrap() as usize;
        let length = if min_length == max_length {
            Length::Fixed(min_length)
        } else {
            Length::Range(min_length, max_length)
        };
        let strand = match yaml["strand"].as_str().unwrap() {
            "pos" => true,
            "neg" => false,
            x => bail!("Invalid strand: {:?}", x),
        };
        Ok(Self { id, name, modality, primer_id, length, strand })
    }
}

#[derive(Debug, Clone)]
pub enum SequenceType {
    Fixed(String),
    OnList(Vec<String>),
    Random,
    Joined,
}

#[derive(Debug, PartialEq, Eq, Clone)]
pub enum RegionType {
    Fastq,
    Barcode,
    UMI,
    CDNA,
    GDNA,
    Other(String),
}

impl From<&str> for RegionType {
    fn from(s: &str) -> Self {
        match s.trim().to_lowercase().as_str() {
            "fastq" => Self::Fastq,
            "barcode" => Self::Barcode,
            "umi" => Self::UMI,
            "cdna" => Self::CDNA,
            "gdna" => Self::GDNA,
            _ => Self::Other(s.to_string()),
        }
    }
}

#[derive(Debug, Copy, Clone)]
pub enum Length {
    Fixed(usize),
    Range(usize, usize),
}

#[derive(Debug, Clone)]
pub struct Region {
    id: String,
    region_type: RegionType,
    name: String,
    sequence_type: SequenceType,
    length: Length,
    sub_regions: Vec<Region>,
}

impl TryFrom<Yaml> for Region {
    type Error = anyhow::Error;

    fn try_from(yaml: Yaml) -> Result<Self> {
        let id = yaml["region_id"].as_str().unwrap().to_string();
        let region_type = yaml["region_type"].as_str().unwrap().into();
        let name = yaml["name"].as_str().unwrap().to_string();
        let sequence_type = match yaml["sequence_type"].as_str().unwrap() {
            "fixed" => SequenceType::Fixed(yaml["sequence"].as_str().unwrap().to_string()),
            "onlist" => SequenceType::OnList(Vec::new()),
            "random" => SequenceType::Random,
            "joined" => SequenceType::Joined,
            _ => bail!("Invalid sequence type: {:?}", yaml["sequence_type"]),
        };
        let min_length = yaml["min_len"].as_i64().unwrap() as usize;
        let max_length = yaml["max_len"].as_i64().unwrap() as usize;
        let length = if min_length == max_length {
            Length::Fixed(min_length)
        } else {
            Length::Range(min_length, max_length)
        };
        let sub_regions = yaml["regions"].clone().into_iter().map(Region::try_from).collect::<Result<Vec<_>>>()?;
        Ok(Self { id, region_type, name, sequence_type, length, sub_regions })
    }
}

#[derive(Debug, Clone, Copy)]
pub struct Range {
    pub start: usize,
    pub end: Option<usize>,
}

impl Region {
    /// Return the list of fastq regions in the region tree.
    pub fn fastqs(&self) -> Vec<&Region> {
        let mut result = Vec::new();
        for x in &self.sub_regions {
            if x.region_type == RegionType::Fastq {
                result.push(x);
            } else {
                let remains = x.fastqs();
                result.extend(remains);
            }
        }
        result
    }

    /// Return the start and end position of each subregion in a fastq region.
    pub fn subregion_range(&self) -> impl Iterator<Item = (RegionType, Range)> + '_ {
        if self.region_type != RegionType::Fastq {
            panic!("must be called on a fastq region");
        }

        let regions = &self.sub_regions;
        let len = regions.len();
        let mut current_pos = 0;
        let mut i = 0;

        std::iter::from_fn(move || {
            if i >= len {
                return None;
            }

            let region = &regions[i];
            let range = match region.length {
                Length::Fixed(n) => {
                    let r = Range {start: current_pos, end: Some(current_pos + n)};
                    current_pos += n;
                    r
                }
                Length::Range(_, _) => {
                    if i != len - 1 {
                        panic!("Variable length region must be the last region in a region list");
                    } else {
                        Range {start: current_pos, end: None}
                    }
                }
            };
            i += 1;
            Some((region.region_type.clone(), range))
        })
    }
}


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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_seqspec() {
        let spec = SeqSpec::from_path("tests/data/spec.yaml").unwrap();
        for fq in spec.modality("rna").unwrap().fastqs() {
            println!("{:?}", fq.subregion_range().collect::<Vec<_>>());
        }
    }
}
