use anyhow::{bail, Context, Result};
use noodles::fastq;
use std::{collections::HashMap, io::{BufRead, BufReader}, path::{Path, PathBuf}};
use yaml_rust::{Yaml, YamlLoader};
use cached_path::Cache;
use itertools::Itertools;

use crate::io::open_file_for_read;

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
    pub sequence_spec: HashMap<Modality, Vec<Read>>,
}

impl TryFrom<&Yaml> for SeqSpec {
    type Error = anyhow::Error;

    fn try_from(yaml: &Yaml) -> Result<Self> {
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
        let sequence_spec = yaml["sequence_spec"].clone().into_iter()
            .map(|read| (read["modality"].as_str().unwrap().to_string(), Read::try_from(&read).unwrap()))
            .sorted_by(|a, b| a.0.cmp(&b.0))
            .chunk_by(|x| x.0.clone())
            .into_iter()
            .map(|(k, g)| (k, g.map(|x| x.1).collect()))
            .collect();
        Ok(Self {
            version, id, name, doi, date, description, lib_struct, library_protocol,
            library_kit, sequence_protocol, sequence_kit, library_spec, sequence_spec,
        })
    }
}

impl TryFrom<Yaml> for SeqSpec {
    type Error = anyhow::Error;

    fn try_from(yaml: Yaml) -> Result<Self> {
        Self::try_from(&yaml)
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

    /// Return an iterator over all regions in the region tree.
    pub fn iter_regions(&self) -> impl Iterator<Item = &Region> {
        self.library_spec.values().flat_map(|x| {
            Box::new(std::iter::once(x).chain(x.iter_regions())) as Box<dyn Iterator<Item = &Region>>
        })
    }

    pub fn get_read_by_primer_id(&self, modality: &str, primer_id: &str) -> Option<&Read> {
        self.sequence_spec.get(modality).unwrap().iter().find(|x| x.primer_id == primer_id)
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Read {
    pub id: Vec<String>,   // The file paths of the sequencing reads from multiple lanes.
    pub name: String,
    pub modality: Modality,
    pub primer_id: String,  // the primer_id should maps to the correct region id off
                        // of which the sequencing read is generated in the library_spec.
    pub length: Length,
    strand: bool,
}

impl Read {
    pub fn read_fastq<P: AsRef<Path>>(&self, base_dir: P) -> fastq::Reader<impl BufRead> {
        let reader = multi_reader::MultiReader::new(
            self.id.clone().into_iter().map(move |x| {
                let mut path = PathBuf::from(x);
                path = if path.is_absolute() {
                    path
                } else {
                    base_dir.as_ref().join(path)
                };
                open_file_for_read(path)
            })
        );
        fastq::Reader::new(BufReader::new(reader))
    }

    pub fn is_reverse(&self) -> bool {
        !self.strand
    }
}

impl TryFrom<&Yaml> for Read {
    type Error = anyhow::Error;

    fn try_from(yaml: &Yaml) -> Result<Self> {
        let id = yaml["read_id"].as_str().unwrap().split(',').map(|x| x.trim().to_string()).collect();
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
    OnList {
        filepath: String,
        md5: Option<String>,
        local: bool,
    },
    Random,
    Joined,
}

impl SequenceType {
    pub fn as_str(&self) -> &'static str {
        match self {
            Self::Fixed(_) => "fixed",
            Self::OnList { .. } => "onlist",
            Self::Random => "random",
            Self::Joined => "joined",
        }
    }

    pub(crate) fn fetch_onlist(&self) -> Result<Vec<String>> {
        match self {
            Self::OnList { filepath, .. } => {
                let cache = Cache::new()?;
                let file = cache.cached_path(&filepath)?;
                let reader = std::io::BufReader::new(open_file_for_read(file));
                Ok(reader.lines().map(|x| x.unwrap()).collect())
            }
            _ => panic!("Not an onlist sequence type"),
        }
    }
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

#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub enum Length {
    Fixed(usize),
    Range(usize, usize),
}

#[derive(Debug, Clone)]
pub struct Region {
    pub id: String,
    pub region_type: RegionType,
    pub name: String,
    pub sequence_type: SequenceType,
    pub length: Length,
    pub sub_regions: Vec<Region>,
}

impl TryFrom<Yaml> for Region {
    type Error = anyhow::Error;

    fn try_from(yaml: Yaml) -> Result<Self> {
        let id = yaml["region_id"].as_str().unwrap().to_string();
        let region_type = yaml["region_type"].as_str().unwrap().into();
        let name = yaml["name"].as_str().unwrap().to_string();
        let sequence_type = match yaml["sequence_type"].as_str().unwrap() {
            "fixed" => SequenceType::Fixed(yaml["sequence"].as_str().unwrap().to_string()),
            "onlist" => {
                let md5 = yaml["onlist"].as_hash().unwrap().get(&Yaml::String("md5".to_string()))
                    .and_then(|x| x.as_str()).map(|x| x.to_string());
                let filepath = yaml["onlist"]["filename"].as_str().unwrap().to_string();
                let local = match yaml["onlist"]["location"].as_str().unwrap() {
                    "local" => true,
                    "remote" => false,
                    x => bail!("Invalid location: {:?}", x),
                };
                SequenceType::OnList { filepath, md5, local }
            },
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

    /// Return an iterator over all regions in the region tree.
    pub fn iter_regions(&self) -> impl Iterator<Item = &Region> {
        self.sub_regions.iter().flat_map(|x| {
            Box::new(std::iter::once(x).chain(x.iter_regions())) as Box<dyn Iterator<Item = &Region>>
        })
    }

    /// Return the start and end position of each subregion in a fastq region.
    pub fn subregion_range(&self) -> Vec<(RegionType, Range)> {
        if self.region_type != RegionType::Fastq {
            panic!("must be called on a fastq region");
        }

        get_region_range(&self.sub_regions).collect()
    }

    /// Return the start and end position of each subregion in a fastq region.
    /// Make adjustments as if the fastq read is reverse complemented.
    pub fn subregion_range_rev(&self) -> Vec<(RegionType, Range)> {
        if self.region_type != RegionType::Fastq {
            panic!("must be called on a fastq region");
        }

        let ranges = get_region_range(&self.sub_regions).collect::<Vec<_>>();
        if ranges.len() <=1 {
            ranges
        } else {
            let total_len = ranges.last().unwrap().1.end.unwrap();
            ranges.into_iter().map(move |(region_type, range)| {
                let start = total_len - range.end.unwrap();
                let end = start + (range.end.unwrap() - range.start);
                (region_type, Range {start, end: Some(end)})
            }).collect()
        }
    }
}

/// Return the start and end position of each subregion in a fastq region.
fn get_region_range(regions: &[Region]) -> impl Iterator<Item = (RegionType, Range)> + '_ {
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_seqspec_io() {
        let spec = SeqSpec::from_path("tests/data/spec.yaml").unwrap();

        for fq in spec.modality("rna").unwrap().fastqs() {
            println!("{:?}", fq.subregion_range());
        }
    }

    #[test]
    fn test_onlist() {
        let spec = SeqSpec::from_path("tests/data/spec.yaml").unwrap();

        spec.iter_regions().for_each(|region| {
            if let SequenceType::OnList { filepath, .. } = &region.sequence_type {
                if region.sequence_type.fetch_onlist().is_err() {
                    panic!("Failed to fetch onlist: {:?}", region.region_type);
                }
            }
        });
    }
}
