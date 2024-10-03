mod read;
mod region;
pub mod utils;

pub use read::{Read, File, UrlType, Strand};
pub use region::{Region, RegionType, SequenceType, Onlist};

use serde::{Deserialize, Deserializer, Serialize};
use serde_yaml::{self, Value};
use std::{fs, ops::Range, str::FromStr};
use anyhow::{bail, Result};
use std::path::Path;

#[derive(Deserialize, Serialize, Debug, Clone, PartialEq)]
pub struct Assay {
    pub seqspec_version: String,
    pub assay_id: String,
    pub name: String,
    pub doi: String,
    pub date: String,
    pub description: String,
    pub modalities: Vec<Modality>,
    pub lib_struct: String,
    pub library_protocol: LibraryProtocol,
    pub library_kit: LibraryKit,
    pub sequence_protocol: SequenceProtocol,
    pub sequence_kit: SequenceKit,
    pub sequence_spec: read::Sequences,
    pub library_spec: Vec<Region>,
}

impl Assay {
    pub fn from_path<P: AsRef<Path>>(path: P) -> Result<Self> {
        let yaml_str = fs::read_to_string(path)?;
        Ok(serde_yaml::from_str(&yaml_str)?)
    }

    /// Update read information. If the read does not exist, it will be created.
    pub fn update_read<P: AsRef<Path>>(
        &mut self,
        read_id: &str,
        modality: Option<&str>,
        primer_id: Option<&str>,
        is_reverse: Option<bool>,
        fastqs: Option<&[P]>,
        min_len: Option<usize>,
        max_len: Option<usize>,
    ) -> Result<()> {
        let all_reads = &mut self.sequence_spec;
        let mut read_exist = false;
        let mut read_buffer = Read::default();
        read_buffer.read_id = read_id.to_string();

        let read = if let Some(r) = all_reads.get_mut(read_id) {
            read_exist = true;
            r
        } else {
            &mut read_buffer
        };

        if !read_exist {
            assert!(modality.is_some(), "modality must be provided for a new read");
            assert!(primer_id.is_some(), "primer_id must be provided for a new read");
            assert!(is_reverse.is_some(), "is_reverse must be provided for a new read");
        }

        if let Some(rev) = is_reverse {
            read.strand = if rev { Strand::Neg } else { Strand::Pos };
        }
        if let Some(modality) = modality {
            read.modality = Modality::from_str(modality)?;
        }
        if let Some(primer_id) = primer_id {
            read.primer_id = primer_id.to_string();
        }
        if let Some(fastq) = fastqs {
            read.files = Some(fastq.iter().map(|path| File::from_fastq(path)).collect::<Result<Vec<File>>>()?);
        }

        if (min_len.is_none() || max_len.is_none()) && read.files.is_some() {
            let len = read.actual_len("./")?;
            read.min_len = min_len.unwrap_or(len) as u32;
            read.max_len = max_len.unwrap_or(len) as u32;
        } else {
            read.min_len = min_len.unwrap() as u32;
            read.max_len = max_len.unwrap() as u32;
        }

        if !read_exist {
            all_reads.insert(read_buffer.read_id.clone(), read_buffer);
        }

        Ok(())
    }

    pub fn delete_read(&mut self, read_id: &str) -> Option<Read> {
        self.sequence_spec.shift_remove(read_id)
    }

    /// Get the index of atomic regions of each read in the sequence spec.
    pub fn get_index_of(&self, modality: Modality) -> impl Iterator<Item = (&Read, Vec<(&Region, Range<u32>)>)> {
        self.sequence_spec.values().filter_map(move |read| if read.modality == modality {
            let region = self.get_region_by_modality(modality)?;
            let index = read.get_index(region)
                .expect(&format!("Region: {} does not have Primer: {}", region.region_id, read.primer_id));
            Some((read, index))
        } else {
            None
        })
    }

    pub fn iter_reads(&self, modality: Modality) -> impl Iterator<Item = &Read> {
        self.sequence_spec.values().filter(move |read| read.modality == modality)
    }

    /// Get the top-level region for a given modality, i.e., the region that contains all other regions.
    pub fn get_region_by_modality(&self, modality: Modality) -> Option<&Region> {
        self.library_spec.iter().find(|region| {
            region.sequence_type == SequenceType::Joined &&
                region.region_type == RegionType::Modality(modality)
        })
    }
}

#[derive(Serialize, Debug, Copy, Clone, PartialEq, Eq, Hash)]
#[serde(rename_all = "lowercase")]
pub enum Modality {
    DNA,
    RNA,
    Tag,
    Protein,
    ATAC,
    Crispr,
}

impl<'de> Deserialize<'de> for Modality {
    fn deserialize<D>(deserializer: D) -> Result<Modality, D::Error>
    where
        D: Deserializer<'de>,
    {
        let value = Value::deserialize(deserializer)?;
        match value {
            Value::String(s) => Modality::from_str(&s).map_err(serde::de::Error::custom),
            _ => Err(serde::de::Error::custom(format!("invalid value: {:?}", value))),
        }
    }
}

impl FromStr for Modality {
    type Err = anyhow::Error;

    fn from_str(s: &str) -> Result<Self> {
        match s.to_lowercase().as_str() {
            "dna" => Ok(Modality::DNA),
            "rna" => Ok(Modality::RNA),
            "tag" => Ok(Modality::Tag),
            "protein" => Ok(Modality::Protein),
            "atac" => Ok(Modality::ATAC),
            "crispr" => Ok(Modality::Crispr),
            _ => bail!("Invalid modality: {}", s),
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub enum LibraryProtocol {
    Standard(String),
    Custom(Vec<ProtocolItem>),
}

#[derive(Deserialize, Serialize, Clone, Debug, PartialEq)]
pub struct ProtocolItem {
    pub protocol_id: String,
    pub name: Option<String>,
    pub modality: Modality,
}

impl<'de> Deserialize<'de> for LibraryProtocol {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        let value = Value::deserialize(deserializer)?;
        match value {
            Value::String(s) => Ok(LibraryProtocol::Standard(s)),
            Value::Sequence(seq) => {
                let items = seq.into_iter().map(|item| serde_yaml::from_value(item).unwrap())
                    .collect::<Vec<ProtocolItem>>();
                Ok(LibraryProtocol::Custom(items))
            }
            _ => Err(serde::de::Error::custom(format!("invalid value: {:?}", value))),
        }
    }
}

impl Serialize for LibraryProtocol {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::ser::Serializer
    {
        match self {
            LibraryProtocol::Standard(s) => s.serialize(serializer),
            LibraryProtocol::Custom(items) => items.serialize(serializer),
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub enum LibraryKit {
    Standard(String),
    Custom(Vec<KitItem>),
}

#[derive(Deserialize, Serialize, Clone, Debug, PartialEq)]
pub struct KitItem {
    pub kit_id: String,
    pub name: Option<String>,
    pub modality: Modality,
}

impl<'de> Deserialize<'de> for LibraryKit {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        let value = Value::deserialize(deserializer)?;
        match value {
            Value::String(s) => Ok(LibraryKit::Standard(s)),
            Value::Sequence(seq) => {
                let items = seq.into_iter().map(|item| serde_yaml::from_value(item).unwrap())
                    .collect::<Vec<KitItem>>();
                Ok(LibraryKit::Custom(items))
            }
            _ => Err(serde::de::Error::custom(format!("invalid value: {:?}", value))),
        }
    }
}

impl Serialize for LibraryKit {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::ser::Serializer
    {
        match self {
            LibraryKit::Standard(s) => s.serialize(serializer),
            LibraryKit::Custom(items) => items.serialize(serializer),
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub enum SequenceProtocol {
    Custom(Vec<ProtocolItem>),
    Standard(String),
}

impl <'de> Deserialize<'de> for SequenceProtocol {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        let value = Value::deserialize(deserializer)?;
        match value {
            Value::String(s) => Ok(SequenceProtocol::Standard(s)),
            Value::Sequence(seq) => {
                let items = seq.into_iter().map(|item| serde_yaml::from_value(item).unwrap())
                    .collect::<Vec<ProtocolItem>>();
                Ok(SequenceProtocol::Custom(items))
            }
            _ => Err(serde::de::Error::custom(format!("invalid value: {:?}", value))),
        }
    }
}

impl Serialize for SequenceProtocol {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::ser::Serializer
    {
        match self {
            SequenceProtocol::Standard(s) => s.serialize(serializer),
            SequenceProtocol::Custom(items) => items.serialize(serializer),
        }
    }
}

#[derive(Debug, PartialEq, Clone)]
pub enum SequenceKit {
    Standard(String),
    Custom(Vec<KitItem>),
}

impl<'de> Deserialize<'de> for SequenceKit {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        let value = Value::deserialize(deserializer)?;
        match value {
            Value::String(s) => Ok(SequenceKit::Standard(s)),
            Value::Sequence(seq) => {
                let items = seq.into_iter().map(|item| serde_yaml::from_value(item).unwrap())
                    .collect::<Vec<KitItem>>();
                Ok(SequenceKit::Custom(items))
            }
            _ => Err(serde::de::Error::custom(format!("invalid value: {:?}", value))),
        }
    }
}

impl Serialize for SequenceKit {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::ser::Serializer
    {
        match self {
            SequenceKit::Standard(s) => s.serialize(serializer),
            SequenceKit::Custom(items) => items.serialize(serializer),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    
    const YAML_FILE: &str = "../seqspec_templates/10x_rna_atac.yaml";

    #[test]
    fn test_parse() {
        let yaml_str = fs::read_to_string(YAML_FILE).expect("Failed to read file");
        let assay: Assay = serde_yaml::from_str(&yaml_str).expect("Failed to parse YAML");

        println!("{:?}", assay);
    }

    #[test]
    fn test_serialize() {
        fn se_de(yaml_str: &str) {
            let assay: Assay = serde_yaml::from_str(&yaml_str).expect("Failed to parse YAML");
            let yaml_str_ = serde_yaml::to_string(&assay).unwrap();
            let assay_ = serde_yaml::from_str(&yaml_str_).expect("Failed to parse YAML");
            assert_eq!(assay, assay_);
        }

        se_de(&fs::read_to_string(YAML_FILE).expect("Failed to read file"));
    }

    #[test]
    fn test_index() {
        let yaml_str = fs::read_to_string(YAML_FILE).expect("Failed to read file");

        let assay: Assay = serde_yaml::from_str(&yaml_str).expect("Failed to parse YAML");
        for (read, regions) in assay.get_index_of(Modality::RNA) {
            println!("{}: {:?}", read.read_id, regions.into_iter().map(|x| (x.0.region_type, x.1)).collect::<Vec<_>>());
        }
        for (read, regions) in assay.get_index_of(Modality::ATAC) {
            println!("{}: {:?}", read.read_id, regions.into_iter().map(|x| (x.0.region_type, x.1)).collect::<Vec<_>>());
        }
        for (read, regions) in assay.get_index_of(Modality::Protein) {
            println!("{}: {:?}", read.read_id, regions.into_iter().map(|x| (x.0.region_type, x.1)).collect::<Vec<_>>());
        }
    }
}