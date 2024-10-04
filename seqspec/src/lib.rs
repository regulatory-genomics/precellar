mod read;
mod region;
pub mod utils;

use log::warn;
use read::{ReadSpan, RegionIndex};
pub use read::{Read, File, UrlType, Strand};
use region::LibSpec;
pub use region::{Region, RegionType, SequenceType, Onlist};

use serde::{Deserialize, Deserializer, Serialize};
use serde_yaml::{self, Value};
use std::{fs, str::FromStr};
use anyhow::{bail, anyhow, Result};
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
    pub sequence_spec: read::SeqSpec,
    pub library_spec: LibSpec,
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
        let mut read_buffer = if let Some(r) = self.sequence_spec.get(read_id) {
            r.clone()
        } else {
            assert!(modality.is_some(), "modality must be provided for a new read");
            assert!(primer_id.is_some(), "primer_id must be provided for a new read");
            assert!(is_reverse.is_some(), "is_reverse must be provided for a new read");
            Read { read_id: read_id.to_string(), ..Default::default() }
        };

        if let Some(rev) = is_reverse {
            read_buffer.strand = if rev { Strand::Neg } else { Strand::Pos };
        }
        if let Some(modality) = modality {
            read_buffer.modality = Modality::from_str(modality)?;
        }
        if let Some(primer_id) = primer_id {
            read_buffer.primer_id = primer_id.to_string();
        }
        if let Some(fastq) = fastqs {
            read_buffer.files = Some(fastq.iter().map(|path| File::from_fastq(path)).collect::<Result<Vec<File>>>()?);
        }

        if (min_len.is_none() || max_len.is_none()) && read_buffer.files.is_some() {
            let len = read_buffer.actual_len("./")?;
            read_buffer.min_len = min_len.unwrap_or(len) as u32;
            read_buffer.max_len = max_len.unwrap_or(len) as u32;
        } else {
            read_buffer.min_len = min_len.unwrap() as u32;
            read_buffer.max_len = max_len.unwrap() as u32;
        }

        self.validate_reads(&read_buffer)?;
        self.sequence_spec.insert(read_id.to_string(), read_buffer);

        Ok(())
    }

    pub fn delete_read(&mut self, read_id: &str) -> Option<Read> {
        self.sequence_spec.shift_remove(read_id)
    }

    /// Get the index of atomic regions of each read in the sequence spec.
    pub fn get_index_by_modality(&self, modality: Modality) -> impl Iterator<Item = (&Read, RegionIndex)> {
        self.sequence_spec.values().filter_map(move |read| if read.modality == modality {
            let index = self.get_index(&read.read_id)
                .expect(&format!("Cannot find index for Read: {}", read.read_id));
            Some((read, index))
        } else {
            None
        })
    }

    /// Get the index of atomic regions of a read in the sequence spec.
    pub fn get_index(&self, read_id: &str) -> Option<RegionIndex> {
        let read = self.sequence_spec.get(read_id)?;
        let region = self.library_spec.get_parent(&read.primer_id)?;
        read.get_index(region)
    }

    pub fn iter_reads(&self, modality: Modality) -> impl Iterator<Item = &Read> {
        self.sequence_spec.values().filter(move |read| read.modality == modality)
    }

    pub fn validate_reads(&self, read: &Read) -> Result<()> {
        let region = self.library_spec.get_parent(&read.primer_id)
            .ok_or_else(|| anyhow!("Primer not found: {}", read.primer_id))?;
        // Check if the primer exists
        if let Some(index) = read.get_index(region) {
            match index.readlen_info {
                ReadSpan::Covered | ReadSpan::Span(_) => {},
                ReadSpan::NotEnough => {
                    warn!("'{}' does not cover the region", read.read_id);
                },
                ReadSpan::MayReadThrough(id) => {
                    warn!("'{}' may read through and contain sequences from: '{}'", read.read_id, id);
                },
                ReadSpan::ReadThrough(id) => {
                    warn!("'{}' length exceeds maximum length of the variable-length region (insertion), \
                    truncating the reads to the maximum length of the region. \
                    Read reads through and contains sequences from: '{}'.
                    If this is not the desired behavior, please adjust the region lengths.", read.read_id, id);
                }
            }
        }
        Ok(())
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
        for (read, index) in assay.get_index_by_modality(Modality::RNA) {
            println!("{}: {:?}", read.read_id, index.index.into_iter().map(|x| (x.1, x.2)).collect::<Vec<_>>());
        }
        for (read, index) in assay.get_index_by_modality(Modality::ATAC) {
            println!("{}: {:?}", read.read_id, index.index.into_iter().map(|x| (x.1, x.2)).collect::<Vec<_>>());
        }
        for (read, index) in assay.get_index_by_modality(Modality::Protein) {
            println!("{}: {:?}", read.read_id, index.index.into_iter().map(|x| (x.1, x.2)).collect::<Vec<_>>());
        }
    }
}