use log::warn;
use serde::{Deserialize, Deserializer};
use serde_yaml::{self, Value};
use std::{fs, ops::Range, str::FromStr};
use anyhow::{bail, Result};
use noodles::fastq;
use std::{io::{BufRead, BufReader}, path::{Path, PathBuf}};
use cached_path::Cache;

use crate::io::open_file_for_read;

#[derive(Deserialize, Debug)]
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
    pub sequence_spec: Option<Vec<Read>>,
    pub library_spec: Option<Vec<Region>>,
}

impl Assay {
    pub fn from_path<P: AsRef<Path>>(path: P) -> Result<Self> {
        let yaml_str = fs::read_to_string(path)?;
        Ok(serde_yaml::from_str(&yaml_str)?)
    }

    /// Get the index of atomic regions of each read in the sequence spec.
    pub fn get_index_of(&self, modality: Modality) -> impl Iterator<Item = (&Read, Vec<(&Region, Range<u32>)>)> {
        self.sequence_spec.as_ref().expect("sequence_spec is empty")
            .iter().filter_map(move |read| if read.modality == modality {
                let region = self.get_region_by_modality(modality)?;
                let index = read.get_index(region).expect("primer not found");
                Some((read, index))
            } else {
                None
            })
    }

    pub fn iter_reads(&self, modality: Modality) -> impl Iterator<Item = &Read> {
        self.sequence_spec.as_ref().expect("sequence_spec is empty").iter()
            .filter(move |read| read.modality == modality)
    }

    /// Get the top-level region for a given modality, i.e., the region that contains all other regions.
    pub fn get_region_by_modality(&self, modality: Modality) -> Option<&Region> {
        self.library_spec.as_ref()?.iter().find(|region| {
            region.sequence_type == SequenceType::Joined &&
                region.region_type == RegionType::Modality(modality)
        })
    }
}

#[derive(Deserialize, Debug, Copy, Clone, PartialEq, Eq, Hash)]
#[serde(rename_all = "lowercase")]
pub enum Modality {
    Dna,
    Rna,
    Tag,
    Protein,
    Atac,
    Crispr,
}

impl FromStr for Modality {
    type Err = anyhow::Error;

    fn from_str(s: &str) -> Result<Self> {
        match s {
            "dna" => Ok(Modality::Dna),
            "rna" => Ok(Modality::Rna),
            "tag" => Ok(Modality::Tag),
            "protein" => Ok(Modality::Protein),
            "atac" => Ok(Modality::Atac),
            "crispr" => Ok(Modality::Crispr),
            _ => bail!("Invalid modality: {}", s),
        }
    }
}

#[derive(Debug)]
pub enum LibraryProtocol {
    Standard(String),
    Custom(Vec<ProtocolItem>),
}

#[derive(Deserialize, Debug)]
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
            _ => Err(serde::de::Error::custom("invalid value")),
        }
    }
}

#[derive(Debug)]
pub enum LibraryKit {
    Standard(String),
    Custom(Vec<KitItem>),
}

#[derive(Deserialize, Debug)]
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
            _ => Err(serde::de::Error::custom("invalid value")),
        }
    }
}

#[derive(Debug)]
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
            _ => Err(serde::de::Error::custom("invalid value")),
        }
    }
}

#[derive(Debug)]
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
            _ => Err(serde::de::Error::custom("invalid value")),
        }
    }
}

#[derive(Deserialize, Debug)]
pub struct Read {
    read_id: String,
    pub name: Option<String>,
    pub modality: Modality,
    pub primer_id: String,
    pub min_len: u32,
    pub max_len: u32,
    pub strand: Strand,
    pub files: Option<Vec<File>>,
}

impl Read {
    pub fn id(&self) -> &str {
        &self.read_id
    }

    pub fn read_fastq<P: AsRef<Path>>(&self, base_dir: P) -> fastq::Reader<impl BufRead> {
        let base_dir = base_dir.as_ref().to_path_buf();
        let reader = multi_reader::MultiReader::new(
            self.files.clone().unwrap().into_iter().map(move |file| file.open(&base_dir))
        );
        fastq::Reader::new(BufReader::new(reader))
    }

    fn get_index<'a>(&'a self, region: &'a Region) -> Option<Vec<(&'a Region, Range<u32>)>> {
        if region.sequence_type != SequenceType::Joined {
            return None;
        }

        let mut found_primer = false;

        let result = if self.is_reverse() {
            self.components(
                region.regions.as_ref().unwrap().iter().rev()
                    .skip_while(|region| {
                        let found = !(region.region_type.is_sequencing_primer() && region.region_id == self.primer_id);
                        if found {
                            found_primer = true;
                        }
                        found
                    }).skip(1)
            )
        } else {
            self.components(
                region.regions.as_ref().unwrap().iter()
                    .skip_while(|region| {
                        let found = !(region.region_type.is_sequencing_primer() && region.region_id == self.primer_id);
                        if found {
                            found_primer = true;
                        }
                        found
                    }).skip(1)
            )
        };
        
        if found_primer {
            Some(result)
        } else {
            None
        }
    }

    fn components<'a, I: Iterator<Item = &'a Region>>(&self, regions: I) -> Vec<(&'a Region, Range<u32>)> {
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
                warn!("Reads ({}) may contain additional bases downstream of the variable-length region, e.g., adapter sequences.", self.read_id);
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

#[derive(Deserialize, Debug)]
#[serde(rename_all = "lowercase")]
pub enum Strand {
    Pos,
    Neg,
}

#[derive(Deserialize, Debug, Clone)]
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
    pub fn open<P: AsRef<Path>>(&self, base_dir: P) -> Box<dyn std::io::Read> {
        match self.urltype {
            UrlType::Local => {
                let mut path = PathBuf::from(&self.url);
                path = if path.is_absolute() {
                    path
                } else {
                    base_dir.as_ref().join(path)
                };
                Box::new(open_file_for_read(path))
            }
            _ => {
                let cache = Cache::new().unwrap();
                let file = cache.cached_path(&self.url).unwrap();
                Box::new(open_file_for_read(file))
            }
        }
    }
}

#[derive(Deserialize, Debug, Copy, Clone)]
#[serde(rename_all = "lowercase")]
pub enum UrlType {
    Local,
    Ftp,
    Http,
    Https,
}

#[derive(Deserialize, Debug, Clone)]
pub struct Region {
    pub region_id: String,
    pub region_type: RegionType,
    pub sequence_type: SequenceType,
    pub sequence: String,
    pub min_len: u32,
    pub max_len: u32,
    pub onlist: Option<Onlist>,
    pub regions: Option<Vec<Region>>,
}

impl Region {
    /// Return an iterator over all regions in the region tree.
    pub fn iter_regions(&self) -> impl Iterator<Item = &Region> {
        self.regions.as_ref().unwrap().iter()
    }
}

#[derive(Deserialize, Debug, Copy, Clone, PartialEq)]
#[serde(rename_all = "lowercase")]
pub enum RegionType {
    Barcode,
    Cdna,
    #[serde(rename = "custom_primer")]
    CustomPrimer,
    Fastq,
    FastqLink,
    Gdna,
    Hic,
    #[serde(rename = "illumina_p5")]
    IlluminaP5,
    #[serde(rename = "illumina_p7")]
    IlluminaP7,
    Index5,
    Index7,
    Linker,
    Me1,
    Me2,
    Methyl,
    Named,
    #[serde(rename = "nextera_read1")]
    NexteraRead1,
    #[serde(rename = "nextera_read2")]
    NexteraRead2,
    #[serde(rename = "poly_a")]
    PolyA,
    #[serde(rename = "poly_g")]
    PolyG,
    #[serde(rename = "poly_t")]
    PolyT,
    #[serde(rename = "poly_c")]
    PolyC,
    S5,
    S7,
    #[serde(rename = "truseq_read1")]
    TruseqRead1,
    #[serde(rename = "truseq_read2")]
    TruseqRead2,
    Umi,
    #[serde(untagged)]
    Modality(Modality),
}

impl RegionType {
    pub fn is_barcode(&self) -> bool {
        match self {
            RegionType::Barcode => true,
            _ => false,
        }
    }

    pub fn is_umi(&self) -> bool {
        match self {
            RegionType::Umi => true,
            _ => false,
        }
    }

    /// Check if the region contains genomic sequences.
    pub fn is_dna(&self) -> bool {
        match self {
            RegionType::Gdna | RegionType::Cdna => true,
            _ => false,
        }
    }

    pub fn is_sequencing_primer(&self) -> bool {
        match self {
            RegionType::CustomPrimer |
                RegionType::TruseqRead1 | RegionType::TruseqRead2 |
                RegionType::NexteraRead1 | RegionType::NexteraRead2 |
                RegionType::IlluminaP5 | RegionType::IlluminaP7 => true,
            _ => false,
        }
    }
}

#[derive(Deserialize, Debug, Copy, Clone, PartialEq)]
#[serde(rename_all = "lowercase")]
pub enum SequenceType {
    Fixed,  // sequence string is known and fixed in length and nucleotide composition
    Random,  // the sequence is not known a-priori
    Onlist,  // the sequence is derived from an onlist
    Joined,  // the sequence is created from nested regions and the regions property must contain Regions
}

#[derive(Deserialize, Debug, Clone)]
pub struct Onlist {
    pub file_id: String,
    pub filename: String,
    pub filetype: String,
    pub filesize: u64,
    url: String,
    pub urltype: UrlType,
    pub location: Option<Location>,
    pub md5: String,
}

impl Onlist {
    pub fn read(&self) -> Result<Vec<String>> {
        let cache = Cache::new()?;
        let file = cache.cached_path(&self.url)?;
        let reader = std::io::BufReader::new(open_file_for_read(file));
        Ok(reader.lines().map(|x| x.unwrap()).collect())
    }
}

#[derive(Deserialize, Debug, Clone, Copy)]
#[serde(rename_all = "lowercase")]
pub enum Location {
    Local,
    Remote,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse() {
        let protocol: LibraryProtocol = serde_yaml::from_str("10X").unwrap();
        println!("{:?}", protocol);

        let protocol: LibraryProtocol = serde_yaml::from_str(
            "- !LibProtocol
            protocol_id: CG000338 Chromium Next GEM Multiome ATAC + Gene Expression Rev. D protocol (10x Genomics)
            name: DogmaSeq-DIG
            modality: rna"
        ).unwrap();
        println!("{:?}", protocol);
    }

    #[test]
    fn test_parse_yaml() {
        let yaml_str = fs::read_to_string("tests/data/spec.yaml").expect("Failed to read file");

        let assay: Assay = serde_yaml::from_str(&yaml_str).expect("Failed to parse YAML");
        for x in assay.get_index_of(Modality::Protein) {
            println!("{:?}", x);
        }
    }
}