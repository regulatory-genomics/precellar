use crate::read::UrlType;
use crate::Modality;

use anyhow::Result;
use cached_path::Cache;
use indexmap::{IndexMap, IndexSet};
use serde::{Deserialize, Serialize};
use itertools::Itertools;
use std::{
    collections::HashMap,
    io::BufRead,
    ops::Deref,
    path::Path,
    sync::{Arc, RwLock},
};

#[derive(Debug, Clone)]
pub struct LibSpec {
    modalities: IndexMap<Modality, Arc<RwLock<Region>>>,
    parent_map: HashMap<String, Arc<RwLock<Region>>>,
    region_map: HashMap<String, Arc<RwLock<Region>>>,
}

impl PartialEq for LibSpec {
    fn eq(&self, other: &Self) -> bool {
        self.modalities.keys().all(|k| {
            self.modalities.get(k).unwrap().read().unwrap().deref()
                == other.modalities.get(k).unwrap().read().unwrap().deref()
        })
    }
}

impl Serialize for LibSpec {
    fn serialize<S: serde::Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        self.modalities
            .values()
            .collect::<Vec<_>>()
            .serialize(serializer)
    }
}

impl<'de> Deserialize<'de> for LibSpec {
    fn deserialize<D: serde::Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        let regions = Vec::<Region>::deserialize(deserializer)?;
        Ok(LibSpec::new(regions).unwrap())
    }
}

impl LibSpec {
    /// Create a new library specification from top-level regions. The only allowed
    /// region types are modality types.
    pub fn new<I: IntoIterator<Item = Region>>(regions: I) -> Result<Self> {
        let mut modalities = IndexMap::new();
        let mut region_map = HashMap::new();
        let mut parent_map = HashMap::new();
        for region in regions {
            if let RegionType::Modality(modality) = region.region_type {
                let region = Arc::new(RwLock::new(region));
                let region_id = region.read().unwrap().region_id.clone();
                if modalities.insert(modality, region.clone()).is_some() {
                    return Err(anyhow::anyhow!("Duplicate modality: {:?}", modality));
                }
                if region_map
                    .insert(region_id.clone(), region.clone())
                    .is_some()
                {
                    return Err(anyhow::anyhow!("Duplicate region id: {}", region_id));
                }
                for subregion in region.read().unwrap().subregions.iter() {
                    let id = subregion.read().unwrap().region_id.clone();
                    if region_map.insert(id.clone(), subregion.clone()).is_some() {
                        return Err(anyhow::anyhow!("Duplicate region id: {}", id));
                    }
                    parent_map.insert(id, region.clone());
                }
            } else {
                return Err(anyhow::anyhow!("Top-level regions must be modalities"));
            };
        }
        Ok(Self {
            modalities,
            region_map,
            parent_map,
        })
    }

    /// Iterate over all regions with modality type in the library.
    pub fn modalities(&self) -> impl Iterator<Item = &Arc<RwLock<Region>>> {
        self.modalities.values()
    }

    /// Iterate over all regions in the library.
    pub fn regions(&self) -> impl Iterator<Item = &Arc<RwLock<Region>>> {
        self.region_map.values()
    }

    pub fn get_modality(&self, modality: &Modality) -> Option<&Arc<RwLock<Region>>> {
        self.modalities.get(modality)
    }

    pub fn get(&self, region_id: &str) -> Option<&Arc<RwLock<Region>>> {
        self.region_map.get(region_id)
    }

    pub fn get_parent(&self, region_id: &str) -> Option<&Arc<RwLock<Region>>> {
        self.parent_map.get(region_id)
    }

    pub fn cat_barcodes(&self, modality: &Modality) -> Option<Vec<Vec<u8>>> {
        let region = self.get_modality(modality)?;
        region.read().unwrap().subregions.iter().filter_map(|r| {
            let r = r.read().unwrap();
            if r.region_type.is_barcode() && r.onlist.is_some() {
                Some(r.onlist.as_ref().unwrap().read().unwrap().into_iter().collect::<Vec<_>>())
            } else {
                None
            }
        }).reduce(|a, b| {
            a.into_iter().cartesian_product(b).map(|(mut a, b)| {
                a.extend(b);
                a
            }).collect()
        })
    }
}

pub type RegionId = String;

#[derive(Deserialize, Serialize, Debug, Clone)]
pub struct Region {
    pub region_id: RegionId,
    pub region_type: RegionType,
    pub name: String,
    pub sequence_type: SequenceType,
    pub sequence: String,
    pub min_len: u32,
    pub max_len: u32,
    pub onlist: Option<Onlist>,
    #[serde(rename = "regions", deserialize_with = "deserialize_regions")]
    pub subregions: Vec<Arc<RwLock<Region>>>,
}

impl PartialEq for Region {
    fn eq(&self, other: &Self) -> bool {
        self.region_id == other.region_id
            && self.region_type == other.region_type
            && self.name == other.name
            && self.sequence_type == other.sequence_type
            && self.sequence == other.sequence
            && self.min_len == other.min_len
            && self.max_len == other.max_len
            && self.onlist == other.onlist
            && self
                .subregions
                .iter()
                .zip(other.subregions.iter())
                .all(|(a, b)| a.read().unwrap().deref() == b.read().unwrap().deref())
    }
}

impl Region {
    pub fn is_fixed_length(&self) -> bool {
        self.min_len == self.max_len
    }

    pub fn len(&self) -> Option<u32> {
        if self.min_len == self.max_len {
            Some(self.min_len)
        } else {
            None
        }
    }

    // https://rust-lang.github.io/rust-clippy/master/index.html#len_without_is_empty
    // It is good custom to have both methods, because for some data structures, asking about the length will be a costly operation, whereas .is_empty() can usually answer in constant time. Also it used to lead to false positives on the len_zero lint – currently that lint will ignore such entities.
    pub fn is_empty(&self) -> bool {
        self.min_len != self.max_len
    }
}

fn deserialize_regions<'de, D>(deserializer: D) -> Result<Vec<Arc<RwLock<Region>>>, D::Error>
where
    D: serde::Deserializer<'de>,
{
    if let Some(regions) = Option::<Vec<Region>>::deserialize(deserializer)? {
        Ok(regions
            .into_iter()
            .map(|x| Arc::new(RwLock::new(x)))
            .collect())
    } else {
        Ok(Vec::new())
    }
}
#[derive(Debug, Copy, Clone, PartialEq, Serialize)]
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
    pub fn is_modality(&self) -> bool {
        matches!(self, RegionType::Modality(_))
    }

    pub fn is_barcode(&self) -> bool {
        matches!(self, RegionType::Barcode)
    }

    pub fn is_umi(&self) -> bool {
        matches!(self, RegionType::Umi)
    }

    /// Check if the region contains region of interest (insertion).
    pub fn is_target(&self) -> bool {
        matches!(self, RegionType::Gdna | RegionType::Cdna)
    }

    pub fn is_poly_nucl(&self) -> Option<u8> {
        match self {
            RegionType::PolyA => Some(b'A'),
            RegionType::PolyG => Some(b'G'),
            RegionType::PolyT => Some(b'T'),
            RegionType::PolyC => Some(b'C'),
            _ => None,
        }
    }

    pub fn is_sequencing_primer(&self) -> bool {
        matches!(
            self,
            RegionType::CustomPrimer
                | RegionType::TruseqRead1
                | RegionType::TruseqRead2
                | RegionType::NexteraRead1
                | RegionType::NexteraRead2
                | RegionType::IlluminaP5
                | RegionType::IlluminaP7
        )
    }

    pub fn from_str_case_insensitive(s: &str) -> Option<Self> {
        match s.to_lowercase().as_str() {
            "barcode" => Some(RegionType::Barcode),
            "cdna" => Some(RegionType::Cdna),
            "custom_primer" | "customprime" => Some(RegionType::CustomPrimer),
            "fastq" => Some(RegionType::Fastq),
            "fastq_link" | "fastqlink" => Some(RegionType::FastqLink),
            "gdna" => Some(RegionType::Gdna),
            "hic" => Some(RegionType::Hic),
            "illumina_p5" | "illuminap5" => Some(RegionType::IlluminaP5),
            "illumina_p7" | "illuminap7" => Some(RegionType::IlluminaP7),
            "index5" => Some(RegionType::Index5),
            "index7" => Some(RegionType::Index7),
            "linker" => Some(RegionType::Linker),
            "me1" => Some(RegionType::Me1),
            "me2" => Some(RegionType::Me2),
            "methyl" => Some(RegionType::Methyl),
            "named" => Some(RegionType::Named),
            "nextera_read1" | "nexteraread1" => Some(RegionType::NexteraRead1),
            "nextera_read2" | "nexteraread2" => Some(RegionType::NexteraRead2),
            "poly_a" | "polya" => Some(RegionType::PolyA),
            "poly_g" | "polyg" => Some(RegionType::PolyG),
            "poly_t" | "polyt" => Some(RegionType::PolyT),
            "poly_c" | "polyc" => Some(RegionType::PolyC),
            "s5" => Some(RegionType::S5),
            "s7" => Some(RegionType::S7),
            "truseq_read1" | "truseqread1" => Some(RegionType::TruseqRead1),
            "truseq_read2" | "truseqread2" => Some(RegionType::TruseqRead2),
            "umi" => Some(RegionType::Umi),
            // Handle Modality separately since it's an untagged enum
            s if Modality::from_str_case_insensitive(s).is_some() => {
                Some(RegionType::Modality(Modality::from_str_case_insensitive(s).unwrap()))
            }
            _ => None,
        }
    }
}

impl<'de> Deserialize<'de> for RegionType {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        let s = String::deserialize(deserializer)?;
        RegionType::from_str_case_insensitive(&s)
            .ok_or_else(|| serde::de::Error::custom(format!("Invalid region type: {}", s)))
    }
}

#[derive(Deserialize, Serialize, Debug, Copy, Clone, PartialEq)]
#[serde(rename_all = "lowercase")]
pub enum SequenceType {
    Fixed,  // sequence string is known and fixed in length and nucleotide composition
    Random, // the sequence is not known a-priori
    Onlist, // the sequence is derived from an onlist
    Joined, // the sequence is created from nested regions and the regions property must contain Regions
}

#[derive(Deserialize, Serialize, Debug, Clone, PartialEq)]
pub struct Onlist {
    pub file_id: String,
    pub filename: String,
    pub filetype: String,
    pub filesize: u64,
    pub url: String,
    pub urltype: UrlType,
    pub location: Option<Location>,
    pub md5: String,
}

impl Onlist {
    pub fn read(&self) -> Result<IndexSet<Vec<u8>>> {
        let mut cache = Cache::new()?;
        cache.dir = home::home_dir().unwrap().join(".cache/seqspec");
        let file = cache.cached_path(&self.url)?;
        let reader = std::io::BufReader::new(crate::utils::open_file(file)?);
        reader.lines().map(|x| Ok(x?.into_bytes())).collect()
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
}

#[derive(Deserialize, Serialize, Debug, Clone, Copy, PartialEq)]
#[serde(rename_all = "lowercase")]
pub enum Location {
    Local,
    Remote,
}
