use crate::Modality;
use crate::read::UrlType;

use cached_path::Cache;
use indexmap::IndexMap;
use serde::{Deserialize, Serialize};
use std::{collections::{HashMap, HashSet}, io::BufRead, ops::Deref, path::Path, sync::{Arc, RwLock}};
use anyhow::Result;

#[derive(Debug, Clone)]
pub struct LibSpec {
    modalities: IndexMap<Modality, Arc<RwLock<Region>>>,
    parent_map: HashMap<String, Arc<RwLock<Region>>>,
    region_map: HashMap<String, Arc<RwLock<Region>>>,
}

impl PartialEq for LibSpec {
    fn eq(&self, other: &Self) -> bool {
        self.modalities.keys().all(|k| {
            self.modalities.get(k).unwrap().read().unwrap().deref() ==
                other.modalities.get(k).unwrap().read().unwrap().deref()
        })
    }
}

impl Serialize for LibSpec {
    fn serialize<S: serde::Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        self.modalities.values().collect::<Vec<_>>().serialize(serializer)
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
                if region_map.insert(region_id.clone(), region.clone()).is_some() {
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
        Ok(Self { modalities, region_map, parent_map })
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
}

#[derive(Deserialize, Serialize, Debug, Clone)]
pub struct Region {
    pub region_id: String,
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
        self.region_id == other.region_id &&
            self.region_type == other.region_type &&
            self.name == other.name &&
            self.sequence_type == other.sequence_type &&
            self.sequence == other.sequence &&
            self.min_len == other.min_len &&
            self.max_len == other.max_len &&
            self.onlist == other.onlist &&
            self.subregions.iter().zip(other.subregions.iter())
                .all(|(a, b)| a.read().unwrap().deref() == b.read().unwrap().deref())
    }
}

impl Region {
    pub fn is_fixed_length(&self) -> bool {
        if self.min_len == self.max_len {
            true
        } else {
            false
        }
    }
    
    pub fn len(&self) -> Option<u32> {
        if self.min_len == self.max_len {
            Some(self.min_len)
        } else {
            None
        }
    }
}

fn deserialize_regions<'de, D>(deserializer: D) -> Result<Vec<Arc<RwLock<Region>>>, D::Error>
where
    D: serde::Deserializer<'de>,
{
    if let Some(regions) = Option::<Vec::<Region>>::deserialize(deserializer)? {
        Ok(regions.into_iter().map(|x| Arc::new(RwLock::new(x))).collect())
    } else {
        Ok(Vec::new())
    }
}

#[derive(Deserialize, Serialize, Debug, Copy, Clone, PartialEq)]
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
        match self {
            RegionType::Modality(_) => true,
            _ => false,
        }
    }

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

    /// Check if the region contains region of interest (insertion).
    pub fn is_target(&self) -> bool {
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

#[derive(Deserialize, Serialize, Debug, Copy, Clone, PartialEq)]
#[serde(rename_all = "lowercase")]
pub enum SequenceType {
    Fixed,  // sequence string is known and fixed in length and nucleotide composition
    Random,  // the sequence is not known a-priori
    Onlist,  // the sequence is derived from an onlist
    Joined,  // the sequence is created from nested regions and the regions property must contain Regions
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
    pub fn read(&self) -> Result<HashSet<Vec<u8>>> {
        let mut cache = Cache::new()?;
        cache.dir = home::home_dir().unwrap().join(".cache/seqspec");
        let file = cache.cached_path(&self.url)?;
        let reader = std::io::BufReader::new(crate::utils::open_file_for_read(file)?);
        reader.lines().map(|x| Ok(x?.into_bytes())).collect()
    }

    pub(crate) fn normalize_path<P: AsRef<Path>>(&mut self, work_dir: P) -> Result<()> {
        if self.urltype == UrlType::Local {
            self.url = crate::utils::normalize_path(work_dir, &mut self.url)?
                .to_str().unwrap().to_owned();
        }
        Ok(())
    }

    pub(crate) fn unnormalize_path<P: AsRef<Path>>(&mut self, work_dir: P) -> Result<()> {
        if self.urltype == UrlType::Local {
            self.url = crate::utils::unnormalize_path(work_dir, &mut self.url)?
                .to_str().unwrap().to_owned();
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
