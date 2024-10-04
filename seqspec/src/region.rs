use crate::Modality;
use crate::read::UrlType;

use cached_path::Cache;
use indexmap::IndexMap;
use serde::{Deserialize, Serialize};
use std::{io::BufRead, sync::Arc};
use anyhow::Result;

#[derive(Debug, Clone, PartialEq)]
pub struct LibSpec {
    modalities: IndexMap<Modality, Arc<Region>>,
    region_map: IndexMap<String, Arc<Region>>,
    parent_map: IndexMap<String, Arc<Region>>,
}

impl Serialize for LibSpec {
    fn serialize<S: serde::Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        self.region_map.values().collect::<Vec<_>>().serialize(serializer)
    }
}

impl<'de> Deserialize<'de> for LibSpec {
    fn deserialize<D: serde::Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        let regions = Vec::<Region>::deserialize(deserializer)?;
        Ok(LibSpec::new(regions))
    }
}

impl LibSpec {
    pub fn new<I: IntoIterator<Item = Region>>(regions: I) -> Self {
        let mut modalities = IndexMap::new();
        let mut region_map = IndexMap::new();
        let mut parent_map = IndexMap::new();
        for region in regions {
            let region = Arc::new(region);
            if region_map.insert(region.region_id.clone(), region.clone()).is_some() {
                panic!("Duplicate region id: {}", region.region_id);
            }
            if let RegionType::Modality(modality) = region.region_type {
                if modalities.insert(modality, region.clone()).is_some() {
                    panic!("Duplicate modality: {:?}", modality);
                }
            }
            for subregion in region.subregions.iter() {
                parent_map.insert(subregion.region_id.clone(), region.clone());
            }
        }
        Self { modalities, region_map, parent_map }
    }

    /// Iterate over all regions with modality type in the library.
    pub fn modalities(&self) -> impl Iterator<Item = &Arc<Region>> {
        self.modalities.values()
    }

    /// Iterate over all regions in the library.
    pub fn regions(&self) -> impl Iterator<Item = &Arc<Region>> {
        self.region_map.values()
    }

    pub fn get_modality(&self, modality: &Modality) -> Option<&Arc<Region>> {
        self.modalities.get(modality)
    }

    pub fn get_modality_mut(&mut self, modality: &Modality) -> Option<&mut Arc<Region>> {
        self.modalities.get_mut(modality)
    }

    pub fn get(&self, region_id: &str) -> Option<&Arc<Region>> {
        self.region_map.get(region_id)
    }

    pub fn get_mut(&mut self, region_id: &str) -> Option<&mut Arc<Region>> {
        self.region_map.get_mut(region_id)
    }

    pub fn get_parent(&self, region_id: &str) -> Option<&Arc<Region>> {
        self.parent_map.get(region_id)
    }

    pub fn get_parent_mut(&mut self, region_id: &str) -> Option<&mut Arc<Region>> {
        self.parent_map.get_mut(region_id)
    }
}

#[derive(Deserialize, Serialize, Debug, Clone, PartialEq)]
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
    pub subregions: Vec<Arc<Region>>,
}

impl Region {
    pub fn is_fixed_length(&self) -> bool {
        if self.min_len == self.max_len {
            true
        } else {
            false
        }
    }
}

fn deserialize_regions<'de, D>(deserializer: D) -> Result<Vec<Arc<Region>>, D::Error>
where
    D: serde::Deserializer<'de>,
{
    if let Some(regions) = Option::<Vec::<Region>>::deserialize(deserializer)? {
        Ok(regions.into_iter().map(Arc::new).collect())
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
    /// Either a barcode, UMI, or DNA/cDNA region.
    pub fn is_target(&self) -> bool {
        self.is_barcode() || self.is_umi() || self.is_dna()
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
    pub fn read(&self) -> Result<Vec<String>> {
        let cache = Cache::new()?;
        let file = cache.cached_path(&self.url)?;
        let reader = std::io::BufReader::new(crate::utils::open_file_for_read(file));
        Ok(reader.lines().map(|x| x.unwrap()).collect())
    }
}

#[derive(Deserialize, Serialize, Debug, Clone, Copy, PartialEq)]
#[serde(rename_all = "lowercase")]
pub enum Location {
    Local,
    Remote,
}