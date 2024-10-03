use crate::Modality;
use crate::read::UrlType;

use cached_path::Cache;
use serde::{Deserialize, Serialize};
use std::io::BufRead;
use anyhow::Result;

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
    pub regions: Option<Vec<Region>>,
}

impl Default for Region {
    fn default() -> Self {
        Self {
            region_id: "".to_string(),
            region_type: RegionType::Named,
            name: "".to_string(),
            sequence_type: SequenceType::Fixed,
            sequence: "".to_string(),
            min_len: 0,
            max_len: 0,
            onlist: None,
            regions: None,
        }
    }
}

impl Region {
    /// Return an iterator over all regions in the region tree.
    pub fn iter_regions(&self) -> impl Iterator<Item = &Region> {
        self.regions.as_ref().unwrap().iter()
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

