use anyhow::{anyhow, Result};
use std::collections::HashMap;
use std::sync::{Arc, RwLock};

pub mod barcode_extractor;
pub mod sequence_aligner;

use crate::barcode::OligoFrequncy;
use seqspec::region::{Region, LibSpec};

/// Long-read barcode extraction result
#[derive(Debug, Clone)]
pub struct LongReadBarcodeResult {
    /// Extracted barcode sequence (standardized from whitelist)
    pub barcode: Option<Vec<u8>>,
    /// Extraction confidence
    pub confidence: f64,
    /// Whether extraction was successful
    pub success: bool,
}

/// Main interface for long-read processing
#[derive(Debug)]
pub struct LongReadProcessor {
    /// Barcode whitelists
    whitelists: HashMap<String, OligoFrequncy>,
}

impl LongReadProcessor {
    /// Create a new long-read processor
    pub fn new(
        whitelists: HashMap<String, OligoFrequncy>,
    ) -> Self {
        Self {
            whitelists,
        }
    }

    /// Extract barcode from FASTQ record
    pub fn extract_barcode(
        &self,
        record: &noodles::fastq::Record,
        lib_spec: &LibSpec,
        modality: &seqspec::Modality,
    ) -> Result<LongReadBarcodeResult> {
        let extractor = barcode_extractor::BarcodeExtractor::new(lib_spec, modality)?;

        extractor.extract_barcode(record, &self.whitelists)
    }
}

/// End type information for describing 5' or 3' end of sequence
#[derive(Debug, Clone, PartialEq)]
pub enum EndType {
    FivePrime,
    ThreePrime,
}

/// Region collection information containing all relevant regions for barcode extraction in one end
#[derive(Debug, Clone)]
pub struct EndRegions {
    /// End type
    pub end_type: EndType,
    /// Region list ordered from outside to inside
    pub regions: Vec<Arc<RwLock<Region>>>,
    /// Whether contains barcode region
    pub has_barcode: bool,
    /// Maximum length (sum of all region max_lens)
    pub max_len: u32,
}

impl EndRegions {
    /// Create new end region collection
    pub fn new(end_type: EndType) -> Self {
        Self {
            end_type,
            regions: Vec::new(),
            has_barcode: false,
            max_len: 0,
        }
    }

    /// Add region to collection
    pub fn add_region(&mut self, region: Arc<RwLock<Region>>) {
        let region_guard = region.read().unwrap();
        if region_guard.region_type.is_barcode() {
            self.has_barcode = true;
    }
        self.max_len += region_guard.max_len;
        drop(region_guard);
        self.regions.push(region);
    }

    /// Calculate cut length (sum of max_len + 15% buffer)
    /// 
    /// extra 15% is for potential indel in long read sequencing.
    pub fn calculate_cut_length(&self) -> usize {
        if self.has_barcode {
            let base_len = self.max_len as f64;
            (base_len * 1.15).ceil() as usize
        } else {
            // If no barcode, only cut 100bp
            100
        }
    }

    /// Get all fixed sequence regions (sequence_type == fixed and length >= 15)
    pub fn get_fixed_regions(&self) -> Vec<Arc<RwLock<Region>>> {
        self.regions
            .iter()
            .filter(|region| {
                let r = region.read().unwrap();
                r.sequence_type.is_fixed() && r.min_len >= 15 // fixed sequence must be at least 15bp
            })
            .cloned()
            .collect()
    }

    /// Get all barcode regions
    pub fn get_barcode_regions(&self) -> Vec<Arc<RwLock<Region>>> {
        self.regions
            .iter()
            .filter(|region| {
                let r = region.read().unwrap();
                r.region_type.is_barcode()
            })
            .cloned()
            .collect()
    }

    /// Build position map: region_id -> position (1-based, outermost is 1)
    pub fn build_position_map(&self) -> HashMap<String, usize> {
        let mut position_map = HashMap::new();
        for (idx, region) in self.regions.iter().enumerate() {
            let region_guard = region.read().unwrap();
            position_map.insert(region_guard.region_id.clone(), idx + 1);
        }
        position_map
    }
}

/// Analyze library specification and find innermost regions at both ends
/// 
/// Trimmed EndRegions will be used for barcode extraction.
pub fn find_innermost_regions(lib_spec: &LibSpec, modality: &seqspec::Modality) -> Result<(EndRegions, EndRegions)> {
    let modality_region = lib_spec
        .get_modality(modality)
        .ok_or_else(|| anyhow!("Cannot find specified modality: {:?}", modality))?;

    let modality_guard = modality_region.read().unwrap();
    let subregions = &modality_guard.subregions;

    if subregions.is_empty() {
        return Err(anyhow!("Modality region has no subregions"));
    }

    let mut five_prime_regions = EndRegions::new(EndType::FivePrime);
    let mut three_prime_regions = EndRegions::new(EndType::ThreePrime);

    // Traverse from 5' end
    for region in subregions.iter() {
        let region_guard = region.read().unwrap();
        let region_type = &region_guard.region_type;

        // Stop if encounter cDNA or gDNA
        if region_type.is_target() {
            break;
        }

        // Add to 5' collection if fixed or barcode type
        if region_guard.sequence_type.is_fixed() || region_type.is_barcode() {
            drop(region_guard);
            five_prime_regions.add_region(region.clone());
        }
    }

    // Traverse from 3' end (reverse)
    for region in subregions.iter().rev() {
        let region_guard = region.read().unwrap();
        let region_type = &region_guard.region_type;

        // Stop if encounter cDNA or gDNA
        if region_type.is_target() {
            break;
        }

        // Add to 3' collection if fixed or barcode type
        if region_guard.sequence_type.is_fixed() || region_type.is_barcode() {
            drop(region_guard);
            three_prime_regions.add_region(region.clone());
        }
    }

    Ok((five_prime_regions, three_prime_regions))
}

#[cfg(test)]
mod tests {
    use super::*;
    use seqspec::region::Region;
    use seqspec::{Modality, RegionType, SequenceType};

    fn create_test_region(
        id: &str,
        region_type: RegionType,
        sequence_type: SequenceType,
        min_len: u32,
        max_len: u32,
        sequence: &str,
    ) -> Region {
        Region {
            region_id: id.to_string(),
            region_type,
            name: id.to_string(),
            sequence_type,
            sequence: sequence.to_string(),
            min_len,
            max_len,
            onlist: None,
            subregions: vec![],
        }
    }

    #[test]
    fn test_end_regions() {
        let mut end_regions = EndRegions::new(EndType::FivePrime);
        
        // Add a barcode region
        let barcode_region = Arc::new(RwLock::new(create_test_region(
            "bc1",
            RegionType::Barcode,
            SequenceType::Onlist,
            8,
            8,
            "",
        )));
        
        end_regions.add_region(barcode_region);
        
        assert!(end_regions.has_barcode);
        assert_eq!(end_regions.max_len, 8);
        assert_eq!(end_regions.calculate_cut_length(), 10); // 8 * 1.15 = 9.2 -> 10
    }

    #[test]
    fn test_find_innermost_regions() {
        // Create test regions
        let barcode1 = Arc::new(RwLock::new(create_test_region(
            "bc1",
            RegionType::Barcode,
            SequenceType::Onlist,
            8,
            8,
            "",
        )));
        
        let linker = Arc::new(RwLock::new(create_test_region(
            "linker",
            RegionType::Linker,
            SequenceType::Fixed,
            20,
            20,
            "AGATCGGAAGAGCGTCGTGT",
        )));
        
        let cdna = Arc::new(RwLock::new(create_test_region(
            "cdna",
            RegionType::Cdna,
            SequenceType::Random,
            100,
            1000,
            "",
        )));
        
        let barcode2 = Arc::new(RwLock::new(create_test_region(
            "bc2",
            RegionType::Barcode,
            SequenceType::Onlist,
            8,
            8,
            "",
        )));

        // Create modality region
        let modality_region = Region {
            region_id: "rna".to_string(),
            region_type: RegionType::Modality(Modality::RNA),
            name: "RNA".to_string(),
            sequence_type: SequenceType::Joined,
            sequence: "".to_string(),
            min_len: 0,
            max_len: 0,
            onlist: None,
            subregions: vec![barcode1, linker, cdna, barcode2],
        };

        let lib_spec = LibSpec::new(vec![modality_region]).unwrap();
        let (five_prime, three_prime) = find_innermost_regions(&lib_spec, &Modality::RNA).unwrap();

        // 5' end should contain barcode1 and linker
        assert_eq!(five_prime.regions.len(), 2);
        assert!(five_prime.has_barcode);
        
        // 3' end should contain barcode2
        assert_eq!(three_prime.regions.len(), 1);
        assert!(three_prime.has_barcode);
    }
}