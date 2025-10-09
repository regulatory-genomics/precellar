use anyhow::{anyhow, bail, Result};
use indexmap::{IndexMap, IndexSet};
use std::collections::HashMap;
use std::sync::{Arc, RwLock};

pub mod barcode_extractor;
pub mod sequence_aligner;

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
    /// Barcode whitelists: String for region_id
    whitelists: IndexMap<String, IndexSet<Vec<u8>>>,
}

impl LongReadProcessor {
    /// Create a new long-read processor
    pub fn new(
        whitelists: IndexMap<String, IndexSet<Vec<u8>>>,
    ) -> Result<Self> {
        // Check that all barcode regions have non-empty whitelists
        for (region_id, whitelist) in &whitelists {
            if whitelist.is_empty() {
                bail!(
                    "Barcode region '{}' does not have a whitelist. Long-read processing requires all barcode regions to have whitelists.", 
                    region_id
                );
            }
        }
        
        Ok(Self {
            whitelists,
        })
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

    #[test]
    fn test_extraction_with_fastq_record() {
        use noodles::fastq::Record;

        // Define fixed region sequences
        let fixed_seq1 = "ACGTACGTACGTACGT"; // 16bp fixed sequence 1
        let fixed_seq2 = "TGCATGCATGCATGCA"; // 16bp fixed sequence 2
        
        // Define barcode whitelist (3 valid barcodes)
        let barcode1 = "AAAAAAAAAA"; // 10bp
        let barcode2 = "CCCCCCCCCC"; // 10bp  
        let barcode3 = "GGGGGGGGGG"; // 10bp
        
        // Create barcode whitelist using IndexSet
        let mut whitelist = IndexSet::new();
        whitelist.insert(barcode1.as_bytes().to_vec());
        whitelist.insert(barcode2.as_bytes().to_vec());
        whitelist.insert(barcode3.as_bytes().to_vec());
        
        let mut whitelists = IndexMap::new();
        whitelists.insert("bc1".to_string(), whitelist);
        
        // Create test library specification
        let fixed_region1 = Arc::new(RwLock::new(create_test_region(
            "fixed1",
            RegionType::Linker,
            SequenceType::Fixed,
            16,
            16,
            fixed_seq1,
        )));
        
        let barcode_region = Arc::new(RwLock::new(create_test_region(
            "bc1",
            RegionType::Barcode,
            SequenceType::Onlist,
            10,
            10,
            "",
        )));
        
        let fixed_region2 = Arc::new(RwLock::new(create_test_region(
            "fixed2", 
            RegionType::Linker,
            SequenceType::Fixed,
            16,
            16,
            fixed_seq2,
        )));
        
        let cdna = Arc::new(RwLock::new(create_test_region(
            "cdna",
            RegionType::Cdna,
            SequenceType::Random,
            100,
            1000,
            "",
        )));

        let modality_region = Region {
            region_id: "rna".to_string(),
            region_type: RegionType::Modality(Modality::RNA),
            name: "RNA".to_string(),
            sequence_type: SequenceType::Joined,
            sequence: "".to_string(),
            min_len: 0,
            max_len: 0,
            onlist: None,
            subregions: vec![fixed_region1, barcode_region, fixed_region2, cdna],
        };

        let lib_spec = LibSpec::new(vec![modality_region]).unwrap();
        let processor = LongReadProcessor::new(whitelists).unwrap();

        // Test sequences with indels and mismatches
        // Record 1: Perfect match with barcode1, 1 mismatch in fixed region
        let seq1 = format!("{}{}{}ATCGATCGATCGATCGATCGATCGATCG", 
                          "ACGTACGTACGTACGG", // 1 mismatch in fixed1 (T->G)
                          barcode1,           // Perfect barcode1
                          fixed_seq2);        // Perfect fixed2
        let record1 = Record::new(
            noodles::fastq::record::Definition::new("read1", ""), 
            seq1.clone(), 
            "I".repeat(seq1.len())
        );

        // Record 2: 1 insertion in fixed region, 1 deletion in barcode2
        let seq2 = format!("{}{}{}ATCGATCGATCGATCGATCGATCGATCG",
                          "ACGTACGTACGTACGTA", // 1 insertion in fixed1 (extra A)
                          "CCCCCCCCC",         // 1 deletion in barcode2 (missing 1 C)
                          fixed_seq2);         // Perfect fixed2
        let record2 = Record::new(
            noodles::fastq::record::Definition::new("read2", ""), 
            seq2.clone(), 
            "I".repeat(seq2.len())
        );

        // Record 3: 2 mismatches in fixed regions, 1 mismatch in barcode3
        let seq3 = format!("{}{}{}ATCGATCGATCGATCGATCGATCGATCG",
                          "ACGTACGTACGTACGA", // 1 mismatch in fixed1 (T->A)
                          "GGGGGGGGTG",       // 1 mismatch in barcode3 (G->T)
                          "TGCATGCATGCATGCT");// 1 mismatch in fixed2 (A->T)
        let record3 = Record::new(
            noodles::fastq::record::Definition::new("read3", ""), 
            seq3.clone(), 
            "I".repeat(seq3.len())
        );

        // Record 4: Too many errors - should not match any barcode
        let seq4 = format!("{}{}{}ATCGATCGATCGATCGATCGATCGATCG",
                          "CCCACGTACGTACGTACGTGG", // Multiple insertion in fixed1
                          "TTTTTTTTTT",            // No match to any barcode
                          fixed_seq2);             // Perfect fixed2
        let record4 = Record::new(
            noodles::fastq::record::Definition::new("read4", ""), 
            seq4.clone(), 
            "I".repeat(seq4.len())
        );

        // Test barcode extraction
        let result1 = processor.extract_barcode(&record1, &lib_spec, &Modality::RNA).unwrap();
        let result2 = processor.extract_barcode(&record2, &lib_spec, &Modality::RNA).unwrap();
        let result3 = processor.extract_barcode(&record3, &lib_spec, &Modality::RNA).unwrap();
        let result4 = processor.extract_barcode(&record4, &lib_spec, &Modality::RNA).unwrap();

        // Verify results
        // Record 1 should successfully extract barcode1
        assert!(result1.success, "Record 1 should successfully extract barcode");
        assert!(result1.barcode.is_some(), "Record 1 should have extracted barcode");
        assert_eq!(result1.barcode.unwrap(), barcode1.as_bytes(), "Record 1 should match barcode1");
        
        // Record 2 should successfully extract barcode2 (despite deletion)
        assert!(result2.success, "Record 2 should successfully extract barcode");
        assert!(result2.barcode.is_some(), "Record 2 should have extracted barcode");
        assert_eq!(result2.barcode.unwrap(), barcode2.as_bytes(), "Record 2 should match barcode2");
        
        // Record 3 should successfully extract barcode3 (despite mismatches)
        assert!(result3.success, "Record 3 should successfully extract barcode");
        assert!(result3.barcode.is_some(), "Record 3 should have extracted barcode");
        assert_eq!(result3.barcode.unwrap(), barcode3.as_bytes(), "Record 3 should match barcode3");
        
        // Record 4 should fail to extract any barcode
        assert!(!result4.success, "Record 4 should fail to extract barcode");
        assert!(result4.barcode.is_none(), "Record 4 should not have extracted barcode");
        
        println!("All barcode extraction tests with indels/mismatches passed!");
    }
}