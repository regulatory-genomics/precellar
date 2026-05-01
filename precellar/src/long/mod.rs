use anyhow::{anyhow, Result};
use std::sync::{Arc, RwLock};

pub mod barcode_extractor;
pub mod barcode_index;
pub mod sequence_aligner;

pub use barcode_extractor::BarcodeExtractor;
pub use barcode_index::BarcodeIndex;
use seqspec::region::{LibSpec, Region};
use sequence_aligner::{CompositePattern, CompositeRegionSpan};

/// End type information for the 5' or 3' side of the designed library structure.
#[derive(Debug, Copy, Clone, PartialEq)]
pub enum EndType {
    FivePrime,
    ThreePrime,
}

/// Region collection information containing all relevant regions for barcode extraction in one end
#[derive(Debug, Clone)]
pub struct EndRegions {
    /// End type
    pub end_type: EndType,
    /// Region list in forward order of the designed library structure (5' -> 3').
    /// In end-segment coordinates, left = 5' side and right = 3' side.
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

    /// Calculate sampled-window length (sum of max_len + 15% buffer)
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

    /// Build a composite pattern from regions between the leftmost and rightmost fixed regions.
    /// Fixed regions contribute their actual sequence.
    /// Non-fixed regions between them contribute N's of average length.
    /// If only one fixed region exists, the composite contains just that fixed sequence.
    pub fn build_composite_pattern(&self) -> CompositePattern {
        // Find the leftmost and rightmost fixed region indices in forward order.
        let first_fixed_idx = self
            .regions
            .iter()
            .position(|r| r.read().unwrap().sequence_type.is_fixed());
        let last_fixed_idx = self
            .regions
            .iter()
            .rposition(|r| r.read().unwrap().sequence_type.is_fixed());

        let (start_idx, end_idx) = match (first_fixed_idx, last_fixed_idx) {
            (Some(first), Some(last)) => (first, last),
            _ => {
                return CompositePattern {
                    pattern: Vec::new(),
                    spans: Vec::new(),
                    pos_to_span: Vec::new(),
                }
            }
        };

        let mut pattern = Vec::new();
        let mut spans = Vec::new();

        for idx in start_idx..=end_idx {
            let region = &self.regions[idx];
            let region_guard = region.read().unwrap();
            let start = pattern.len();
            let is_fixed = region_guard.sequence_type.is_fixed();

            if is_fixed {
                pattern.extend_from_slice(region_guard.sequence.as_bytes());
            } else {
                let avg_len = ((region_guard.min_len + region_guard.max_len) / 2) as usize;
                let spacer_len = avg_len.max(if region_guard.max_len > 0 { 1 } else { 0 });
                pattern.extend(std::iter::repeat(b'N').take(spacer_len));
            }

            let end = pattern.len();
            spans.push(CompositeRegionSpan {
                region: region.clone(),
                pattern_start: start,
                pattern_end: end,
                is_spacer: !is_fixed,
            });
        }

        // Build pos_to_span lookup
        let mut pos_to_span = vec![0usize; pattern.len()];
        for (span_idx, span) in spans.iter().enumerate() {
            for p in span.pattern_start..span.pattern_end {
                pos_to_span[p] = span_idx;
            }
        }

        CompositePattern {
            pattern,
            spans,
            pos_to_span,
        }
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
}

/// Collect the 5' and 3' end-region windows used for long-read barcode extraction.
///
/// Both returned collections stay in forward order of the designed library structure (5' -> 3').
/// The 5' collection includes all regions from the library start through the last fixed/barcode
/// before the target. The 3' collection includes all regions from the first fixed/barcode after
/// the target through the library end.
pub fn collect_end_regions(
    lib_spec: &LibSpec,
    modality: &seqspec::Modality,
) -> Result<(EndRegions, EndRegions)> {
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

    // Traverse from the 5' side and keep all regions through the last fixed/barcode before the target.
    let mut last_5p_anchor_idx = None;
    for (idx, region) in subregions.iter().enumerate() {
        let region_guard = region.read().unwrap();
        let region_type = &region_guard.region_type;

        // Stop if encounter cDNA or gDNA
        if region_type.is_target() {
            break;
        }

        // Track the last fixed or barcode region
        if region_guard.sequence_type.is_fixed() || region_type.is_barcode() {
            last_5p_anchor_idx = Some(idx);
        }
    }

    // Keep all regions from the 5' side through the last fixed/barcode before the target.
    if let Some(last_idx) = last_5p_anchor_idx {
        for (idx, region) in subregions.iter().enumerate() {
            let region_guard = region.read().unwrap();
            let region_type = &region_guard.region_type;

            // Stop if encounter cDNA or gDNA
            if region_type.is_target() {
                break;
            }

            // Add all regions up to and including the last fixed/barcode before the target.
            if idx <= last_idx {
                drop(region_guard);
                five_prime_regions.add_region(region.clone());
            }
        }
    }

    // Traverse in forward order, then keep all regions from the first fixed/barcode after the
    // target through the 3' side of the library structure.
    let mut seen_target = false;
    let mut first_3p_idx = None;
    for (idx, region) in subregions.iter().enumerate() {
        let region_guard = region.read().unwrap();
        let region_type = &region_guard.region_type;

        if region_type.is_target() {
            seen_target = true;
            continue;
        }

        if seen_target && (region_guard.sequence_type.is_fixed() || region_type.is_barcode()) {
            first_3p_idx = Some(idx);
            break;
        }
    }

    if let Some(first_idx) = first_3p_idx {
        for (idx, region) in subregions.iter().enumerate() {
            if idx >= first_idx {
                three_prime_regions.add_region(region.clone());
            }
        }
    }

    Ok((five_prime_regions, three_prime_regions))
}

/// Collect the flanking non-target regions for 5' and 3' ends, used for read trimming.
pub(crate) fn collect_target_flanks(
    lib_spec: &LibSpec,
    modality: &seqspec::Modality,
) -> Result<(Vec<Arc<RwLock<Region>>>, Vec<Arc<RwLock<Region>>>)> {
    let modality_region = lib_spec
        .get_modality(modality)
        .ok_or_else(|| anyhow!("Cannot find specified modality: {:?}", modality))?;

    let modality_guard = modality_region.read().unwrap();
    let subregions = &modality_guard.subregions;

    if subregions.is_empty() {
        return Err(anyhow!("Modality region has no subregions"));
    }

    let mut five_prime_flank = Vec::new();
    let mut three_prime_flank = Vec::new();
    let mut seen_target = false;

    for region in subregions {
        let region_type = region.read().unwrap().region_type.clone();
        if region_type.is_target() {
            seen_target = true;
            continue;
        }
        if seen_target {
            three_prime_flank.push(region.clone());
        } else {
            five_prime_flank.push(region.clone());
        }
    }

    Ok((five_prime_flank, three_prime_flank))
}

#[cfg(test)]
mod tests {
    use super::*;
    use indexmap::{IndexMap, IndexSet};
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
    fn test_build_composite_pattern() {
        let mut end_regions = EndRegions::new(EndType::FivePrime);

        // fixed1(16bp) -> barcode(10bp) -> umi(4~8bp) -> fixed2(16bp)
        let fixed1 = Arc::new(RwLock::new(create_test_region(
            "fixed1",
            RegionType::Linker,
            SequenceType::Fixed,
            16,
            16,
            "ACGTACGTACGTACGT",
        )));
        let barcode = Arc::new(RwLock::new(create_test_region(
            "bc1",
            RegionType::Barcode,
            SequenceType::Onlist,
            10,
            10,
            "",
        )));
        let umi = Arc::new(RwLock::new(create_test_region(
            "umi1",
            RegionType::Umi,
            SequenceType::Random,
            4,
            8,
            "",
        )));
        let fixed2 = Arc::new(RwLock::new(create_test_region(
            "fixed2",
            RegionType::Linker,
            SequenceType::Fixed,
            16,
            16,
            "TGCATGCATGCATGCA",
        )));

        end_regions.add_region(fixed1);
        end_regions.add_region(barcode);
        end_regions.add_region(umi);
        end_regions.add_region(fixed2);

        let composite = end_regions.build_composite_pattern();
        assert_eq!(composite.pattern.len(), 48); // 16 + 10 + 6 + 16
        assert_eq!(composite.spans.len(), 4);
        assert!(!composite.spans[0].is_spacer);
        assert!(composite.spans[1].is_spacer);
        assert!(composite.spans[2].is_spacer);
        assert!(!composite.spans[3].is_spacer);

        // Verify N's in spacer region
        assert!(composite.pattern[16..32].iter().all(|&b| b == b'N'));

        // Verify fixed sequences
        assert_eq!(&composite.pattern[0..16], b"ACGTACGTACGTACGT");
        assert_eq!(&composite.pattern[32..48], b"TGCATGCATGCATGCA");
    }

    #[test]
    fn test_collect_end_regions() {
        // Test structure:
        // primer (random) -> barcode1 (onlist) -> linker (fixed) -> cdna (target)
        // -> umi3 (random) -> barcode2 (onlist) -> adapter3 (fixed)

        let primer = Arc::new(RwLock::new(create_test_region(
            "primer",
            RegionType::IlluminaP5,
            SequenceType::Random,
            30,
            60,
            "",
        )));

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

        let umi3 = Arc::new(RwLock::new(create_test_region(
            "umi3",
            RegionType::Umi,
            SequenceType::Random,
            10,
            12,
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

        let adapter3 = Arc::new(RwLock::new(create_test_region(
            "adapter3",
            RegionType::Linker,
            SequenceType::Fixed,
            18,
            18,
            "TTAACCGGTTAACCGGTT",
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
            subregions: vec![primer, barcode1, linker, cdna, umi3, barcode2, adapter3],
        };

        let lib_spec = LibSpec::new(vec![modality_region]).unwrap();
        let (five_prime, three_prime) = collect_end_regions(&lib_spec, &Modality::RNA).unwrap();

        // 5' end should contain primer, barcode1, and linker (all regions through the
        // last fixed/barcode before the target).
        assert_eq!(five_prime.regions.len(), 3);
        assert!(five_prime.has_barcode);

        // Check that the sampled-window length includes the primer
        // max_len sum: 60 + 8 + 20 = 88
        // With 15% buffer: 88 * 1.15 = 101.2 -> 102
        assert_eq!(five_prime.max_len, 88);
        assert_eq!(five_prime.calculate_cut_length(), 102);

        // 3' end should skip the random region right after target, then keep all regions
        // from the first fixed/barcode through the library 3' side in forward order.
        let three_prime_region_ids: Vec<_> = three_prime
            .regions
            .iter()
            .map(|r| r.read().unwrap().region_id.clone())
            .collect();
        assert_eq!(three_prime_region_ids, vec!["bc2", "adapter3"]);
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
    
        let extractor = BarcodeExtractor::new(&lib_spec, &Modality::RNA, whitelists).unwrap();

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
        let result1 = extractor.extract_barcode(&record1).unwrap();
        let result2 = extractor.extract_barcode(&record2).unwrap();
        let result3 = extractor.extract_barcode(&record3).unwrap();
        let result4 = extractor.extract_barcode(&record4).unwrap();

        // Verify results
        // Record 1 should successfully extract barcode1
        assert!(result1.is_success(), "Record 1 should successfully extract barcode");
        assert!(result1.barcode.is_some(), "Record 1 should have extracted barcode");
        assert_eq!(result1.barcode.unwrap(), barcode1.as_bytes(), "Record 1 should match barcode1");
        
        // Record 2 should successfully extract barcode2 (despite deletion)
        assert!(result2.is_success(), "Record 2 should successfully extract barcode");
        assert!(result2.barcode.is_some(), "Record 2 should have extracted barcode");
        assert_eq!(result2.barcode.unwrap(), barcode2.as_bytes(), "Record 2 should match barcode2");
        
        // Record 3 should successfully extract barcode3 (despite mismatches)
        assert!(result3.is_success(), "Record 3 should successfully extract barcode");
        assert!(result3.barcode.is_some(), "Record 3 should have extracted barcode");
        assert_eq!(result3.barcode.unwrap(), barcode3.as_bytes(), "Record 3 should match barcode3");
        
        // Record 4 should fail to extract any barcode
        assert!(!result4.is_success(), "Record 4 should fail to extract barcode");
        assert!(result4.barcode.is_none(), "Record 4 should not have extracted barcode");
        
        println!("All barcode extraction tests with indels/mismatches passed!");
    }
}