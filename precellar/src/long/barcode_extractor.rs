use anyhow::Result;
use log::warn;
use std::collections::HashMap;
use std::sync::{Arc, RwLock};

use crate::barcode::OligoFrequncy;
use seqspec::region::{LibSpec, Region};
use seqspec::Modality;

use super::{
    find_innermost_regions, EndRegions, LongReadBarcodeResult,
    sequence_aligner::{FittingAligner, FixedSequenceAlignment, find_non_overlapping_alignments, fitting_alignment_distance}
};

/// Successfully extracted barcode with metadata
#[derive(Debug, Clone)]
pub struct ExtractedBarcode {
    pub region_id: String,
    pub barcode: Vec<u8>,  // Standardized barcode sequence from whitelist
    pub confidence: f64,   // confidence = 1.0 - (min_edit_distance / barcode_length)
}

/// Anchor with position information
#[derive(Debug, Clone)]
pub struct AnchorWithPosition {
    pub position: usize,  // Position in EndRegions (1-based)
    pub alignment: FixedSequenceAlignment,
}

/// Find and validate fixed sequence alignments
pub struct AnchorFinder {
    aligner: FittingAligner,
}

impl AnchorFinder {
    /// Create a new anchor finder
    pub fn new() -> Result<Self> {
        Ok(Self {
            aligner: FittingAligner::new()?,
        })
    }

    /// Find anchor positions using fixed sequence alignment
    pub fn find_anchors(
        &mut self,
        sequence: &[u8],
        fixed_regions: &[Arc<RwLock<Region>>],
        end_regions: &EndRegions,
    ) -> Result<Vec<AnchorWithPosition>> {
        if fixed_regions.is_empty() {
            return Ok(Vec::new());
        }

        let alignments = self.aligner.align_fixed_sequences(sequence, fixed_regions)?;
        // Greedy algorithm: save high score fixed alignments and discard lower overlapping ones
        let non_overlapping = find_non_overlapping_alignments(alignments);
        
        // Build position map from EndRegions
        let position_map = end_regions.build_position_map();
        
        // Convert alignments to AnchorWithPosition
        let mut anchors_with_position = Vec::new();
        for alignment in non_overlapping {
            let region_guard = alignment.region.read().unwrap();
            if let Some(&position) = position_map.get(&region_guard.region_id) {
                drop(region_guard);
                anchors_with_position.push(AnchorWithPosition {
                    position,
                    alignment,
                });
            }
        }
        
        Ok(anchors_with_position)
    }

    /// Validate that fixed region alignments are in correct order (outside to inside)
    pub fn validate_order(
        &self,
        anchors: &[AnchorWithPosition],
    ) -> Result<bool> {
        if anchors.len() <= 1 {
            return Ok(true); // Single or no alignment is always valid
        }

        // Sort anchors by their position (outside to inside)
        let mut sorted_anchors = anchors.to_vec();
        sorted_anchors.sort_by_key(|anchor| anchor.position);

        // Check that alignments are in correct positional order
        for i in 0..sorted_anchors.len() - 1 {
            let current = &sorted_anchors[i].alignment;
            let next = &sorted_anchors[i + 1].alignment;
            
            // For outside-to-inside order, the end position of outer region
            // should be less than or equal to start position of inner region
            if current.query_end > next.query_start {
                let current_region = current.region.read().unwrap();
                let next_region = next.region.read().unwrap();
                warn!(
                    "Fixed region order violation: {} (end: {}) overlaps with {} (start: {})",
                    current_region.region_id, current.query_end,
                    next_region.region_id, next.query_start
                );
                return Ok(false);
            }
        }

        Ok(true)
    }
}

/// Barcode locator component responsible for determining barcode extraction ranges
pub struct BarcodeLocator;

impl BarcodeLocator {
    /// Create a new barcode locator
    pub fn new() -> Self {
        Self
    }

    /// Locate the potential barcode window using position-based logic
    pub fn locate(
        &self,
        barcode_region: &Arc<RwLock<Region>>,
        anchors: &[AnchorWithPosition],
        end_regions: &EndRegions,
        sequence_length: usize,
    ) -> Option<(usize, usize)> {
        let region_guard = barcode_region.read().unwrap();
        let barcode_length = region_guard.max_len as usize; // For barcode regions, the max_len should be equal to the min_len
        let region_id = region_guard.region_id.clone();
        drop(region_guard);

        // Get barcode position using position map
        let position_map = end_regions.build_position_map();
        let barcode_pos = position_map.get(&region_id).copied()?;

        // Find closest outer and inner anchors in one pass
        let mut closest_outer: Option<&AnchorWithPosition> = None;
        let mut closest_inner: Option<&AnchorWithPosition> = None;
        
        for anchor in anchors {
            if anchor.position < barcode_pos {
                // Outer anchor: find the one with highest position (closest to barcode)
                if closest_outer.is_none() || anchor.position > closest_outer.unwrap().position {
                    closest_outer = Some(anchor);
                }
            } else if anchor.position > barcode_pos {
                // Inner anchor: find the one with lowest position (closest to barcode)
                if closest_inner.is_none() || anchor.position < closest_inner.unwrap().position {
                    closest_inner = Some(anchor);
                }
            }
        }
        
        // Determine extraction range based on adjacency and anchor positions
        self.determine_extraction_range(
            barcode_pos,
            barcode_length,
            closest_outer,
            closest_inner,
            sequence_length,
        )
    }


    /// Determine extraction range based on adjacency and anchor positions
    fn determine_extraction_range(
        &self,
        barcode_pos: usize,
        barcode_length: usize,
        closest_outer: Option<&AnchorWithPosition>,
        closest_inner: Option<&AnchorWithPosition>,
        sequence_length: usize,
    ) -> Option<(usize, usize)> {

        // Check if inner anchor is adjacent (position = barcode_pos + 1)
        if let Some(inner_anchor) = closest_inner {
            if inner_anchor.position == barcode_pos + 1 {
                // Inner anchor is adjacent - extract before it
                return self.extract_with_adjacent_inner(
                    closest_outer,
                    inner_anchor,
                    barcode_length,
                    sequence_length,
                );
            }
        }

        // Check if outer anchor is adjacent (position = barcode_pos - 1)
        if let Some(outer_anchor) = closest_outer {
            if outer_anchor.position == barcode_pos - 1 {
                // Outer anchor is adjacent - extract after it
                return self.extract_with_adjacent_outer(
                    outer_anchor,
                    closest_inner,
                    barcode_length,
                    sequence_length,
                );
            }
        }

        // Neither is adjacent - extract between nearest anchors
        self.extract_between_anchors(
            closest_outer,
            closest_inner,
            barcode_length,
            sequence_length,
        )
    }

    /// Extract barcode with adjacent inner anchor
    fn extract_with_adjacent_inner(
        &self,
        closest_outer: Option<&AnchorWithPosition>,
        inner_anchor: &AnchorWithPosition,
        barcode_length: usize,
        sequence_length: usize,
    ) -> Option<(usize, usize)> {
        let target_length = (barcode_length as f64 * 1.2).ceil() as usize;
        let min_length = (barcode_length as f64 * 0.8).floor() as usize;
        
        let end = inner_anchor.alignment.query_start;
        let mut start = if end >= target_length {
            end - target_length
        } else {
            0
        };

        // Adjust start to avoid overlap with outer anchor
        if let Some(outer_anchor) = closest_outer {
            start = start.max(outer_anchor.alignment.query_end);
        }

        let actual_length = end.saturating_sub(start); // avoid negative length
        if actual_length >= min_length && end <= sequence_length {
            Some((start, end))
        } else {
            None
        }
    }

    /// Extract barcode with adjacent outer anchor  
    fn extract_with_adjacent_outer(
        &self,
        outer_anchor: &AnchorWithPosition,
        closest_inner: Option<&AnchorWithPosition>,
        barcode_length: usize,
        sequence_length: usize,
    ) -> Option<(usize, usize)> {
        let target_length = (barcode_length as f64 * 1.2).ceil() as usize;
        let min_length = (barcode_length as f64 * 0.8).floor() as usize;
        
        let start = outer_anchor.alignment.query_end;
        let mut end = start + target_length;

        // Adjust end to avoid overlap with inner anchor or sequence end
        if let Some(inner_anchor) = closest_inner {
            end = end.min(inner_anchor.alignment.query_start);
        }
        end = end.min(sequence_length);

        let actual_length = end.saturating_sub(start);
        if actual_length >= min_length {
            Some((start, end))
        } else {
            None
        }
    }

    /// Extract barcode between nearest anchors when neither is adjacent
    fn extract_between_anchors(
        &self,
        closest_outer: Option<&AnchorWithPosition>,
        closest_inner: Option<&AnchorWithPosition>,
        barcode_length: usize,
        sequence_length: usize,
    ) -> Option<(usize, usize)> {
        let min_length = (barcode_length as f64 * 0.8).floor() as usize;
        
        let start = closest_outer
            .map(|anchor| anchor.alignment.query_end)
            .unwrap_or(0);
        
        let end = closest_inner
            .map(|anchor| anchor.alignment.query_start)
            .unwrap_or(sequence_length);

        let available_length = end.saturating_sub(start);
        
        if available_length >= min_length {
            Some((start, end))
        } else {
            None
        }
    }

}

/// Find the best barcode match using fitting alignment distance
fn find_best_barcode_match(
    candidate_seq: &[u8],
    barcode_region: &Arc<RwLock<Region>>,
    whitelists: &HashMap<String, OligoFrequncy>,
) -> Option<ExtractedBarcode> {
    let region_guard = barcode_region.read().unwrap();
    let region_id = region_guard.region_id.clone();
    drop(region_guard);

    // Get whitelist for this region
    let whitelist = whitelists.get(&region_id)?;
    
    // Find best match using fitting alignment distance
    let (matched_barcode, score) = find_best_fitting_match(candidate_seq, whitelist)?;

    Some(ExtractedBarcode {
        region_id,
        barcode: matched_barcode,
        confidence: score,
    })
}

/// Find best barcode match using fitting alignment distance
fn find_best_fitting_match(
    candidate_seq: &[u8],
    whitelist: &OligoFrequncy,
) -> Option<(Vec<u8>, f64)> {
    if whitelist.is_empty() {
        return None;
    }

    let mut best_barcode = None;
    let mut min_distance = usize::MAX;

    for barcode in whitelist.keys() {
        // Use fitting alignment distance: short sequence (barcode) vs long sequence (candidate)
        let distance = fitting_alignment_distance(barcode, candidate_seq);
        
        if distance < min_distance {
            min_distance = distance;
            best_barcode = Some(barcode.clone());
        }
    }

    if let Some(barcode) = best_barcode {
        // Convert edit distance to confidence score
        let barcode_length = barcode.len().max(1);
        let confidence = 1.0 - (min_distance as f64 / barcode_length as f64);
        
        // Only return matches with reasonable confidence
        if confidence >= 0.7 {
            Some((barcode, confidence))
        } else {
            None
        }
    } else {
        None
    }
}

/// Main barcode extractor for long reads.
/// Acts as an orchestrator, delegating tasks to specialized components.
pub struct BarcodeExtractor {
    five_prime_regions: EndRegions,
    three_prime_regions: EndRegions,
    barcode_locator: BarcodeLocator,
}

impl BarcodeExtractor {
    /// Create a new barcode extractor for the given library_spec and modality.
    pub fn new(lib_spec: &LibSpec, modality: &Modality) -> Result<Self> {
        let (five_prime_regions, three_prime_regions) = find_innermost_regions(lib_spec, modality)?;
        
        Ok(Self {
            five_prime_regions,
            three_prime_regions,
            barcode_locator: BarcodeLocator::new(),
        })
    }

    /// Extract barcode from a FASTQ record.
    pub fn extract_barcode(
        &self,
        record: &noodles::fastq::Record,
        whitelists: &HashMap<String, OligoFrequncy>,
    ) -> Result<LongReadBarcodeResult> {
        let sequence = record.sequence();
        let quality = record.quality_scores();

        // Step 1: Cut segments from both ends
        let five_prime_segment = self.cut_end_segment(sequence, quality, &self.five_prime_regions, true)?;
        let three_prime_segment = self.cut_end_segment(sequence, quality, &self.three_prime_regions, false)?;

        let mut extracted_barcodes = Vec::new();

        // Step 2 & 3: Process 5' end
        if let Some((seq, _qual)) = five_prime_segment {
            let barcode_results = self
                .process_end_for_barcode(&seq, &self.five_prime_regions, whitelists)?;
            extracted_barcodes.extend(barcode_results);
        }

        // Step 2 & 3: Process 3' end  
        if let Some((seq, _qual)) = three_prime_segment {
            let barcode_results = self
                .process_end_for_barcode(&seq, &self.three_prime_regions, whitelists)?;
            extracted_barcodes.extend(barcode_results);
        }

        // Step 4: Combine multiple barcodes if present
        let final_result = self.combine_barcodes(extracted_barcodes)?;
        Ok(final_result)
    }

    /// Process one end to extract all valid barcodes.
    /// Returns a vector of all successfully extracted barcodes from this end.
    fn process_end_for_barcode(
        &self,
        segment_seq: &[u8],
        end_regions: &EndRegions,
        whitelists: &HashMap<String, OligoFrequncy>,
    ) -> Result<Vec<ExtractedBarcode>> {
        // Return early if no barcode region is present on this end
        if !end_regions.has_barcode {
            return Ok(Vec::new());
        }

        let fixed_regions = end_regions.get_fixed_regions();
        let barcode_regions = end_regions.get_barcode_regions();

        if barcode_regions.is_empty() {
            return Ok(Vec::new());
        }

        // Step 1: Find anchor points using the AnchorFinder
        let mut anchor_finder = AnchorFinder::new()?;
        let anchor_positions = anchor_finder.find_anchors(segment_seq, &fixed_regions, end_regions)?;
        
        // Step 2: Validate alignment order
        if !anchor_positions.is_empty() {
            if !anchor_finder.validate_order(&anchor_positions)? {
                warn!("Fixed region alignment order is incorrect, skipping this read");
                return Ok(Vec::new());
            }
        }

        // Step 3: Extract all valid barcodes from all barcode regions
        let mut extracted_barcodes = Vec::new();

        for barcode_region in barcode_regions {
            // Step 3a: Locate the potential barcode sequence range
            if let Some(extraction_range) = self.barcode_locator.locate(
                &barcode_region,
                &anchor_positions,
                end_regions,
                segment_seq.len(),
            ) {
                // BarcodeLocator already validates minimum length requirement via `locate` method
                let candidate_seq = &segment_seq[extraction_range.0..extraction_range.1];

                // Step 3b: Match the candidate against the whitelist
                if let Some(matched) = find_best_barcode_match(
                    candidate_seq,
                    &barcode_region,
                    whitelists,
                ) {
                    extracted_barcodes.push(matched);
                }
            }
        }
        
        Ok(extracted_barcodes)
    }

    /// Cut end segments from one fastq record.
    fn cut_end_segment(
        &self,
        sequence: &[u8],
        quality: &[u8],
        end_regions: &EndRegions,
        is_five_prime: bool,
    ) -> Result<Option<(Vec<u8>, Vec<u8>)>> {
        let cut_length = end_regions.calculate_cut_length();
        
        if sequence.len() < cut_length {
            // If sequence is shorter than cut length, report error and skip this fastq record
            warn!(
                "Sequence length ({}) is shorter than required cut length ({}), skipping record",
                sequence.len(), cut_length
            );
            return Ok(None);
        }

        if sequence.len() > 300 {
            // If segment is too long (>300bp), skip this record
            warn!("Sequence length ({}) is longer than 300bp, skipping record", sequence.len());
            return Ok(None);
        }

        let (seq_segment, qual_segment) = if is_five_prime {
            // Cut from the beginning (5' end)
            (&sequence[..cut_length], &quality[..cut_length])
        } else {
            // Cut from the end (3' end)
            let start = sequence.len() - cut_length;
            (&sequence[start..], &quality[start..])
        };

        Ok(Some((seq_segment.to_vec(), qual_segment.to_vec())))
    }

    /// Combine multiple extracted barcodes into final result.
    fn combine_barcodes(&self, barcodes: Vec<ExtractedBarcode>) -> Result<LongReadBarcodeResult> {
        if barcodes.is_empty() {
            return Ok(LongReadBarcodeResult {
                barcode: None,
                confidence: 0.0,
                success: false,
            });
        }

        if barcodes.len() == 1 {
            let barcode = &barcodes[0];
            return Ok(LongReadBarcodeResult {
                barcode: Some(barcode.barcode.clone()),
                confidence: barcode.confidence,
                success: true,
            });
        }

        // Multiple barcodes - concatenate them
        let mut combined_barcode = Vec::new();
        let mut total_confidence = 0.0;

        for barcode in &barcodes {
            combined_barcode.extend(&barcode.barcode);
            total_confidence += barcode.confidence;
        }

        let average_confidence = total_confidence / barcodes.len() as f64;

        Ok(LongReadBarcodeResult {
            barcode: Some(combined_barcode),
            confidence: average_confidence,
            success: true,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::sequence_aligner::fitting_alignment_distance;

    #[test]
    fn test_fitting_alignment_distance() {
        // Test with exact matches
        assert_eq!(fitting_alignment_distance(b"ATCG", b"ATCG"), 0);
        
        // Test with partial matches (short sequence found in long sequence)
        assert_eq!(fitting_alignment_distance(b"ATC", b"ATCA"), 0);
        assert_eq!(fitting_alignment_distance(b"ATCG", b"GGATCGCC"), 0);
        
        // Test with empty sequences
        assert_eq!(fitting_alignment_distance(b"", b"ATCG"), 0);
        assert_eq!(fitting_alignment_distance(b"ATCG", b""), 4);
    }

    #[test]
    fn test_barcode_locator_integration() {
        use seqspec::{RegionType, SequenceType};
        use seqspec::region::Region;
        use super::super::EndRegions;
        
        let locator = BarcodeLocator::new();
        
        // Create test barcode region
        let barcode_region = Arc::new(RwLock::new(Region {
            region_id: "bc1".to_string(),
            region_type: RegionType::Barcode,
            name: "bc1".to_string(),
            sequence_type: SequenceType::Onlist,
            sequence: "".to_string(),
            min_len: 8,
            max_len: 8,
            onlist: None,
            subregions: vec![],
        }));
        
        // Create test outer anchor region  
        let outer_region = Arc::new(RwLock::new(Region {
            region_id: "outer".to_string(),
            region_type: RegionType::Linker,
            name: "outer".to_string(),
            sequence_type: SequenceType::Fixed,
            sequence: "ATCG".to_string(),
            min_len: 4,
            max_len: 4,
            onlist: None,
            subregions: vec![],
        }));
        
        // Create test inner anchor region
        let inner_region = Arc::new(RwLock::new(Region {
            region_id: "inner".to_string(),
            region_type: RegionType::Linker,
            name: "inner".to_string(),
            sequence_type: SequenceType::Fixed,
            sequence: "GCTA".to_string(),
            min_len: 4,
            max_len: 4,
            onlist: None,
            subregions: vec![],
        }));
        
        // Create EndRegions with proper order: outer -> barcode -> inner
        let mut end_regions = EndRegions::new(super::super::EndType::FivePrime);
        end_regions.add_region(outer_region.clone());
        end_regions.add_region(barcode_region.clone());
        end_regions.add_region(inner_region.clone());
        
        // Create alignments of fixed sequences
        let alignment1 = FixedSequenceAlignment {
            region: outer_region,
            query_start: 2,
            query_end: 6,
            score: 1.0,
            matches: 4,
            alignment_length: 4,
        };
        
        let alignment2 = FixedSequenceAlignment {
            region: inner_region,
            query_start: 20,
            query_end: 24,
            score: 1.0,
            matches: 4,
            alignment_length: 4,
        };
        
        let anchor1 = AnchorWithPosition {
            position: 1, // outer position
            alignment: alignment1,
        };
        
        let anchor2 = AnchorWithPosition {
            position: 3, // inner position
            alignment: alignment2,
        };
        
        let anchors = vec![anchor1, anchor2];
        
        // Test locate method
        let result = locator.locate(&barcode_region, &anchors, &end_regions, 100);
        
        // Should find a range between the outer anchor end (6) and inner anchor start (20)
        assert!(result.is_some());
        let (start, end) = result.unwrap();
        assert!(start >= 6); // After outer anchor
        assert!(end <= 20);  // Before inner anchor
        assert!(end > start); // Valid range

        println!("start: {}, end: {}", start, end);
    }
}