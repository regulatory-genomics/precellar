use anyhow::Result;
use indexmap::{IndexMap, IndexSet};
use log::warn;
use std::sync::{Arc, RwLock};

use seqspec::region::{LibSpec, Region};
use seqspec::Modality;

use super::{
    find_innermost_regions, EndRegions,
    sequence_aligner::{FittingAligner, FixedSequenceAlignment, find_non_overlapping_alignments, fitting_alignment_distance}
};

/// Final long-read barcode extraction result
#[derive(Debug, Clone)]
pub struct LongReadBarcodeResult {
    /// Extracted barcode sequence (standardized from whitelist)
    pub barcode: Option<Vec<u8>>,
    /// Extraction confidence (0.0 if no barcode extracted)
    pub confidence: f64,
    /// Whether the end segments were reverse-complemented during extraction
    pub is_reverse_complemented: bool,
    /// Whether some anchors were not detected in both orientations (for QC tracking)
    pub has_missing_anchors: bool,
}

impl LongReadBarcodeResult {
    /// Whether extraction was successful
    pub fn is_success(&self) -> bool {
        self.barcode.is_some()
    }
}

/// Successfully extracted barcode with metadata (single barcode)
#[derive(Debug, Clone)]
pub struct ExtractedBarcode {
    pub region_id: String,
    pub barcode: Vec<u8>,  // Standardized barcode sequence from whitelist
    pub confidence: f64,   // confidence = 1.0 - (min_edit_distance / barcode_length)
}

/// Evidence for read orientation determination, including found anchors for reuse
#[derive(Debug, Clone)]
pub struct OrientationEvidence {
    /// Number of anchors found
    pub num_anchors: usize,
    /// Whether anchors are in valid order (outside-to-inside)
    pub valid_order: bool,
    /// Average confidence of anchor alignments
    pub avg_confidence: f64,
}

impl OrientationEvidence {
    /// Check if evidence meets conservative thresholds for accepting orientation
    /// Conservative thresholds: ≥2 anchors, valid order, ≥0.8 confidence
    pub fn meets_threshold(&self) -> bool {
        self.num_anchors >= 2 && self.valid_order && self.avg_confidence >= 0.8
    }

    /// Create from anchor list (FixedSequenceAlignment with position set)
    pub fn from_anchors(anchors: &[FixedSequenceAlignment], valid_order: bool) -> Self {
        let num_anchors = anchors.len();
        let avg_confidence = if num_anchors > 0 {
            anchors.iter().map(|a| a.score).sum::<f64>() / num_anchors as f64
        } else {
            0.0
        };

        Self {
            num_anchors,
            valid_order,
            avg_confidence,
        }
    }
}

/// End segment with precomputed anchors for reuse
#[derive(Debug)]
struct EndSegmentWithAnchors {
    /// The cut segment sequence (empty if sequence was too short)
    sequence: Vec<u8>,
    /// The cut segment quality scores
    quality: Vec<u8>,
    /// Precomputed anchors from orientation detection (FixedSequenceAlignment with position set)
    anchors: Vec<FixedSequenceAlignment>,
    /// Whether anchors are in valid order
    valid_order: bool,
    /// Number of expected fixed regions (for tracking alignment failures)
    expected_fixed_count: usize,
}

impl EndSegmentWithAnchors {
    /// Create an empty segment (used when source sequence is too short)
    fn empty() -> Self {
        Self {
            sequence: Vec::new(),
            quality: Vec::new(),
            anchors: Vec::new(),
            valid_order: false,
            expected_fixed_count: 0,
        }
    }

    /// Check if this segment is empty (no sequence was extracted)
    fn is_empty(&self) -> bool {
        self.sequence.is_empty()
    }

    /// Check if all expected fixed regions were found
    fn all_fixed_found(&self) -> bool {
        self.anchors.len() >= self.expected_fixed_count
    }
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

    /// Find anchor positions using fixed sequence alignment.
    /// Returns FixedSequenceAlignment with position field set.
    pub fn find_anchors(
        &mut self,
        sequence: &[u8],
        fixed_regions: &[Arc<RwLock<Region>>],
        end_regions: &EndRegions,
    ) -> Result<Vec<FixedSequenceAlignment>> {
        if fixed_regions.is_empty() {
            return Ok(Vec::new());
        }

        let alignments = self.aligner.align_fixed_sequences(sequence, fixed_regions)?;
        // Greedy algorithm: save high score fixed alignments and discard lower overlapping ones
        let non_overlapping = find_non_overlapping_alignments(alignments);

        // Build position map from EndRegions
        let position_map = end_regions.build_position_map();

        // Set position on each alignment
        let mut anchors = Vec::new();
        for mut alignment in non_overlapping {
            let region_guard = alignment.region.read().unwrap();
            if let Some(&position) = position_map.get(&region_guard.region_id) {
                drop(region_guard);
                alignment.position = Some(position);
                anchors.push(alignment);
            }
        }

        Ok(anchors)
    }

    /// Validate that fixed region alignments are in correct order (outside to inside)
    pub fn validate_order(
        &self,
        anchors: &[FixedSequenceAlignment],
    ) -> Result<bool> {
        if anchors.len() <= 1 {
            return Ok(true); // Single or no alignment is always valid
        }

        // Sort anchors by their position (outside to inside)
        let mut sorted_anchors = anchors.to_vec();
        sorted_anchors.sort_by_key(|anchor| anchor.position);

        // Check that alignments are in correct positional order
        for i in 0..sorted_anchors.len() - 1 {
            let current = &sorted_anchors[i];
            let next = &sorted_anchors[i + 1];

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


/// Find the best barcode match using fitting alignment distance
fn find_best_barcode_match(
    candidate_seq: &[u8],
    barcode_region: &Arc<RwLock<Region>>,
    whitelists: &IndexMap<String, IndexSet<Vec<u8>>>,
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
    whitelist: &IndexSet<Vec<u8>>,
) -> Option<(Vec<u8>, f64)> {
    if whitelist.is_empty() {
        return None;
    }

    let mut best_barcode = None;
    let mut min_distance = usize::MAX;

    for barcode in whitelist.iter() {
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

/// Locate the potential barcode window using position-based logic
fn locate_barcode(
    barcode_region: &Arc<RwLock<Region>>,
    anchors: &[FixedSequenceAlignment],
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
    let mut closest_outer: Option<&FixedSequenceAlignment> = None;
    let mut closest_inner: Option<&FixedSequenceAlignment> = None;

    for anchor in anchors {
        let anchor_pos = anchor.position.unwrap_or(0);
        if anchor_pos < barcode_pos {
            // Outer anchor: find the one with highest position (closest to barcode)
            if closest_outer.is_none() || anchor_pos > closest_outer.unwrap().position.unwrap_or(0) {
                closest_outer = Some(anchor);
            }
        } else if anchor_pos > barcode_pos {
            // Inner anchor: find the one with lowest position (closest to barcode)
            if closest_inner.is_none() || anchor_pos < closest_inner.unwrap().position.unwrap_or(usize::MAX) {
                closest_inner = Some(anchor);
            }
        }
    }

    // Determine extraction range based on adjacency and anchor positions
    determine_extraction_range(
        barcode_pos,
        barcode_length,
        closest_outer,
        closest_inner,
        sequence_length,
    )
}

/// Determine extraction range based on adjacency and anchor positions
fn determine_extraction_range(
    barcode_pos: usize,
    barcode_length: usize,
    closest_outer: Option<&FixedSequenceAlignment>,
    closest_inner: Option<&FixedSequenceAlignment>,
    sequence_length: usize,
) -> Option<(usize, usize)> {
    // Check if inner anchor is adjacent (position = barcode_pos + 1)
    if let Some(inner_anchor) = closest_inner {
        if inner_anchor.position == Some(barcode_pos + 1) {
            // Inner anchor is adjacent - extract before it
            return extract_with_adjacent_inner(
                closest_outer,
                inner_anchor,
                barcode_length,
                sequence_length,
            );
        }
    }

    // Check if outer anchor is adjacent (position = barcode_pos - 1)
    if let Some(outer_anchor) = closest_outer {
        if outer_anchor.position == Some(barcode_pos - 1) {
            // Outer anchor is adjacent - extract after it
            return extract_with_adjacent_outer(
                outer_anchor,
                closest_inner,
                barcode_length,
                sequence_length,
            );
        }
    }

    // Neither is adjacent - extract between nearest anchors
    extract_between_anchors(
        closest_outer,
        closest_inner,
        barcode_length,
        sequence_length,
    )
}

/// Extract barcode with adjacent inner anchor
fn extract_with_adjacent_inner(
    closest_outer: Option<&FixedSequenceAlignment>,
    inner_anchor: &FixedSequenceAlignment,
    barcode_length: usize,
    sequence_length: usize,
) -> Option<(usize, usize)> {
    let target_length = (barcode_length as f64 * 1.2).ceil() as usize;
    let min_length = (barcode_length as f64 * 0.8).floor() as usize;

    let end = inner_anchor.query_start;
    let mut start = if end >= target_length {
        end - target_length
    } else {
        0
    };

    // Adjust start to avoid overlap with outer anchor
    if let Some(outer_anchor) = closest_outer {
        start = start.max(outer_anchor.query_end);
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
    outer_anchor: &FixedSequenceAlignment,
    closest_inner: Option<&FixedSequenceAlignment>,
    barcode_length: usize,
    sequence_length: usize,
) -> Option<(usize, usize)> {
    let target_length = (barcode_length as f64 * 1.2).ceil() as usize;
    let min_length = (barcode_length as f64 * 0.8).floor() as usize;

    let start = outer_anchor.query_end;
    let mut end = start + target_length;

    // Adjust end to avoid overlap with inner anchor or sequence end
    if let Some(inner_anchor) = closest_inner {
        end = end.min(inner_anchor.query_start);
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
    closest_outer: Option<&FixedSequenceAlignment>,
    closest_inner: Option<&FixedSequenceAlignment>,
    barcode_length: usize,
    sequence_length: usize,
) -> Option<(usize, usize)> {
    let min_length = (barcode_length as f64 * 0.8).floor() as usize;

    let start = closest_outer
        .map(|anchor| anchor.query_end)
        .unwrap_or(0);

    let end = closest_inner
        .map(|anchor| anchor.query_start)
        .unwrap_or(sequence_length);

    let available_length = end.saturating_sub(start);

    if available_length >= min_length {
        Some((start, end))
    } else {
        None
    }
}

/// Main barcode extractor for long reads.
/// Acts as an orchestrator, delegating tasks to specialized components.
#[derive(Debug)]
pub struct BarcodeExtractor {
    five_prime_regions: EndRegions,
    three_prime_regions: EndRegions,
    /// Barcode whitelists: String for region_id
    whitelists: IndexMap<String, IndexSet<Vec<u8>>>,
}

impl BarcodeExtractor {
    /// Create a new barcode extractor for the given library_spec, modality, and whitelists.
    pub fn new(
        lib_spec: &LibSpec, 
        modality: &Modality,
        whitelists: IndexMap<String, IndexSet<Vec<u8>>>,
    ) -> Result<Self> {
        // Check that all barcode regions have non-empty whitelists
        for (region_id, whitelist) in &whitelists {
            if whitelist.is_empty() {
                anyhow::bail!(
                    "Barcode region '{}' does not have a whitelist. Long-read processing requires all barcode regions to have whitelists.", 
                    region_id
                );
            }
        }

        let (five_prime_regions, three_prime_regions) = find_innermost_regions(lib_spec, modality)?;

        Ok(Self {
            five_prime_regions,
            three_prime_regions,
            whitelists,
        })
    }

    /// Extract barcode from a FASTQ record with automatic orientation detection.
    ///
    /// This method uses an adaptive forward-first strategy with anchor reuse:
    /// 1. Cut end segments and find anchors for forward orientation
    /// 2. If forward meets threshold, extract barcodes using cached anchors (fast path)
    /// 3. Otherwise, cut and RC end segments, find anchors for reverse orientation
    /// 4. Compare evidence, choose better orientation, extract barcodes using cached anchors
    ///
    /// IMPORTANT: Only end segments are processed and RC'd, not the full read.
    /// Anchors found during orientation detection are reused for barcode extraction.
    pub fn extract_barcode(
        &self,
        record: &noodles::fastq::Record,
    ) -> Result<LongReadBarcodeResult> {
        let sequence = record.sequence();
        let quality = record.quality_scores();

        // Step 1: Cut end segments and find anchors for forward orientation
        let forward_5p = self.cut_and_analyze_segment(sequence, quality, &self.five_prime_regions, true, false)?;
        let forward_3p = self.cut_and_analyze_segment(sequence, quality, &self.three_prime_regions, false, false)?;

        let forward_evidence = Self::combine_evidence(&forward_5p, &forward_3p);

        // Step 2: If forward meets threshold, use it directly (fast path)
        if forward_evidence.meets_threshold() {
            return self.extract_barcodes_from_segments(forward_5p, forward_3p, false);
        }

        // Step 3: Forward didn't meet threshold, try reverse orientation
        // Only RC end segments, not full read
        let reverse_5p = self.cut_and_analyze_segment(sequence, quality, &self.five_prime_regions, true, true)?;
        let reverse_3p = self.cut_and_analyze_segment(sequence, quality, &self.three_prime_regions, false, true)?;

        let reverse_evidence = Self::combine_evidence(&reverse_5p, &reverse_3p);

        // Step 4: Compare and choose better orientation
        let should_rc = self.compare_evidence(&forward_evidence, &reverse_evidence);

        // Check if both orientations failed to find all fixed regions (for QC tracking)
        let forward_all_found = forward_5p.all_fixed_found() && forward_3p.all_fixed_found();
        let reverse_all_found = reverse_5p.all_fixed_found() && reverse_3p.all_fixed_found();
        let has_missing_anchors = !forward_all_found && !reverse_all_found;

        #[cfg(debug_assertions)]
        if has_missing_anchors {
            let read_name = std::str::from_utf8(record.name()).unwrap_or("<invalid>");
            eprintln!(
                "[DEBUG] Read '{}': fixed region alignment failed in both orientations \
                (forward: {}/{} 5', {}/{} 3'; reverse: {}/{} 5', {}/{} 3')",
                read_name,
                forward_5p.anchors.len(), forward_5p.expected_fixed_count,
                forward_3p.anchors.len(), forward_3p.expected_fixed_count,
                reverse_5p.anchors.len(), reverse_5p.expected_fixed_count,
                reverse_3p.anchors.len(), reverse_3p.expected_fixed_count,
            );
            eprintln!(
                "[DEBUG] Failed FASTQ record:\n@{}\n{}\n+\n{}",
                read_name,
                std::str::from_utf8(sequence).unwrap_or("<invalid seq>"),
                std::str::from_utf8(quality).unwrap_or("<invalid qual>"),
            );
        }

        let mut result = if should_rc {
            self.extract_barcodes_from_segments(reverse_5p, reverse_3p, true)?
        } else {
            self.extract_barcodes_from_segments(forward_5p, forward_3p, false)?
        };
        result.has_missing_anchors = has_missing_anchors;
        Ok(result)
    }

    /// Cut end segment and find anchors in one pass.
    /// Returns the segment with precomputed anchors for reuse.
    /// Returns an empty segment if the FASTQ sequence is too short.
    fn cut_and_analyze_segment(
        &self,
        sequence: &[u8],
        quality: &[u8],
        end_regions: &EndRegions,
        is_five_prime: bool,
        should_rc: bool,
    ) -> Result<EndSegmentWithAnchors> {
        let cut_length = end_regions.calculate_cut_length();

        if sequence.len() < cut_length {
            return Ok(EndSegmentWithAnchors::empty());
        }

        // Determine which end to cut from based on orientation
        let cut_from_beginning = should_rc != is_five_prime;

        let (seq_segment, qual_segment) = if cut_from_beginning {
            (&sequence[..cut_length], &quality[..cut_length])
        } else {
            let start = sequence.len() - cut_length;
            (&sequence[start..], &quality[start..])
        };

        // Apply RC to the segment if needed
        let (final_seq, final_qual) = if should_rc {
            let rc_seq = seqspec::utils::rev_compl(seq_segment);
            let rc_qual: Vec<u8> = qual_segment.iter().copied().rev().collect();
            (rc_seq, rc_qual)
        } else {
            (seq_segment.to_vec(), qual_segment.to_vec())
        };

        // Find anchors
        let fixed_regions = end_regions.get_fixed_regions();
        let expected_fixed_count = fixed_regions.len();
        let mut anchor_finder = AnchorFinder::new()?;
        let anchors = anchor_finder.find_anchors(&final_seq, &fixed_regions, end_regions)?;

        let valid_order = if !anchors.is_empty() {
            anchor_finder.validate_order(&anchors)?
        } else {
            false
        };

        Ok(EndSegmentWithAnchors {
            sequence: final_seq,
            quality: final_qual,
            anchors,
            valid_order,
            expected_fixed_count,
        })
    }

    /// Combine evidence from both ends
    fn combine_evidence(
        five_prime: &EndSegmentWithAnchors,
        three_prime: &EndSegmentWithAnchors,
    ) -> OrientationEvidence {
        let total_anchors = five_prime.anchors.len() + three_prime.anchors.len();

        let avg_confidence = if total_anchors > 0 {
            let five_conf: f64 = five_prime.anchors.iter().map(|a| a.score).sum();
            let three_conf: f64 = three_prime.anchors.iter().map(|a| a.score).sum();
            (five_conf + three_conf) / total_anchors as f64
        } else {
            0.0
        };

        OrientationEvidence {
            num_anchors: total_anchors,
            valid_order: five_prime.valid_order && three_prime.valid_order,
            avg_confidence,
        }
    }

    /// Extract barcodes from precomputed segments with cached anchors
    fn extract_barcodes_from_segments(
        &self,
        five_prime: EndSegmentWithAnchors,
        three_prime: EndSegmentWithAnchors,
        is_reverse_complemented: bool,
    ) -> Result<LongReadBarcodeResult> {
        let mut extracted_barcodes = Vec::new();

        // Process 5' end using cached anchors
        if !five_prime.is_empty() {
            let barcode_results = self.extract_barcodes_with_cached_anchors(
                &five_prime,
                &self.five_prime_regions,
            )?;
            extracted_barcodes.extend(barcode_results);
        }

        // Process 3' end using cached anchors
        if !three_prime.is_empty() {
            let barcode_results = self.extract_barcodes_with_cached_anchors(
                &three_prime,
                &self.three_prime_regions,
            )?;
            extracted_barcodes.extend(barcode_results);
        }

        let mut result = self.combine_barcodes(extracted_barcodes)?;
        result.is_reverse_complemented = is_reverse_complemented;
        Ok(result)
    }

    /// Extract barcodes using precomputed anchors
    fn extract_barcodes_with_cached_anchors(
        &self,
        segment: &EndSegmentWithAnchors,
        end_regions: &EndRegions,
    ) -> Result<Vec<ExtractedBarcode>> {
        if !end_regions.has_barcode {
            return Ok(Vec::new());
        }

        // Skip if anchor order is invalid
        if !segment.anchors.is_empty() && !segment.valid_order {
            warn!("Fixed region alignment order is incorrect, skipping this end");
            return Ok(Vec::new());
        }

        let barcode_regions = end_regions.get_barcode_regions();
        if barcode_regions.is_empty() {
            return Ok(Vec::new());
        }

        let mut extracted_barcodes = Vec::new();

        for barcode_region in barcode_regions {
            if let Some(extraction_range) = locate_barcode(
                &barcode_region,
                &segment.anchors,
                end_regions,
                segment.sequence.len(),
            ) {
                let candidate_seq = &segment.sequence[extraction_range.0..extraction_range.1];

                if let Some(matched) = find_best_barcode_match(
                    candidate_seq,
                    &barcode_region,
                    &self.whitelists,
                ) {
                    extracted_barcodes.push(matched);
                }
            }
        }

        Ok(extracted_barcodes)
    }

    /// Get whitelists reference for external access
    pub fn whitelists(&self) -> &IndexMap<String, IndexSet<Vec<u8>>> {
        &self.whitelists
    }

    /// Get 5' end regions for external access
    pub fn five_prime_regions(&self) -> &EndRegions {
        &self.five_prime_regions
    }

    /// Get 3' end regions for external access
    pub fn three_prime_regions(&self) -> &EndRegions {
        &self.three_prime_regions
    }

    /// Compare evidence and choose which orientation is better
    /// Returns true if reverse-complement should be used
    fn compare_evidence(&self, forward: &OrientationEvidence, reverse: &OrientationEvidence) -> bool {
        // Choose reverse if: more anchors, OR same anchors but better quality
        if reverse.num_anchors > forward.num_anchors {
            true
        } else if reverse.num_anchors == forward.num_anchors {
            // Same number of anchors, use tie-breaker logic
            if reverse.valid_order && !forward.valid_order {
                true
            } else if forward.valid_order && !reverse.valid_order {
                false
            } else {
                // Both have same order status, choose higher confidence
                reverse.avg_confidence > forward.avg_confidence
            }
        } else {
            false
        }
    }

    /// Combine multiple extracted barcodes into final result.
    fn combine_barcodes(&self, barcodes: Vec<ExtractedBarcode>) -> Result<LongReadBarcodeResult> {
        if barcodes.is_empty() {
            return Ok(LongReadBarcodeResult {
                barcode: None,
                confidence: 0.0,
                is_reverse_complemented: false,
                has_missing_anchors: false,
            });
        }

        if barcodes.len() == 1 {
            let barcode = &barcodes[0];
            return Ok(LongReadBarcodeResult {
                barcode: Some(barcode.barcode.clone()),
                confidence: barcode.confidence,
                is_reverse_complemented: false,
                has_missing_anchors: false,
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
            is_reverse_complemented: false,
            has_missing_anchors: false,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::sequence_aligner::fitting_alignment_distance;

    #[test]
    fn test_orientation_evidence_thresholds() {
        // Case 1: Meets all thresholds
        let evidence = OrientationEvidence {
            num_anchors: 2,
            valid_order: true,
            avg_confidence: 0.90,
        };
        assert!(evidence.meets_threshold());

        // Case 2: Below confidence threshold (0.8)
        let evidence = OrientationEvidence {
            num_anchors: 2,
            valid_order: true,
            avg_confidence: 0.79,
        };
        assert!(!evidence.meets_threshold());

        // Case 3: Below anchor count
        let evidence = OrientationEvidence {
            num_anchors: 1,
            valid_order: true,
            avg_confidence: 0.90,
        };
        assert!(!evidence.meets_threshold());

        // Case 4: Invalid order
        let evidence = OrientationEvidence {
            num_anchors: 2,
            valid_order: false,
            avg_confidence: 0.90,
        };
        assert!(!evidence.meets_threshold());

        // Case 5: Edge case - exactly at threshold
        let evidence = OrientationEvidence {
            num_anchors: 2,
            valid_order: true,
            avg_confidence: 0.80,
        };
        assert!(evidence.meets_threshold());
    }

    #[test]
    fn test_orientation_evidence_from_anchors_empty() {
        let anchors: Vec<FixedSequenceAlignment> = vec![];
        let evidence = OrientationEvidence::from_anchors(&anchors, false);

        assert_eq!(evidence.num_anchors, 0);
        assert!(!evidence.valid_order);
        assert_eq!(evidence.avg_confidence, 0.0);
    }

    #[test]
    fn test_orientation_detection_with_alignment_scores() {
        use seqspec::{Modality, RegionType, SequenceType};
        use seqspec::region::{LibSpec, Region};
        use indexmap::{IndexMap, IndexSet};

        // Create regions with one fixed sequence
        // lib structure: fixed(16bp "AAAAAAAAAAAAAAAA") -> barcode(10bp) -> cdna
        let fixed_seq = "AAAAAAAAAAAAAAAA"; // 16bp of A's

        let fixed_5p = Arc::new(RwLock::new(Region {
            region_id: "fixed_5p".to_string(),
            region_type: RegionType::Linker,
            name: "fixed_5p".to_string(),
            sequence_type: SequenceType::Fixed,
            sequence: fixed_seq.to_string(),
            min_len: 16,
            max_len: 16,
            onlist: None,
            subregions: vec![],
        }));

        let barcode = Arc::new(RwLock::new(Region {
            region_id: "bc1".to_string(),
            region_type: RegionType::Barcode,
            name: "bc1".to_string(),
            sequence_type: SequenceType::Onlist,
            sequence: "".to_string(),
            min_len: 10,
            max_len: 10,
            onlist: None,
            subregions: vec![],
        }));

        let cdna = Arc::new(RwLock::new(Region {
            region_id: "cdna".to_string(),
            region_type: RegionType::Cdna,
            name: "cdna".to_string(),
            sequence_type: SequenceType::Random,
            sequence: "".to_string(),
            min_len: 100,
            max_len: 1000,
            onlist: None,
            subregions: vec![],
        }));

        let modality_region = Region {
            region_id: "rna".to_string(),
            region_type: RegionType::Modality(Modality::RNA),
            name: "RNA".to_string(),
            sequence_type: SequenceType::Joined,
            sequence: "".to_string(),
            min_len: 0,
            max_len: 0,
            onlist: None,
            subregions: vec![fixed_5p, barcode, cdna],
        };

        let lib_spec = LibSpec::new(vec![modality_region]).unwrap();
        let whitelists: IndexMap<String, IndexSet<Vec<u8>>> = IndexMap::new();
        let extractor = BarcodeExtractor::new(&lib_spec, &Modality::RNA, whitelists).unwrap();

        let cut_length = extractor.five_prime_regions().calculate_cut_length();

        // Create a test sequence that indicates forward orientation:
        let sequence_with_error = format!(
            "{}{}{}",
            "AAAAACAAAAACAAAAA",                                // 17bp: 15 A's + 2 insertion (C)
            "C".repeat(cut_length + 10),                        // C's in the middle
            "TTTTCTTTTCTTTTCCC"                                 // 17bp: 12 T's + 5 insertion (C)
        );
        let test_sequence = sequence_with_error.as_bytes();
        let test_quality = "I".repeat(test_sequence.len());

        // Test forward orientation
        let forward_5p = extractor.cut_and_analyze_segment(
            test_sequence,
            test_quality.as_bytes(),
            extractor.five_prime_regions(),
            true,  // is_five_prime
            false, // should_rc = false (forward)
        ).unwrap();

        // Test reverse orientation - cut from end and RC
        let reverse_5p = extractor.cut_and_analyze_segment(
            test_sequence,
            test_quality.as_bytes(),
            extractor.five_prime_regions(),
            true, // is_five_prime
            true, // should_rc = true (reverse)
        ).unwrap();

        // Forward should have better anchor match
        // Reverse would have RC of T's at end = A's
        println!("Forward anchors: {:?}", forward_5p.anchors.len());
        println!("Reverse anchors: {:?}", reverse_5p.anchors.len());

        // Combine evidence and compare
        let forward_evidence = BarcodeExtractor::combine_evidence(&forward_5p, &EndSegmentWithAnchors::empty());
        let reverse_evidence = BarcodeExtractor::combine_evidence(&reverse_5p, &EndSegmentWithAnchors::empty());

        println!("Forward evidence: anchors={}, valid_order={}, confidence={:.3}",
            forward_evidence.num_anchors, forward_evidence.valid_order, forward_evidence.avg_confidence);
        println!("Reverse evidence: anchors={}, valid_order={}, confidence={:.3}",
            reverse_evidence.num_anchors, reverse_evidence.valid_order, reverse_evidence.avg_confidence);

        let should_rc = extractor.compare_evidence(&forward_evidence, &reverse_evidence);

        if forward_evidence.num_anchors > 0 || reverse_evidence.num_anchors > 0 {
            println!("should_rc = {}", should_rc);
        }
    }

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

        // Create alignments of fixed sequences with position set
        let anchor1 = FixedSequenceAlignment {
            region: outer_region,
            query_start: 2,
            query_end: 6,
            score: 1.0,
            matches: 4,
            alignment_length: 4,
            position: Some(1), // outer position
        };

        let anchor2 = FixedSequenceAlignment {
            region: inner_region,
            query_start: 20,
            query_end: 24,
            score: 1.0,
            matches: 4,
            alignment_length: 4,
            position: Some(3), // inner position
        };

        let anchors = vec![anchor1, anchor2];

        // Test locate_barcode method
        let result = locate_barcode(&barcode_region, &anchors, &end_regions, 100);

        // Should find a range between the outer anchor end (6) and inner anchor start (20)
        assert!(result.is_some());
        let (start, end) = result.unwrap();
        assert!(start >= 6); // After outer anchor
        assert!(end <= 20);  // Before inner anchor
        assert!(end > start); // Valid range

        println!("start: {}, end: {}", start, end);
    }
}