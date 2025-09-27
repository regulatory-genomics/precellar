use anyhow::Result;
use log::warn;
use std::collections::HashMap;
use std::sync::{Arc, RwLock};

use crate::barcode::{BarcodeCorrector, OligoFrequncy};
use seqspec::region::{LibSpec, Region};
use seqspec::Modality;

use super::{
    find_innermost_regions, EndRegions, LongReadBarcodeResult,
    sequence_aligner::{LongReadSequenceAligner, FixedSequenceAlignment, find_non_overlapping_alignments}
};

/// Main barcode extractor for long reads
pub struct BarcodeExtractor {
    /// 5' end regions
    five_prime_regions: EndRegions,
    /// 3' end regions  
    three_prime_regions: EndRegions,
}

impl BarcodeExtractor {
    /// Create a new barcode extractor for the given library specification and modality
    pub fn new(lib_spec: &LibSpec, modality: &Modality) -> Result<Self> {
        let (five_prime_regions, three_prime_regions) = find_innermost_regions(lib_spec, modality)?;
        
        Ok(Self {
            five_prime_regions,
            three_prime_regions,
        })
    }

    /// Extract barcode from a FASTQ record
    pub fn extract_barcode(
        &self,
        record: &noodles::fastq::Record,
        whitelists: &HashMap<String, OligoFrequncy>,
        corrector: Option<&BarcodeCorrector>,
    ) -> Result<LongReadBarcodeResult> {
        let sequence = record.sequence();
        let quality = record.quality_scores();

        // Step 1: Cut sequences from both ends
        let five_prime_segment = self.cut_end_segment(sequence, quality, &self.five_prime_regions, true)?;
        let three_prime_segment = self.cut_end_segment(sequence, quality, &self.three_prime_regions, false)?;

        let mut extracted_barcodes = Vec::new();

        // Step 2 & 3: Process 5' end
        if let Some((seq, qual)) = five_prime_segment {
            if let Some(barcode_result) = self
            .process_end_for_barcode(&seq, &qual, &self.five_prime_regions, whitelists, corrector)? 
            {
                extracted_barcodes.push(barcode_result);
            }
        }

        // Step 2 & 3: Process 3' end  
        if let Some((seq, qual)) = three_prime_segment {
            if let Some(barcode_result) = self
            .process_end_for_barcode(&seq, &qual, &self.three_prime_regions, whitelists, corrector)? 
            {
                extracted_barcodes.push(barcode_result);
            }
        }

        // Step 4: Combine multiple barcodes if present
        let final_result = self.combine_barcodes(extracted_barcodes)?;

        Ok(final_result)
    }

    /// Cut end sequence from one fastq record
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

    /// Process one end to extract barcode
    fn process_end_for_barcode(
        &self,
        segment_seq: &[u8],
        segment_qual: &[u8],
        end_regions: &EndRegions,
        whitelists: &HashMap<String, OligoFrequncy>,
        corrector: Option<&BarcodeCorrector>,
    ) -> Result<Option<ExtractedBarcode>> {

        // If one end has no barcode region, return None
        if !end_regions.has_barcode {
            return Ok(None);
        }

        let fixed_regions = end_regions.get_fixed_regions();
        let barcode_regions = end_regions.get_barcode_regions();

        if barcode_regions.is_empty() {
            return Ok(None);
        }

        // Step1: Align fixed sequences to find anchor points
        let anchor_positions = if !fixed_regions.is_empty() {
            self.find_anchor_positions(segment_seq, &fixed_regions)?
        } else {
            Vec::new()
        };

        // Extract barcodes based on anchor positions or sliding window
        let barcode_candidates = self.extract_barcode_candidates(
            segment_seq,
            segment_qual,
            &barcode_regions,
            &anchor_positions,
        )?;

        // Find best barcode match using whitelist
        let best_barcode = self.find_best_barcode_match(
            barcode_candidates,
            whitelists,
            corrector,
        )?;

        Ok(best_barcode)
    }

    /// Find anchor positions using fixed sequence alignment
    fn find_anchor_positions(
        &self,
        sequence: &[u8],
        fixed_regions: &[Arc<RwLock<Region>>],
    ) -> Result<Vec<FixedSequenceAlignment>> {
        let mut aligner = LongReadSequenceAligner::new()?;
        let alignments = aligner.align_fixed_sequences(sequence, fixed_regions)?;
        
        // Filter alignments with good scores (e.g., > 0.7)
        let good_alignments: Vec<_> = alignments
            .into_iter()
            .filter(|aln| aln.score > 0.7)
            .collect();

        // Find non-overlapping alignments
        Ok(find_non_overlapping_alignments(good_alignments))
    }

    /// Extract barcode candidates from the sequence
    fn extract_barcode_candidates(
        &self,
        sequence: &[u8],
        quality: &[u8],
        barcode_regions: &[Arc<RwLock<Region>>],
        anchor_positions: &[FixedSequenceAlignment],
    ) -> Result<Vec<BarcodeCandidate>> {
        let mut candidates = Vec::new();

        for barcode_region in barcode_regions {
            let region_guard = barcode_region.read().unwrap();
            let barcode_length = region_guard.max_len as usize;
            let region_id = region_guard.region_id.clone();
            drop(region_guard);

            // Determine potential barcode windows based on anchor positions
            let windows = if anchor_positions.is_empty() {
                // No anchors - use sliding window across entire sequence
                self.generate_sliding_windows(sequence.len(), barcode_length)
            } else {
                // Use anchor positions to determine barcode windows
                self.generate_anchor_based_windows(anchor_positions, barcode_length)
            };

            // Extract candidates from each window
            for (start, end) in windows {
                if end <= sequence.len() && start < end {
                    let candidate_seq = &sequence[start..end];
                    let candidate_qual = &quality[start..end];
                    
                    candidates.push(BarcodeCandidate {
                        region_id: region_id.clone(),
                        sequence: candidate_seq.to_vec(),
                        quality: candidate_qual.to_vec(),
                        start_pos: start,
                        end_pos: end,
                    });
                }
            }
        }

        Ok(candidates)
    }

    /// Generate sliding windows across the sequence
    fn generate_sliding_windows(&self, seq_length: usize, barcode_length: usize) -> Vec<(usize, usize)> {
        let mut windows = Vec::new();
        
        if barcode_length > seq_length {
            return windows;
        }

        // Slide window across sequence with step size of 1
        for start in 0..=(seq_length - barcode_length) {
            windows.push((start, start + barcode_length));
        }

        windows
    }

    /// Generate barcode windows based on anchor positions
    fn generate_anchor_based_windows(
        &self,
        anchors: &[FixedSequenceAlignment],
        barcode_length: usize,
    ) -> Vec<(usize, usize)> {
        let mut windows = Vec::new();

        if anchors.is_empty() {
            return windows;
        }

        // For each pair of adjacent anchors, the barcode could be between them
        for i in 0..anchors.len() {
            let current_anchor = &anchors[i];
            
            // Barcode could be before this anchor
            if current_anchor.query_start >= barcode_length {
                let start = current_anchor.query_start - barcode_length;
                windows.push((start, current_anchor.query_start));
            }
            
            // Barcode could be after this anchor
            let start = current_anchor.query_end;
            windows.push((start, start + barcode_length));
            
            // If there's a next anchor, barcode could be between them
            if i + 1 < anchors.len() {
                let next_anchor = &anchors[i + 1];
                let space_between = next_anchor.query_start.saturating_sub(current_anchor.query_end);
                
                if space_between >= barcode_length {
                    // Multiple possible positions between anchors
                    for pos in current_anchor.query_end..=(next_anchor.query_start - barcode_length) {
                        windows.push((pos, pos + barcode_length));
                    }
                }
            }
        }

        windows
    }

    /// Find the best barcode match using whitelists and error correction
    fn find_best_barcode_match(
        &self,
        candidates: Vec<BarcodeCandidate>,
        whitelists: &HashMap<String, OligoFrequncy>,
        corrector: Option<&BarcodeCorrector>,
    ) -> Result<Option<ExtractedBarcode>> {
        let mut best_candidate = None;
        let mut best_score = 0.0;

        for candidate in candidates {
            if let Some(whitelist) = whitelists.get(&candidate.region_id) {
                let score = if let Some(corrector) = corrector {
                    // Use corrector to find best match
                    match corrector.correct(whitelist, &candidate.sequence, &candidate.quality) {
                        Ok(_) => {
                            // High confidence correction
                            1.0
                        }
                        Err(_) => {
                            // Try direct match
                            if whitelist.contains_key(&candidate.sequence) {
                                0.8
                            } else {
                                0.0
                            }
                        }
                    }
                } else {
                    // Direct whitelist lookup
                    if whitelist.contains_key(&candidate.sequence) {
                        0.9
                    } else {
                        // Calculate edit distance to closest whitelist entry
                        self.calculate_whitelist_similarity(whitelist, &candidate.sequence)
                    }
                };

                if score > best_score {
                    best_score = score;
                    best_candidate = Some(ExtractedBarcode {
                        region_id: candidate.region_id.clone(),
                        sequence: candidate.sequence.clone(),
                        quality: candidate.quality.clone(),
                        confidence: score,
                        start_pos: candidate.start_pos,
                        end_pos: candidate.end_pos,
                    });
                }
            }
        }

        Ok(best_candidate)
    }

    /// Calculate similarity to whitelist entries using edit distance
    fn calculate_whitelist_similarity(&self, whitelist: &OligoFrequncy, query: &[u8]) -> f64 {
        if whitelist.is_empty() {
            return 0.0;
        }

        let mut min_distance = usize::MAX;
        
        for barcode in whitelist.keys() {
            let distance = edit_distance(query, barcode);
            min_distance = min_distance.min(distance);
        }

        // Convert edit distance to similarity score
        let max_distance = query.len().max(1);
        1.0 - (min_distance as f64 / max_distance as f64)
    }

    /// Combine multiple extracted barcodes into final result
    fn combine_barcodes(&self, barcodes: Vec<ExtractedBarcode>) -> Result<LongReadBarcodeResult> {
        if barcodes.is_empty() {
            return Ok(LongReadBarcodeResult {
                barcode: None,
                quality: None,
                confidence: 0.0,
                success: false,
            });
        }

        if barcodes.len() == 1 {
            let barcode = &barcodes[0];
            return Ok(LongReadBarcodeResult {
                barcode: Some(barcode.sequence.clone()),
                quality: Some(barcode.quality.clone()),
                confidence: barcode.confidence,
                success: true,
            });
        }

        // Multiple barcodes - concatenate them
        let mut combined_sequence = Vec::new();
        let mut combined_quality = Vec::new();
        let mut total_confidence = 0.0;

        for barcode in &barcodes {
            combined_sequence.extend(&barcode.sequence);
            combined_quality.extend(&barcode.quality);
            total_confidence += barcode.confidence;
        }

        let average_confidence = total_confidence / barcodes.len() as f64;

        Ok(LongReadBarcodeResult {
            barcode: Some(combined_sequence),
            quality: Some(combined_quality),
            confidence: average_confidence,
            success: true,
        })
    }
}

/// Barcode candidate extracted from sequence
#[derive(Debug, Clone)]
struct BarcodeCandidate {
    region_id: String,
    sequence: Vec<u8>,
    quality: Vec<u8>,
    start_pos: usize,
    end_pos: usize,
}

/// Successfully extracted barcode with metadata
#[derive(Debug, Clone)]
struct ExtractedBarcode {
    #[allow(dead_code)]
    region_id: String,
    sequence: Vec<u8>,
    quality: Vec<u8>,
    confidence: f64,
    #[allow(dead_code)]
    start_pos: usize,
    #[allow(dead_code)]
    end_pos: usize,
}

/// Calculate edit distance between two sequences
fn edit_distance(a: &[u8], b: &[u8]) -> usize {
    let m = a.len();
    let n = b.len();
    
    if m == 0 {
        return n;
    }
    if n == 0 {
        return m;
    }

    let mut dp = vec![vec![0; n + 1]; m + 1];
    
    // Initialize base cases
    for i in 0..=m {
        dp[i][0] = i;
    }
    for j in 0..=n {
        dp[0][j] = j;
    }
    
    // Fill DP table
    for i in 1..=m {
        for j in 1..=n {
            let cost = if a[i-1] == b[j-1] { 0 } else { 1 };
            dp[i][j] = (dp[i-1][j] + 1)
                .min(dp[i][j-1] + 1)
                .min(dp[i-1][j-1] + cost);
        }
    }
    
    dp[m][n]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_edit_distance() {
        assert_eq!(edit_distance(b"ATCG", b"ATCG"), 0);
        assert_eq!(edit_distance(b"ATCG", b"ATCA"), 1);
        assert_eq!(edit_distance(b"ATCG", b"GCTA"), 4);
        assert_eq!(edit_distance(b"", b"ATCG"), 4);
        assert_eq!(edit_distance(b"ATCG", b""), 4);
    }

    #[test]
    fn test_generate_sliding_windows() {
        use super::super::EndType;
        let extractor = BarcodeExtractor {
            five_prime_regions: EndRegions::new(EndType::FivePrime),
            three_prime_regions: EndRegions::new(EndType::ThreePrime),
        };
        
        let windows = extractor.generate_sliding_windows(10, 4);
        assert_eq!(windows.len(), 7); // 10 - 4 + 1 = 7
        assert_eq!(windows[0], (0, 4));
        assert_eq!(windows[6], (6, 10));
    }
}
