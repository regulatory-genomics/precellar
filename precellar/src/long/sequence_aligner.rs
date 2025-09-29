use anyhow::Result;
use log::warn;
use std::sync::{Arc, RwLock};
use bio::alignment::pairwise::{Aligner, Scoring, MIN_SCORE};
use bio::alignment::AlignmentOperation;

use seqspec::region::Region;

/// Alignment result for a fixed sequence region
#[derive(Debug, Clone)]
pub struct FixedSequenceAlignment {
    /// Region that was aligned
    pub region: Arc<RwLock<Region>>,
    /// Start position in the trimmed end sequence
    pub query_start: usize,
    /// End position in the trimmed end sequence
    pub query_end: usize,
    /// Match rate (matches / alignment_length)
    pub score: f64,
    /// Number of matches
    pub matches: usize,
    /// Total alignment length
    pub alignment_length: usize,
}


/// Sequence aligner for finding fixed regions in trimmed EndRegions from long reads
/// Uses rust-bio's pairwise aligner with fitting alignment mode
pub struct FittingAligner {
    /// Minimum match score threshold
    min_score: f64,
    /// Scoring configuration for fitting alignment
    scoring: Scoring<fn(u8, u8) -> i32>,
}

// Define a match function for scoring initialization
fn match_func(a: u8, b: u8) -> i32 {
    if a == b { 2i32 } else { -1i32 }
}

impl FittingAligner {
    /// Create a new sequence aligner with fitting alignment configuration
    pub fn new() -> Result<Self> {
        // Configure scoring for fitting alignment:
        // match = +2, mismatch = -1, gap_penalty = -1
        // For fitting alignment: 
        // - xclip_prefix = MIN_SCORE (no gaps allowed at start of pattern)
        // - yclip_prefix = 0 (gaps allowed at start of text, no penalty)
        // - xclip_suffix = MIN_SCORE (no gaps allowed at end of pattern)  
        // - yclip_suffix = 0 (gaps allowed at end of text, no penalty)
        let scoring = Scoring {
            gap_open: -1,
            gap_extend: -1,
            match_fn: match_func as fn(u8, u8) -> i32,
            match_scores: Some((2, -1)), // (match, mismatch)
            xclip_prefix: MIN_SCORE,  // No gaps allowed at start of pattern
            xclip_suffix: MIN_SCORE,  // No gaps allowed at end of pattern
            yclip_prefix: 0,          // Allow gaps at start of text, no penalty
            yclip_suffix: 0,          // Allow gaps at end of text, no penalty
        };
        
        Ok(Self {
            min_score: 0.7,
            scoring,
        })
    }

    /// Align fixed sequences to find their positions in the end sequence using fitting alignment
    pub fn align_fixed_sequences(
        &mut self,
        end_sequence: &[u8],
        fixed_regions: &[Arc<RwLock<Region>>],
    ) -> Result<Vec<FixedSequenceAlignment>> {
        let mut results = Vec::new();

        for fixed_region in fixed_regions {
            let fixed_sequence = {
                let region_guard = fixed_region.read().unwrap();
                let seq = region_guard.sequence.as_bytes().to_vec();
                
                // // Skip if fixed sequence is too short
                // if seq.len() < 15 {
                //     continue;
                // }
                seq
            };

            let mut aligner = Aligner::with_capacity_and_scoring(
                fixed_sequence.len(), 
                end_sequence.len(), 
                self.scoring
            );
            
            // Customize with fitting alignment parameters
            let alignment = aligner.custom(&fixed_sequence, end_sequence);
            
            // Calculate number of matches and alignment length
            let (matches, alignment_length) = self.calculate_alignment_stats(&alignment);

            // println!("Alignment score from rust-bio: {}", alignment.score);

            let score = if alignment_length > 0 {
                matches as f64 / alignment_length as f64
            } else {
                0.0
            };

            if score >= self.min_score {

                // ystart and yend are the start and end positions in the end sequence
                let (query_start, query_end) = (alignment.ystart, alignment.yend);
                
                results.push(FixedSequenceAlignment {
                    region: fixed_region.clone(),
                    query_start,
                    query_end,
                    score,
                    matches,
                    alignment_length,
                });
            } else {
                // Warn about low-quality alignment
                let region_guard = fixed_region.read().unwrap();
                warn!(
                    "Fixed region '{}' (sequence: '{}') alignment score ({:.3}) below threshold ({:.3}), skipping",
                    region_guard.region_id, region_guard.sequence, score, self.min_score
                );
            }
        }

        Ok(results)
    }

    /// Calculate alignment statistics from bio::alignment operations
    fn calculate_alignment_stats(&self, alignment: &bio::alignment::Alignment) -> (usize, usize) {
        let mut matches = 0;
        let mut total_length = 0;
        
        for operation in &alignment.operations {
            match operation {
                AlignmentOperation::Match => {
                    matches += 1;
                    total_length += 1;
                }
                AlignmentOperation::Subst => {
                    total_length += 1;
                }
                AlignmentOperation::Del | AlignmentOperation::Ins => {
                    total_length += 1;
                }
                _ => {}
            }
        }
        
        (matches, total_length)
    }

}

/// Selects the best non-overlapping alignments using a greedy algorithm.
///
/// Sorts all candidates by score in descending order, then iteratively picks the
/// next best-scoring alignment that does not conflict with those already selected.
/// The final result is sorted by the alignment start position.
pub fn find_non_overlapping_alignments(
    alignments: Vec<FixedSequenceAlignment>
) -> Vec<FixedSequenceAlignment> {
    if alignments.is_empty() {
        return Vec::new();
    }

    let mut sorted_alignments = alignments;
    // Sort by score (descending)
    sorted_alignments.sort_by(|a, b| b.score.partial_cmp(&a.score).unwrap_or(std::cmp::Ordering::Equal));

    let mut selected = Vec::new();
    
    for alignment in sorted_alignments {
        // Check if this alignment overlaps with any already selected
        let overlaps = selected.iter().any(|selected_aln: &FixedSequenceAlignment| {
            !(alignment.query_end <= selected_aln.query_start || 
              alignment.query_start >= selected_aln.query_end)
        });
        
        // Greedy algorithm: save high score fixed alignments and discard lower overlapping ones
        if !overlaps {
            selected.push(alignment);
        }
    }
    
    // Sort selected alignments by position
    selected.sort_by_key(|a| a.query_start);
    
    selected
}

/// Calculate fitting alignment distance between short sequence and long sequence
pub fn fitting_alignment_distance(short_seq: &[u8], long_seq: &[u8]) -> usize {
    let m = short_seq.len();
    let n = long_seq.len();
    
    if m == 0 {
        return 0; // Empty short sequence can always match
    }
    if n == 0 {
        return m; // Cannot fit non-empty short sequence in empty long sequence
    }

    let mut dp = vec![vec![0; n + 1]; m + 1];
    
    // First row: no penalty for gaps at the start of long sequence (fitting alignment)
    for j in 0..=n {
        dp[0][j] = 0;
    }
    // First column: penalty for gaps in short sequence (cannot skip characters in short sequence)
    for i in 1..=m {
        dp[i][0] = i;
    }
    
    // Fill DP table
    for i in 1..=m {
        for j in 1..=n {
            let cost = if short_seq[i-1] == long_seq[j-1] { 0 } else { 1 };
            dp[i][j] = (dp[i-1][j] + 1)
                .min(dp[i][j-1] + 1)
                .min(dp[i-1][j-1] + cost);
        }
    }
    
    // Return minimum value in the last row (fitting alignment)
    dp[m].iter().min().copied().unwrap_or(m)
}

#[cfg(test)]
mod tests {
    use super::*;
    use seqspec::{RegionType, SequenceType};

    // Create test fixed sequence region
    fn create_test_region(
        id: &str,
        sequence: &str,
        min_len: u32,
        max_len: u32,
    ) -> Arc<RwLock<Region>> {
        Arc::new(RwLock::new(Region {
            region_id: id.to_string(),
            region_type: RegionType::Linker,
            name: id.to_string(),
            sequence_type: SequenceType::Fixed,
            sequence: sequence.to_string(),
            min_len,
            max_len,
            onlist: None,
            subregions: vec![],
        }))
    }

    #[test]
    fn test_find_non_overlapping_alignments() {
        let region1 = create_test_region("r1", "ATCG", 4, 4);
        let region2 = create_test_region("r2", "GCTA", 4, 4);
        let region3 = create_test_region("r3", "TTTT", 4, 4);

        let alignments = vec![
            FixedSequenceAlignment {
                region: region1,
                query_start: 10,
                query_end: 14,
                score: 0.9,
                matches: 4,
                alignment_length: 4,
            },
            FixedSequenceAlignment {
                region: region2,
                query_start: 12, // Overlaps with first
                query_end: 16,
                score: 0.8,
                matches: 4,
                alignment_length: 4,
            },
            FixedSequenceAlignment {
                region: region3,
                query_start: 20, // No overlap
                query_end: 24,
                score: 0.7,
                matches: 4,
                alignment_length: 4,
            },
        ];

        let non_overlapping = find_non_overlapping_alignments(alignments);
        
        // Should select the first (highest score) and third (no overlap)
        assert_eq!(non_overlapping.len(), 2);
        assert_eq!(non_overlapping[0].query_start, 10);
        assert_eq!(non_overlapping[1].query_start, 20);
    }

    #[test]
    fn test_fitting_alignment_basic() {
        // Test the specific example: short sequence v="AT", long sequence w="GCATG"
        // Scoring: match = +2, mismatch = -1, gap_penalty = -1
        // This test verifies that fitting alignment works correctly
        let mut aligner = FittingAligner::new().unwrap();
        let long_seq = b"GCATG";  // text to search in
        
        // Create a test fixed sequence region
        let region = create_test_region("test", "AT", 2, 2);
        let fixed_regions = vec![region];
        
        let alignments = aligner.align_fixed_sequences(long_seq, &fixed_regions).unwrap();
        
        
        // Should find the "AT" at position 2-4 in "GCATG"
        assert!(!alignments.is_empty(), "Should find at least one alignment");
        let best_alignment = &alignments[0];
        
        // The "AT" should be found at positions 2-4 in "GCATG"
        assert_eq!(best_alignment.query_start, 2);
        assert_eq!(best_alignment.query_end, 4);
        assert_eq!(best_alignment.matches, 2);  // Both A and T match
        assert_eq!(best_alignment.alignment_length, 2);
        assert_eq!(best_alignment.score, 1.0);  // Perfect match
    }
    
    #[test]
    fn test_fitting_alignment_properties() {
        // This test verifies the key properties of fitting alignment:
        // 1. First row initialization to 0 (no penalty for gaps at start of long sequence)
        // 2. Finding maximum in last row (no penalty for gaps at end of long sequence)
        
        let mut aligner = FittingAligner::new().unwrap();
        
        // Test case where the pattern appears at the very beginning
        let text_start = b"GCATTTTT";
        let region_start = create_test_region("start", "GC", 2, 2);
        let alignments_start = aligner.align_fixed_sequences(text_start, &vec![region_start]).unwrap();
        
        assert!(!alignments_start.is_empty());
        assert_eq!(alignments_start[0].query_start, 0);  // Found at start
        assert_eq!(alignments_start[0].query_end, 2);
        
        // Test case where the pattern appears at the very end
        let text_end = b"GGGGGGTT";
        let region_end = create_test_region("end", "TT", 2, 2);
        let alignments_end = aligner.align_fixed_sequences(text_end, &vec![region_end]).unwrap();
        
        assert!(!alignments_end.is_empty());
        assert_eq!(alignments_end[0].query_start, 6);  // Found at end
        assert_eq!(alignments_end[0].query_end, 8);
        
        // Both should have perfect scores since they're exact matches
        assert_eq!(alignments_start[0].score, 1.0);
        assert_eq!(alignments_end[0].score, 1.0);
    }
}