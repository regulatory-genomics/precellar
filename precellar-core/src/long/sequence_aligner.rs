use anyhow::Result;
use bio::alignment::pairwise::{Aligner, Scoring, MIN_SCORE};
use bio::alignment::AlignmentOperation;
use std::ops::Range;
use std::sync::{Arc, RwLock};

use seqspec::region::Region;

/// A composite alignment pattern built from regions between the leftmost and rightmost fixed
/// regions in forward order of the designed library structure.
#[derive(Debug, Clone)]
pub struct CompositePattern {
    /// The concatenated pattern bytes (fixed sequences + N-spacers)
    pub pattern: Vec<u8>,
    /// Region spans ordered by position in the pattern
    pub spans: Vec<CompositeRegionSpan>,
    /// Lookup: pattern position -> span index (precomputed, same length as pattern)
    pub pos_to_span: Vec<usize>,
}

impl CompositePattern {
    /// Total length of fixed (non-spacer) sequences in the pattern.
    pub fn total_fixed_len(&self) -> usize {
        self.spans
            .iter()
            .filter(|s| !s.is_spacer)
            .map(|s| s.pattern_end - s.pattern_start)
            .sum()
    }
}

/// Describes one region's span within the CompositePattern.
#[derive(Debug, Clone)]
pub struct CompositeRegionSpan {
    /// The region this span corresponds to
    pub region: Arc<RwLock<Region>>,
    /// Start offset in the composite pattern (0-based, inclusive)
    pub pattern_start: usize,
    /// End offset in the composite pattern (0-based, exclusive)
    pub pattern_end: usize,
    /// Whether this region used N-spacers (non-fixed) or actual sequence (fixed)
    pub is_spacer: bool,
}

/// Result of a composite fitting alignment, with per-region mapping.
#[derive(Debug, Clone)]
pub struct CompositeAlignmentResult {
    /// Raw alignment score from rust-bio
    pub score: i32,
    /// Per-region boundaries in the read.
    pub region_mappings: Vec<RegionMapping>,
}

/// Per-region mapping result from a composite alignment.
#[derive(Debug, Clone)]
pub struct RegionMapping {
    /// The region this mapping corresponds to
    pub region: Arc<RwLock<Region>>,
    /// Start position in the read (inclusive)
    pub read_start: usize,
    /// End position in the read (exclusive)
    pub read_end: usize,
    /// Whether this region's pattern segment was a spacer (N's)
    pub is_spacer: bool,
    /// Match rate (matches / alignment_length) for this region segment
    pub match_rate: f64,
}

/// Match function for composite alignment.
/// N in the pattern (spacer) matches any base with score +1.
/// Otherwise, standard match (+2) / mismatch (-1).
fn composite_match_func(a: u8, b: u8) -> i32 {
    if a == b'N' {
        1 // Spacer position: any base scores +1
    } else if a == b {
        2 // Fixed position match
    } else {
        -1 // Fixed position mismatch
    }
}

/// Sequence aligner for composite fitting alignment of long reads.
/// Aligns a composite pattern (fixed sequences + N-spacers) against a read segment.
pub struct FittingAligner;

impl FittingAligner {
    /// Perform a single composite fitting alignment of the pattern against the read segment.
    /// Returns a CompositeAlignmentResult with per-region read boundaries.
    pub fn align_composite(
        end_sequence: &[u8],
        composite: &CompositePattern,
    ) -> Result<CompositeAlignmentResult> {
        if composite.pattern.is_empty() {
            return Ok(CompositeAlignmentResult {
                score: 0,
                region_mappings: Vec::new(),
            });
        }

        // Build scoring with composite match function.
        // match_scores must be None so rust-bio uses match_fn instead of SIMD fast path.
        let scoring = Scoring {
            gap_open: -1, // same as gap_extend due to indels in long_read sequencing
            gap_extend: -1,
            match_fn: composite_match_func as fn(u8, u8) -> i32,
            match_scores: None,
            xclip_prefix: MIN_SCORE, // Must align entire pattern (no prefix clip)
            xclip_suffix: MIN_SCORE, // Must align entire pattern (no suffix clip)
            yclip_prefix: 0,         // Free read clipping at start (fitting)
            yclip_suffix: 0,         // Free read clipping at end (fitting)
        };

        // Avoid unnecessary memory allocations with expeced size hints.
        let mut aligner = Aligner::with_capacity_and_scoring(
            composite.pattern.len(),
            end_sequence.len(),
            scoring,
        );

        let alignment = aligner.custom(&composite.pattern, end_sequence);

        // Extract per-region mappings by walking the alignment operations
        let region_mappings = Self::extract_region_mappings(&alignment, composite);

        Ok(CompositeAlignmentResult {
            score: alignment.score,
            region_mappings,
        })
    }

    /// Walk alignment operations to extract per-region read boundaries.
    ///
    /// In Custom mode, operations start from (x=0, y=0) with Yclip/Xclip as
    /// explicit operations in the list.
    fn extract_region_mappings(
        alignment: &bio::alignment::Alignment,
        composite: &CompositePattern,
    ) -> Vec<RegionMapping> {
        if composite.spans.is_empty() {
            return Vec::new();
        }

        let num_spans = composite.spans.len();
        let mut span_read_start: Vec<Option<usize>> = vec![None; num_spans];
        let mut span_read_end: Vec<usize> = vec![0; num_spans];
        let mut span_matches: Vec<usize> = vec![0; num_spans];
        let mut span_total: Vec<usize> = vec![0; num_spans];

        // In Custom mode, operations start from (0, 0)
        let mut x_pos: usize = 0;
        let mut y_pos: usize = 0;

        for op in &alignment.operations {
            match op {
                AlignmentOperation::Match | AlignmentOperation::Subst => {
                    if x_pos < composite.pattern.len() {
                        let span_idx = composite.pos_to_span[x_pos];
                        if span_read_start[span_idx].is_none() {
                            span_read_start[span_idx] = Some(y_pos);
                        }
                        span_read_end[span_idx] = y_pos + 1;
                        span_total[span_idx] += 1;
                        if matches!(op, AlignmentOperation::Match) {
                            span_matches[span_idx] += 1;
                        }
                    }
                    x_pos += 1;
                    y_pos += 1;
                }
                AlignmentOperation::Ins => {
                    // Gap in read (y): consumes 1 from pattern, 0 from read
                    if x_pos < composite.pattern.len() {
                        let span_idx = composite.pos_to_span[x_pos];
                        if span_read_start[span_idx].is_none() {
                            span_read_start[span_idx] = Some(y_pos);
                        }
                        // Keep span_read_end at least y_pos so fully-deleted spans
                        // get a zero-width region at the correct read position
                        span_read_end[span_idx] = span_read_end[span_idx].max(y_pos);
                        span_total[span_idx] += 1;
                    }
                    x_pos += 1;
                }
                AlignmentOperation::Del => {
                    // Gap in pattern (x): consumes 0 from pattern, 1 from read
                    if x_pos < composite.pattern.len() {
                        let span_idx = composite.pos_to_span[x_pos];
                        if span_read_start[span_idx].is_none() {
                            span_read_start[span_idx] = Some(y_pos);
                        }
                        span_read_end[span_idx] = y_pos + 1;
                        span_total[span_idx] += 1;
                    } else if num_spans > 0 {
                        // Deletion past end of pattern: attribute to last span
                        let span_idx = num_spans - 1;
                        span_read_end[span_idx] = y_pos + 1;
                        span_total[span_idx] += 1;
                    }
                    y_pos += 1;
                }
                AlignmentOperation::Yclip(len) => {
                    y_pos += len;
                }
                AlignmentOperation::Xclip(len) => {
                    x_pos += len;
                }
            }
        }

        // Build RegionMapping for each span
        composite
            .spans
            .iter()
            .enumerate()
            .map(|(idx, span)| {
                let read_start = span_read_start[idx].unwrap_or(0);
                let read_end = span_read_end[idx];
                let match_rate = if span_total[idx] > 0 {
                    span_matches[idx] as f64 / span_total[idx] as f64
                } else {
                    0.0
                };

                RegionMapping {
                    region: span.region.clone(),
                    read_start,
                    read_end,
                    is_spacer: span.is_spacer,
                    match_rate,
                }
            })
            .collect()
    }
}

/// Calculate fitting alignment distance between short sequence and long sequence
/// Used for comparing a barcode candidate against a whitelist entry.
pub fn fitting_alignment_distance(short_seq: &[u8], long_seq: &[u8]) -> usize {
    let m = short_seq.len();

    if m == 0 {
        return 0; // Empty short sequence can always match
    }
    if long_seq.is_empty() {
        return m; // Cannot fit non-empty short sequence in empty long sequence
    }

    let dp = fitting_alignment_matrix(short_seq, long_seq);

    // Return minimum value in the last row (fitting alignment)
    dp[m].iter().min().copied().unwrap_or(m)
}

/// Calculate fitting-alignment distance and the observed span in `long_seq`.
///
/// Barcode candidate filtering should keep using [`fitting_alignment_distance`];
/// this function is intended only for the first tied-best candidate after
/// filtering has completed. The span is the `ystart..yend` returned by
/// [`Aligner::semiglobal`], which aligns all of `short_seq` to a substring of
/// `long_seq`.
pub fn fitting_alignment_with_span(short_seq: &[u8], long_seq: &[u8]) -> (usize, Range<usize>) {
    let m = short_seq.len();
    let n = long_seq.len();

    if m == 0 {
        return (0, 0..0);
    }
    if n == 0 {
        return (m, 0..0);
    }

    // With this linear scoring, score == -Levenshtein distance:
    // mismatch = -1 and a gap of length k costs -1 + -1 * (k - 1) = -k.
    let scoring = Scoring::from_scores(-1, -1, 0, -1);
    let mut aligner = Aligner::with_capacity_and_scoring(m, n, scoring);
    let alignment = aligner.semiglobal(short_seq, long_seq);
    debug_assert!(alignment.score <= 0);

    (
        (-alignment.score) as usize,
        alignment.ystart..alignment.yend,
    )
}

fn fitting_alignment_matrix(short_seq: &[u8], long_seq: &[u8]) -> Vec<Vec<usize>> {
    let m = short_seq.len();
    let n = long_seq.len();
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
            let cost = if short_seq[i - 1] == long_seq[j - 1] {
                0
            } else {
                1
            };
            dp[i][j] = (dp[i - 1][j] + 1)
                .min(dp[i][j - 1] + 1)
                .min(dp[i - 1][j - 1] + cost);
        }
    }

    dp
}

#[cfg(test)]
mod tests {
    use super::*;
    use seqspec::{RegionType, SequenceType};

    fn create_fixed_region(id: &str, sequence: &str) -> Arc<RwLock<Region>> {
        let len = sequence.len() as u32;
        Arc::new(RwLock::new(Region {
            region_id: id.to_string(),
            region_type: RegionType::Linker,
            name: id.to_string(),
            sequence_type: SequenceType::Fixed,
            sequence: sequence.to_string(),
            min_len: len,
            max_len: len,
            onlist: None,
            subregions: vec![],
        }))
    }

    fn create_barcode_region(id: &str, len: u32) -> Arc<RwLock<Region>> {
        Arc::new(RwLock::new(Region {
            region_id: id.to_string(),
            region_type: RegionType::Barcode,
            name: id.to_string(),
            sequence_type: SequenceType::Onlist,
            sequence: "N".repeat(len as usize),
            min_len: len,
            max_len: len,
            onlist: None,
            subregions: vec![],
        }))
    }

    fn create_test_composite() -> CompositePattern {
        // Pattern: ACGTACGTACGTACGT (16bp fixed) + NNNNNNNNNN (10bp spacer) + TGCATGCATGCATGCA (16bp fixed)
        let pattern = b"ACGTACGTACGTACGTNNNNNNNNNNTGCATGCATGCATGCA".to_vec();
        let spans = vec![
            CompositeRegionSpan {
                region: create_fixed_region("r1", "ACGTACGTACGTACGT"),
                pattern_start: 0,
                pattern_end: 16,
                is_spacer: false,
            },
            CompositeRegionSpan {
                region: create_barcode_region("bc1", 10),
                pattern_start: 16,
                pattern_end: 26,
                is_spacer: true,
            },
            CompositeRegionSpan {
                region: create_fixed_region("r2", "TGCATGCATGCATGCA"),
                pattern_start: 26,
                pattern_end: 42,
                is_spacer: false,
            },
        ];
        let mut pos_to_span = vec![0usize; pattern.len()];
        for (idx, span) in spans.iter().enumerate() {
            for p in span.pattern_start..span.pattern_end {
                pos_to_span[p] = idx;
            }
        }
        CompositePattern {
            pattern,
            spans,
            pos_to_span,
        }
    }

    #[test]
    fn test_composite_alignment_perfect_match() {
        let composite = create_test_composite();

        let read = b"ACGTACGTACGTACGTAAAAAAAAAATGCATGCATGCATGCA";
        let result = FittingAligner::align_composite(read, &composite).unwrap();

        assert_eq!(result.region_mappings.len(), 3);
        // Fixed region 1
        assert_eq!(result.region_mappings[0].read_start, 0);
        assert_eq!(result.region_mappings[0].read_end, 16);
        assert!((result.region_mappings[0].match_rate - 1.0).abs() < 0.01);
        // Barcode spacer
        assert_eq!(result.region_mappings[1].read_start, 16);
        assert_eq!(result.region_mappings[1].read_end, 26);
        // Fixed region 2
        assert_eq!(result.region_mappings[2].read_start, 26);
        assert_eq!(result.region_mappings[2].read_end, 42);
        assert!((result.region_mappings[2].match_rate - 1.0).abs() < 0.01);
    }

    #[test]
    fn test_composite_alignment_fitting_in_longer_read() {
        // Composite pattern should fit within a longer read (Yclip at both ends)
        let composite = create_test_composite();

        // Read with flanking sequences
        let read = b"GGGGGACGTACGTACGTACGTCCCCCCCCCCGTGCATGCATGCATGCATTTTT";
        let result = FittingAligner::align_composite(read, &composite).unwrap();

        assert_eq!(result.region_mappings.len(), 3);
        // Fixed region 1 should be found after the flanking G's
        assert_eq!(result.region_mappings[0].read_start, 5);
        assert_eq!(result.region_mappings[0].read_end, 21);
        // Fixed region 2 should end before the trailing T's
        assert_eq!(result.region_mappings[2].read_end, 48);
        // Both fixed regions should have high match rates
        assert!(result.region_mappings[0].match_rate >= 0.9);
        assert!(result.region_mappings[2].match_rate >= 0.9);
    }

    #[test]
    fn test_composite_alignment_with_indels() {
        // Test with 1 mismatch in fixed region
        let composite = create_test_composite();

        // 1 mismatch in fixed1 (T->G at end), perfect barcode + 2 del in fixed2
        let read = b"ACGTACGTACGTACGGAAAAAAAAAATGCATGCATCATCAGGGG";
        let result = FittingAligner::align_composite(read, &composite).unwrap();

        assert_eq!(result.region_mappings.len(), 3);
        assert_eq!(result.region_mappings[0].read_start, 0);
        assert_eq!(result.region_mappings[0].read_end, 16);
        // Match rate should be 15/16 = 0.9375
        assert!(
            (result.region_mappings[0].match_rate - 15.0 / 16.0).abs() < 0.01,
            "Expected ~0.9375, got {}",
            result.region_mappings[0].match_rate
        );
        // Barcode region
        assert_eq!(result.region_mappings[1].read_start, 16);
        assert_eq!(result.region_mappings[1].read_end, 26);

        assert_eq!(result.region_mappings[2].read_start, 26);
        assert_eq!(result.region_mappings[2].read_end, 40);
    }

    #[test]
    fn test_fitting_alignment_span() {
        assert_eq!(fitting_alignment_with_span(b"ATCG", b"GGATCGCC"), (0, 2..6));
        assert_eq!(
            fitting_alignment_with_span(b"ATCG", b"GGATCCGCC"),
            (1, 2..5)
        );
        assert_eq!(fitting_alignment_with_span(b"ATCG", b"GGACGCC"), (1, 2..5));
    }

    #[test]
    fn test_fitting_alignment_span_uses_bio_tied_coordinates() {
        assert_eq!(fitting_alignment_with_span(b"AT", b"ATAT"), (0, 2..4));
    }
}
