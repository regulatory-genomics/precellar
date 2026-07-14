//! Seed-and-extend extraction of short insertions from sequencing reads.
//!
//! [`InsertionExtractor`] is intended for reads with the following layout:
//!
//! ```text
//! 5' fixed flank | variable insertion | 3' fixed flank
//! ```
//!
//! Construction precomputes every `kmer_size`-base substring of both fixed
//! flanks and stores them in an Aho-Corasick automaton. Extraction then scans a
//! read once for exact seed matches. Each seed identifies one junction; the
//! other junction is searched fuzzily in a small window around the expected
//! insertion length. This avoids running dynamic programming across reads that
//! contain unrelated genomic DNA.
//!
//! Coordinates in this module are byte offsets. This is appropriate for DNA
//! sequences represented by ASCII bytes and keeps extraction allocation-free
//! until a valid insertion is returned.

use aho_corasick::{AhoCorasick, AhoCorasickBuilder};

/// Maximum displacement searched around the expected insertion boundary.
const SEARCH_RADIUS: usize = 5;
/// Minimum aggregate exact matches required across both flanks.
const MIN_MATCHING_FLANK_BASES: usize = 20;

/// Identifies which fixed flank supplied an automaton seed.
#[derive(Clone, Copy)]
enum Flank {
    FivePrime,
    ThreePrime,
}

/// Metadata associated with one k-mer pattern in the automaton.
///
/// `offset` is normalized to point from the seed start to the insertion
/// junction for a 5' seed, and from the junction to the seed start for a 3'
/// seed. Keeping this direction-specific coordinate here makes extraction use
/// the same `match.start()` anchor while preserving the correct arithmetic for
/// either orientation.
#[derive(Clone, Copy)]
struct Seed {
    flank: Flank,
    offset: usize,
}

/// Extracts a short variable insertion between two fixed flanking sequences.
///
/// The extractor is reusable across reads and should be constructed once per
/// pair of flanks. Its automaton contains all k-mers from the 5' and 3' flanks,
/// so construction performs the relatively expensive precomputation while
/// [`extract`](Self::extract) performs only a linear seed scan plus a bounded
/// local comparison.
///
/// The flanks are supplied as `&str` for convenient construction, but reads are
/// processed as `&[u8]`. The extractor does not validate alphabet characters or
/// normalize case; matching is exact and therefore `b'A'` and `b'a'` differ.
pub struct InsertionExtractor {
    automaton: Option<AhoCorasick>,
    seeds: Vec<Seed>,
    flank_5: Vec<u8>,
    flank_3: Vec<u8>,
    expected_len: usize,
    mismatch_tolerance: usize,
}

impl InsertionExtractor {
    /// Builds an extractor for a pair of fixed flanks.
    ///
    /// `flank_5` and `flank_3` are written in read orientation. `kmer_size`
    /// controls the exact seed length, `expected_len` is the anticipated
    /// insertion length, and `mismatch_tolerance` is the maximum number of
    /// mismatches among all sequenced flank bases used for validation.
    ///
    /// Every possible k-mer in each flank is inserted into the automaton. A
    /// zero `kmer_size`, or a k-mer size longer than both flanks, creates an
    /// extractor that simply returns `None` for every read. Empty flanks are
    /// accepted; their available boundary is empty and consequently contributes
    /// no mismatches.
    ///
    /// # Example
    ///
    /// ```
    /// use precellar::utils::insertion_extractor::InsertionExtractor;
    ///
    /// let extractor = InsertionExtractor::new(
    ///     "ACGTCAGTGGCA",
    ///     "TTGGAACCTTGG",
    ///     12,
    ///     4,
    ///     1,
    /// );
    /// assert_eq!(
    ///     extractor.extract(b"ACGTCAGTGGCAACGTTTGGAACCTTAG"),
    ///     Some(b"ACGT".to_vec()),
    /// );
    /// ```
    pub fn new(
        flank_5: &str,
        flank_3: &str,
        kmer_size: usize,
        expected_len: usize,
        mismatch_tolerance: usize,
    ) -> Self {
        let flank_5 = flank_5.as_bytes();
        let flank_3 = flank_3.as_bytes();
        let mut patterns = Vec::new();
        let mut seeds = Vec::new();

        if kmer_size > 0 {
            for (flank, sequence) in [(Flank::FivePrime, flank_5), (Flank::ThreePrime, flank_3)] {
                if sequence.len() >= kmer_size {
                    for offset in 0..=sequence.len() - kmer_size {
                        patterns.push(sequence[offset..offset + kmer_size].to_vec());
                        let junction_offset = match flank {
                            Flank::FivePrime => sequence.len() - offset,
                            Flank::ThreePrime => offset,
                        };
                        seeds.push(Seed {
                            flank,
                            offset: junction_offset,
                        });
                    }
                }
            }
        }

        let automaton = if patterns.is_empty() {
            None
        } else {
            Some(
                AhoCorasickBuilder::new()
                    .build(patterns)
                    .expect("valid k-mer patterns should build an Aho-Corasick automaton"),
            )
        };

        Self {
            automaton,
            seeds,
            flank_5: flank_5.to_vec(),
            flank_3: flank_3.to_vec(),
            expected_len,
            mismatch_tolerance,
        }
    }

    /// Extracts the insertion from `read_seq`, if a valid structure is found.
    ///
    /// The Aho-Corasick matches are tried in scan order until one produces a
    /// valid extraction. This is important when a flank contains repeated
    /// k-mers: an early match can describe a different physical occurrence
    /// than the one bordering the insertion. A 5' seed fixes the insertion
    /// start and searches for the 3' junction. A 3' seed fixes the insertion
    /// end and searches backward for the 5' junction.
    ///
    /// The fuzzy search examines junctions from
    /// `expected_len - 5` through `expected_len + 5`, clamped at zero for
    /// underflow. For each candidate, it compares all available bases adjacent
    /// to both junctions. A candidate must have at least 20 aggregate exact
    /// flank matches, regardless of which flank supplies those matches.
    /// Candidates are ranked by matching bases, then mismatch count, then
    /// distance from the expected junction. The returned vector contains only
    /// the insertion bytes; flanking bases are not included.
    ///
    /// Returns `None` when no seed is present, coordinates fall outside the
    /// read, a required boundary exceeds the read, the mismatch tolerance is
    /// exceeded, or the resulting coordinates would be reversed.
    pub fn extract(&self, read_seq: &[u8]) -> Option<Vec<u8>> {
        let automaton = self.automaton.as_ref()?;
        automaton
            .find_iter(read_seq)
            .find_map(|matched| self.extract_from_match(read_seq, matched))
    }

    fn extract_from_match(&self, read_seq: &[u8], matched: aho_corasick::Match) -> Option<Vec<u8>> {
        let seed = *self.seeds.get(matched.pattern().as_usize())?;

        match seed.flank {
            Flank::FivePrime => {
                // `offset` is the distance from the seed start to the 5' junction.
                // Thus match.start() + offset is the exclusive end of the fixed
                // 5' flank and the inclusive start of the insertion.
                let insert_start = matched.start().checked_add(seed.offset)?;
                let expected_end = insert_start.checked_add(self.expected_len)?;
                let insert_end = self.find_boundary(read_seq, insert_start, expected_end)?;
                if insert_end < insert_start {
                    return None;
                }
                Some(read_seq[insert_start..insert_end].to_vec())
            }
            Flank::ThreePrime => {
                // The seed starts `offset` bases into the 3' flank, so its junction
                // is match.start - offset. This is the insertion's exclusive end.
                let insert_end = matched.start().checked_sub(seed.offset)?;
                let expected_start = insert_end.checked_sub(self.expected_len)?;
                let insert_start =
                    self.find_boundary_reverse(read_seq, expected_start, insert_end)?;
                if insert_start > insert_end {
                    return None;
                }
                Some(read_seq[insert_start..insert_end].to_vec())
            }
        }
    }

    /// Searches forward-facing 3' junction candidates and scores both flanks.
    fn find_boundary(&self, read: &[u8], insert_start: usize, expected: usize) -> Option<usize> {
        let first = expected.saturating_sub(SEARCH_RADIUS);
        let last = expected.checked_add(SEARCH_RADIUS)?;
        (first..=last)
            .filter_map(|junction| {
                let (matches, mismatches) = self.score_boundaries(read, insert_start, junction);
                (matches >= MIN_MATCHING_FLANK_BASES && mismatches <= self.mismatch_tolerance)
                    .then_some((matches, mismatches, junction.abs_diff(expected), junction))
            })
            .max_by_key(|&(matches, mismatches, distance, junction)| {
                (
                    matches,
                    usize::MAX - mismatches,
                    usize::MAX - distance,
                    junction,
                )
            })
            .map(|(_, _, _, junction)| junction)
    }

    /// Searches backward-facing 5' junction candidates and scores both flanks.
    fn find_boundary_reverse(
        &self,
        read: &[u8],
        expected: usize,
        insert_end: usize,
    ) -> Option<usize> {
        let first = expected.saturating_sub(SEARCH_RADIUS);
        let last = expected.checked_add(SEARCH_RADIUS)?;
        (first..=last)
            .filter_map(|junction| {
                let (matches, mismatches) = self.score_boundaries(read, junction, insert_end);
                (matches >= MIN_MATCHING_FLANK_BASES && mismatches <= self.mismatch_tolerance)
                    .then_some((matches, mismatches, junction.abs_diff(expected), junction))
            })
            .max_by_key(|&(matches, mismatches, distance, junction)| {
                (
                    matches,
                    usize::MAX - mismatches,
                    usize::MAX - distance,
                    junction,
                )
            })
            .map(|(_, _, _, junction)| junction)
    }

    /// Counts matches and mismatches in all available flank bases adjacent to
    /// the proposed insertion. Missing sequence at either read end is ignored.
    fn score_boundaries(
        &self,
        read: &[u8],
        insert_start: usize,
        insert_end: usize,
    ) -> (usize, usize) {
        if insert_start > insert_end || insert_end > read.len() {
            return (0, usize::MAX);
        }

        let left_len = insert_start.min(self.flank_5.len());
        let left_read_start = insert_start - left_len;
        let left_expected_start = self.flank_5.len() - left_len;
        let (left_matches, left_mismatches) = match_slices(
            &read[left_read_start..insert_start],
            &self.flank_5[left_expected_start..],
        );

        let right_len = read
            .len()
            .saturating_sub(insert_end)
            .min(self.flank_3.len());
        let (right_matches, right_mismatches) = match_slices(
            &read[insert_end..insert_end + right_len],
            &self.flank_3[..right_len],
        );

        (
            left_matches + right_matches,
            left_mismatches + right_mismatches,
        )
    }
}

/// Counts exact matches and substitutions between equally sized byte slices.
///
/// The caller determines how many bases are available before comparison.
fn match_slices(left: &[u8], right: &[u8]) -> (usize, usize) {
    left.iter()
        .zip(right)
        .fold((0, 0), |(matches, mismatches), (left, right)| {
            if left == right {
                (matches + 1, mismatches)
            } else {
                (matches, mismatches + 1)
            }
        })
}

#[cfg(test)]
mod tests {
    use super::InsertionExtractor;

    #[test]
    fn extracts_from_five_prime_anchor() {
        let extractor = InsertionExtractor::new("ACGTCAGTGGCA", "TTGGAACCTTGG", 12, 4, 0);
        assert_eq!(
            extractor.extract(b"ACGTCAGTGGCAACGTTTGGAACCTTGG"),
            Some(b"ACGT".to_vec())
        );
    }

    #[test]
    fn extracts_from_three_prime_anchor_with_mismatch() {
        let extractor = InsertionExtractor::new("ACGTCAGTGGCA", "TTGGAACCTTGG", 12, 4, 1);
        assert_eq!(
            extractor.extract(b"ACGTCAGTGGCAACGTTTGGAACCTTAG"),
            Some(b"ACGT".to_vec())
        );
    }

    #[test]
    fn extracts_from_three_prime_anchor_with_repeat() {
        let extractor = InsertionExtractor::new(
            "ACGTCAGTGGCAGGGGGGGGGGGGGACGTCAGTGGCA",
            "TTGGAACCTTGG",
            12,
            4,
            0,
        );
        assert_eq!(
            extractor.extract(b"ACGTCAGTGGCAACGTTTGGAACCTTGG"),
            Some(b"ACGT".to_vec())
        );
    }

    #[test]
    fn accepts_when_only_one_flank_is_sequenced() {
        let extractor =
            InsertionExtractor::new("ACGTCAGTGGCACTGATCGTACGA", "TTGGAACCTTGG", 12, 4, 0);
        assert_eq!(
            extractor.extract(b"ACGTCAGTGGCACTGATCGTACGAACGA"),
            Some(b"ACGA".to_vec())
        );
    }

    #[test]
    fn rejects_missing_or_invalid_boundary() {
        let extractor = InsertionExtractor::new("ACGTCAGTGGCA", "TTGGAACCTTGG", 12, 4, 0);
        assert_eq!(extractor.extract(b"GATTGATTGATT"), None);
        assert_eq!(extractor.extract(b"ACGTCAGTGGCAGATCAACCTTGG"), None);
    }
}
