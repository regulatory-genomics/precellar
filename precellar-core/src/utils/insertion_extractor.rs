//! Seeded extraction of insertions from reads containing a fixed-sequence template.
//!
//! A template consists of ordered fixed sequences and insertions between them:
//!
//! ```text
//! fixed_0 | insertion_0 | fixed_1 | insertion_1 | fixed_2
//! ```
//!
//! The extractor scans a read for k-mer seeds from the fixed sequences. Each
//! seed is extended through the template in both directions. The first seed
//! that produces a valid extension wins; later seeds are not considered after a
//! successful extension. Substitutions and gaps are allowed while aligning
//! fixed sequences, but insertion intervals are returned as contiguous read
//! ranges.

use aho_corasick::{AhoCorasick, AhoCorasickBuilder};
use anyhow::{bail, Result};
use std::collections::HashMap;
use std::ops::Range;

/// Metadata for one k-mer occurrence in a fixed sequence.
#[derive(Clone, Copy)]
struct Seed {
    fixed_index: usize,
    offset: usize,
}

/// A placement of an observed complete or read-edge partial fixed sequence.
#[derive(Clone, Copy)]
struct Placement {
    fixed_index: usize,
    start: usize,
    end: usize,
    edits: usize,
    observed: usize,
}

/// A fixed-sequence alignment anchored at one read coordinate.
#[derive(Clone, Copy)]
struct Span {
    consumed: usize,
    edits: usize,
    observed: usize,
}

/// Extracts insertions from an ordered template of fixed sequences.
pub struct InsertionExtractor {
    automaton: Option<AhoCorasick>,
    seeds: Vec<Vec<Seed>>,
    flanks: Vec<Vec<u8>>,
    expected_lens: Vec<usize>,
    kmer_size: usize,
    max_fixed_edit_rate: f64,
    max_alignment_edits: usize,
}

impl InsertionExtractor {
    /// Builds an extractor for an ordered fixed-sequence template.
    ///
    /// `flanks` contains the fixed sequences in template order. There must be
    /// one entry in `expected_lens` for every gap between adjacent fixed
    /// sequences. The insertion lengths are exact. Fixed-sequence alignment
    /// permits substitutions and gaps. The global edit rate is calculated as
    /// edits divided by observed fixed-sequence bases.
    /// Insertions themselves are extracted as contiguous read ranges.
    pub fn new(
        flanks: Vec<String>,
        expected_lens: Vec<usize>,
        kmer_size: usize,
        max_fixed_edit_rate: f64,
    ) -> Result<Self> {
        if flanks.len() < 2 {
            bail!("at least two fixed sequences are required");
        }
        if expected_lens.len() != flanks.len() - 1 {
            bail!(
                "expected_lens must contain one length per insertion gap (expected {}, got {})",
                flanks.len() - 1,
                expected_lens.len()
            );
        }
        if kmer_size == 0 {
            bail!("kmer_size must be greater than zero");
        }
        if !max_fixed_edit_rate.is_finite() || !(0.0..=1.0).contains(&max_fixed_edit_rate) {
            bail!("max_fixed_edit_rate must be finite and between zero and one");
        }

        let flanks: Vec<Vec<u8>> = flanks.into_iter().map(|flank| flank.into_bytes()).collect();
        if flanks.iter().all(|flank| flank.len() < kmer_size) {
            bail!("kmer_size must not exceed every fixed-sequence length");
        }
        let total_fixed_bases: usize = flanks.iter().map(Vec::len).sum();
        let max_alignment_edits = (max_fixed_edit_rate * total_fixed_bases as f64).floor() as usize;
        let mut patterns = Vec::new();
        let mut seeds_by_pattern: Vec<Vec<Seed>> = Vec::new();
        let mut pattern_ids = HashMap::<Vec<u8>, usize>::new();

        for (fixed_index, sequence) in flanks.iter().enumerate() {
            if sequence.len() < kmer_size {
                continue;
            }
            for offset in 0..=sequence.len() - kmer_size {
                let pattern = sequence[offset..offset + kmer_size].to_vec();
                let pattern_id = if let Some(&pattern_id) = pattern_ids.get(&pattern) {
                    pattern_id
                } else {
                    let pattern_id = patterns.len();
                    pattern_ids.insert(pattern.clone(), pattern_id);
                    patterns.push(pattern);
                    seeds_by_pattern.push(Vec::new());
                    pattern_id
                };
                seeds_by_pattern[pattern_id].push(Seed {
                    fixed_index,
                    offset,
                });
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

        Ok(Self {
            automaton,
            seeds: seeds_by_pattern,
            flanks,
            expected_lens,
            kmer_size,
            max_fixed_edit_rate,
            max_alignment_edits,
        })
    }

    /// Returns the number of insertion gaps in the template.
    pub fn num_gaps(&self) -> usize {
        self.expected_lens.len()
    }

    /// Returns the configured insertion lengths in template order.
    pub fn expected_lens(&self) -> &[usize] {
        &self.expected_lens
    }

    /// Extracts all insertion sequences from the first successful template
    /// extension.
    pub fn extract(&self, read_seq: &[u8]) -> Option<Vec<(usize, Vec<u8>)>> {
        self.extract_ranges(read_seq).map(|ranges| {
            ranges
                .into_iter()
                .map(|(index, range)| (index, read_seq[range].to_vec()))
                .collect()
        })
    }

    /// Extracts indexed insertion ranges from the first successful template
    /// extension. A failed seed is skipped and the next seed is tried. Once a
    /// seed extends successfully through the template, its result is returned
    /// immediately and no later seed is considered.
    pub fn extract_ranges(&self, read_seq: &[u8]) -> Option<Vec<(usize, Range<usize>)>> {
        let automaton = self.automaton.as_ref()?;

        for matched in automaton.find_overlapping_iter(read_seq) {
            let pattern_seeds = self.seeds.get(matched.pattern().as_usize())?;
            for &seed in pattern_seeds {
                if let Some(placements) = self.extend_from_seed(read_seq, matched, seed) {
                    let ranges = self.insertion_ranges(read_seq.len(), &placements);
                    if !ranges.is_empty() {
                        return Some(ranges);
                    }
                }
            }
        }
        None
    }

    fn extend_from_seed(
        &self,
        read: &[u8],
        matched: aho_corasick::Match,
        seed: Seed,
    ) -> Option<Vec<Placement>> {
        let fixed = self.flanks.get(seed.fixed_index)?;
        let seed_end = matched.end();
        let seed_start = matched.start();
        let seed_end_in_fixed = seed.offset.checked_add(self.kmer_size)?;
        if seed_end_in_fixed > fixed.len() {
            return None;
        }

        let prefix = align_backward_at_edge(
            &fixed[..seed.offset],
            &read[..seed_start],
            self.max_alignment_edits,
        )?;
        let suffix = align_forward_at_edge(
            &fixed[seed_end_in_fixed..],
            &read[seed_end..],
            self.max_alignment_edits,
        )?;
        let edits = prefix.edits.checked_add(suffix.edits)?;
        let observed = prefix.observed + self.kmer_size + suffix.observed;

        let placement = Placement {
            fixed_index: seed.fixed_index,
            start: seed_start.checked_sub(prefix.consumed)?,
            end: seed_end.checked_add(suffix.consumed)?,
            edits,
            observed,
        };

        let mut backward = Vec::with_capacity(seed.fixed_index);
        let mut current = placement;
        let mut total_edits = edits;
        let mut total_observed = observed;

        while current.fixed_index > 0 {
            let gap_index = current.fixed_index - 1;
            let Some(expected_end) = current.start.checked_sub(self.expected_lens[gap_index])
            else {
                break;
            };
            if expected_end == 0 {
                break;
            }

            let fixed_index = current.fixed_index - 1;
            let fixed = &self.flanks[fixed_index];
            let previous = if expected_end < fixed.len() {
                let span =
                    align_partial_backward(fixed, &read[..expected_end], self.max_alignment_edits)?;
                Placement {
                    fixed_index,
                    start: 0,
                    end: expected_end,
                    edits: span.edits,
                    observed: span.observed,
                }
            } else {
                self.backward_placement(read, fixed_index, expected_end)?
            };
            total_edits = total_edits.checked_add(previous.edits)?;
            total_observed = total_observed.checked_add(previous.observed)?;
            backward.push(previous);
            current = previous;

            if previous.start == 0 {
                break;
            }
        }

        backward.reverse();
        backward.push(placement);
        current = placement;

        while current.fixed_index + 1 < self.flanks.len() {
            let gap_index = current.fixed_index;
            let expected_start = current.end.checked_add(self.expected_lens[gap_index])?;
            if expected_start >= read.len() {
                break;
            }

            let fixed_index = current.fixed_index + 1;
            let fixed = &self.flanks[fixed_index];
            let next = if read.len() - expected_start < fixed.len() {
                let span = align_partial_forward(
                    fixed,
                    &read[expected_start..],
                    self.max_alignment_edits,
                )?;
                Placement {
                    fixed_index,
                    start: expected_start,
                    end: read.len(),
                    edits: span.edits,
                    observed: span.observed,
                }
            } else {
                self.forward_placement(read, fixed_index, expected_start)?
            };
            total_edits = total_edits.checked_add(next.edits)?;
            total_observed = total_observed.checked_add(next.observed)?;
            backward.push(next);
            current = next;

            if current.end == read.len() {
                break;
            }
        }

        if total_observed == 0
            || (total_edits as f64) > self.max_fixed_edit_rate * total_observed as f64
        {
            return None;
        }
        Some(backward)
    }

    fn forward_placement(
        &self,
        read: &[u8],
        fixed_index: usize,
        start: usize,
    ) -> Option<Placement> {
        let fixed = &self.flanks[fixed_index];
        let span = align_forward(fixed, &read[start..], self.max_alignment_edits)?;
        Some(Placement {
            fixed_index,
            start,
            end: start.checked_add(span.consumed)?,
            edits: span.edits,
            observed: span.observed,
        })
    }

    fn backward_placement(&self, read: &[u8], fixed_index: usize, end: usize) -> Option<Placement> {
        let fixed = &self.flanks[fixed_index];
        let span = align_backward(fixed, &read[..end], self.max_alignment_edits)?;
        Some(Placement {
            fixed_index,
            start: end.checked_sub(span.consumed)?,
            end,
            edits: span.edits,
            observed: span.observed,
        })
    }

    fn insertion_ranges(
        &self,
        read_len: usize,
        placements: &[Placement],
    ) -> Vec<(usize, Range<usize>)> {
        let first = placements.first().expect("seed creates a placement");
        let last = placements.last().expect("seed creates a placement");
        let mut ranges = Vec::with_capacity(placements.len());

        if first.fixed_index > 0 {
            let index = first.fixed_index - 1;
            let end = first.start;
            let start = end.saturating_sub(self.expected_lens[index]);
            if end - start == self.expected_lens[index] && start < end {
                ranges.push((index, start..end));
            }
        }

        for pair in placements.windows(2) {
            let left = pair[0];
            let right = pair[1];
            if left.end < right.start
                && right.start - left.end == self.expected_lens[left.fixed_index]
            {
                ranges.push((left.fixed_index, left.end..right.start));
            }
        }

        if last.fixed_index + 1 < self.flanks.len() {
            let index = last.fixed_index;
            let start = last.end;
            let end = start
                .saturating_add(self.expected_lens[index])
                .min(read_len);
            if end - start == self.expected_lens[index] && start < end {
                ranges.push((index, start..end));
            }
        }

        ranges
    }
}

/// Aligns a complete fixed sequence to the beginning of `read`.
fn align_forward(fixed: &[u8], read: &[u8], max_edits: usize) -> Option<Span> {
    edit_span(fixed, read, max_edits)
}

/// Aligns a complete fixed sequence to the end of `read`.
fn align_backward(fixed: &[u8], read: &[u8], max_edits: usize) -> Option<Span> {
    let max_read = fixed.len().saturating_add(max_edits);
    let read = &read[read.len().saturating_sub(max_read)..];
    let fixed: Vec<u8> = fixed.iter().rev().copied().collect();
    let read: Vec<u8> = read.iter().rev().copied().collect();
    edit_span(&fixed, &read, max_edits)
}

/// Aligns a fixed-sequence suffix at the beginning of a read. Fixed bases that
/// precede the read are clipped, while edits in all observed read bases count.
fn align_partial_backward(fixed: &[u8], read: &[u8], max_edits: usize) -> Option<Span> {
    let fixed: Vec<u8> = fixed.iter().rev().copied().collect();
    let read: Vec<u8> = read.iter().rev().copied().collect();
    let aligned = edit_span(&read, &fixed, max_edits)?;
    Some(Span {
        consumed: read.len(),
        edits: aligned.edits,
        observed: aligned.observed,
    })
}

/// Aligns a fixed-sequence prefix at the end of a read. Fixed bases beyond the
/// read are clipped, while edits in all observed read bases count.
fn align_partial_forward(fixed: &[u8], read: &[u8], max_edits: usize) -> Option<Span> {
    let aligned = edit_span(read, fixed, max_edits)?;
    Some(Span {
        consumed: read.len(),
        edits: aligned.edits,
        observed: aligned.observed,
    })
}

fn align_backward_at_edge(fixed: &[u8], read: &[u8], max_edits: usize) -> Option<Span> {
    if read.len() < fixed.len() {
        align_partial_backward(fixed, read, max_edits)
    } else {
        align_backward(fixed, read, max_edits)
    }
}

fn align_forward_at_edge(fixed: &[u8], read: &[u8], max_edits: usize) -> Option<Span> {
    if read.len() < fixed.len() {
        align_partial_forward(fixed, read, max_edits)
    } else {
        align_forward(fixed, read, max_edits)
    }
}

/// Finds the lowest-edit alignment of all of `fixed` against a prefix of
/// `read`. The consumed read length is allowed to differ from the fixed length
/// by the edit budget, which permits gaps in fixed regions.
fn edit_span(fixed: &[u8], read: &[u8], max_edits: usize) -> Option<Span> {
    if fixed.is_empty() {
        return Some(Span {
            consumed: 0,
            edits: 0,
            observed: 0,
        });
    }
    if max_edits == 0 {
        return (read.len() >= fixed.len() && read[..fixed.len()] == *fixed).then_some(Span {
            consumed: fixed.len(),
            edits: 0,
            observed: fixed.len(),
        });
    }

    let max_read = read.len().min(fixed.len().saturating_add(max_edits));
    let columns = max_read + 1;
    #[derive(Clone, Copy)]
    struct Cell {
        edits: usize,
        observed: usize,
    }

    let mut previous = vec![
        Cell {
            edits: 0,
            observed: 0
        };
        columns
    ];
    let mut current = vec![
        Cell {
            edits: 0,
            observed: 0
        };
        columns
    ];
    for (column, value) in previous.iter_mut().enumerate() {
        *value = Cell {
            edits: column,
            observed: 0,
        };
    }

    for (row, expected) in fixed.iter().enumerate() {
        current.fill(Cell {
            edits: 0,
            observed: 0,
        });
        current[0] = Cell {
            edits: row + 1,
            observed: 0,
        };
        for column in 1..columns {
            let diagonal = Cell {
                edits: previous[column - 1].edits + usize::from(*expected != read[column - 1]),
                observed: previous[column - 1].observed + 1,
            };
            let deletion = Cell {
                edits: previous[column].edits + 1,
                observed: previous[column].observed,
            };
            let insertion = Cell {
                edits: current[column - 1].edits + 1,
                observed: current[column - 1].observed,
            };
            current[column] = [diagonal, deletion, insertion]
                .into_iter()
                .min_by_key(|cell| (cell.edits, usize::MAX - cell.observed))
                .expect("alignment transition candidates are non-empty");
        }
        std::mem::swap(&mut previous, &mut current);
    }

    (0..=max_read)
        .filter(|&consumed| previous[consumed].edits <= max_edits)
        .min_by_key(|&consumed| {
            (
                previous[consumed].edits,
                consumed.abs_diff(fixed.len()),
                usize::MAX - consumed,
            )
        })
        .map(|consumed| Span {
            consumed,
            edits: previous[consumed].edits,
            observed: previous[consumed].observed,
        })
}

#[cfg(test)]
mod tests {
    use super::InsertionExtractor;

    fn extractor() -> InsertionExtractor {
        InsertionExtractor::new(
            vec!["ACGTCAGTGGCA".into(), "TTGGAACCTTGG".into()],
            vec![4],
            12,
            0.1,
        )
        .unwrap()
    }

    #[test]
    fn extracts_from_first_fixed_sequence() {
        let extractor = extractor();
        assert_eq!(
            extractor.extract_ranges(b"ACGTCAGTGGCAACGTTTGGAACCTTGG"),
            Some(vec![(0, 12..16)])
        );
        assert_eq!(
            extractor.extract(b"ACGTCAGTGGCAACGTTTGGAACCTTGG"),
            Some(vec![(0, b"ACGT".to_vec())])
        );
        assert_eq!(
            extractor.extract(b"ACGTCAGTGGCAACGT"),
            Some(vec![(0, b"ACGT".to_vec())])
        );
        assert_eq!(extractor.extract(b"ACGTCAGTGGCAAC"), None);
        assert_eq!(extractor.extract(b"ACTTGGAACCTTGG"), None);
        assert_eq!(
            extractor.extract(b"ACGTCAGTGGCAACGTTTGGAA"),
            Some(vec![(0, b"ACGT".to_vec())])
        );
        assert_eq!(
            extractor.extract(b"ACGTCAGTGGCAACGTTTGGAT"),
            Some(vec![(0, b"ACGT".to_vec())])
        );
        assert_eq!(extractor.extract(b"ACGTCAGTGGCAACGTTTGGTT"), None);
    }

    #[test]
    fn counts_edits_in_a_partial_fixed_sequence_at_read_start() {
        let extractor = extractor();
        assert_eq!(
            extractor.extract(b"GTGGCAACGTTTGGAACCTTGG"),
            Some(vec![(0, b"ACGT".to_vec())])
        );
        assert_eq!(
            extractor.extract(b"GTGGCTACGTTTGGAACCTTGG"),
            Some(vec![(0, b"ACGT".to_vec())])
        );
        assert_eq!(extractor.extract(b"GTGGTTACGTTTGGAACCTTGG"), None);
    }

    #[test]
    fn extracts_multiple_insertions() {
        let extractor = InsertionExtractor::new(
            vec!["AAAA".into(), "CCCC".into(), "GGGG".into()],
            vec![2, 3],
            3,
            0.0,
        )
        .unwrap();
        assert_eq!(
            extractor.extract_ranges(b"AAAATTCCCCTTTGGGG"),
            Some(vec![(0, 4..6), (1, 10..13)])
        );
    }

    #[test]
    fn omits_truncated_edge_insertions_but_keeps_full_insertions() {
        let extractor = InsertionExtractor::new(
            vec!["AAAA".into(), "CCCC".into(), "GGGG".into()],
            vec![2, 3],
            3,
            0.0,
        )
        .unwrap();
        assert_eq!(
            extractor.extract_ranges(b"AAAATTCCCCA"),
            Some(vec![(0, 4..6)])
        );
        assert_eq!(
            extractor.extract_ranges(b"TCCCCAAAGGGG"),
            Some(vec![(1, 5..8)])
        );
    }

    #[test]
    fn a_failed_seed_does_not_block_a_later_seed() {
        let extractor =
            InsertionExtractor::new(vec!["AAAA".into(), "CCCC".into()], vec![2], 3, 0.0).unwrap();
        assert_eq!(
            extractor.extract_ranges(b"AAATTTAAAATTCCCC"),
            Some(vec![(0, 10..12)])
        );
    }

    #[test]
    fn missing_template_ends_are_allowed() {
        let extractor = InsertionExtractor::new(
            vec!["AAAA".into(), "CCCC".into(), "GGGG".into()],
            vec![2, 3],
            3,
            0.0,
        )
        .unwrap();
        assert_eq!(
            extractor.extract_ranges(b"AAAATTCCCC"),
            Some(vec![(0, 4..6)])
        );
    }

    #[test]
    fn infers_an_insertion_from_one_observed_boundary() {
        let extractor =
            InsertionExtractor::new(vec!["AAAA".into(), "CCCC".into()], vec![4], 3, 0.0).unwrap();
        assert_eq!(extractor.extract_ranges(b"AAAAGGGG"), Some(vec![(0, 4..8)]));
        assert_eq!(extractor.extract_ranges(b"GGGGCCCC"), Some(vec![(0, 0..4)]));
    }

    #[test]
    fn rejects_a_missing_internal_fixed_sequence() {
        let extractor = InsertionExtractor::new(
            vec!["AAAA".into(), "CCCC".into(), "GGGG".into()],
            vec![2, 2],
            3,
            0.0,
        )
        .unwrap();
        assert_eq!(extractor.extract_ranges(b"AAAATTAAGGGG"), None);
    }

    #[test]
    fn permits_gaps_in_fixed_sequence_alignment() {
        let extractor =
            InsertionExtractor::new(vec!["AAAA".into(), "CCCC".into()], vec![2], 3, 0.2).unwrap();
        assert_eq!(
            extractor.extract_ranges(b"AAAATTCCACC"),
            Some(vec![(0, 4..6)])
        );
    }

    #[test]
    fn applies_global_edit_rate_across_both_directions() {
        let extractor = InsertionExtractor::new(
            vec!["AAAA".into(), "CCCC".into(), "GGGG".into()],
            vec![2, 2],
            3,
            0.1,
        )
        .unwrap();
        assert_eq!(
            extractor.extract_ranges(b"AAAATTCCCCAAGGGT"),
            Some(vec![(0, 4..6), (1, 10..12)])
        );
        assert_eq!(extractor.extract_ranges(b"AAATTTCCCCAAGGGT"), None);
    }

    #[test]
    fn rejects_invalid_template_configuration() {
        assert!(InsertionExtractor::new(vec!["AAAA".into()], vec![], 3, 0.0).is_err());
        assert!(
            InsertionExtractor::new(vec!["AAAA".into(), "CCCC".into()], vec![], 3, 0.0).is_err()
        );
        assert!(
            InsertionExtractor::new(vec!["AAAA".into(), "CCCC".into()], vec![2], 5, 0.0).is_err()
        );
        assert!(
            InsertionExtractor::new(vec!["AAAA".into(), "CCCC".into()], vec![2], 3, -0.1).is_err()
        );
        assert!(
            InsertionExtractor::new(vec!["AAAA".into(), "CCCC".into()], vec![2], 3, 1.1).is_err()
        );
        assert!(
            InsertionExtractor::new(vec!["AAAA".into(), "CCCC".into()], vec![2], 3, f64::NAN)
                .is_err()
        );
    }
}
