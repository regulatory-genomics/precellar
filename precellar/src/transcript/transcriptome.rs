use anyhow::{bail, ensure, Result};
use noodles::sam::alignment::record::cigar::op::Kind;
use noodles::sam::alignment::record::cigar::Op;
use noodles::sam::alignment::record_buf::{Cigar, RecordBuf};
use noodles::sam::record::data::field::value::base_modifications::group::Strand;
use std::cmp;
use log::debug;


/// 0-based, half-open
#[derive(Debug, Clone)]
pub struct Transcript {
    pub id: String,
    pub chrom: String,
    pub start: u64,
    pub end: u64, // exclusive
    pub strand: Strand,
    pub gene: Gene,
    pub exons: Exons, // make public for test purpose
    pub introns: Introns, // make public for test purpose
}

impl TryFrom<star_aligner::transcript::Transcript> for Transcript {
    type Error = anyhow::Error;

    fn try_from(transcript: star_aligner::transcript::Transcript) -> Result<Self> {
        let start = transcript.start;
        let end = transcript.end;
        let strand = match transcript.strand {
            star_aligner::transcript::Strand::Forward => Strand::Forward,
            star_aligner::transcript::Strand::Reverse => Strand::Reverse,
            _ => bail!("Strand must be Forward or Reverse"),
        };
        let exons = Exons::new(transcript.exons.iter().map(|exon| {
            assert!(exon.start < exon.end, "Exon start must be less than exon end");
            assert!(exon.start >= start, "Exon start must be greater than transcript start");
            assert!(exon.end <= end, "Exon end must be less than transcript end");
            (exon.start, exon.end)
        }))?;
        Ok(Self {
            id: transcript.id,
            chrom: transcript.chrom,
            start,
            end,
            strand,
            gene: Gene {
                id: transcript.gene_id,
                name: transcript.gene_name,
            },
            exons,
            introns: Introns::new(std::iter::empty()).unwrap(),
        })
    }
}

impl Transcript {
    /// Transcript length is the sum of the lengths of the exons.
    pub fn len(&self) -> u64 {
        self.exons().iter().map(|x| x.len()).sum()
    }

    pub fn exons(&self) -> &[Exon] {
        self.exons.as_ref()
    }

    pub fn introns(&self) -> &[Intron] {
        self.introns.as_ref()
    }

    /// Convert a coordinate in the genome to a coordinate in the exons/transcript.
    /// The coordinate starts at the beginning of the transcript.
    pub fn get_offset(&self, coord: u64) -> Option<i64> {
        let mut cum_len = 0;
        for exon in &self.exons.0 {
            if coord < exon.end {
                return Some((cum_len + coord) as i64 - exon.start as i64);
            }
            cum_len += exon.len();
        }
        None
    }

    pub fn make_intron_by_exons(&mut self) {
        if self.exons().len() < 2 {
            self.introns = Introns::new(std::iter::empty()).unwrap();
            return;
        }
        let introns = (0..self.exons().len() - 1)
            .map(|i| (self.exons()[i].end+1, self.exons()[i + 1].start-1));
        self.introns = Introns::new(introns).unwrap();
    }
    

        /// Validate a specific intron by index
    pub fn validate_intron(&mut self, index: usize) -> bool {
        if let Some(intron) = self.introns.0.get_mut(index) {
            intron.set_validated(true);
            true
        } else {
            false
        }
    }

    pub fn is_intron_validated(&self, index: usize) -> bool {
        self.introns.0.get(index)
            .map(|intron| intron.validated())
            .unwrap_or(false)
    }


    /// Get count of validated introns
    pub fn validated_intron_count(&self) -> usize {
        self.introns.0.iter()
            .filter(|intron| intron.validated())
            .count()
    }

    /// Get total intron count
    pub fn total_intron_count(&self) -> usize {
        self.introns.0.len()
    }
}
#[derive(Debug, Clone)]
pub struct Exons(Vec<Exon>);

impl Exons {
    pub fn new(iter: impl IntoIterator<Item = (u64, u64)>) -> Result<Self> {
        let mut prev_end = None;
        let exon = iter
            .into_iter()
            .map(|(start, end)| {
                ensure!(
                    prev_end.is_none() || start >= prev_end.unwrap(),
                    "Exons must be non-overlapping and in order"
                );
                ensure!(
                    end > start,
                    "End coordinate must be greater than start coordinate"
                );
                prev_end = Some(end);
                Ok(Exon { start, end })
            })
            .collect::<Result<Vec<_>>>()?;
        Ok(Self(exon))
    }
}

impl AsRef<[Exon]> for Exons {
    fn as_ref(&self) -> &[Exon] {
        &self.0
    }
}

/// Exon coordinates are 0-based, half-open.
#[derive(Eq, PartialEq, Debug, Clone, Ord, PartialOrd)]
pub struct Exon {
    start: u64,
    end: u64,
}

impl Exon {
    pub fn len(&self) -> u64 {
        self.end - self.start // Why not self.end - self.start + 1???
    }
    pub fn start(&self) -> u64 {
        self.start
    }
    pub fn end(&self) -> u64 {
        self.end
    }
}


/// I try to make intron struct to store the validated intron information
/// I found Intron maybe can be removed later for simplicity
#[derive(Eq, PartialEq, Debug, Clone, Ord, PartialOrd)]
pub struct Intron {
    start: u64,
    end: u64,
    validated: bool,
}

impl Intron {
    pub fn len(&self) -> u64 {
        self.end - self.start
    }
    pub fn start(&self) -> u64 {
        self.start
    }

    pub fn end(&self) -> u64 {
        self.end
    }

    pub fn validated(&self) -> bool {
        self.validated
    }

    pub fn set_validated(&mut self, validated: bool) {
        self.validated = validated;
    }
}

#[derive(Debug, Clone)]
pub struct Introns(Vec<Intron>);

impl Introns {
    pub fn new(iter: impl IntoIterator<Item = (u64, u64)>) -> Result<Self> {
        let introns = iter.into_iter().map(|(start, end)| Intron { start, end, validated: false }).collect();
        Ok(Self(introns))
    }
}
impl AsRef<[Intron]> for Introns {
    fn as_ref(&self) -> &[Intron] {
        &self.0
    }
}


#[derive(Hash, Eq, PartialEq, Debug, Clone, Ord, PartialOrd)]
pub struct Gene {
    pub id: String,
    pub name: String,
}

/// SpliceSegment represents a contiguous block of cigar operations not containing
/// any "Skip" sections. The SpliceSegment is 0-based, half-open with respect to the reference.
#[derive(Debug)]
struct SpliceSegment {
    start: u64,
    end: u64,
    cigar: Cigar,
}

/// SpliceSegments is used to represent the alignment of a read to a transcript.
/// It consists of the left and right clipping operations, and a list of SpliceSegments.
pub struct SpliceSegments {
    left_clip: Cigar,
    right_clip: Cigar,
    segments: Vec<SpliceSegment>,
}

impl SpliceSegments {
    /// The leftmost position of all segments.
    pub fn start(&self) -> u64 {
        self.segments.first().map_or(0, |segment| segment.start)
    }

    /// The rightmost position of all segments.
    pub fn end(&self) -> u64 {
        self.segments.last().map_or(0, |segment| segment.end)
    }

    /// Determine if the read aligns to exonic regions of a transcript. A read is considered exonic
    /// if it aligns to all exons with at least `min_overlap_frac` overlap.
    pub fn is_exonic(&self, transcript: &Transcript, min_overlap_frac: f64,is_any : bool) -> bool {
        if is_any {
            self.segments.iter().any(|segment| {
                // find first exon that ends to the right of the segment start
                let idx = transcript
                    .exons()
                    .binary_search_by_key(&segment.start, |ex| ex.end - 1)
                    .unwrap_or_else(std::convert::identity);
                transcript.exons().get(idx).map_or(false, |exon| {
                    get_overlap(segment.start, segment.end, exon.start, exon.end) >= min_overlap_frac
                })
            })
        } else {
            self.segments.iter().all(|segment| {
                // find first exon that ends to the right of the segment start
                let idx = transcript
                    .exons()
                    .binary_search_by_key(&segment.start, |ex| ex.end - 1)
                    .unwrap_or_else(std::convert::identity);
                transcript.exons().get(idx).map_or(false, |exon| {
                    get_overlap(segment.start, segment.end, exon.start, exon.end) >= min_overlap_frac
                })
            })
        }
    }

    /// Align to a transcript. Returns the aligned cigar and the number of aligned bases.
    pub fn align_junctions(
        &self,
        transcript: &Transcript,
        tolerance: u64,
        intergenic_trim_bases: u64,
        intronic_trim_bases: u64,
    ) -> Option<(Cigar, u64)> {  // (cigar, aligned_bases)
        let (ex_start, ex_end) = find_exons(
            &transcript.exons(),
            self.start(),
            self.end(),
            intergenic_trim_bases,
            intronic_trim_bases,
            false, // I suppose we don't need validate intergenic & intronic trim in intron validation
        )?;
        let (cigar, aligned_bases) = self._align_junctions_helper(
            &transcript.exons()[ex_start..=ex_end],
            tolerance,
        )?;
        return Some((cigar, aligned_bases));
    }

    pub fn annotate_splice(&self,
        transcript: &Transcript
    ) -> Option<(bool, Vec<usize>)> {

        let mut has_spanning = false;
        let mut intron_mapped = Vec::new();
        
        let introns = transcript.introns();

        for (seg_idx, current_segment) in self.segments.iter().enumerate() {
            // Handle the case where find_exons returns None gracefully
            let (ex_start, ex_end) = match find_exons(
                &transcript.exons(),
                current_segment.start,
                current_segment.end,
                0,
                0,
                true,
            ) {
                Some(result) => result,
                None => {
                    // Skip this segment if find_exons fails (e.g., due to overhang or no overlap)
                    continue;
                }
            };
            
            if introns.is_empty() {
                continue;
            }
            
            let start_idx = ex_start.saturating_sub(1);
            let end_idx = std::cmp::min(ex_end + 1, introns.len() - 1);
            
            if start_idx <= end_idx {
                let start_idx = ex_start.saturating_sub(1);
                let end_idx = std::cmp::min(ex_end + 1, introns.len().saturating_sub(1));
                let intron_selected = &introns[start_idx..=end_idx];
                
                for (intron_idx, intron) in intron_selected.iter().enumerate() {
                    let intron_start = intron.start;
                    let intron_end = intron.end;
                    let global_intron_idx = intron_idx + start_idx as usize;
                    
                    // Check if segment spans intron boundaries
                    if (intron_start > current_segment.start && intron_start < current_segment.end) || 
                       (intron_end > current_segment.start && intron_end < current_segment.end) {
                        has_spanning = true;
                    }
                    
                    // Check for overlap using efficient comparison
                    if intron_start < current_segment.end && intron_end > current_segment.start {
                        intron_mapped.push(global_intron_idx);
                    }
                }
            }
        }
        
        Some((has_spanning, intron_mapped))
    }

    /// Align the read to the exons. Returns the aligned cigar and the number of aligned bases.
    fn _align_junctions_helper(&self,
         exons: &[Exon], 
         tolerance: u64,
        ) -> Option<(Cigar, u64)> {
        // check if the number of segments matches the number of exons
        if self.segments.len() != exons.len() {
            return None;
        }

        let mut full_cigar = self.left_clip.clone();
        let mut aligned_bases = 0;
        for i in 0..self.segments.len() {
            let curr_segment = &self.segments[i];
            let curr_exon = &exons[i];
            aligned_bases += curr_exon.len();
            let mut tmp_cigar = curr_segment.cigar.clone();

            // align the start
            let start_diff = curr_exon.start as i64 - curr_segment.start as i64;
            //mark the overhang as validated intron

            if i == 0 {
                // first segment
                if start_diff > 0 {
                    // overhang -> softclip
                    tmp_cigar = mask_read_bases(
                        &mut tmp_cigar,
                        Op::new(Kind::SoftClip, start_diff as usize),
                        false,
                    );


                } else if start_diff < 0 {
                    // underhang -> decrement aligned bases
                    aligned_bases -= start_diff.unsigned_abs();
                }
            } else if start_diff.unsigned_abs() > tolerance {
                return None; // can't align properly
            } else if start_diff > 0 {
                // overhang -> mark as insertion
                tmp_cigar = mask_read_bases(
                    &mut tmp_cigar,
                    Op::new(Kind::Insertion, start_diff as usize),
                    false,
                );
            } else if start_diff < 0 {
                // underhang -> mark as deletion
                tmp_cigar = mark_deleted_ref_bases(
                    &mut tmp_cigar,
                    start_diff.unsigned_abs().try_into().unwrap(),
                    false,
                );
            }

            // align the end
            let end_diff = curr_segment.end as i64 - curr_exon.end as i64 - 1;
            if i == self.segments.len() - 1 {
                // last segment
                if end_diff > 0 {
                    // overhang -> softclip
                    tmp_cigar = mask_read_bases(
                        &mut tmp_cigar,
                        Op::new(Kind::SoftClip, end_diff as usize),
                        true,
                    );
                } else if end_diff < 0 {
                    // underhang -> decrement aligned bases
                    aligned_bases -= end_diff.unsigned_abs();
                }
            } else if end_diff.unsigned_abs() > tolerance {
                return None; // can't align properly
            } else if end_diff > 0 {
                // overhang -> mark as insertion
                tmp_cigar = mask_read_bases(
                    &mut tmp_cigar,
                    Op::new(Kind::Insertion, end_diff as usize),
                    true,
                );
            } else if end_diff < 0 {
                // underhang -> mark as deletion
                tmp_cigar = mark_deleted_ref_bases(
                    &mut tmp_cigar,
                    end_diff.unsigned_abs().try_into().unwrap(),
                    true,
                );
            }

            // extend
            full_cigar.extend(Vec::from(tmp_cigar).into_iter());
        }
        full_cigar.extend(self.right_clip.as_ref().iter().copied());

        Some((full_cigar, aligned_bases))
    }
}

impl From<&RecordBuf> for SpliceSegments {
    fn from(read: &RecordBuf) -> Self {
        let cigar = read.cigar();
        let alignment_start = read.alignment_start().unwrap().get();

        let mut left_clip: Vec<Op> = Vec::new();
        let mut right_clip: Vec<Op> = Vec::new();
        let mut splice_segments: Vec<SpliceSegment> = Vec::new();
        let mut seen_nonclips = false; // whether we've seen non-clip bases yet
        let mut curr_segment = SpliceSegment {
            start: alignment_start as u64,
            end: alignment_start as u64,
            cigar: Vec::new().into(),
        };

        for c in cigar.as_ref() {
            match c.kind() {
                Kind::HardClip | Kind::SoftClip => {
                    if seen_nonclips {
                        right_clip.push(*c);
                    } else {
                        left_clip.push(*c);
                    }
                }
                Kind::Skip => {
                    seen_nonclips = true;
                    let next_start = curr_segment.end + c.len() as u64;
                    splice_segments.push(curr_segment);
                    curr_segment = SpliceSegment {
                        start: next_start,
                        end: next_start,
                        cigar: Vec::new().into(),
                    };
                }
                Kind::Insertion => {
                    seen_nonclips = true;
                    curr_segment.cigar.as_mut().push(*c);
                }
                Kind::Match | Kind::Deletion | Kind::SequenceMatch | Kind::SequenceMismatch => {
                    seen_nonclips = true;
                    curr_segment.end += c.len() as u64;
                    curr_segment.cigar.as_mut().push(*c);
                }
                Kind::Pad => unreachable!(),
            }
        }
        splice_segments.push(curr_segment);

        Self {
            left_clip: left_clip.into(),
            right_clip: right_clip.into(),
            segments: splice_segments,
        }
    }
}

/// Fraction of read interval covered by ref interval
fn get_overlap(read_start: u64, read_end: u64, ref_start: u64, ref_end: u64) -> f64 {
    let mut overlap_bases =
        cmp::min(ref_end, read_end) as f64 - cmp::max(ref_start, read_start) as f64;
    if overlap_bases < 0.0 {
        overlap_bases = 0.0;
    }
    overlap_bases / ((read_end - read_start) as f64)
}

/// Find the exons that the read aligns to. Returns the indices of the first and last exons.
fn find_exons(
    exon_info: &[Exon],
    read_start: u64,
    read_end: u64, // inclusive
    intergenic_trim_bases: u64,
    intronic_trim_bases: u64,
    skip_trim_validation: bool, // I suppose we don't need validate intergenic & intronic trim in intron validation
) -> Option<(usize, usize)> {
    // find first exon that ends to the right of the read start
    let ex_start = exon_info
        .binary_search_by_key(&read_start, |ex| ex.end - 1)
        .map_or_else(|i| i, |i| i);
    // find first exon that starts to the left of the read end
    let ex_end = exon_info
        .binary_search_by_key(&read_end, |ex| ex.start)
        .map_or_else(|i| if i > 0 { Some(i - 1) } else { None }, |i| Some(i))?;
    if ex_start >= exon_info.len() {
        return None;
    }
    if skip_trim_validation {
        return Some((ex_start, ex_end));
    }
    let starting_exon = &exon_info[ex_start];
    let ending_exon = &exon_info[ex_end];

    if read_start < starting_exon.start {
        // read overhangs exon on the left
        let overhang = starting_exon.start - read_start;
        let trim_bases = if ex_start == 0 {
            intergenic_trim_bases
        } else {
            intronic_trim_bases
        };
        if overhang > trim_bases {
            // too much overhang
            return None;
        };
    }

    if read_end > ending_exon.end {
        // read overhangs exon on the right
        let overhang = read_end - ending_exon.end;
        let trim_bases = if ex_end >= exon_info.len() {
            intergenic_trim_bases
        } else {
            intronic_trim_bases
        };
        if overhang > trim_bases {
            // too much overhang
            return None;
        };
    }

    Some((ex_start, ex_end))
}

fn mask_read_bases(cigar: &mut Cigar, mask: Op, reverse: bool) -> Cigar {
    // NOTE: this assumes that refskips have been removed
    let mut new_cigar = Vec::new();
    let mask_len = mask.len();
    let mut consumed_bases = 0;
    new_cigar.push(mask);
    if reverse {
        cigar.as_mut().reverse();
    }
    for c in cigar.as_ref() {
        if consumed_bases < mask_len {
            // this op should be masked
            let read_bases = match c.kind() {
                Kind::Deletion => 0, // deletions don't consume read bases
                _ => c.len(),
            };
            if consumed_bases + read_bases >= mask_len {
                let truncated = Op::new(c.kind(), read_bases + consumed_bases - mask_len);
                new_cigar.push(truncated);
            };
            consumed_bases += read_bases;
        } else {
            // just copy the op
            new_cigar.push(*c);
        };
    }
    if reverse {
        new_cigar.reverse();
    }
    new_cigar.into()
}

fn mark_deleted_ref_bases(cigar: &mut Cigar, del_len: usize, reverse: bool) -> Cigar {
    let del = Op::new(Kind::Deletion, del_len);
    if reverse {
        let mut new_cigar: Cigar = vec![del].into();
        new_cigar.extend(cigar.as_ref().iter().copied());
        new_cigar
    } else {
        let mut new_cigar = cigar.clone();
        new_cigar.as_mut().push(del);
        new_cigar
    }
}



#[cfg(test)]
mod tests {
    use super::*;
    use crate::transcript::{Transcript, Gene, Exons, Introns};
    use noodles::sam::record::data::field::value::base_modifications::group::Strand;

    fn create_test_gene() -> Gene {
        Gene {
            id: "GENE001".to_string(),
            name: "TestGene".to_string(),
        }
    }

    fn create_test_transcript_with_gaps() -> Transcript {
        // Create a transcript with 3 exons that have gaps between them (will create introns)
        let exons = Exons::new(vec![
            (100, 200),  // Exon 1: 100-200
            (300, 400),  // Exon 2: 300-400 (gap: 200-300)
            (500, 600),  // Exon 3: 500-600 (gap: 400-500)
        ]).unwrap();

        let introns = Introns::new(std::iter::empty()).unwrap(); // Start with empty introns

        Transcript {
            id: "TRANSCRIPT001".to_string(),
            chrom: "chr1".to_string(),
            start: 100,
            end: 600,
            strand: Strand::Forward,
            gene: create_test_gene(),
            exons,
            introns,
        }
    }


    #[test]
    fn test_introns_new_method() {
        // Test creating introns from coordinate pairs
        let intron_coords = vec![
            (200, 300),
            (400, 500),
        ];

        let introns = Introns::new(intron_coords).unwrap();
        let intron_slice = introns.as_ref();

        assert_eq!(intron_slice.len(), 2);

        // Check first intron
        assert_eq!(intron_slice[0].start(), 200);
        assert_eq!(intron_slice[0].end(), 300);
        assert_eq!(intron_slice[0].len(), 100);
        assert!(!intron_slice[0].validated());

        // Check second intron
        assert_eq!(intron_slice[1].start(), 400);
        assert_eq!(intron_slice[1].end(), 500);
        assert_eq!(intron_slice[1].len(), 100);
        assert!(!intron_slice[1].validated());
    }

    #[test]
    fn test_transcript_make_intron_by_exons_with_gaps() {
        let mut transcript = create_test_transcript_with_gaps();

        // Initially, transcript should have no introns
        assert_eq!(transcript.introns().len(), 0);

        // Generate introns from exons
        transcript.make_intron_by_exons();

        // Should now have 2 introns
        let introns = transcript.introns();
        assert_eq!(introns.len(), 2);

        // Check first intron (between exon 1 and exon 2)
        assert_eq!(introns[0].start(), 201);  // End of first exon
        assert_eq!(introns[0].end(), 299);    // Start of second exon
        assert_eq!(introns[0].len(), 98);
        assert!(!introns[0].validated());

        // Check second intron (between exon 2 and exon 3)
        assert_eq!(introns[1].start(), 401);  // End of second exon
        assert_eq!(introns[1].end(), 499);    // Start of third exon
        assert_eq!(introns[1].len(), 98);
        assert!(!introns[1].validated());
    }



    #[test]
    fn test_transcript_introns_accessor() {
        let mut transcript = create_test_transcript_with_gaps();
        transcript.make_intron_by_exons();

        // Test that we can access introns through the accessor method
        let introns = transcript.introns();
        assert_eq!(introns.len(), 2);

        // Verify we get the same data through the accessor
        assert_eq!(introns[0].start(), 201);
        assert_eq!(introns[0].end(), 299);
        assert_eq!(introns[1].start(), 401);
        assert_eq!(introns[1].end(), 499);
    }

    #[test]
    fn test_single_exon_transcript() {
        // Test transcript with only one exon (should have no introns)
        let exons = Exons::new(vec![(100, 200)]).unwrap();
        let introns = Introns::new(std::iter::empty()).unwrap();

        let mut transcript = Transcript {
            id: "TRANSCRIPT_SINGLE".to_string(),
            chrom: "chr1".to_string(),
            start: 100,
            end: 200,
            strand: Strand::Forward,
            gene: create_test_gene(),
            exons,
            introns,
        };

        transcript.make_intron_by_exons();

        // Should have no introns
        assert_eq!(transcript.introns().len(), 0);
    }

    #[test]
    fn test_empty_introns() {
        let introns = Introns::new(vec![]).unwrap();
        let intron_slice = introns.as_ref();

        assert_eq!(intron_slice.len(), 0);
        assert!(intron_slice.is_empty());
    }

    #[test]
    fn test_multiple_transcripts_with_introns() {
        // Test that multiple transcripts can have their own introns
        let mut transcript1 = create_test_transcript_with_gaps();

        transcript1.make_intron_by_exons();

        // Transcript 1 should have 2 introns
        assert_eq!(transcript1.introns().len(), 2);


        // Verify they don't interfere with each other
        assert_eq!(transcript1.introns()[0].start(), 201);
    }


    #[test]
    fn test_intron_length_calculation() {
        let intron_coords = vec![
            (100, 200),   // Length: 100
            (500, 1000),  // Length: 500
            (2000, 2001), // Length: 1
        ];

        let introns = Introns::new(intron_coords).unwrap();
        let intron_slice = introns.as_ref();

        assert_eq!(intron_slice[0].len(), 100);
        assert_eq!(intron_slice[1].len(), 500);
        assert_eq!(intron_slice[2].len(), 1);
    }
}