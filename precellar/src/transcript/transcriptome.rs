use anyhow::{bail, ensure, Result};
use noodles::sam::alignment::record::cigar::op::Kind;
use noodles::sam::alignment::record::cigar::Op;
use noodles::sam::alignment::record_buf::Cigar;
use noodles::sam::alignment::record_buf::RecordBuf;
use noodles::sam::record::data::field::value::base_modifications::group::Strand;
use std::cmp;

/// 0-based, half-open
#[derive(Debug, Clone)]
pub struct Transcript {
    pub id: String,
    pub chrom: String,
    pub start: u64,
    pub end: u64, // exclusive
    pub strand: Strand,
    pub gene: Gene,
    exons: Exons,
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
        self.end - self.start
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
    pub fn is_exonic(&self, transcript: &Transcript, min_overlap_frac: f64) -> bool {
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

    /// Align to a transcript. Returns the aligned cigar and the number of aligned bases.
    pub fn align_junctions(
        &self,
        transcript: &Transcript,
        tolerance: u64,
        intergenic_trim_bases: u64,
        intronic_trim_bases: u64,
    ) -> Option<(Cigar, u64)> {
        let (ex_start, ex_end) = find_exons(
            &transcript.exons(),
            self.start(),
            self.end(),
            intergenic_trim_bases,
            intronic_trim_bases,
        )?;
        self._align_junctions_helper(&transcript.exons()[ex_start..=ex_end], tolerance)
    }

    /// Align the read to the exons. Returns the aligned cigar and the number of aligned bases.
    fn _align_junctions_helper(&self, exons: &[Exon], tolerance: u64) -> Option<(Cigar, u64)> {
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