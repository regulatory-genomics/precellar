use anyhow::Result;
use noodles::sam;
use noodles::sam::alignment::{
    record::cigar::{op::Kind, Op},
    record_buf::{Cigar, RecordBuf},
};
use noodles::sam::record::data::field::value::base_modifications::group::Strand;
use std::cmp;

use crate::fragment::Fragment;
use crate::transcriptome::{Exon, Transcript};

#[derive(Debug, Copy, Clone)]
pub enum ChemistryStrandness {
    Forward,
    Reverse,
    Unstranded,
}

#[derive(Debug, Clone)]
pub struct JunctionAlignOptions {
    /// Minimum overlap fraction required for a region to be considered exonic.
    pub region_min_overlap: f64,
    /// The maximum number of bases that can be misaligned at each junction.
    pub junction_trim_bases: u64,
    /// The number of bases to allow overhangs into intergenic regions.
    pub intergenic_trim_bases: u64,
    /// The number of bases to allow overhangs into intronic regions.
    pub intronic_trim_bases: u64,
    /// Strandness of the chemistry.
    pub chemistry_strandedness: Option<ChemistryStrandness>,
}

impl Default for JunctionAlignOptions {
    fn default() -> Self {
        Self {
            region_min_overlap: 0.5,
            junction_trim_bases: 0,
            intergenic_trim_bases: 0,
            intronic_trim_bases: 0,
            chemistry_strandedness: None,
        }
    }
}

#[derive(Eq, PartialEq, Debug, Clone)]
pub struct TranscriptAlignment {
    pub gene_id: String,
    pub transcript_id: String,
    pub strand: Strand,
    pub exon_align: Option<Cigar>,
}

impl TranscriptAlignment {
    pub fn is_exonic(&self) -> bool {
        self.exon_align.is_some()
    }

    pub fn is_intronic(&self) -> bool {
        !self.is_exonic()
    }
}

/// Segment represents a contiguous block of cigar operations not containing
/// any "Skip" sections. The SpliceSegment is 0-based, half-open with respect to the reference.
#[derive(Debug, Clone)]
struct Segment {
    start: u64,
    end: u64,
    cigar: Cigar,
}

impl bincode::Encode for Segment {
    fn encode<E: bincode::enc::Encoder>(
        &self,
        encoder: &mut E,
    ) -> core::result::Result<(), bincode::error::EncodeError> {
        bincode::Encode::encode(&self.start, encoder)?;
        bincode::Encode::encode(&self.end, encoder)?;
        bincode::Encode::encode(&encode_cigar(&self.cigar).unwrap(), encoder)?;
        Ok(())
    }
}

impl<Context> bincode::Decode<Context> for Segment {
    fn decode<D: bincode::de::Decoder<Context = Context>>(
        decoder: &mut D,
    ) -> core::result::Result<Self, bincode::error::DecodeError> {
        let start = bincode::Decode::decode(decoder)?;
        let end = bincode::Decode::decode(decoder)?;
        let cigar_bytes: Vec<u8> = bincode::Decode::decode(decoder)?;
        Ok(Self {
            start,
            end,
            cigar: decode_cigar(&cigar_bytes).unwrap(),
        })
    }
}

impl<'de, Context> bincode::BorrowDecode<'de, Context> for Segment {
    fn borrow_decode<D: bincode::de::BorrowDecoder<'de, Context = Context>>(
        decoder: &mut D,
    ) -> core::result::Result<Self, bincode::error::DecodeError> {
        let start = bincode::BorrowDecode::borrow_decode(decoder)?;
        let end = bincode::BorrowDecode::borrow_decode(decoder)?;
        let cigar_bytes: Vec<u8> = bincode::BorrowDecode::borrow_decode(decoder)?;
        Ok(Self {
            start,
            end,
            cigar: decode_cigar(&cigar_bytes).unwrap(),
        })
    }
}

fn encode_cigar(cigar: &Cigar) -> Result<Vec<u8>> {
    let mut buf = Vec::new();
    sam::io::writer::record::write_cigar(&mut buf, &cigar)?;
    Ok(buf)
}

fn decode_cigar(cigar: &[u8]) -> Result<Cigar> {
    if cigar[0] == b'*' {
        Ok(Vec::new().into())
    } else {
        let cigar = noodles::sam::record::Cigar::new(cigar);
        Ok(Cigar::try_from(cigar)?)
    }
}

/// SplicedRecord represents segments of a larger whole that have been cut and joined together
/// to form a new, continuous sequence.
/// It consists of the left and right clipping operations, and a list of components.
/// A SplicedRecord can be constructed from a BAM record by splitting the CIGAR string
/// at each "Skip" operation.
#[derive(Debug, Clone)]
pub struct SplicedRecord {
    pub chrom: String,
    left_clip: Cigar,
    right_clip: Cigar,
    segments: Vec<Segment>,
    pub strand: bed_utils::bed::Strand, // Strand information of the read.
}

impl bincode::Encode for SplicedRecord {
    fn encode<E: bincode::enc::Encoder>(
        &self,
        encoder: &mut E,
    ) -> core::result::Result<(), bincode::error::EncodeError> {
        bincode::Encode::encode(&self.chrom, encoder)?;
        bincode::Encode::encode(&encode_cigar(&self.left_clip).unwrap(), encoder)?;
        bincode::Encode::encode(&encode_cigar(&self.right_clip).unwrap(), encoder)?;
        bincode::Encode::encode(&self.segments, encoder)?;
        bincode::Encode::encode(&self.strand, encoder)?;
        Ok(())
    }
}

impl<Context> bincode::Decode<Context> for SplicedRecord {
    fn decode<D: bincode::de::Decoder<Context = Context>>(
        decoder: &mut D,
    ) -> core::result::Result<Self, bincode::error::DecodeError> {
        let chrom = bincode::Decode::decode(decoder)?;
        let left_cigar_bytes: Vec<u8> = bincode::Decode::decode(decoder)?;
        let right_cigar_bytes: Vec<u8> = bincode::Decode::decode(decoder)?;
        let segments = bincode::Decode::decode(decoder)?;
        let strand = bincode::Decode::decode(decoder)?;
        Ok(Self {
            chrom,
            left_clip: decode_cigar(&left_cigar_bytes).unwrap(),
            right_clip: decode_cigar(&right_cigar_bytes).unwrap(),
            segments,
            strand,
        })
    }
}

impl<'de, Context> bincode::BorrowDecode<'de, Context> for SplicedRecord {
    fn borrow_decode<D: bincode::de::BorrowDecoder<'de, Context = Context>>(
        decoder: &mut D,
    ) -> core::result::Result<Self, bincode::error::DecodeError> {
        let chrom = bincode::BorrowDecode::borrow_decode(decoder)?;
        let left_cigar_bytes: Vec<u8> = bincode::BorrowDecode::borrow_decode(decoder)?;
        let right_cigar_bytes: Vec<u8> = bincode::BorrowDecode::borrow_decode(decoder)?;
        let segments = bincode::BorrowDecode::borrow_decode(decoder)?;
        let strand = bincode::Decode::decode(decoder)?;
        Ok(Self {
            chrom,
            left_clip: decode_cigar(&left_cigar_bytes).unwrap(),
            right_clip: decode_cigar(&right_cigar_bytes).unwrap(),
            segments,
            strand,
        })
    }
}

impl SplicedRecord {
    pub fn new(read: &RecordBuf, header: &sam::Header) -> Result<Option<Self>> {
        if read.flags().is_unmapped() {
            return Ok(None);
        }

        let chrom = read.reference_sequence(header).unwrap()?.0;
        let chrom = std::str::from_utf8(chrom)?.to_string();
        let cigar = read.cigar();
        let alignment_start = read.alignment_start().unwrap().get();

        let mut left_clip: Vec<Op> = Vec::new();
        let mut right_clip: Vec<Op> = Vec::new();
        let mut splice_segments: Vec<Segment> = Vec::new();
        let mut seen_nonclips = false; // whether we've seen non-clip bases yet
        let mut curr_segment = Segment {
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
                    curr_segment = Segment {
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

        let flag = read.flags();
        let strand = if flag.is_reverse_complemented() {
            if flag.is_segmented() && flag.is_last_segment() {
                bed_utils::bed::Strand::Forward
            } else {
                bed_utils::bed::Strand::Reverse
            }
        } else {
            if flag.is_segmented() && flag.is_last_segment() {
                bed_utils::bed::Strand::Reverse
            } else {
                bed_utils::bed::Strand::Forward
            }
        };
        Ok(Some(Self {
            chrom,
            left_clip: left_clip.into(),
            right_clip: right_clip.into(),
            segments: splice_segments,
            strand,
        }))
    }

    /// The leftmost position of the alignment.
    pub fn start(&self) -> u64 {
        self.segments.first().map_or(0, |segment| segment.start)
    }

    /// The rightmost position of the alignment (exclusive).
    pub fn end(&self) -> u64 {
        self.segments.last().map_or(0, |segment| segment.end)
    }

    /// Determine if the alignment belongs to exonic regions of a transcript.
    /// An alignment is considered exonic if each of its components (segments)
    /// overlaps with exons with at least `min_overlap_frac`.
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

    /// Aligns a read to a transcript and determines the region type (exonic, intronic, or intergenic).
    pub fn align_transcript(
        &self,
        transcript: &Transcript,
        opts: &JunctionAlignOptions,
    ) -> Option<TranscriptAlignment> {
        // figure out coordinates
        let tx_start = transcript.start;
        let tx_end = transcript.end;
        let genomic_start = self.start();
        let genomic_end = self.end();

        let is_exonic = self.is_exonic(transcript, opts.region_min_overlap);
        if is_exonic || get_overlap(genomic_start, genomic_end, tx_start, tx_end) >= 1.0 {
            let is_antisense = match opts.chemistry_strandedness {
                Some(ChemistryStrandness::Forward) => transcript.strand != self.strand,
                Some(ChemistryStrandness::Reverse) => transcript.strand == self.strand,
                _ => false, // if unstranded, always sense
            };
            let tx_strand = if is_antisense {
                Strand::Reverse
            } else {
                Strand::Forward
            };

            let mut alignment = TranscriptAlignment {
                gene_id: transcript.gene_id.clone(),
                transcript_id: transcript.id.clone(),
                strand: tx_strand,
                exon_align: None,
            };

            if is_exonic {
                // compute offsets
                let mut tx_offset = transcript.get_offset(genomic_start).unwrap().max(0) as u64;
                let tx_length = transcript.exon_len();

                // align the read to the exons
                if let Some((mut tx_cigar, tx_aligned_bases)) =
                    self.align_junctions(transcript, opts)
                {
                    // flip reverse strand
                    if transcript.strand == bed_utils::bed::Strand::Reverse {
                        tx_offset = tx_length - (tx_offset + tx_aligned_bases);
                        tx_cigar.as_mut().reverse();
                    };
                    alignment.exon_align = Some(tx_cigar);
                }
            }
            Some(alignment)
        } else {
            None
        }
    }

    pub fn to_fragments(&self) -> impl Iterator<Item = Fragment> + '_ {
        self.segments.iter().map(move |segment| Fragment {
            chrom: self.chrom.clone(),
            start: segment.start, 
            end: segment.end,
            barcode: None,
            count: 1,
            strand: Some(self.strand),
            extended: None,
        })
    }

    /// Determine the orientation of the read with respect to a transcript.
    pub fn orientation(
        &self,
        transcript: &Transcript,
        min_overlap_frac: f64,
    ) -> Option<Strand> {
        // figure out coordinates
        let tx_start = transcript.start;
        let tx_end = transcript.end;
        let genomic_start = self.start();
        let genomic_end = self.end();

        let is_exonic = self.is_exonic(transcript, min_overlap_frac);
        if is_exonic || get_overlap(genomic_start, genomic_end, tx_start, tx_end) >= 1.0 {
            if transcript.strand == self.strand {
                return Some(Strand::Forward);  // sense
            } else {
                return Some(Strand::Reverse);  // antisense
            }
        } else {
            None
        }
    }

    /// Align the read to exons of a transcript. Returns the aligned cigar and the number of aligned bases.
    ///
    /// Returns None if the alignment cannot be aligned to the transcript.
    fn align_junctions(
        &self,
        transcript: &Transcript,
        opts: &JunctionAlignOptions,
    ) -> Option<(Cigar, u64)> {
        let exons = find_exons(
            &transcript.exons(),
            self.start(),
            self.end(),
            opts.intergenic_trim_bases,
            opts.intronic_trim_bases,
        )?;

        // The number of segments should match the number of exons
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
            } else if start_diff.unsigned_abs() > opts.junction_trim_bases {
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
            } else if end_diff.unsigned_abs() > opts.junction_trim_bases {
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
) -> Option<&[Exon]> {
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

    Some(&exon_info[ex_start..=ex_end])
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