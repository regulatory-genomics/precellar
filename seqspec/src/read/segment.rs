use std::ops::Range;

use editdistancek::edit_distance_bounded;
use noodles::fastq::{self, record::Definition};

use crate::{utils::hamming_distance, Region, RegionType, SequenceType};

/// Result of anchor search when indels are allowed: position, actual match length, and edit distance.
#[derive(Debug, Clone, Copy)]
pub struct AnchorMatch {
    pub start_pos: usize,
    pub match_len: usize,
    pub edit_distance: usize,
}

#[derive(Debug, Clone)]
pub enum SegmentType<'a> {
    Single {
        region_id: &'a str,
        region_type: RegionType,
    },
    Joined(Vec<RegionType>),
}

impl core::fmt::Display for SegmentType<'_> {
    fn fmt(&self, f: &mut core::fmt::Formatter) -> core::fmt::Result {
        match self {
            SegmentType::Single { region_type, .. } => write!(f, "{}", region_type),
            SegmentType::Joined(region_types) => {
                let joined = region_types
                    .iter()
                    .map(|r| r.to_string())
                    .collect::<Vec<String>>()
                    .join("");
                write!(f, "{}", joined)
            }
        }
    }
}

/// A segment of a read, consisting of a region type, sequence, and quality scores.
#[derive(Debug, Clone)]
pub struct Segment<'a> {
    pub region_type: SegmentType<'a>,
    pub seq: &'a [u8],
    pub qual: &'a [u8],
}

impl core::fmt::Display for Segment<'_> {
    fn fmt(&self, f: &mut core::fmt::Formatter) -> core::fmt::Result {
        write!(
            f,
            "[{}]{}",
            self.region_type,
            std::str::from_utf8(self.seq).unwrap()
        )
    }
}

impl<'a> Segment<'a> {
    pub fn new(region_type: SegmentType<'a>, seq: &'a [u8], qual: &'a [u8]) -> Self {
        Self {
            region_type,
            seq,
            qual,
        }
    }

    pub fn region_id(&self) -> &str {
        match &self.region_type {
            SegmentType::Single { region_id, .. } => region_id,
            SegmentType::Joined(_) => panic!("Joined segments have no region id"),
        }
    }

    pub fn into_fq(&self, defi: &Definition) -> fastq::Record {
        fastq::Record::new(defi.clone(), self.seq, self.qual)
    }

    pub fn is_barcode(&self) -> bool {
        match &self.region_type {
            SegmentType::Single { region_type, .. } => region_type.is_barcode(),
            _ => false,
        }
    }

    pub fn is_umi(&self) -> bool {
        match &self.region_type {
            SegmentType::Single { region_type, .. } => region_type.is_umi(),
            _ => false,
        }
    }

    pub fn contains_target(&self) -> bool {
        match &self.region_type {
            SegmentType::Single { region_type, .. } => region_type.is_target(),
            SegmentType::Joined(region_types) => region_types.iter().any(|r| r.is_target()),
        }
    }
}

#[derive(Debug, Clone)]
pub struct SegmentInfo {
    elems: Vec<SegmentInfoElem>,
    is_reverse: bool,
}

#[derive(Debug, Clone)]
pub struct SplitResult<'a> {
    pub segments: Vec<Segment<'a>>,
    pub mismatches: usize,
    pub is_reverse: bool,
}

/// Result of trying both orientations when parsing an unstranded read.
/// Used for two-phase strand selection (e.g. with barcode validation).
#[derive(Debug)]
pub enum SplitAutoRcBothResult<'a> {
    /// Only one orientation parsed successfully.
    One(SplitResult<'a>),
    /// Both forward and reverse parsed; caller can choose by mismatch count or barcode.
    Both {
        forward: SplitResult<'a>,
        reverse: SplitResult<'a>,
    },
}

#[derive(Debug, Copy, Clone)]
pub enum SplitError {
    AnchorNotFound,
    PatternMismatch(usize),
}

impl std::fmt::Display for SplitError {
    fn fmt(&self, f: &mut core::fmt::Formatter) -> core::fmt::Result {
        match self {
            SplitError::AnchorNotFound => write!(f, "Anchor not found"),
            SplitError::PatternMismatch(d) => write!(f, "Pattern mismatch ({} mismatches)", d),
        }
    }
}

impl std::error::Error for SplitError {}

fn make_rc_record(read: &fastq::Record) -> fastq::Record {
    let rc_seq = crate::utils::rev_compl(read.sequence());
    let mut rc_qual = read.quality_scores().to_vec();
    rc_qual.reverse();
    fastq::Record::new(
        fastq::record::Definition::new("", ""),
        rc_seq,
        rc_qual,
    )
}

/// Internal: mapping from RC segment coordinates to original read (for building SplitResult).
struct SegmentMapping<'a> {
    region_id: &'a str,
    region_type_variant: RegionType,
    is_joined: bool,
    joined_types: Vec<RegionType>,
    orig_start: usize,
    orig_end: usize,
}

impl SegmentInfo {
    pub fn new<T, S>(iter: T, is_reverse: bool) -> Self
    where
        T: IntoIterator<Item = S>,
        S: Into<SegmentInfoElem>,
    {
        let elems = iter.into_iter().map(|x| x.into()).collect::<Vec<_>>();
        Self { elems, is_reverse }
    }

    pub fn len(&self) -> usize {
        self.elems.len()
    }

    pub fn is_reverse(&self) -> bool {
        self.is_reverse
    }

    pub fn iter<'a>(&'a self) -> impl Iterator<Item = &'a SegmentInfoElem> {
        self.elems.iter()
    }

    pub fn into_iter(self) -> impl Iterator<Item = SegmentInfoElem> {
        self.elems.into_iter()
    }

    /// Truncate the segment information to the given length L, such that L is larger
    /// than the sum of the maximum lengths of the n-1 segments + minimum length of the last segment.
    pub fn truncate_max(self, len: usize) -> Self {
        let is_reverse = self.is_reverse;
        let mut sum = 0;
        let mut result = Vec::new();
        let mut segments = self.into_iter();
        while let Some(elem) = segments.next() {
            if sum + elem.max_len() > len {
                if sum + elem.min_len() <= len {
                    result.push(elem);
                }
                break;
            } else {
                sum += elem.max_len();
                result.push(elem);
            }
        }
        Self {
            elems: result,
            is_reverse: is_reverse,
        }
    }

    /// Split a read into segments according to the segment information.
    /// Uses default tolerances (1.0, 0.2) and no anchor indels (ratio 0).
    pub fn split<'a>(&'a self, read: &'a fastq::Record) -> Result<Vec<Segment<'a>>, SplitError> {
        self.split_with_tolerance(read, 1.0, 0.2, 0.0)
            .map(|res| res.segments)
    }

    /// Split a read into segments according to the segment information.
    ///
    /// # Arguments
    ///
    /// * `read` - The read to split
    /// * `linker_tolerance` - The fraction of mismatches allowed for the linker sequence
    /// * `anchor_tolerance` - The fraction of mismatches allowed for the anchor sequence
    /// * `anchor_max_indels_ratio` - Ratio of anchor length for max indels (e.g. 0.1 → cap 1 for 8 bp). 0 = Hamming only.
    ///
    /// # Returns
    /// `SplitResult` with segments, mismatch count, and orientation.
    pub fn split_with_tolerance<'a>(
        &'a self,
        read: &'a fastq::Record,
        linker_tolerance: f64,
        anchor_tolerance: f64,
        anchor_max_indels_ratio: f64,
    ) -> Result<SplitResult<'a>, SplitError> {
        fn consume_buf<'a>(
            buf: &mut Vec<&'a SegmentInfoElem>,
            mut seq: &'a [u8],
            mut qual: &'a [u8],
            result: &mut Vec<Segment<'a>>,
        ) {
            let n1 = result.len();
            while let Some(elem) = buf.pop() {
                if !elem.is_fixed_len() {
                    buf.push(elem);
                    break;
                }

                let l = seq.len();
                let i = l - elem.max_len();
                let s = &seq[i..l];
                let q = &qual[i..l];
                seq = &seq[..i];
                qual = &qual[..i];
                result.push(Segment::new(
                    SegmentType::Single {
                        region_id: elem.region_id.as_str(),
                        region_type: elem.region_type,
                    },
                    s,
                    q,
                ));
            }

            if !buf.is_empty() {
                let segment_type = if buf.len() == 1 {
                    SegmentType::Single {
                        region_id: buf[0].region_id.as_str(),
                        region_type: buf[0].region_type,
                    }
                } else {
                    SegmentType::Joined(buf.iter().map(|s| s.region_type).collect())
                };
                buf.clear();
                result.push(Segment::new(segment_type, seq, qual));
            }

            let n2 = result.len();
            if n2 > n1 {
                result[n1..n2].reverse();
            }
        }

        let mut result = Vec::new();
        let mut min_offset = 0;
        let mut max_offset = 0;
        let mut total_mismatches = 0usize;
        let mut buffer: Vec<&SegmentInfoElem> = Vec::new();
        let len = read.sequence().len();

        for segment in self.iter() {
            let min_len = segment.min_len();
            let max_len = segment.max_len();

            // Check if the segment is within the read bounds.
            // If it's a fixed segment (anchor), allow it to be shorter by up to cap_from_ratio (deletions).
            let anchor_cap = if segment.sequence_type.is_fixed() && anchor_max_indels_ratio > 0.0 {
                (min_len as f64 * anchor_max_indels_ratio).ceil() as usize
            } else {
                0
            };
            let allowed_min_len = if segment.sequence_type.is_fixed() {
                min_len.saturating_sub(anchor_cap)
            } else {
                min_len
            };

            if max_offset >= len || min_offset + allowed_min_len > len {
                break;
            }

            if min_len == max_len {
                if !buffer.is_empty() {
                    // if the sequcence contains fixed patterns, we can try to use them as anchors
                    // to find the location of the next segment
                    if segment.sequence_type.is_fixed() {
                        let search_start = min_offset;
                        let needle = if self.is_reverse {
                            segment.rc_sequence.as_slice()
                        } else {
                            segment.sequence.as_bytes()
                        };
                        // Cap from ratio (e.g. 0.1 * 8 = 0.8 → 1)
                        let cap_from_ratio = if anchor_max_indels_ratio > 0.0 {
                            (needle.len() as f64 * anchor_max_indels_ratio).ceil() as usize
                        } else {
                            0
                        };
                        // Allow the window to expand to catch insertions
                        let window_expansion = cap_from_ratio;
                        let search_end = len.min(max_offset + max_len + window_expansion);
                        let seq = read
                            .sequence()
                            .get(search_start..search_end)
                            .unwrap();

                        // Dual-gate: k_bound = min(ratio-based, cap from ratio); k_bound == 0 => Hamming only
                        let allowed_by_ratio =
                            (anchor_tolerance * needle.len() as f64).floor() as usize;
                        let k_bound = allowed_by_ratio.min(cap_from_ratio);

                        let anchor_match = if k_bound == 0 {
                            let (pos_opt, mis) = find_best_pattern_match(seq, needle);
                            pos_opt.map(|pos| AnchorMatch {
                                start_pos: pos,
                                match_len: needle.len(),
                                edit_distance: mis,
                            })
                        } else {
                            find_best_indel_pattern_match(seq, needle, k_bound)
                        };

                        if let Some(match_result) = anchor_match {
                            let max_allowed_by_ratio =
                                (anchor_tolerance * needle.len() as f64).floor() as usize;
                            if match_result.edit_distance <= max_allowed_by_ratio {
                                total_mismatches += match_result.edit_distance;
                                let offset_left =
                                    min_offset - buffer.iter().map(|s| s.min_len()).sum::<usize>();
                                let offset_right = min_offset + match_result.start_pos;

                                consume_buf(
                                    &mut buffer,
                                    read.sequence().get(offset_left..offset_right).unwrap(),
                                    read.quality_scores()
                                        .get(offset_left..offset_right)
                                        .unwrap(),
                                    &mut result,
                                );

                                min_offset = offset_right + match_result.match_len;
                                max_offset = min_offset;

                                let anchor_seq = read
                                    .sequence()
                                    .get(offset_right..offset_right + match_result.match_len)
                                    .unwrap();
                                let anchor_qual = read
                                    .quality_scores()
                                    .get(offset_right..offset_right + match_result.match_len)
                                    .unwrap();

                                result.push(Segment::new(
                                    SegmentType::Single {
                                        region_id: segment.region_id.as_str(),
                                        region_type: segment.region_type,
                                    },
                                    anchor_seq,
                                    anchor_qual,
                                ));
                                // Offsets already advanced by match_len; skip loop-tail increment
                                continue;
                            } else if max_offset + max_len > len {
                                // Read too short to contain anchor; break and flush buffer.
                                break;
                            } else {
                                return Err(SplitError::PatternMismatch(
                                    match_result.edit_distance,
                                ));
                            }
                        } else {
                            // No anchor match; if read is too short, break and flush buffer.
                            if max_offset + max_len > len {
                                break;
                            } else {
                                return Err(SplitError::AnchorNotFound);
                            }
                        }
                    } else {
                        buffer.push(segment);
                    }
                } else {
                    let seq = read
                        .sequence()
                        .get(max_offset..max_offset + max_len)
                        .unwrap();
                    let qual = read
                        .quality_scores()
                        .get(max_offset..max_offset + max_len)
                        .unwrap();

                    if linker_tolerance < 1.0 && segment.sequence_type.is_fixed() {
                        let d = if self.is_reverse {
                            hamming_distance(seq, segment.rc_sequence.as_slice()).unwrap()
                        } else {
                            hamming_distance(seq, segment.sequence.as_bytes()).unwrap()
                        };

                        if d as f64 > linker_tolerance * seq.len() as f64 {
                            return Err(SplitError::PatternMismatch(d as usize));
                        }
                        total_mismatches += d as usize;
                    }

                    result.push(Segment::new(
                        SegmentType::Single {
                            region_id: segment.region_id.as_str(),
                            region_type: segment.region_type,
                        },
                        seq,
                        qual,
                    ))
                }
            } else {
                buffer.push(segment);
            }

            min_offset += min_len;
            max_offset += max_len;
        }

        if !buffer.is_empty() {
            let offset_left = min_offset - buffer.iter().map(|s| s.min_len()).sum::<usize>();
            consume_buf(
                &mut buffer,
                read.sequence()
                    .get(offset_left..max_offset.min(len))
                    .unwrap(),
                read.quality_scores()
                    .get(offset_left..max_offset.min(len))
                    .unwrap(),
                &mut result,
            );
        }

        Ok(SplitResult {
            segments: result,
            mismatches: total_mismatches,
            is_reverse: self.is_reverse,
        })
    }

    /// Automatically try both orientations using sequence RC and coordinate mapping.
    ///
    /// Delegates to `split_auto_rc_both`. When both parse, prefers forward (two-phase callers use barcode to choose).
    pub fn split_auto_rc<'a>(
        fwd_info: &'a SegmentInfo,
        read: &'a fastq::Record,
        linker_tolerance: f64,
        anchor_tolerance: f64,
        anchor_max_indels_ratio: f64,
    ) -> Result<SplitResult<'a>, SplitError> {
        let both_res =
            Self::split_auto_rc_both(fwd_info, read, linker_tolerance, anchor_tolerance, anchor_max_indels_ratio)?;
        match both_res {
            SplitAutoRcBothResult::One(res) => Ok(res),
            SplitAutoRcBothResult::Both { forward, .. } => Ok(forward),
        }
    }

    /// Try both forward and reverse orientations; return both results when both parse.
    /// Used for two-phase strand selection (e.g. with barcode validation).
    pub fn split_auto_rc_both<'a>(
        fwd_info: &'a SegmentInfo,
        read: &'a fastq::Record,
        linker_tolerance: f64,
        anchor_tolerance: f64,
        anchor_max_indels_ratio: f64,
    ) -> Result<SplitAutoRcBothResult<'a>, SplitError> {
        let res_fwd =
            fwd_info.split_with_tolerance(read, linker_tolerance, anchor_tolerance, anchor_max_indels_ratio);
        let rc_read = make_rc_record(read);
        let res_rev = fwd_info.split_with_tolerance(
            &rc_read,
            linker_tolerance,
            anchor_tolerance,
            anchor_max_indels_ratio,
        );
        let read_len = read.sequence().len();

        match (res_fwd, res_rev) {
            (Ok(fwd), Ok(rev)) => {
                let (mismatches, mappings) = Self::rev_to_mappings(&rev, &rc_read, read_len, fwd_info);
                let reverse = Self::build_reverse_mapped_result(read, fwd_info, mappings, mismatches);
                Ok(SplitAutoRcBothResult::Both { forward: fwd, reverse })
            }
            (Ok(fwd), Err(_)) => Ok(SplitAutoRcBothResult::One(fwd)),
            (Err(_), Ok(rev)) => {
                let (mismatches, mappings) = Self::rev_to_mappings(&rev, &rc_read, read_len, fwd_info);
                let reverse = Self::build_reverse_mapped_result(read, fwd_info, mappings, mismatches);
                Ok(SplitAutoRcBothResult::One(reverse))
            }
            (Err(e), Err(_)) => Err(e.clone()),
        }
    }

    /// Build (mismatches, Vec<SegmentMapping>) from reverse split result and RC read.
    /// Uses fwd_info to get region_id slices so mappings borrow from fwd_info, not rc_read.
    fn rev_to_mappings<'a>(
        rev: &SplitResult,
        rc_read: &fastq::Record,
        read_len: usize,
        fwd_info: &'a SegmentInfo,
    ) -> (usize, Vec<SegmentMapping<'a>>) {
        let rc_seq_ptr = rc_read.sequence().as_ptr() as usize;
        let mut mappings = Vec::with_capacity(rev.segments.len());
        for seg in &rev.segments {
            let seg_ptr = seg.seq.as_ptr() as usize;
            let start = seg_ptr - rc_seq_ptr;
            let end = start + seg.seq.len();
            let orig_start = read_len - end;
            let orig_end = read_len - start;
            match &seg.region_type {
                SegmentType::Single { region_id, region_type } => {
                    let region_id: &'a str = fwd_info
                        .iter()
                        .find(|e| e.region_id == *region_id)
                        .map(|e| e.region_id.as_str())
                        .expect("segment region_id should exist in fwd_info");
                    mappings.push(SegmentMapping {
                        region_id,
                        region_type_variant: *region_type,
                        is_joined: false,
                        joined_types: Vec::new(),
                        orig_start,
                        orig_end,
                    });
                }
                SegmentType::Joined(types) => {
                    mappings.push(SegmentMapping {
                        region_id: "",
                        region_type_variant: RegionType::Linker,
                        is_joined: true,
                        joined_types: types.clone(),
                        orig_start,
                        orig_end,
                    });
                }
            }
        }
        (rev.mismatches, mappings)
    }

    /// Build SplitResult (reverse orientation) from original read and mappings.
    fn build_reverse_mapped_result<'a>(
        read: &'a fastq::Record,
        fwd_info: &'a SegmentInfo,
        mappings: Vec<SegmentMapping<'a>>,
        mismatches: usize,
    ) -> SplitResult<'a> {
        let mut mapping_map: std::collections::HashMap<&'a str, Vec<SegmentMapping<'a>>> =
            std::collections::HashMap::new();
        let mut joined_mappings = Vec::new();
        for mapping in mappings {
            if mapping.is_joined {
                joined_mappings.push(mapping);
            } else if !mapping.region_id.is_empty() {
                mapping_map
                    .entry(mapping.region_id)
                    .or_insert_with(Vec::new)
                    .push(mapping);
            }
        }
        let mut mapped_segments = Vec::with_capacity(fwd_info.len() + joined_mappings.len());
        for elem in fwd_info.iter() {
            if let Some(mapping_vec) = mapping_map.get_mut(elem.region_id.as_str()) {
                if let Some(mapping) = mapping_vec.pop() {
                    let region_type = SegmentType::Single {
                        region_id: elem.region_id.as_str(),
                        region_type: mapping.region_type_variant,
                    };
                    let seq_slice = &read.sequence()[mapping.orig_start..mapping.orig_end];
                    let qual_slice = &read.quality_scores()[mapping.orig_start..mapping.orig_end];
                    mapped_segments.push(Segment::new(region_type, seq_slice, qual_slice));
                }
            }
        }
        for mapping in joined_mappings {
            let region_type = SegmentType::Joined(mapping.joined_types);
            let seq_slice = &read.sequence()[mapping.orig_start..mapping.orig_end];
            let qual_slice = &read.quality_scores()[mapping.orig_start..mapping.orig_end];
            mapped_segments.push(Segment::new(region_type, seq_slice, qual_slice));
        }
        SplitResult {
            segments: mapped_segments,
            mismatches,
            is_reverse: true,
        }
    }
}

#[derive(Debug, Clone)]
pub struct SegmentInfoElem {
    pub region_id: String,
    pub region_type: RegionType,
    pub sequence_type: SequenceType,
    pub sequence: String,
    /// Precomputed reverse complement of sequence; avoids allocs in split_with_tolerance when is_reverse.
    pub rc_sequence: Vec<u8>,
    pub len: Range<usize>, // length of the segment
}

impl From<Region> for SegmentInfoElem {
    fn from(region: Region) -> Self {
        let rc_sequence = crate::utils::rev_compl(region.sequence.as_bytes());
        Self {
            region_id: region.region_id,
            region_type: region.region_type,
            sequence_type: region.sequence_type,
            sequence: region.sequence,
            rc_sequence,
            len: region.min_len as usize..region.max_len as usize,
        }
    }
}

impl From<&Region> for SegmentInfoElem {
    fn from(region: &Region) -> Self {
        region.clone().into()
    }
}

impl SegmentInfoElem {
    pub fn new(
        region_id: &str,
        region_type: RegionType,
        sequence_type: SequenceType,
        sequence: &str,
        min_len: usize,
        max_len: usize,
    ) -> Self {
        let sequence = sequence.to_string();
        let rc_sequence = crate::utils::rev_compl(sequence.as_bytes());
        Self {
            region_id: region_id.to_string(),
            region_type,
            sequence_type,
            sequence,
            rc_sequence,
            len: min_len..max_len,
        }
    }

    pub fn min_len(&self) -> usize {
        self.len.start
    }

    pub fn max_len(&self) -> usize {
        self.len.end
    }

    pub fn is_fixed_len(&self) -> bool {
        self.len.start == self.len.end
    }

    pub fn is_barcode(&self) -> bool {
        self.region_type.is_barcode()
    }

    pub fn is_umi(&self) -> bool {
        self.region_type.is_umi()
    }
}

/// Zero-allocation pattern search. For short DNA sequences, slice windows are faster than KMP.
/// I have benchmarked these two methods: this implementation is 2x faster than KMP.
fn find_pattern(haystack: &[u8], needle: &[u8]) -> Option<usize> {
    if needle.is_empty() || haystack.is_empty() || needle.len() > haystack.len() {
        return None;
    }
    haystack.windows(needle.len()).position(|w| w == needle)
}

/// Find the best match for a pattern within a haystack, always returning the position
/// with the minimal number of mismatches. Highly optimized for short patterns (< 10 characters).
///
/// # Arguments
///
/// Sliding-window k-bounded Levenshtein search for anchor matching with indels.
/// Returns best match by smallest edit distance, then by length closest to needle.len().
fn find_best_indel_pattern_match(
    haystack: &[u8],
    needle: &[u8],
    max_indels: usize,
) -> Option<AnchorMatch> {
    if needle.is_empty() || haystack.is_empty() {
        return None;
    }
    let m = needle.len();
    let n = haystack.len();
    let min_l = m.saturating_sub(max_indels);
    let max_l = m.saturating_add(max_indels);
    if min_l > max_l || n < min_l {
        return None;
    }

    let mut best_match: Option<AnchorMatch> = None;
    for i in 0..=n.saturating_sub(min_l) {
        let max_window_len = (n - i).min(max_l);
        for l in min_l..=max_window_len {
            let window = &haystack[i..i + l];
            if let Some(d) = edit_distance_bounded(needle, window, max_indels) {
                let is_better = match &best_match {
                    None => true,
                    Some(cur) => {
                        if d < cur.edit_distance {
                            true // 1. Lowest edit distance always wins
                        } else if d == cur.edit_distance {
                            if i < cur.start_pos {
                                true // 2. Prefer leftmost match (closest to expected biological coordinate)
                            } else if i == cur.start_pos {
                                // 3. Same start: prefer length closest to expected needle length
                                let new_len_diff = l.abs_diff(m);
                                let cur_len_diff = cur.match_len.abs_diff(m);
                                if new_len_diff < cur_len_diff {
                                    true
                                } else if new_len_diff == cur_len_diff {
                                    // 4. Ultimate tie-breaker: prefer shorter consumption
                                    l < cur.match_len
                                } else {
                                    false
                                }
                            } else {
                                false // i > cur.start_pos: deeper in read (false positive risk)
                            }
                        } else {
                            false
                        }
                    }
                };
                if is_better {
                    best_match = Some(AnchorMatch {
                        start_pos: i,
                        match_len: l,
                        edit_distance: d,
                    });
                }
            }
        }
    }
    best_match
}

/// * `haystack` - The sequence to search in
/// * `needle` - The pattern to search for
///
/// # Returns
///
/// A tuple of (Option<position>, mismatch_count) with the best match position and its mismatch count
fn find_best_pattern_match(haystack: &[u8], needle: &[u8]) -> (Option<usize>, usize) {
    if needle.is_empty() || haystack.is_empty() || needle.len() > haystack.len() {
        return (None, usize::MAX);
    }
    // Try exact matching first - it's very efficient
    if let Some(pos) = find_pattern(haystack, needle) {
        return (Some(pos), 0);
    }
    let n = haystack.len();
    let m = needle.len();

    // Ultra-optimized special case for 1-character patterns
    if m == 1 {
        // For single character patterns, we either have an exact match or nothing
        // We already checked for exact matches above
        return (None, 1);
    }
    let mut min_mismatches = usize::MAX;
    let mut best_pos = None;
    // Simple and efficient approach for very short patterns
    for i in 0..=n.saturating_sub(m) {
        let mut mismatches = 0;
        // Count mismatches at this position
        for j in 0..m {
            if haystack[i + j] != needle[j] {
                mismatches += 1;
                // Early termination if we exceed current best
                if mismatches >= min_mismatches && min_mismatches != usize::MAX {
                    break;
                }
            }
        }
        // Update if this is a better match
        if mismatches < min_mismatches {
            min_mismatches = mismatches;
            best_pos = Some(i);
        }
    }
    (best_pos, min_mismatches)
}

#[cfg(test)]
mod tests {
    use itertools::Itertools;

    use super::*;

    fn new_fq(seq: &str) -> fastq::Record {
        fastq::Record::new(
            fastq::record::Definition::default(),
            seq.chars()
                .filter(|c| !c.is_whitespace())
                .collect::<String>()
                .as_bytes(),
            '@'.to_string().repeat(seq.len()).as_bytes(),
        )
    }

    fn split_seq_by<I: IntoIterator<Item = SegmentInfoElem>>(seq: &str, elems: I) -> String {
        let segment_info = SegmentInfo::new(elems, false);
        let fq = new_fq(seq);
        segment_info
            .split(&fq)
            .unwrap()
            .into_iter()
            .map(|s| s.to_string())
            .join(",")
    }

    fn bc(n: usize) -> SegmentInfoElem {
        SegmentInfoElem::new("bc", RegionType::Barcode, SequenceType::Random, "N", n, n)
    }

    fn umi(n: usize) -> SegmentInfoElem {
        SegmentInfoElem::new("umi", RegionType::Umi, SequenceType::Random, "N", n, n)
    }

    fn cdna(min: usize, max: usize) -> SegmentInfoElem {
        SegmentInfoElem::new(
            "cdna",
            RegionType::Cdna,
            SequenceType::Random,
            "N",
            min,
            max,
        )
    }

    fn poly(min: usize, max: usize) -> SegmentInfoElem {
        SegmentInfoElem::new(
            "polya",
            RegionType::PolyA,
            SequenceType::Random,
            "N",
            min,
            max,
        )
    }

    fn linker(seq: &str) -> SegmentInfoElem {
        SegmentInfoElem::new(
            "linker",
            RegionType::Linker,
            SequenceType::Fixed,
            seq,
            seq.len(),
            seq.len(),
        )
    }

    fn var(min: usize, max: usize) -> SegmentInfoElem {
        SegmentInfoElem::new(
            "var",
            RegionType::Linker,
            SequenceType::Random,
            "N",
            min,
            max,
        )
    }

    #[test]
    fn test_fixed() {
        assert_eq!(
            split_seq_by("AAAAAA", [cdna(1, 5), linker("ATGT")]),
            "[T]AAAAA",
        );

        assert_eq!(
            split_seq_by("AAAA TT CCC", [bc(4), umi(2), cdna(2, 4), poly(2, 4)]),
            "[B]AAAA,[U]TT,[T]CCC",
        );

        assert_eq!(
            split_seq_by("AATT", [bc(4), umi(2), cdna(2, 4), poly(2, 4)]),
            "[B]AATT",
        );

        assert_eq!(
            split_seq_by("AATT GG T", [bc(4), umi(2), cdna(2, 4), poly(2, 4)]),
            "[B]AATT,[U]GG",
        );

        assert_eq!(
            split_seq_by("AAAA TT CCG GG", [bc(4), umi(2), cdna(2, 4), poly(2, 4)]),
            "[B]AAAA,[U]TT,[TO]CCGGG",
        );

        assert_eq!(
            split_seq_by(
                "AAAA TT CCCCCCC GGGGGGG AAAA TT",
                [bc(4), umi(2), cdna(2, 10), poly(2, 10), bc(4), umi(2)]
            ),
            "[B]AAAA,[U]TT,[TO]CCCCCCCGGGGGGGAAAATT",
        );
    }

    #[test]
    fn test_variable() {
        assert_eq!(
            split_seq_by(
                "AAAA A TT ACTG CCC",
                [bc(4), var(1, 3), umi(2), linker("ACTG")]
            ),
            "[B]AAAA,[O]A,[U]TT,[O]ACTG",
        );

        assert_eq!(
            split_seq_by(
                "AAAA A TT ACCGG CCCACTGG",
                [bc(4), var(1, 3), umi(2), linker("ACTGG")]
            ),
            "[B]AAAA,[O]A,[U]TT,[O]ACCGG",
        );

        assert_eq!(
            split_seq_by(
                "AAAA A TT ACCGG CCCACTG",
                [bc(4), var(1, 3), umi(2), linker("ACTGG"), cdna(10, 100)]
            ),
            "[B]AAAA,[O]A,[U]TT,[O]ACCGG",
        );

        assert_eq!(
            split_seq_by(
                "AAAA A TT ACCGG CCCACTG",
                [bc(4), var(1, 3), umi(2), linker("ACTGG"), cdna(1, 100)]
            ),
            "[B]AAAA,[O]A,[U]TT,[O]ACCGG,[T]CCCACTG",
        );

        assert_eq!(
            split_seq_by(
                "AAAA A TT ACCGG CC CA GGGG TTTT AT GGGG ACTGGGGGACTG",
                [
                    bc(4),
                    var(1, 3),
                    umi(2),
                    linker("ACTGG"),
                    var(2, 4),
                    bc(2),
                    linker("GGGG"),
                    var(2, 4),
                    bc(2),
                    linker("GGGG"),
                    cdna(1, 100),
                ]
            ),
            "[B]AAAA,[O]A,[U]TT,[O]ACCGG,[O]CC,[B]CA,[O]GGGG,[O]TTTT,[B]AT,[O]GGGG,[T]ACTGGGGGACTG",
        );
    }

    #[test]
    fn test_truncate() {
        let info = [
            bc(4),
            var(2, 4),
            bc(2),
            linker("GGGG"),
            cdna(1, 1000),
        ];
        let info = SegmentInfo::new(info, false);
        println!("{:?}", info.truncate_max(50));
    }

    #[test]
    fn test_anchor_indel_tolerance() {
        // Layout: bc(4) + var(1,2) + linker("ACGT"). Read 9 bp: "AAAA" + "A" + "ACGX" (1 sub in linker).
        let elems = [bc(4), var(1, 2), linker("ACGT")];
        let info = SegmentInfo::new(elems, false);
        let read = new_fq("AAAAAACGX");
        // anchor_tolerance 0.3 -> max_allowed_by_ratio = 1 so 1 edit is accepted
        // ratio 0.1 → cap 1 for 4 bp needle (0.1 * 4 = 0.4 → ceil = 1)
        let res_indel = info.split_with_tolerance(&read, 0.3, 0.3, 0.1);
        assert!(
            res_indel.is_ok(),
            "split with anchor_max_indels_ratio=0.1 (cap 1) should succeed: {:?}",
            res_indel
        );
        let split = res_indel.unwrap();
        assert_eq!(split.segments.len(), 3, "barcode + var + linker");
        let anchor_seg = split.segments.last().unwrap();
        // Coordinate-priority tie-breaker: same i, prefer length closest to needle (4 bp). So we take "ACGX" (4 bp, 1 sub).
        assert_eq!(anchor_seg.seq, b"ACGX", "anchor segment (1 edit accepted; length closest to needle)");
    }

    /// True insertion: linker grows by 1 bp; search window expands and downstream cDNA is still correct.
    #[test]
    fn test_anchor_insertion_tolerance() {
        let elems = [bc(4), var(1, 2), linker("ACGT"), cdna(3, 10)];
        let info = SegmentInfo::new(elems, false);
        // "AAAA" (bc) + "A" (var) + "ACXGT" (linker with inserted X) + "CCCC" (cDNA) = 14 bp
        let read = new_fq("AAAAAACXGTCCCC");

        // ratio 0.1 → cap 1 for 4 bp linker
        let res = info.split_with_tolerance(&read, 1.0, 0.3, 0.1);
        assert!(res.is_ok(), "Insertion should be accepted: {:?}", res);

        let split = res.unwrap();
        assert_eq!(split.segments.len(), 4);

        assert_eq!(split.segments[2].seq, b"ACXGT", "Anchor should be 5 bp due to insertion");
        assert_eq!(split.segments[3].seq, b"CCCC", "Downstream cDNA still correct");
    }

    /// Dual-gate rejection: ratio limits allowed edits; read with 2 substitutions is rejected even if absolute cap is high.
    /// Uses var(0, 4) so the anchor branch is taken. Read length must be >= max_offset + max_len (12) so we don't
    /// trigger the "short read" graceful break (which would flush buffer instead of returning an error).
    #[test]
    fn test_anchor_dual_gate_rejection() {
        let elems = [bc(4), var(0, 4), linker("ACGT")];
        let info = SegmentInfo::new(elems, false);
        // Read: "AAAA" + "ACXX" + "YYY" so len=12 -> max_offset+max_len (12) not > len; we return error, not break.
        let read = new_fq("AAAAACXXXYYY");

        // anchor_tolerance 0.3 -> k_bound = 1; ratio 1.25 → cap 5 for 4 bp; we reject because ratio limits edits to 1.
        let res = info.split_with_tolerance(&read, 1.0, 0.3, 1.25);

        assert!(
            res.is_err(),
            "Should reject read because ratio limits edits to 1, even if cap from ratio is 5"
        );

        match res {
            Err(SplitError::AnchorNotFound) => {}
            Err(SplitError::PatternMismatch(dist)) => assert_eq!(dist, 2, "Should report 2 mismatches"),
            Ok(_) => panic!("Expected an error"),
        }
    }

    /// scNanopore-like layout; Read 2 has one error per anchor (seq1=+1bp, seq2=1 sub, seq3=1 del).
    /// Verifies that with tolerance and anchor_max_indels_ratio the split finds all anchors.
    #[test]
    fn test_scnanopore_anchor_one_error_each() {
        use crate::RegionType;
        let leading = SegmentInfoElem::new(
            "leading_primer",
            RegionType::Linker,
            SequenceType::Fixed,
            "TGTACTTCGT",
            10,
            10,
        );
        let seq1 = linker("AATGATACGGCGACCACCGAGATCT");
        let seq2 = linker("TCGTCGGCAGCGTC");
        let seq3 = linker("AGATGTGTATAAGAGACAG");
        let gdna = SegmentInfoElem::new("gdna", RegionType::Gdna, SequenceType::Random, "N", 1, 100_000);
        let elems = [
            leading,
            var(0, 80),
            seq1,
            var(22, 26), // pcr_barcode
            seq2,
            var(22, 26), // tn5_barcode
            seq3,
            gdna,
        ];
        let info = SegmentInfo::new(elems, false);
        // Read 2: seq1 has extra C (ins), seq2 has A instead of G (sub), seq3 has one T missing (del).
        let read2 = "TGTACTTCGTTCAGTTACGTACTAATGATACGGCGACCACCGAGATCTCCACAGTGTCAACTAGAGCCTCTCTCGTCAGCAGCGTCAAGCGTTGAAACCTTTGTCCTCTCAGATGTATAAGAGACAGGTTGAATACACACAAC";
        let record = new_fq(read2);
        // Relaxed tolerance so one edit per anchor is accepted.
        let res = info.split_with_tolerance(&record, 0.1, 0.2, 0.1);
        assert!(
            res.is_ok(),
            "Read 2 (seq1=+1bp, seq2=1sub, seq3=1del) should parse: {:?}",
            res
        );
        let split = res.unwrap();
        assert_eq!(
            split.segments.len(),
            8,
            "expected 8 segments: leading, unfixed, seq1, pcr_barcode, seq2, tn5_barcode, seq3, gdna"
        );
        let gdna_seg = split.segments.last().unwrap();
        assert!(!gdna_seg.seq.is_empty(), "gdna segment should be non-empty");
    }

   
}
