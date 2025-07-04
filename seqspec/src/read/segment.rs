use std::ops::Range;

use noodles::fastq::{self, record::Definition};

use crate::{utils::rev_compl, Region, RegionType, SequenceType};

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
    /// This function uses a default tolerance of 1.0 for the linker sequence and 0.2 for the anchor sequence,
    /// which means that we do not check the sequence of the linker and allow up to 20% mismatches for the anchor sequence.
    pub fn split<'a>(&'a self, read: &'a fastq::Record) -> Result<Vec<Segment<'a>>, SplitError> {
        self.split_with_tolerance(read, 1.0, 0.2)
    }

    /// Split a read into segments according to the segment information.
    ///
    /// # Arguments
    ///
    /// * `read` - The read to split
    /// * `linker_tolerance` - The fraction of mismatches allowed for the linker sequence
    /// * `anchor_tolerance` - The fraction of mismatches allowed for the anchor sequence
    pub fn split_with_tolerance<'a>(
        &'a self,
        read: &'a fastq::Record,
        linker_tolerance: f64,
        anchor_tolerance: f64,
    ) -> Result<Vec<Segment<'a>>, SplitError> {
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
        let mut buffer: Vec<&SegmentInfoElem> = Vec::new();
        let len = read.sequence().len();

        for segment in self.iter() {
            let min_len = segment.min_len();
            let max_len = segment.max_len();

            // Check if the segment is within the read bounds
            if max_offset >= len || min_offset + min_len > len {
                //if max_offset + min_len > len {
                break;
            }

            if min_len == max_len {
                if !buffer.is_empty() {
                    // if the sequcence contains fixed patterns, we can try to use them as anchors
                    // to find the location of the next segment
                    if segment.sequence_type.is_fixed() {
                        // by definition, the pattern can only be found within the window defined by [min_offset, max_offset+max_len]
                        let seq = read
                            .sequence()
                            .get(min_offset..len.min(max_offset + max_len))
                            .unwrap();
                        let match_result = if self.is_reverse {
                            find_best_pattern_match(&rev_compl(seq), segment.sequence.as_bytes())
                        } else {
                            find_best_pattern_match(seq, segment.sequence.as_bytes())
                        };
                        if let (Some(pos), mis) = match_result
                        {
                            if mis as f64 <= anchor_tolerance * segment.sequence.len() as f64 {
                                // [offset_left, offset_right] is the region of the read that
                                // contains segments prior to the matched pattern
                                let offset_left =
                                    min_offset - buffer.iter().map(|s| s.min_len()).sum::<usize>();
                                let offset_right = min_offset + pos;
                                consume_buf(
                                    &mut buffer,
                                    read.sequence().get(offset_left..offset_right).unwrap(),
                                    read.quality_scores()
                                        .get(offset_left..offset_right)
                                        .unwrap(),
                                    &mut result,
                                );

                                // as the lengths of the previous segments are identified, we can now
                                // update the min_offset and max_offset
                                min_offset = offset_right;
                                max_offset = offset_right;

                                let seq = read
                                    .sequence()
                                    .get(offset_right..offset_right + max_len)
                                    .unwrap();
                                let qual = read
                                    .quality_scores()
                                    .get(offset_right..offset_right + max_len)
                                    .unwrap();
                                result.push(Segment::new(
                                    SegmentType::Single {
                                        region_id: segment.region_id.as_str(),
                                        region_type: segment.region_type,
                                    },
                                    seq,
                                    qual,
                                ))
                            } else if max_offset + max_len > len {
                                // if the pattern is not found and we are at the end of the read
                                break;
                            } else {
                                return Err(SplitError::PatternMismatch(mis));
                            }
                        } else {
                            return Err(SplitError::AnchorNotFound);
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
                            hamming::distance(&rev_compl(seq), segment.sequence.as_bytes())
                        } else {
                            hamming::distance(&seq, segment.sequence.as_bytes())
                        };
                        if d as f64 > linker_tolerance * seq.len() as f64 {
                            return Err(SplitError::PatternMismatch(d as usize));
                        }
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

        Ok(result)
    }
}

#[derive(Debug, Clone)]
pub struct SegmentInfoElem {
    pub region_id: String,
    pub region_type: RegionType,
    pub sequence_type: SequenceType,
    pub sequence: String,
    pub len: Range<usize>, // length of the segment
}

impl From<Region> for SegmentInfoElem {
    fn from(region: Region) -> Self {
        Self {
            region_id: region.region_id,
            region_type: region.region_type,
            sequence_type: region.sequence_type,
            sequence: region.sequence,
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
        Self {
            region_id: region_id.to_string(),
            region_type,
            sequence_type,
            sequence: sequence.to_string(),
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
}

/// Find a pattern with the Knuth-Morris-Pratt (KMP) algorithm - more efficient for large patterns
fn find_pattern(haystack: &[u8], needle: &[u8]) -> Option<usize> {
    if needle.is_empty() || haystack.is_empty() || needle.len() > haystack.len() {
        return None;
    }

    // Compute the KMP failure function
    let mut lps = vec![0; needle.len()];
    let mut len = 0;
    let mut i = 1;

    while i < needle.len() {
        if needle[i] == needle[len] {
            len += 1;
            lps[i] = len;
            i += 1;
        } else if len != 0 {
            len = lps[len - 1];
        } else {
            lps[i] = 0;
            i += 1;
        }
    }

    // Search for the pattern
    let mut i = 0; // index for haystack
    let mut j = 0; // index for needle

    while i < haystack.len() {
        if needle[j] == haystack[i] {
            i += 1;
            j += 1;
        }

        if j == needle.len() {
            return Some(i - j);
        } else if i < haystack.len() && needle[j] != haystack[i] {
            if j != 0 {
                j = lps[j - 1];
            } else {
                i += 1;
            }
        }
    }

    None
}

/// Find the best match for a pattern within a haystack, always returning the position
/// with the minimal number of mismatches. Highly optimized for short patterns (< 10 characters).
///
/// # Arguments
///
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
}
