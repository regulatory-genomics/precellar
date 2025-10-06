use anyhow::Result;
use bincode::{Decode, Encode};
use noodles::sam;
use noodles::sam::alignment::{record::cigar::op::Kind, record_buf::RecordBuf};

use bed_utils::bed::map::GIntervalMap;
use bed_utils::bed::GenomicRange;
use itertools::Itertools;
use log::info;
use rayon::iter::ParallelIterator;
use rayon::slice::ParallelSlice;

use crate::align::MultiMapR;
use crate::fragment::Fragment;
use crate::transcriptome::{Exon, Transcript};

#[derive(Debug, Copy, Clone)]
pub enum ChemistryStrandness {
    Forward,
    Reverse,
    Unstranded,
}

#[derive(Debug, Copy, Clone, PartialEq, Eq, Encode, Decode)]
pub enum Orientation {
    Sense,
    Antisense,
}

/// Indicates whether a transcript alignment is exonic, intronic, or spans exon-intron boundaries.
#[derive(Eq, PartialEq, Debug, Clone, Encode, Decode)]
pub enum AlignType {
    Exonic,
    Intronic,
    Spanning,
    Discordant,
}

#[derive(Eq, PartialEq, Debug, Clone, Encode, Decode)]
pub struct TxAlignResult {
    pub gene_id: String,
    pub transcript_id: String,
    pub strand: Orientation,
    align_type: AlignType,
}

impl TxAlignResult {
    pub fn is_exonic(&self) -> bool {
        self.align_type == AlignType::Exonic
    }

    pub fn is_intronic(&self) -> bool {
        self.align_type == AlignType::Intronic
    }

    pub fn is_spanning(&self) -> bool {
        self.align_type == AlignType::Spanning
    }

    pub fn is_discordant(&self) -> bool {
        self.align_type == AlignType::Discordant
    }

    pub fn is_transcriptomic(&self) -> bool {
        !self.is_discordant() && self.strand == Orientation::Sense
    }
}

/// Segment represents a contiguous block of cigar operations not containing
/// any "Skip" sections. The SpliceSegment is 0-based, half-open with respect to the reference.
#[derive(Debug, Clone, Encode, Decode)]
struct Segment {
    start: u64,
    end: u64,
}

/// SplicedRecord represents segments of a larger whole that have been cut and joined together
/// to form a new, continuous sequence.
/// It consists of the left and right clipping operations, and a list of components.
/// A SplicedRecord can be constructed from a BAM record by splitting the CIGAR string
/// at each "Skip" operation.
#[derive(Debug, Clone, Encode, Decode)]
pub struct SplicedRecord {
    chrom: String,
    segments: Vec<Segment>,
    strand: bed_utils::bed::Strand, // Strand information of the read.
}

impl SplicedRecord {
    fn new(read: &RecordBuf, header: &sam::Header) -> Result<Option<Self>> {
        if read.flags().is_unmapped() {
            return Ok(None);
        }

        let chrom = read.reference_sequence(header).unwrap()?.0;
        let chrom = std::str::from_utf8(chrom)?.to_string();
        let cigar = read.cigar();
        let alignment_start = read.alignment_start().unwrap().get() - 1;

        let mut splice_segments: Vec<Segment> = Vec::new();
        let mut curr_segment = Segment {
            start: alignment_start as u64,
            end: alignment_start as u64,
        };

        for c in cigar.as_ref() {
            match c.kind() {
                Kind::HardClip | Kind::SoftClip | Kind::Insertion => {}
                Kind::Skip => {
                    let next_start = curr_segment.end + c.len() as u64;
                    splice_segments.push(curr_segment);
                    curr_segment = Segment {
                        start: next_start,
                        end: next_start,
                    };
                }
                Kind::Match | Kind::Deletion | Kind::SequenceMatch | Kind::SequenceMismatch => {
                    curr_segment.end += c.len() as u64;
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
            segments: splice_segments,
            strand,
        }))
    }

    /// The leftmost position of the alignment.
    fn start(&self) -> u64 {
        self.segments.first().map_or(0, |segment| segment.start)
    }

    /// The rightmost position of the alignment (exclusive).
    fn end(&self) -> u64 {
        self.segments.last().map_or(0, |segment| segment.end)
    }

    /// Aligns a read to a transcript and determines the region type (exonic, intronic, or intergenic).
    fn align_transcript(
        &self,
        transcript: &Transcript,
        chemistry_strandness: ChemistryStrandness,
    ) -> TxAlignResult {
        let is_antisense = match chemistry_strandness {
            ChemistryStrandness::Forward => transcript.strand != self.strand,
            ChemistryStrandness::Reverse => transcript.strand == self.strand,
            ChemistryStrandness::Unstranded => false, // if unstranded, always sense
        };
        let strand = if is_antisense {
            Orientation::Antisense
        } else {
            Orientation::Sense
        };

        TxAlignResult {
            gene_id: transcript.gene_id.clone(),
            transcript_id: transcript.id.clone(),
            strand,
            align_type: self.align_helper(transcript),
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
    fn orientation(&self, transcript: &Transcript) -> Option<Orientation> {
        if self.start() >= transcript.start && self.end() <= transcript.end {
            if transcript.strand == self.strand {
                Some(Orientation::Sense)
            } else {
                Some(Orientation::Antisense)
            }
        } else {
            None
        }
    }

    /// Align segments to exons and determine if the alignment is exonic or is intronic or
    /// spans exon-intron boundaries.
    ///
    /// 1. If the read does not overlap the transcript, it is considered unaligned.
    ///
    /// Transcript:  |||||||----|||||-------|||||---||
    /// Exonic:        XXXXX----XXXX--------XXXXX
    /// Intronic:            XX
    /// Spanning:        XXXXXXXXXXXX-------XXXXX
    fn align_helper(&self, transcript: &Transcript) -> AlignType {
        if self.end() <= transcript.start || self.start() >= transcript.end {
            panic!("Read does not overlap transcript");
        }

        let exons = find_exons(transcript.exons(), self.start(), self.end() - 1);
        if exons.is_empty() {
            // No exon overlaps the read, but the read overlaps the transcript.
            return AlignType::Intronic;
        }

        let segments = &self.segments;
        let n_ex = exons.len();
        let n_segment = segments.len();
        if n_ex == n_segment {
            if n_ex == 1 {
                if segments[0].start >= exons[0].start && segments[0].end <= exons[0].end {
                    // Fully contained within the exon.
                    AlignType::Exonic
                } else {
                    AlignType::Spanning
                }
            } else if segments[1..n_segment - 1]
                .iter()
                .zip(exons[1..n_ex - 1].iter())
                .all(|(s, e)| s.start == e.start && s.end == e.end)
            {
                // All internal segments (except the 1st and last) match exactly.
                if segments[0].end == exons[0].end
                    && segments[n_segment - 1].start == exons[n_ex - 1].start
                {
                    // first and last segments also match.
                    if segments[0].start >= exons[0].start
                        && segments[n_segment - 1].end <= exons[n_ex - 1].end
                    {
                        AlignType::Exonic
                    } else {
                        AlignType::Spanning
                    }
                } else {
                    AlignType::Discordant
                }
            } else {
                AlignType::Discordant
            }
        } else if n_ex > n_segment {
            // There are more exons than segments. Possibly spanning.
            if exons[n_ex - 1].start >= segments[n_segment - 1].start
                && exons[0].end > segments[0].start
                && exons[1..n_ex - 1].iter().all(|_| true)
            // FIXME: this is a bit hacky
            {
                AlignType::Spanning
            } else {
                AlignType::Discordant
            }
        } else {
            // There are more segments than exons. Cannot be concordant.
            AlignType::Discordant
        }
    }
}

/// Returns the exons that overlap the interval.
/// ||||||-------|||||||
///    ^^^^^^^^^^^^^
fn find_exons(
    exon_info: &[Exon],
    read_start: u64,
    read_end: u64, // inclusive
) -> &[Exon] {
    // find first exon that ends to the right of the read start
    let ex_start = exon_info
        .binary_search_by_key(&read_start, |ex| ex.end - 1)
        .map_or_else(|i| i, |i| i);
    // find first exon that starts to the left of the read end
    if let Some(ex_end) = exon_info
        .binary_search_by_key(&read_end, |ex| ex.start)
        .map_or_else(|i| if i > 0 { Some(i - 1) } else { None }, |i| Some(i))
    {
        if ex_start >= exon_info.len() {
            &[]
        } else {
            &exon_info[ex_start..=ex_end]
        }
    } else {
        &[]
    }
}

#[derive(Debug, Encode, Decode)]
pub enum TxAlignment {
    Intergenic,
    Discordant,
    Multimapped,
    Antisense,
    SeAligned {
        read: SplicedRecord,
        alignment: Vec<TxAlignResult>,
        barcode: String,
        umi: Option<String>,
    },
}

impl TxAlignment {
    pub fn alignments(&self) -> Box<dyn Iterator<Item = &TxAlignResult> + '_> {
        match self {
            TxAlignment::SeAligned { alignment, .. } => Box::new(alignment.iter()),
            _ => Box::new(std::iter::empty()),
        }
    }

    /// Returns the gene if the read is confidently mapped to a single gene.
    pub fn uniquely_mapped_gene(&self) -> Option<&str> {
        match self {
            TxAlignment::SeAligned { alignment, .. } => Some(alignment[0].gene_id.as_str()),
            _ => None,
        }
    }

    pub fn barcode(&self) -> Option<&str> {
        match self {
            TxAlignment::SeAligned { barcode, .. } => Some(barcode.as_str()),
            _ => None,
        }
    }

    pub fn umi(&self) -> Option<&str> {
        match self {
            TxAlignment::SeAligned { umi, .. } => umi.as_deref(),
            _ => None,
        }
    }

    pub fn is_exonic_only(&self) -> bool {
        match self {
            TxAlignment::SeAligned { alignment, .. } => alignment.iter().all(|x| x.is_exonic()),
            _ => false,
        }
    }

    pub fn is_intronic_only(&self) -> bool {
        match self {
            TxAlignment::SeAligned { alignment, .. } => alignment.iter().all(|x| x.is_intronic()),
            _ => false,
        }
    }

    /// A alignment is spanning if it spans the exon-intron boundary at least in
    /// one transcript, and does NOT align to exons in any other transcripts.
    pub fn is_spanning(&self) -> bool {
        match self {
            TxAlignment::SeAligned { alignment, .. } => {
                alignment.iter().any(|x| x.is_spanning())
                    && alignment.iter().all(|x| !x.is_exonic())
            }
            _ => false,
        }
    }

    pub fn to_fragments(&self) -> Box<dyn Iterator<Item = Fragment> + '_> {
        match self {
            TxAlignment::SeAligned { read, .. } => Box::new(read.to_fragments()),
            _ => Box::new(std::iter::empty()),
        }
    }
}

/// Manages the annotation of alignments using transcriptome data.
#[derive(Debug, Clone)]
pub struct TxAligner {
    /// Map of genomic intervals to transcripts.
    transcripts: GIntervalMap<Transcript>,
    chemistry_strandness: Option<ChemistryStrandness>,
    header: sam::Header,
}

impl TxAligner {
    /// Creates a new `AlignmentAnnotator` with the provided transcripts.
    pub fn new(
        transcripts: impl IntoIterator<Item = Transcript>,
        header: sam::Header,
        chemistry_strandness: Option<ChemistryStrandness>,
    ) -> Self {
        let transcripts = transcripts
            .into_iter()
            .map(|x| (GenomicRange::new(&x.chrom, x.start, x.end), x))
            .collect();
        Self {
            transcripts,
            header,
            chemistry_strandness,
        }
    }

    pub fn transcripts(&self) -> impl Iterator<Item = &Transcript> {
        self.transcripts.iter().map(|(_, t)| t)
    }

    pub fn align<'a>(
        &'a self,
        data: impl IntoIterator<Item = Vec<(Option<MultiMapR>, Option<MultiMapR>)>> + 'a,
    ) -> impl Iterator<Item = Option<TxAlignment>> + 'a {
        fn helper<'a>(
            tx_aligner: &'a TxAligner,
            records: &'a [(Option<MultiMapR>, Option<MultiMapR>)],
            strandness: ChemistryStrandness,
        ) -> impl Iterator<Item = Option<TxAlignment>> + 'a {
            records.into_iter().map(move |(read1, read2)| {
                if read1.is_some() && read2.is_some() {
                    tx_aligner._align_pe(read1.as_ref().unwrap(), read2.as_ref().unwrap())
                } else {
                    let read = if read1.is_some() {
                        read1.as_ref().unwrap()
                    } else {
                        read2.as_ref().unwrap()
                    };
                    tx_aligner._align_se(&read, strandness)
                }
            })
        }

        let mut chemistry_strandness = self.chemistry_strandness;
        data.into_iter().flat_map(move |recs| {
            if chemistry_strandness.is_none() {
                chemistry_strandness = Some(self.detect_strandness(&recs));
            }
            recs.par_chunks(4096)
                .flat_map_iter(|chunk| helper(self, chunk, chemistry_strandness.unwrap()))
                .collect::<Vec<_>>()
        })
    }

    /// Annotate the alignments by mapping them to the transcriptome. If multiple
    /// alignments are present, we will try to find the confident ones and promote
    /// them to primary. A read may align to multiple transcripts and genes, but
    /// it is only considered confidently mapped to the transcriptome if it is
    /// mapped to a single gene. The confident alignment will be returned if found.
    ///
    /// # Arguments
    /// * `header` - Reference to the SAM header.
    /// * `multi_map` - Vector of single-end alignment records.
    ///
    /// Returns `Some(AnnotatedAlignment)` if a confident alignment is found, otherwise `None`.
    fn _align_se(&self, multi_map: &MultiMapR, strandness: ChemistryStrandness) -> Option<TxAlignment> {
        let barcode = multi_map.barcode().unwrap()?.clone();
        let aln: Vec<_> = multi_map
            .iter()
            .flat_map(|rec| self._align_one(rec, strandness))
            .collect();
        if aln.is_empty() {
            None
        } else if aln.iter().all(|(_, items)| items.is_empty()) {
            Some(TxAlignment::Intergenic)
        } else {
            let num_unique_genes: usize = aln
                .iter()
                .flat_map(|(_, items)| {
                    items.iter().filter_map(|x| {
                        if x.is_transcriptomic() {
                            Some(x.gene_id.clone())
                        } else {
                            None
                        }
                    })
                })
                .unique()
                .count();

            if num_unique_genes > 1 {
                Some(TxAlignment::Multimapped)
            } else if num_unique_genes == 0 {
                if aln.len() > 1 {
                    Some(TxAlignment::Discordant)
                } else {
                    let has_antisense = aln
                        .iter()
                        .flat_map(|(_, items)| {
                            items.iter().filter_map(|x| {
                                if x.is_discordant() {
                                    None
                                } else {
                                    Some(x.strand == Orientation::Antisense)
                                }
                            })
                        })
                        .any(|x| x);
                    if has_antisense {
                        Some(TxAlignment::Antisense)
                    } else {
                        Some(TxAlignment::Discordant)
                    }
                }
            } else {
                // Now we know there's exactly one gene, we choose the first non-discordant alignment
                let (read, alignment) = aln
                    .into_iter()
                    .map(|(rec, x)| {
                        let result: Vec<_> =
                            x.into_iter().filter(|x| x.is_transcriptomic()).collect();
                        (rec, result)
                    })
                    .find(|(_, x)| !x.is_empty())
                    .unwrap();
                Some(TxAlignment::SeAligned {
                    read,
                    alignment,
                    barcode,
                    umi: multi_map.umi().unwrap().clone(),
                })
            }
        }
    }

    fn _align_pe(&self, multi_map1: &MultiMapR, multi_map2: &MultiMapR) -> Option<TxAlignment> {
        todo!()
    }

    fn _align_one(&self, rec: &RecordBuf, strandness: ChemistryStrandness) -> Option<(SplicedRecord, Vec<TxAlignResult>)> {
        let spliced_rec = SplicedRecord::new(rec, &self.header).unwrap()?;
        let region = GenomicRange::new(
            &spliced_rec.chrom,
            spliced_rec.start() as u64,
            spliced_rec.end() as u64,
        );
        let alignments = self
            .transcripts
            .find(&region)
            .map(|(_, transcript)| spliced_rec.align_transcript(transcript, strandness))
            .collect();
        Some((spliced_rec, alignments))
    }

    /// Detects the strandness of the RNA-seq data based on the alignments.
    fn detect_strandness(
        &self,
        alignments: &[(Option<MultiMapR>, Option<MultiMapR>)],
    ) -> ChemistryStrandness {
        let mut num_sense = 0;
        let mut num_antisense = 0;
        alignments
            .into_iter()
            .filter(|(rec1, rec2)| {
                rec1.as_ref().map_or(true, |x| x.is_confidently_mapped())
                    && rec2.as_ref().map_or(true, |x| x.is_confidently_mapped())
            })
            .for_each(|(rec1, rec2)| {
                let mut is_sense = false;
                let mut is_antisense = false;

                rec1.iter().chain(rec2.iter()).for_each(|aln| {
                    if let Some(rec) = SplicedRecord::new(&aln.primary, &self.header).unwrap() {
                        let region =
                            GenomicRange::new(&rec.chrom, rec.start() as u64, rec.end() as u64);
                        self.transcripts.find(&region).for_each(|(_, transcript)| {
                            match rec.orientation(transcript) {
                                Some(Orientation::Sense) => is_sense = true,
                                Some(Orientation::Antisense) => is_antisense = true,
                                None => {}
                            }
                        });
                        if is_sense {
                            num_sense += 1;
                        } else if is_antisense {
                            num_antisense += 1;
                        }
                    }
                });
            });

        let total = num_sense + num_antisense;
        let percent_sense = num_sense as f64 / total as f64 * 100.0;
        let percent_antisense = num_antisense as f64 / total as f64 * 100.0;
        info!(
            "Found {:.3}% sense and {:.3}% antisense alignments",
            percent_sense, percent_antisense
        );
        if percent_sense > 75.0 {
            info!("Chemistry strandness is set to: Forward");
            ChemistryStrandness::Forward
        } else if percent_antisense > 75.0 {
            info!("Chemistry strandness is set to: Reverse");
            ChemistryStrandness::Reverse
        } else {
            info!("Chemistry strandness is set to: Unstranded");
            ChemistryStrandness::Unstranded
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::transcriptome::Exons;

    use super::*;

    #[test]
    fn test_align() {
        let transcript = Transcript {
            id: "tx1".to_string(),
            chrom: "chr1".to_string(),
            start: 5,
            end: 8000,
            strand: bed_utils::bed::Strand::Forward,
            gene_id: "gene1".to_string(),
            gene_name: "GeneName1".to_string(),
            exons: Exons(vec![
                Exon { start: 10, end: 20 },
                Exon { start: 30, end: 40 },
                Exon { start: 50, end: 60 },
                Exon {
                    start: 500,
                    end: 700,
                },
                Exon {
                    start: 1500,
                    end: 1550,
                },
                Exon {
                    start: 5000,
                    end: 6000,
                },
            ]),
        };

        let mut record = SplicedRecord {
            chrom: "chr1".to_string(),
            segments: vec![],
            strand: bed_utils::bed::Strand::Forward,
        };

        let mut test_fun = |segments: Vec<(u64, u64)>, expected: AlignType| {
            record.segments = segments
                .into_iter()
                .map(|(s, e)| Segment { start: s, end: e })
                .collect();
            assert_eq!(record.align_helper(&transcript), expected);
        };

        // Incompatible
        test_fun(vec![(12, 18), (32, 38)], AlignType::Discordant);
        test_fun(vec![(21, 29), (32, 38)], AlignType::Discordant);
        test_fun(
            vec![(5, 40), (50, 60), (500, 700), (1400, 1501)],
            AlignType::Discordant,
        );

        // Exonic
        test_fun(vec![(12, 20), (30, 38)], AlignType::Exonic);

        // Intronic
        test_fun(vec![(6, 8)], AlignType::Intronic);

        // Spanning
        test_fun(vec![(29, 38)], AlignType::Spanning);
        test_fun(
            vec![(5, 40), (50, 60), (500, 700), (1500, 1501)],
            AlignType::Spanning,
        );
    }
}
