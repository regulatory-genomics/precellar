use std::collections::HashMap;

use anyhow::Result;
use bitcode::{Decode, Encode};
use noodles::sam;
use noodles::sam::alignment::{record::cigar::op::Kind, record_buf::RecordBuf};

use bed_utils::bed::{map::GIntervalMap, GenomicRange};
use itertools::Itertools;
use log::{info, warn};
use rayon::iter::ParallelIterator;
use rayon::slice::ParallelSlice;
use seqspec::ChemistryStrandedness;

use crate::align::MultiMapR;
use crate::fragment::Fragment;
use crate::transcriptome::{Exon, Transcript};

#[derive(Debug, Copy, Clone, PartialEq, Eq, Encode, Decode)]
pub enum Orientation {
    Sense,
    Antisense,
}

/// Indicates whether a transcript alignment is exonic, intronic, or spans exon-intron boundaries.
#[derive(Eq, PartialEq, Debug, Clone, Encode, Decode)]
pub enum AlignType {
    Exonic,     // Read maps entirely within exons
    Intronic,   // Read maps entirely within introns
    Spanning,   // Read spans exon-intron boundaries
    Discordant, // Read maps discordantly to the transcript (splice junctions do not match)
}

impl AlignType {
    pub fn is_exonic(&self) -> bool {
        matches!(self, AlignType::Exonic)
    }

    pub fn is_intronic(&self) -> bool {
        matches!(self, AlignType::Intronic)
    }

    pub fn is_spanning(&self) -> bool {
        matches!(self, AlignType::Spanning)
    }

    pub fn is_discordant(&self) -> bool {
        matches!(self, AlignType::Discordant)
    }
}

#[derive(Eq, PartialEq, Debug, Clone, Encode, Decode)]
pub struct TxAlignResult {
    pub gene_id: String,
    pub transcript_id: String,
    pub strand: Orientation,
    align_type: AlignType,
}

impl TxAlignResult {
    /// Returns true if the alignment is to the transcriptome (i.e., not discordant and sense strand).
    pub fn is_transcriptomic(&self) -> bool {
        !self.align_type.is_discordant() && self.strand == Orientation::Sense
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
    fn new(read: &RecordBuf, header: &sam::Header, is_read2: bool) -> Result<Option<Self>> {
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
                Kind::HardClip | Kind::SoftClip | Kind::Insertion => {}
                Kind::Pad => unreachable!(),
            }
        }
        splice_segments.push(curr_segment);

        let flag = read.flags();
        let strand = if flag.is_reverse_complemented() {
            if is_read2 {
                bed_utils::bed::Strand::Forward
            } else {
                bed_utils::bed::Strand::Reverse
            }
        } else {
            if is_read2 {
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
        strandedness: ChemistryStrandedness,
    ) -> TxAlignResult {
        let is_antisense = match strandedness {
            ChemistryStrandedness::Forward => transcript.strand != self.strand,
            ChemistryStrandedness::Reverse => transcript.strand == self.strand,
            ChemistryStrandedness::Unstranded => false, // if unstranded, always sense
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

/// Represents the result of aligning a read to the transcriptome.
#[derive(Debug, Encode, Decode)]
pub enum TxAlignment {
    Intergenic,
    /// Read maps to intergenic region
    Discordant,
    /// Read maps discordantly
    Multimapped,
    /// Read maps to multiple genes
    Antisense,
    /// Read maps antisense to a gene
    SeAligned {
        read: SplicedRecord,
        alignment: Vec<(String, AlignType)>, // transcritopmic alignments represented as (transcript_id, AlignType)
        gene_id: String,
        barcode: String,
        umi: Option<String>,
    },
    PeAligned {
        read1: SplicedRecord,
        read2: SplicedRecord,
        alignment: Vec<(String, AlignType)>, // transcritopmic alignments represented as (transcript_id, AlignType)
        gene_id: String,
        barcode: String,
        umi: Option<String>,
    },
}

impl TxAlignment {
    pub fn alignments(&self) -> Box<dyn Iterator<Item = &(String, AlignType)> + '_> {
        match self {
            TxAlignment::SeAligned { alignment, .. } => Box::new(alignment.iter()),
            TxAlignment::PeAligned { alignment, .. } => Box::new(alignment.iter()),
            _ => Box::new(std::iter::empty()),
        }
    }

    /// Returns the gene if the read is confidently mapped to a single gene.
    pub fn uniquely_mapped_gene(&self) -> Option<&str> {
        match self {
            TxAlignment::SeAligned { gene_id, .. } => Some(gene_id.as_str()),
            TxAlignment::PeAligned { gene_id, .. } => Some(gene_id.as_str()),
            _ => None,
        }
    }

    pub fn barcode(&self) -> Option<&str> {
        match self {
            TxAlignment::SeAligned { barcode, .. } => Some(barcode.as_str()),
            TxAlignment::PeAligned { barcode, .. } => Some(barcode.as_str()),
            _ => None,
        }
    }

    pub fn umi(&self) -> Option<&str> {
        match self {
            TxAlignment::SeAligned { umi, .. } => umi.as_deref(),
            TxAlignment::PeAligned { umi, .. } => umi.as_deref(),
            _ => None,
        }
    }

    pub fn is_confidently_mapped(&self) -> bool {
        match self {
            TxAlignment::SeAligned { .. } | TxAlignment::PeAligned { .. } => true,
            _ => false,
        }
    }

    pub fn is_exonic_only(&self) -> bool {
        match self {
            TxAlignment::SeAligned { alignment, .. } => alignment.iter().all(|x| x.1.is_exonic()),
            TxAlignment::PeAligned { alignment, .. } => alignment.iter().all(|x| x.1.is_exonic()),
            _ => false,
        }
    }

    pub fn is_intronic_only(&self) -> bool {
        match self {
            TxAlignment::SeAligned { alignment, .. } => alignment.iter().all(|x| x.1.is_intronic()),
            TxAlignment::PeAligned { alignment, .. } => alignment.iter().all(|x| x.1.is_intronic()),
            _ => false,
        }
    }

    /// A alignment is spanning if it spans the exon-intron boundary at least in
    /// one transcript, and does NOT align to exons in any other transcripts.
    pub fn is_spanning(&self) -> bool {
        match self {
            TxAlignment::SeAligned { alignment, .. } => {
                alignment.iter().any(|x| x.1.is_spanning())
                    && alignment.iter().all(|x| !x.1.is_exonic())
            }
            TxAlignment::PeAligned { alignment, .. } => {
                alignment.iter().any(|x| x.1.is_spanning())
                    && alignment.iter().all(|x| !x.1.is_exonic())
            }
            _ => false,
        }
    }

    pub fn to_fragments(&self) -> Box<dyn Iterator<Item = Fragment> + '_> {
        match self {
            TxAlignment::SeAligned { read, .. } => Box::new(read.to_fragments()),
            TxAlignment::PeAligned { read1, read2, .. } => Box::new(read1.to_fragments().chain(read2.to_fragments())),
            _ => Box::new(std::iter::empty()),
        }
    }
}

/// Manages the annotation of alignments using transcriptome data.
#[derive(Debug, Clone)]
pub struct TxAligner {
    /// Map of genomic intervals to transcripts.
    transcripts: GIntervalMap<Transcript>,
    strandedness: Option<ChemistryStrandedness>,
    header: sam::Header,
}

impl TxAligner {
    /// Creates a new `AlignmentAnnotator` with the provided transcripts.
    pub fn new(
        transcripts: impl IntoIterator<Item = Transcript>,
        header: sam::Header,
        strandedness: Option<ChemistryStrandedness>,
    ) -> Self {
        let transcripts = transcripts
            .into_iter()
            .map(|x| (GenomicRange::new(&x.chrom, x.start, x.end), x))
            .collect();
        Self {
            transcripts,
            header,
            strandedness,
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
            strandedness: ChemistryStrandedness,
        ) -> impl Iterator<Item = Option<TxAlignment>> + 'a {
            records.into_iter().map(move |(read1, read2)| {
                if read1.is_some() && read2.is_some() {
                    tx_aligner._align_pe(
                        read1.as_ref().unwrap(),
                        read2.as_ref().unwrap(),
                        strandedness,
                    )
                } else if read1.is_some() {
                    tx_aligner._align_se(read1.as_ref().unwrap(), false, strandedness)
                } else {
                    tx_aligner._align_se(read2.as_ref().unwrap(), true, strandedness)
                }
            })
        }

        let mut strandedness = self.strandedness;
        data.into_iter().flat_map(move |recs| {
            if strandedness.is_none() {
                warn!("Strandedness not provided, detecting from data ...");
                strandedness = Some(self.detect_strandedness(&recs));
            }
            recs.par_chunks(4096)
                .flat_map_iter(|chunk| helper(self, chunk, strandedness.unwrap()))
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
    /// * `multi_map` - Vector of single-end alignment records.
    /// * `is_read2` - Indicates if the read is the second read in a paired-end alignment.
    /// * `strandedness` - The strandedness of the RNA-seq data.
    ///
    /// Returns `Some(AnnotatedAlignment)` if a confident alignment is found, otherwise `None`.
    fn _align_se(
        &self,
        multi_map: &MultiMapR,
        is_read2: bool,
        strandedness: ChemistryStrandedness,
    ) -> Option<TxAlignment> {
        let barcode = multi_map.barcode().unwrap()?.clone();
        let aln: Vec<_> = multi_map
            .iter()
            .flat_map(|rec| {
                let rec = SplicedRecord::new(rec, &self.header, is_read2).unwrap()?;
                let aln = self._align_one(&rec, strandedness);
                Some((rec, aln))
            })
            .collect();
        if aln.is_empty() {
            None
        } else if aln.iter().all(|(_, items)| items.is_empty()) {
            // All alignments are intergenic
            Some(TxAlignment::Intergenic)
        } else {
            // Count the number of unique genes that the read maps to transcriptomically (i.e., sense and not discordant)
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
                                if x.align_type.is_discordant() {
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
                    gene_id: alignment[0].gene_id.clone(),
                    alignment: alignment
                        .into_iter()
                        .map(|x| (x.transcript_id, x.align_type))
                        .collect(),
                    barcode,
                    umi: multi_map.umi().unwrap().clone(),
                })
            }
        }
    }

    fn _align_pe(
        &self,
        multi_map1: &MultiMapR,
        multi_map2: &MultiMapR,
        strandedness: ChemistryStrandedness,
    ) -> Option<TxAlignment> {
        let barcode = multi_map1.barcode().unwrap()?.clone();
        let mut has_antisense = false;
        let mut has_discordant = false;
        let mut multimapped = false;
        let mut gene_id = None;
        let aln: Vec<_> = multi_map1
            .iter()
            .zip(multi_map2.iter())
            .flat_map(|(rec1, rec2)| {
                let mut alignments: Vec<(String, AlignType)> = Vec::new();
                let rec1 = SplicedRecord::new(rec1, &self.header, false).unwrap()?;
                let rec2 = SplicedRecord::new(rec2, &self.header, true).unwrap()?;
                let aln1: HashMap<_, _> = self
                    ._align_one(&rec1, strandedness)
                    .into_iter()
                    .map(|x| (x.transcript_id.clone(), x))
                    .collect();
                self._align_one(&rec2, strandedness).into_iter().for_each(|a2| {
                    if let Some(a1) = aln1.get(&a2.transcript_id) {
                        if a1.strand == a2.strand {
                            if a1.strand == Orientation::Antisense {
                                has_antisense = true;
                            } else {
                                let ty = match (&a1.align_type, a2.align_type) {
                                    (AlignType::Intronic, AlignType::Intronic) => AlignType::Intronic,
                                    (AlignType::Exonic, AlignType::Exonic) => AlignType::Exonic,
                                    (AlignType::Discordant, _) | (_, AlignType::Discordant) => AlignType::Discordant,
                                    _ => AlignType::Spanning,
                                };
                                if ty == AlignType::Discordant {
                                    has_discordant = true;
                                } else {
                                    if let Some(gid) = &gene_id {
                                        if gid != &a1.gene_id {
                                            multimapped = true;
                                        }
                                    } else {
                                        gene_id = Some(a1.gene_id.clone());
                                    }
                                    alignments.push((a1.transcript_id.clone(), ty));
                                }
                            }
                        }
                    }
                });
                Some((rec1, rec2, alignments))
            }).collect();

        if aln.is_empty() {
            None
        } else if multimapped {
            Some(TxAlignment::Multimapped)
        } else if gene_id.is_none() {
            if has_discordant {
                Some(TxAlignment::Discordant)
            } else if has_antisense {
                Some(TxAlignment::Antisense)
            } else {
                Some(TxAlignment::Intergenic)
            }
        } else {
            // Now we know there's exactly one gene, we choose the first non-discordant alignment
            let (read1, read2, alignment) = aln.into_iter().next().unwrap();
            Some(TxAlignment::PeAligned {
                read1,
                read2,
                gene_id: gene_id.unwrap(),
                alignment,
                barcode,
                umi: multi_map1.umi().unwrap().clone(),
            })
        }
    }

    fn _align_one(
        &self,
        rec: &SplicedRecord,
        strandedness: ChemistryStrandedness,
    ) -> Vec<TxAlignResult> {
        let region = GenomicRange::new(&rec.chrom, rec.start() as u64, rec.end() as u64);
        self.transcripts
            .find(&region)
            .map(|(_, transcript)| rec.align_transcript(transcript, strandedness))
            .collect()
    }

    /// Detects the strandedness of the RNA-seq data based on the alignments.
    fn detect_strandedness(
        &self,
        alignments: &[(Option<MultiMapR>, Option<MultiMapR>)],
    ) -> ChemistryStrandedness {
        let mut num_sense = 0;
        let mut num_antisense = 0;

        alignments
            .into_iter()
            .flat_map(|(rec1, rec2)| {
                let rec1 = rec1.as_ref().and_then(|x| {
                    if x.is_confidently_mapped() {
                        SplicedRecord::new(&x.primary, &self.header, false).unwrap()
                    } else {
                        None
                    }
                });
                let rec2 = rec2.as_ref().and_then(|x| {
                    if x.is_confidently_mapped() {
                        SplicedRecord::new(&x.primary, &self.header, true).unwrap()
                    } else {
                        None
                    }
                });
                [rec1, rec2]
            })
            .flatten()
            .for_each(|aln| {
                let mut is_sense = false;
                let mut is_antisense = false;

                let region = GenomicRange::new(&aln.chrom, aln.start() as u64, aln.end() as u64);
                self.transcripts.find(&region).for_each(|(_, transcript)| {
                    match aln.orientation(transcript) {
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
            });

        let total = num_sense + num_antisense;
        let percent_sense = num_sense as f64 / total as f64 * 100.0;
        let percent_antisense = num_antisense as f64 / total as f64 * 100.0;
        info!(
            "Found {:.3}% sense and {:.3}% antisense alignments",
            percent_sense, percent_antisense
        );
        if percent_sense > 67.0 {
            info!("Chemistry strandedness is set to: Forward");
            ChemistryStrandedness::Forward
        } else if percent_antisense > 67.0 {
            info!("Chemistry strandedness is set to: Reverse");
            ChemistryStrandedness::Reverse
        } else {
            info!("Chemistry strandedness is set to: Unstranded");
            ChemistryStrandedness::Unstranded
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
