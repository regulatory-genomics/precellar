use std::collections::HashMap;

use anyhow::Result;
use bitcode::{Decode, Encode};
use noodles_sam as sam;
use noodles_sam::alignment::{record::cigar::op::Kind, record_buf::RecordBuf};

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

/// Tolerances applied when matching read structure against a transcript model.
#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub struct TxAlignOpts {
    /// Annotated introns at or below this length are treated as alignment
    /// artifacts rather than retained introns: aligners routinely fail to
    /// report short skips, so a read covering one looks like intron retention
    /// when it is not. Zero-length introns (abutting exons) are always
    /// exempted, since they contain no intronic sequence at all. IsoQuant
    /// uses 50.
    pub micro_intron_length: u64,
    /// Tolerance for matching observed junctions against annotated exon
    /// boundaries. Long-read indel noise routinely shifts a junction by a
    /// few bases, which must not be read as a structural mismatch. Only the
    /// inner edges of the alignment are affected; transcript ends are
    /// governed by the unconditional exemption above. IsoQuant uses 4 for
    /// PacBio and 6 for ONT; `0` requires exact equality.
    pub junction_tolerance: u64,
}

impl Default for TxAlignOpts {
    fn default() -> Self {
        Self {
            micro_intron_length: 50,
            junction_tolerance: 0,
        }
    }
}

impl TxAlignOpts {
    /// Options that disable every tolerance, restoring exact structural
    /// matching. Used to reproduce historical classifications. Zero-length
    /// introns remain exempted, as they hold no intronic sequence.
    pub fn exact() -> Self {
        Self {
            micro_intron_length: 0,
            junction_tolerance: 0,
        }
    }

    /// Whether an observed junction coordinate matches an annotated exon
    /// boundary within tolerance.
    fn junction_matches(&self, observed: u64, annotated: u64) -> bool {
        observed.abs_diff(annotated) <= self.junction_tolerance
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

    /// A record with the given segments, for tests that exercise logic keyed on
    /// the transcript alignment list rather than on read coordinates.
    #[cfg(test)]
    pub(crate) fn for_test(chrom: &str, segments: Vec<(u64, u64)>) -> Self {
        Self {
            chrom: chrom.to_string(),
            segments: segments
                .into_iter()
                .map(|(start, end)| Segment { start, end })
                .collect(),
            strand: bed_utils::bed::Strand::Forward,
        }
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
        opts: &TxAlignOpts,
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
            align_type: self.align_helper(transcript, opts),
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
    ///
    /// Outer transcript edges are exempt from the containment checks: TSS/polyA
    /// annotations (also long-read end trimming) are unreliable, so a read overhanging
    /// a terminal exon says nothing about splicing. Inner edges (junctions) are still
    /// checked, within `opts.junction_tolerance`. The exemption is keyed on the
    /// transcript's terminal exons, not the first overlapped exon, so that a read
    /// starting in intron 2 still counts as intron retention rather than an overhang.
    fn align_helper(&self, transcript: &Transcript, opts: &TxAlignOpts) -> AlignType {
        if self.end() <= transcript.start || self.start() >= transcript.end {
            panic!("Read does not overlap transcript");
        }

        let tx_exons = transcript.exons();
        let exons = find_exons(tx_exons, self.start(), self.end() - 1);
        if exons.is_empty() {
            // No exon overlaps the read, but the read overlaps the transcript.
            return AlignType::Intronic;
        }

        let segments = &self.segments;
        let n_ex = exons.len();
        let n_segment = segments.len();

        // The left outer edge is exempt only at the transcript's start, and the
        // right outer edge only at the transcript's end. Exon boundaries are
        // unique within a transcript (exons are non-overlapping and in order),
        // so coordinate equality identifies the terminal exons.
        let at_tx_start = exons[0].start == tx_exons[0].start;
        let at_tx_end = exons[n_ex - 1].end == tx_exons[tx_exons.len() - 1].end;
        if n_ex == n_segment {
            if n_ex == 1 {
                if (at_tx_start || segments[0].start >= exons[0].start)
                    && (at_tx_end || segments[0].end <= exons[0].end)
                {
                    // Fully contained within the exon (terminal edges exempt).
                    AlignType::Exonic
                } else {
                    AlignType::Spanning
                }
            } else if segments[1..n_segment - 1]
                .iter()
                .zip(exons[1..n_ex - 1].iter())
                .all(|(s, e)| {
                    opts.junction_matches(s.start, e.start) && opts.junction_matches(s.end, e.end)
                })
            {
                // All internal segments (except the 1st and last) match within
                // tolerance.
                if opts.junction_matches(segments[0].end, exons[0].end)
                    && opts.junction_matches(segments[n_segment - 1].start, exons[n_ex - 1].start)
                {
                    // first and last segments also match.
                    if (at_tx_start || segments[0].start >= exons[0].start)
                        && (at_tx_end || segments[n_segment - 1].end <= exons[n_ex - 1].end)
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
            // There are more exons than segments. This is only compatible when
            // every observed splice junction lands on adjacent exon boundaries.
            let is_compatible = segments.windows(2).all(|junction| {
                let left_segment = &junction[0];
                let right_segment = &junction[1];

                let left_exons = find_exons(exons, left_segment.start, left_segment.end - 1);
                let right_exons = find_exons(exons, right_segment.start, right_segment.end - 1);
                let (Some(left_exon), Some(right_exon)) = (left_exons.last(), right_exons.first())
                else {
                    return false;
                };

                let Some(left_idx) = exons.iter().position(|exon| exon == left_exon) else {
                    return false;
                };
                let Some(right_idx) = exons.iter().position(|exon| exon == right_exon) else {
                    return false;
                };

                opts.junction_matches(left_segment.end, left_exon.end)
                    && opts.junction_matches(right_segment.start, right_exon.start)
                    && right_idx == left_idx + 1
            });

            if !is_compatible {
                AlignType::Discordant
            } else {
                // Every observed junction lands on adjacent exon boundaries, so
                // the read differs from a plain exonic match only by the introns
                // bridged within single segments. An intron short enough that
                // aligners routinely fail to report it as a skip is an
                // alignment artifact, not retained-intron evidence; if the read
                // bridges no genuine intron, it is exonic. Zero-length introns
                // (abutting exons) never count: they hold no intronic sequence.
                let only_artifacts = exons
                    .windows(2)
                    .filter(|pair| {
                        let (left, right) = (&pair[0], &pair[1]);
                        // The intron is bridged by a single segment that aligns
                        // straight through it without reporting a junction. A
                        // segment that only reaches partway into the intron is
                        // genuine retained-intron sequence, not a missed skip.
                        segments
                            .iter()
                            .any(|segment| segment.start <= left.end && segment.end >= right.start)
                    })
                    .all(|pair| pair[1].start - pair[0].end <= opts.micro_intron_length);
                if only_artifacts {
                    AlignType::Exonic
                } else {
                    AlignType::Spanning
                }
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
            TxAlignment::PeAligned { read1, read2, .. } => {
                Box::new(read1.to_fragments().chain(read2.to_fragments()))
            }
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
    opts: TxAlignOpts,
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
            opts: TxAlignOpts::default(),
        }
    }

    /// Applies custom tolerance options to the annotator.
    pub fn with_opts(mut self, opts: TxAlignOpts) -> Self {
        self.opts = opts;
        self
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
                self._align_one(&rec2, strandedness)
                    .into_iter()
                    .for_each(|a2| {
                        if let Some(a1) = aln1.get(&a2.transcript_id) {
                            if a1.strand == a2.strand {
                                if a1.strand == Orientation::Antisense {
                                    has_antisense = true;
                                } else {
                                    let ty = match (&a1.align_type, a2.align_type) {
                                        (AlignType::Intronic, AlignType::Intronic) => {
                                            AlignType::Intronic
                                        }
                                        (AlignType::Exonic, AlignType::Exonic) => AlignType::Exonic,
                                        (AlignType::Discordant, _) | (_, AlignType::Discordant) => {
                                            AlignType::Discordant
                                        }
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
            })
            .collect();

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
            // Now we know there's exactly one gene, we choose the first non-discordant
            // alignment. `gene_id` is set by any record pair that produced an
            // alignment, which need not be the first pair, so the first pair may
            // carry none.
            let (read1, read2, alignment) = aln
                .into_iter()
                .find(|(_, _, alignment)| !alignment.is_empty())
                .unwrap();
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
            .map(|(_, transcript)| rec.align_transcript(transcript, strandedness, &self.opts))
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
                        Some(Orientation::Sense) => is_sense = true, // R1 orientation is identical to original mRNA
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

        if total == 0 {
            warn!("No alignments found for strandedness detection, defaulting to Unstranded");
            return ChemistryStrandedness::Unstranded;
        }

        let percent_sense = num_sense as f64 / total as f64 * 100.0;
        let percent_antisense = num_antisense as f64 / total as f64 * 100.0;
        info!(
            "Found {:.3}% sense and {:.3}% antisense alignments (n={})",
            percent_sense, percent_antisense, total
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
    use crate::align::MultiMap;
    use crate::transcriptome::Exons;
    use noodles_core::Position;
    use noodles_sam::alignment::record::cigar::Op;
    use noodles_sam::alignment::record::data::field::tag::Tag;
    use noodles_sam::alignment::record::Flags;
    use noodles_sam::alignment::record_buf::data::field::value::Value;
    use noodles_sam::header::record::value::{map::ReferenceSequence, Map};
    use std::num::NonZeroUsize;

    use super::*;

    #[test]
    fn test_align_pe() {
        let header = sam::Header::builder()
            .add_reference_sequence(
                "chr1",
                Map::<ReferenceSequence>::new(NonZeroUsize::new(10_000).unwrap()),
            )
            .build();

        let transcript = Transcript {
            id: "tx1".to_string(),
            chrom: "chr1".to_string(),
            start: 100,
            end: 600,
            strand: bed_utils::bed::Strand::Forward,
            gene_id: "gene1".to_string(),
            gene_name: "GeneName1".to_string(),
            exons: Exons(vec![
                Exon {
                    start: 100,
                    end: 200,
                },
                Exon {
                    start: 300,
                    end: 400,
                },
                Exon {
                    start: 500,
                    end: 600,
                },
            ]),
        };

        let record = |start: u64, cigar: Vec<(Kind, usize)>, is_read2: bool| {
            let mut rec = RecordBuf::default();
            *rec.reference_sequence_id_mut() = Some(0);
            *rec.alignment_start_mut() = Position::new(start as usize + 1);
            *rec.flags_mut() = if is_read2 {
                Flags::SEGMENTED | Flags::REVERSE_COMPLEMENTED
            } else {
                Flags::SEGMENTED
            };
            *rec.cigar_mut() = cigar
                .into_iter()
                .map(|(kind, len)| Op::new(kind, len))
                .collect();
            rec.data_mut()
                .insert(Tag::CELL_BARCODE_ID, Value::String("ACGT".into()));
            rec
        };

        let discordant = |is_read2: bool| {
            record(
                120,
                vec![(Kind::Match, 60), (Kind::Skip, 140), (Kind::Match, 60)],
                is_read2,
            )
        };
        let intronic = |is_read2: bool| record(220, vec![(Kind::Match, 30)], is_read2);

        let read1 = MultiMap::new(discordant(false), Some(vec![intronic(false)]));
        let read2 = MultiMap::new(discordant(true), Some(vec![intronic(true)]));

        let aligner = TxAligner::new(
            vec![transcript],
            header,
            Some(ChemistryStrandedness::Forward),
        );
        let result = aligner
            ._align_pe(&read1, &read2, ChemistryStrandedness::Forward)
            .expect("the secondary pair aligns, so the read is not dropped");

        let TxAlignment::PeAligned {
            read1: reported1,
            read2: reported2,
            alignment,
            gene_id,
            ..
        } = &result
        else {
            panic!("expected PeAligned, got {result:?}");
        };

        assert_eq!(gene_id, "gene1");
        assert_eq!(
            alignment.as_slice(),
            &[("tx1".to_string(), AlignType::Intronic)]
        );
        assert_eq!(reported1.start(), 220);
        assert_eq!(reported2.start(), 220);

        assert!(!result.is_exonic_only());
        assert!(result.is_intronic_only());
        assert!(!result.is_spanning());
    }

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

        let mut test_fun = |segments: Vec<(u64, u64)>, opts: &TxAlignOpts, expected: AlignType| {
            record.segments = segments
                .into_iter()
                .map(|(s, e)| Segment { start: s, end: e })
                .collect();
            assert_eq!(record.align_helper(&transcript, opts), expected);
        };
        let exact = TxAlignOpts::exact();
        let default = TxAlignOpts::default();

        // Incompatible
        test_fun(vec![(12, 18), (32, 38)], &exact, AlignType::Discordant);
        test_fun(vec![(12, 18), (32, 38)], &default, AlignType::Discordant);
        test_fun(vec![(21, 29), (32, 38)], &exact, AlignType::Discordant);
        test_fun(vec![(21, 29), (32, 38)], &default, AlignType::Discordant);
        test_fun(
            vec![(5, 40), (50, 60), (500, 700), (1400, 1501)],
            &exact,
            AlignType::Discordant,
        );
        test_fun(
            vec![(5, 40), (50, 60), (500, 700), (1400, 1501)],
            &default,
            AlignType::Discordant,
        );

        // Exonic
        test_fun(vec![(12, 20), (30, 38)], &exact, AlignType::Exonic);
        test_fun(vec![(12, 20), (30, 38)], &default, AlignType::Exonic);

        // Intronic
        test_fun(vec![(6, 8)], &exact, AlignType::Intronic);
        test_fun(vec![(6, 8)], &default, AlignType::Intronic);

        // Reaching left of a middle exon is genuine intron retention.
        test_fun(vec![(29, 38)], &exact, AlignType::Spanning);
        test_fun(vec![(29, 38)], &default, AlignType::Spanning);

        // The first segment bridges a 10bp intron (exons 10-20 and 30-40):
        // an aligner artifact under the default micro-intron exemption,
        // retained-intron evidence under exact matching.
        test_fun(
            vec![(5, 40), (50, 60), (500, 700), (1500, 1501)],
            &exact,
            AlignType::Spanning,
        );
        test_fun(
            vec![(5, 40), (50, 60), (500, 700), (1500, 1501)],
            &default,
            AlignType::Exonic,
        );
    }

    /// The outer edge of the transcript's first and last exon is TSS/polyA
    /// territory, not splicing evidence, and is exempt unconditionally. The
    /// same edges on interior exons face introns and stay checked.
    #[test]
    fn test_terminal_exemption() {
        let transcript = Transcript {
            id: "tx1".to_string(),
            chrom: "chr1".to_string(),
            start: 100,
            end: 1100,
            strand: bed_utils::bed::Strand::Forward,
            gene_id: "gene1".to_string(),
            gene_name: "GeneName1".to_string(),
            exons: Exons(vec![
                Exon {
                    start: 100,
                    end: 200,
                },
                Exon {
                    start: 500,
                    end: 600,
                },
                Exon {
                    start: 1000,
                    end: 1100,
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
            // The exemption is unconditional: it holds under exact matching too.
            assert_eq!(
                record.align_helper(&transcript, &TxAlignOpts::exact()),
                expected
            );
            assert_eq!(
                record.align_helper(&transcript, &TxAlignOpts::default()),
                expected
            );
        };

        // Reaching past the transcript start / end.
        test_fun(vec![(80, 150)], AlignType::Exonic);
        test_fun(vec![(1050, 1150)], AlignType::Exonic);
        test_fun(vec![(80, 200), (500, 600), (1000, 1150)], AlignType::Exonic);

        // Extending into intron 1 from inside the first exon: the right edge
        // of the first exon is a junction, checked exactly.
        test_fun(vec![(120, 300)], AlignType::Spanning);
        // Read starts in intron 1 and reaches left of exon 2: exon 2 is not
        // the transcript's first exon, so this is genuine intron retention.
        test_fun(vec![(450, 590)], AlignType::Spanning);
        // Multi-segment read reaching left of a middle exon.
        test_fun(vec![(450, 600), (1000, 1100)], AlignType::Spanning);

        // A single-exon transcript: both edges of its only exon are terminal.
        let single = Transcript {
            id: "tx2".to_string(),
            chrom: "chr1".to_string(),
            start: 100,
            end: 200,
            strand: bed_utils::bed::Strand::Forward,
            gene_id: "gene1".to_string(),
            gene_name: "GeneName1".to_string(),
            exons: Exons(vec![Exon {
                start: 100,
                end: 200,
            }]),
        };
        record.segments = vec![Segment {
            start: 80,
            end: 220,
        }];
        assert_eq!(
            record.align_helper(&single, &TxAlignOpts::exact()),
            AlignType::Exonic
        );
        assert_eq!(
            record.align_helper(&single, &TxAlignOpts::default()),
            AlignType::Exonic
        );
    }

    /// An intron short enough that aligners routinely miss the skip is an
    /// alignment artifact, not retained-intron evidence. A segment that only
    /// reaches partway into an intron holds genuine intronic sequence and is
    /// never exempt.
    #[test]
    fn test_micro_intron() {
        // 10bp intron.
        let micro = Transcript {
            id: "tx1".to_string(),
            chrom: "chr1".to_string(),
            start: 10,
            end: 40,
            strand: bed_utils::bed::Strand::Forward,
            gene_id: "gene1".to_string(),
            gene_name: "GeneName1".to_string(),
            exons: Exons(vec![
                Exon { start: 10, end: 20 },
                Exon { start: 30, end: 40 },
            ]),
        };
        // 40bp and 40bp introns: each within the default threshold, but their
        // sum (80bp) exceeds it, so the length check must be per-gap.
        let micro_chain = Transcript {
            id: "tx2".to_string(),
            chrom: "chr1".to_string(),
            start: 10,
            end: 120,
            strand: bed_utils::bed::Strand::Forward,
            gene_id: "gene1".to_string(),
            gene_name: "GeneName1".to_string(),
            exons: Exons(vec![
                Exon { start: 10, end: 20 },
                Exon { start: 60, end: 70 },
                Exon {
                    start: 110,
                    end: 120,
                },
            ]),
        };
        // 50bp intron, exactly at the default threshold.
        let boundary = Transcript {
            id: "tx3".to_string(),
            chrom: "chr1".to_string(),
            start: 10,
            end: 80,
            strand: bed_utils::bed::Strand::Forward,
            gene_id: "gene1".to_string(),
            gene_name: "GeneName1".to_string(),
            exons: Exons(vec![
                Exon { start: 10, end: 20 },
                Exon { start: 70, end: 80 },
            ]),
        };
        // 200bp intron: genuine retention evidence.
        let real = Transcript {
            id: "tx4".to_string(),
            chrom: "chr1".to_string(),
            start: 10,
            end: 230,
            strand: bed_utils::bed::Strand::Forward,
            gene_id: "gene1".to_string(),
            gene_name: "GeneName1".to_string(),
            exons: Exons(vec![
                Exon { start: 10, end: 20 },
                Exon {
                    start: 220,
                    end: 230,
                },
            ]),
        };
        // Zero-length intron (abutting exons): holds no intronic sequence.
        let abutting = Transcript {
            id: "tx5".to_string(),
            chrom: "chr1".to_string(),
            start: 10,
            end: 30,
            strand: bed_utils::bed::Strand::Forward,
            gene_id: "gene1".to_string(),
            gene_name: "GeneName1".to_string(),
            exons: Exons(vec![
                Exon { start: 10, end: 20 },
                Exon { start: 20, end: 30 },
            ]),
        };

        let mut record = SplicedRecord {
            chrom: "chr1".to_string(),
            segments: vec![],
            strand: bed_utils::bed::Strand::Forward,
        };
        let mut test_fun = |transcript: &Transcript,
                            segments: Vec<(u64, u64)>,
                            opts: &TxAlignOpts,
                            expected: AlignType| {
            record.segments = segments
                .into_iter()
                .map(|(s, e)| Segment { start: s, end: e })
                .collect();
            assert_eq!(record.align_helper(transcript, opts), expected);
        };
        let exact = TxAlignOpts::exact();
        let default = TxAlignOpts::default();

        // Single segment through a 10bp intron: artifact at default, retained
        // at exact.
        test_fun(&micro, vec![(12, 38)], &default, AlignType::Exonic);
        test_fun(&micro, vec![(12, 38)], &exact, AlignType::Spanning);
        // Multi-segment read whose first segment bridges a 40bp micro-intron;
        // the observed junction between exons 2 and 3 is not bridging evidence.
        test_fun(
            &micro_chain,
            vec![(12, 70), (110, 115)],
            &default,
            AlignType::Exonic,
        );
        test_fun(
            &micro_chain,
            vec![(12, 70), (110, 115)],
            &exact,
            AlignType::Spanning,
        );
        // 50bp intron sits exactly at the default threshold (`<=`).
        test_fun(&boundary, vec![(12, 75)], &default, AlignType::Exonic);
        test_fun(&boundary, vec![(12, 75)], &exact, AlignType::Spanning);
        // 200bp intron: never exempt.
        test_fun(&real, vec![(12, 225)], &default, AlignType::Spanning);
        test_fun(&real, vec![(12, 225)], &exact, AlignType::Spanning);
        // Zero-length intron: exempt even under exact matching.
        test_fun(&abutting, vec![(12, 28)], &default, AlignType::Exonic);
        test_fun(&abutting, vec![(12, 28)], &exact, AlignType::Exonic);
        // Reaching partway into a 10bp intron is genuine intronic sequence.
        test_fun(&micro, vec![(12, 25)], &default, AlignType::Spanning);
        test_fun(&micro, vec![(12, 25)], &exact, AlignType::Spanning);
    }

    /// Splice sites carry a GT-AG motif, so a junction that misses the
    /// annotated boundary is evidence against the isoform -- but long-read
    /// indel noise shifts junctions by a few bases, and treating that as a
    /// structural mismatch loses the read entirely (`Discordant` is filtered
    /// out before counting). The tolerance applies only to inner edges.
    #[test]
    fn test_junction_tolerance() {
        let transcript = Transcript {
            id: "tx1".to_string(),
            chrom: "chr1".to_string(),
            start: 100,
            end: 1100,
            strand: bed_utils::bed::Strand::Forward,
            gene_id: "gene1".to_string(),
            gene_name: "GeneName1".to_string(),
            exons: Exons(vec![
                Exon {
                    start: 100,
                    end: 200,
                },
                Exon {
                    start: 500,
                    end: 600,
                },
                Exon {
                    start: 1000,
                    end: 1100,
                },
            ]),
        };

        let mut record = SplicedRecord {
            chrom: "chr1".to_string(),
            segments: vec![],
            strand: bed_utils::bed::Strand::Forward,
        };
        let mut test_fun = |segments: Vec<(u64, u64)>, tolerance: u64, expected: AlignType| {
            record.segments = segments
                .into_iter()
                .map(|(s, e)| Segment { start: s, end: e })
                .collect();
            let opts = TxAlignOpts {
                junction_tolerance: tolerance,
                ..TxAlignOpts::exact()
            };
            assert_eq!(
                record.align_helper(&transcript, &opts),
                expected,
                "segments={:?} tolerance={}",
                record.segments,
                tolerance
            );
        };

        // An exact three-exon match is unaffected by tolerance.
        for tol in [0, 4, 6] {
            test_fun(
                vec![(120, 200), (500, 600), (1000, 1080)],
                tol,
                AlignType::Exonic,
            );
        }

        // Inner edge of the first segment off by 3bp: recovered at 4 and 6.
        test_fun(
            vec![(120, 197), (500, 600), (1000, 1080)],
            0,
            AlignType::Discordant,
        );
        test_fun(
            vec![(120, 197), (500, 600), (1000, 1080)],
            4,
            AlignType::Exonic,
        );

        // Off by 5bp: still discordant at 4, recovered at 6 (PacBio vs ONT).
        test_fun(
            vec![(120, 195), (500, 600), (1000, 1080)],
            4,
            AlignType::Discordant,
        );
        test_fun(
            vec![(120, 195), (500, 600), (1000, 1080)],
            6,
            AlignType::Exonic,
        );

        // Off by 7bp exceeds every supported tolerance.
        for tol in [0, 4, 6] {
            test_fun(
                vec![(120, 193), (500, 600), (1000, 1080)],
                tol,
                AlignType::Discordant,
            );
        }

        // Both edges of the internal segment off by 3bp.
        test_fun(
            vec![(120, 200), (503, 597), (1000, 1080)],
            0,
            AlignType::Discordant,
        );
        test_fun(
            vec![(120, 200), (503, 597), (1000, 1080)],
            4,
            AlignType::Exonic,
        );

        // Inner edge of the last segment off by 3bp.
        test_fun(
            vec![(120, 200), (500, 600), (1003, 1080)],
            0,
            AlignType::Discordant,
        );
        test_fun(
            vec![(120, 200), (500, 600), (1003, 1080)],
            4,
            AlignType::Exonic,
        );

        // Tolerance must not rescue a junction that no boundary explains: 300
        // is 100bp past the end of exon 1.
        for tol in [0, 4, 6] {
            test_fun(vec![(120, 300), (500, 600)], tol, AlignType::Discordant);
        }

        // Outer edges are governed by the terminal exemption, not tolerance:
        // overhanging the transcript start by 20bp stays exonic at tolerance 0.
        test_fun(
            vec![(80, 200), (500, 600), (1000, 1120)],
            0,
            AlignType::Exonic,
        );
    }
}
