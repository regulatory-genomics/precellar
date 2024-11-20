use crate::transcript::{Gene, SpliceSegments, Transcript};

use anyhow::{ensure, Result};
use bed_utils::bed::map::GIntervalMap;
use bed_utils::bed::GenomicRange;
use noodles::gtf::record::strand::Strand;
use noodles::sam;
use noodles::sam::alignment::record::Flags;
use noodles::sam::alignment::record_buf::Cigar;
use noodles::sam::alignment::record_buf::RecordBuf;
use std::cmp;
use std::collections::{BTreeMap, HashSet};

/// BAM record with annotation
pub enum AnnotatedAlignment {
    Unmapped(RecordBuf, Option<RecordBuf>),
    SeMapped(RecordBuf, Annotation),
    PeMapped(
        RecordBuf,
        Annotation,
        RecordBuf,
        Annotation,
        PairAnnotationData,
    ),
}

impl AnnotatedAlignment {
    pub fn rec(&self) -> (&RecordBuf, Option<&RecordBuf>) {
        match self {
            AnnotatedAlignment::Unmapped(ref rec, ref rec2) => (rec, rec2.as_ref()),
            AnnotatedAlignment::SeMapped(ref rec, _) => (rec, None),
            AnnotatedAlignment::PeMapped(ref rec1, _, ref rec2, _, _) => (rec1, Some(rec2)),
        }
    }

    pub fn mut_rec(&mut self) -> (&mut RecordBuf, Option<&mut RecordBuf>) {
        match self {
            AnnotatedAlignment::Unmapped(ref mut rec1, ref mut rec2) => (rec1, rec2.as_mut()),
            AnnotatedAlignment::SeMapped(ref mut rec, _) => (rec, None),
            AnnotatedAlignment::PeMapped(ref mut rec1, _, ref mut rec2, _, _) => (rec1, Some(rec2)),
        }
    }

    pub fn annotation(&self) -> (Option<&Annotation>, Option<&Annotation>) {
        match self {
            AnnotatedAlignment::SeMapped(_, ref anno) => (Some(anno), None),
            AnnotatedAlignment::PeMapped(_, ref anno1, _, ref anno2, _) => {
                (Some(anno1), Some(anno2))
            }
            AnnotatedAlignment::Unmapped(_, _) => (None, None),
        }
    }

    fn set_rescued(&mut self) {
        match self {
            AnnotatedAlignment::SeMapped(_, ref mut anno) => anno.rescued = true,
            AnnotatedAlignment::PeMapped(_, ref mut anno1, _, ref mut anno2, _) => {
                anno1.rescued = true;
                anno2.rescued = true;
            }
            AnnotatedAlignment::Unmapped(_, _) => panic!("Unmapped annotations cannot be rescued"),
        }
    }
}

#[derive(Debug, Clone)]
pub struct AnnotationParams {
    pub chemistry_strandedness: Strand,
    pub chemistry_endedness: WhichEnd,
    pub intergenic_trim_bases: i64,
    pub intronic_trim_bases: i64,
    pub junction_trim_bases: i64,
    pub region_min_overlap: f64,
}

#[derive(Debug)]
pub struct AlignmentAnnotator {
    params: AnnotationParams,
    transcripts: GIntervalMap<Transcript>,
}

impl AlignmentAnnotator {
    /// Annotate the alignments by mapping them to the transcriptome. If multiple
    /// alignments are present, we will try to find the confident ones and promote
    /// them to primary. A read may align to multiple transcripts and genes, but
    /// it is only considered confidently mapped to the transcriptome if it is
    /// mapped to a single gene.
    pub fn annotate_alignments_se(
        &self,
        header: &sam::Header,
        rec: Vec<RecordBuf>,
    ) -> Vec<AnnotatedAlignment> {
        let mut results = rec
            .into_iter()
            .map(|rec| self.annotate_alignment_se(header, rec))
            .collect::<Vec<_>>();
        rescue_alignments_se(&mut results);
        results
    }

    pub fn annotate_alignments_pe(
        &self,
        header: &sam::Header,
        rec1: Vec<RecordBuf>,
        rec2: Vec<RecordBuf>,
    ) -> Vec<AnnotatedAlignment> {
        assert!(
            !rec2.is_empty(),
            "There should be at least one paired-end alignment"
        );

        let pair_improper = rec1.len() != rec2.len();
        let mut result: Vec<_> = rec1
            .into_iter()
            .zip(rec2)
            .map(|(r1, r2)| self.annotate_alignment_pe(header, r1, r2))
            .collect();
        rescue_alignments_pe(result.as_mut_slice());
        result
    }

    /// Create annotation for a single alignment
    fn annotate_alignment_se(&self, header: &sam::Header, rec: RecordBuf) -> AnnotatedAlignment {
        if rec.flags().is_unmapped() {
            AnnotatedAlignment::Unmapped(rec, None)
        } else {
            let anno = self.annotate_alignment(header, &rec).unwrap();
            AnnotatedAlignment::SeMapped(rec, anno)
        }
    }

    /// Create annotation for a pair of alignments
    fn annotate_alignment_pe(
        &self,
        header: &sam::Header,
        rec1: RecordBuf,
        rec2: RecordBuf,
    ) -> AnnotatedAlignment {
        // STAR _shouldn't_ return pairs where only a single end is mapped,
        //   but if it does, consider the pair unmapped
        if rec1.flags().is_unmapped() || rec2.flags().is_unmapped() {
            AnnotatedAlignment::Unmapped(rec1, Some(rec2))
        } else {
            let anno1 = self.annotate_alignment(header, &rec1).unwrap();
            let anno2 = self.annotate_alignment(header, &rec2).unwrap();
            let annop = PairAnnotationData::from_pair(&anno1, &anno2);
            AnnotatedAlignment::PeMapped(rec1, anno1, rec2, anno2, annop)
        }
    }

    fn annotate_alignment(&self, header: &sam::Header, read: &RecordBuf) -> Result<Annotation> {
        ensure!(
            !read.flags().is_unmapped(),
            "Unmapped alignments cannot be annotated"
        );
        let chrom = read.reference_sequence(header).unwrap()?.0;
        let chrom = std::str::from_utf8(chrom)?.to_string();

        let region = GenomicRange::new(
            &chrom,
            read.alignment_start().unwrap().get() as u64,
            read.alignment_end().unwrap().get() as u64,
        );
        let alignments = self
            .transcripts
            .find(&region)
            .flat_map(|(_, transcript)| self.align_to_transcript(read, transcript))
            .collect::<Vec<_>>();

        let mut seen_genes = HashSet::new();
        let mut transcripts = BTreeMap::new();
        let mut antisense = BTreeMap::new();
        let mut annotation_data = Annotation::new(chrom);
        if alignments.is_empty() {
            annotation_data.region = AnnotationRegion::Intergenic;
        } else if alignments.iter().any(|x| x.is_exonic()) {
            annotation_data.region = AnnotationRegion::Exonic;
            // Check if there are transcriptome compatible alignments
            alignments.into_iter().rev().for_each(|aln| {
                if let Some(tx_align) = &aln.exon_align {
                    match aln.strand {
                        Strand::Forward => {
                            // Transcript sense alignment
                            seen_genes.insert(aln.gene.clone());
                            transcripts.insert(tx_align.id.clone(), aln);
                        }
                        Strand::Reverse => {
                            // Transcript anti-sense alignment
                            antisense.insert(tx_align.id.clone(), aln);
                        }
                    }
                }
            });
        } else {
            annotation_data.region = AnnotationRegion::Intronic;
            alignments
                .into_iter()
                .rev()
                .for_each(|aln| match aln.strand {
                    Strand::Forward => {
                        seen_genes.insert(aln.gene.clone());
                        transcripts.insert(aln.gene.id.clone(), aln);
                    }
                    Strand::Reverse => {
                        antisense.insert(aln.gene.id.clone(), aln);
                    }
                });
        }

        annotation_data.aln_sense = transcripts.into_values().collect::<Vec<_>>();
        annotation_data.aln_antisense = antisense.into_values().collect::<Vec<_>>();
        annotation_data.genes = seen_genes.into_iter().collect::<Vec<Gene>>();
        // Sorting this makes life easier later.
        annotation_data.genes.sort_unstable();

        Ok(annotation_data)
    }

    fn align_to_transcript(
        &self,
        read: &RecordBuf,
        transcript: &Transcript,
    ) -> Option<TranscriptAlignment> {
        // figure out coordinates
        let tx_start = transcript.start;
        let tx_end = transcript.end;
        let genomic_start = read.alignment_start().unwrap().get() as i64;
        let genomic_end = read.alignment_end().unwrap().get() as i64;

        // compute region
        let splice_segments = SpliceSegments::from(read);
        let is_exonic = splice_segments.is_exonic(transcript, self.params.region_min_overlap);
        let is_intronic =
            !is_exonic && get_overlap(genomic_start, genomic_end, tx_start, tx_end) >= 1.0;
        if !(is_exonic || is_intronic) {
            // If neither exonic nor intronic => intergenic.
            return None;
        }
        // compute strand
        let tx_reverse_strand = transcript.strand == Strand::Reverse;
        let flags = read.flags();
        let mut read_reverse_strand = flags.is_reverse_complemented();
        if flags.is_segmented() && flags.is_last_segment() {
            read_reverse_strand = !read_reverse_strand;
        };
        let is_antisense = match self.params.chemistry_strandedness {
            Strand::Forward => tx_reverse_strand != read_reverse_strand,
            Strand::Reverse => tx_reverse_strand == read_reverse_strand,
        };

        let tx_strand = if is_antisense {
            Strand::Reverse
        } else {
            Strand::Forward
        };

        // Prior to aligning to the transcript we have a gene alignment
        let gene = transcript.gene.clone();
        let gene_aln = TranscriptAlignment {
            gene: gene.clone(),
            strand: tx_strand,
            exon_align: None,
        };
        // if intronic return region and strand
        if is_intronic {
            return Some(gene_aln);
        };

        // must be exonic to continue
        assert!(is_exonic);

        // compute offsets
        let mut tx_offset = transcript.get_offset(genomic_start).unwrap();
        let tx_length = transcript.len();

        // align the read to the exons
        let Some((mut tx_cigar, tx_aligned_bases)) = splice_segments.align_junctions(
            transcript,
            self.params.junction_trim_bases,
            self.params.intergenic_trim_bases,
            self.params.intronic_trim_bases,
        ) else {
            return Some(gene_aln);
        };

        // flip reverse strand
        if tx_reverse_strand {
            tx_offset = tx_length - (tx_offset + tx_aligned_bases);
            tx_cigar.as_mut().reverse();
        };

        // Apparent insert size, assuming a SE read coming from this transcript
        let se_insert_size = match self.params.chemistry_endedness {
            WhichEnd::FivePrime => tx_offset + tx_aligned_bases,
            WhichEnd::ThreePrime => tx_length - tx_offset,
        } as u32;

        let alignment = TranscriptAlignment {
            gene,
            strand: tx_strand,
            exon_align: Some(TxAlignProperties {
                id: transcript.id.clone(),
                pos: tx_offset,
                cigar: tx_cigar,
                alen: tx_aligned_bases,
                se_insert_size,
            }),
        };

        Some(alignment)
    }
}

#[derive(Eq, PartialEq, Debug)]
pub struct PairAnnotationData {
    pub genes: Vec<Gene>,
}

impl PairAnnotationData {
    /// Annotate a pair of alignments
    /// Take the intersection of the non-empty gene sets of the mates
    pub fn from_pair(anno1: &Annotation, anno2: &Annotation) -> PairAnnotationData {
        let genes = match (!anno1.genes.is_empty(), !anno2.genes.is_empty()) {
            (true, false) => anno1.genes.clone(),
            (false, true) => anno2.genes.clone(),
            _ if anno1.chrom == anno2.chrom => anno1
                .genes
                .iter()
                .collect::<HashSet<_>>()
                .intersection(&anno2.genes.iter().collect::<HashSet<_>>())
                .map(|x| (*x).clone())
                .collect(),
            _ => vec![],
        };
        PairAnnotationData { genes }
    }
}

/// Use transcriptome alignments to promote a single genome alignment
/// when none are confidently mapped to the genome.
/// Returns true if rescue took place.
fn rescue_alignments_se(recs: &mut [AnnotatedAlignment]) -> bool {
    let mut rescued = false;

    // If there are multiple alignments, we may need to rescue
    if recs.len() > 1 {
        let mut promote_index: Option<usize> = None;
        let mut ori_primary: Option<usize> = None;
        let mut seen_genes = HashSet::new();

        for (i, rec) in recs.iter_mut().enumerate() {
            match rec {
                AnnotatedAlignment::SeMapped(aln, anno) => {
                    if !aln.flags().is_secondary() {
                        // Track the original primary alignment
                        ori_primary = Some(i);
                        // Make sure all mapped alignments are secondary
                        aln.flags_mut().set(Flags::SECONDARY, true);
                    }

                    // Only consider transcriptomic alignments for rescue
                    if anno.aln_sense.iter().any(|x| x.is_exonic()) {
                        let genes = &anno.genes;
                        // Track which record/record-pair we should promote;
                        // Take the first record/pair with 1 gene
                        if genes.len() == 1 {
                            promote_index = promote_index.or(Some(i));
                        }

                        // Track number of distinct genes we're aligned to
                        seen_genes.extend(genes);
                    }
                }
                _ => unimplemented!("Only single-end alignments can be rescued"),
            }
        }

        // The alignment can be rescued if there is only one uniquely mapped gene
        if seen_genes.len() == 1 && promote_index.is_some() {
            rescued = true
        } else {
            promote_index = ori_primary;
        }

        // Promote a single alignment
        let promote_index = promote_index.expect("No primary alignment found!");
        let rec = &mut recs[promote_index];
        rec.set_rescued();
        rec.mut_rec().0.flags_mut().set(Flags::SECONDARY, false);
    }

    rescued
}

/// Use transcriptome alignments to promote a single genome alignment
/// when none are confidently mapped to the genome.
/// Returns true if rescue took place.
fn rescue_alignments_pe(pairs: &mut [AnnotatedAlignment]) -> bool {
    let mut rescued = false;

    if pairs.len() > 1 {
        // Check if rescue is appropriate and determine which record to promote
        let mut seen_genes = HashSet::new();
        let mut promote_index: Option<usize> = None;
        let mut ori_primary: Option<usize> = None;

        for (i, pair) in pairs.iter_mut().enumerate() {
            match pair {
                AnnotatedAlignment::PeMapped(aln1, _, aln2, _, anno) => {
                    if !aln1.flags().is_secondary() {
                        ori_primary = Some(i);
                        // Make sure all mapped alignments are secondary
                        aln1.flags_mut().set(Flags::SECONDARY, true);
                        aln2.flags_mut().set(Flags::SECONDARY, true);
                    }

                    let genes = &anno.genes;

                    // Track which record/record-pair we should promote;
                    // Take the first record/pair with 1 gene
                    if genes.len() == 1 {
                        promote_index = promote_index.or(Some(i));
                    }

                    // Track number of distinct genes we're aligned to
                    seen_genes.extend(genes);
                }
                _ => unimplemented!("Only paired-end alignments can be rescued"),
            }
        }

        // The alignment can be rescued if there is only one uniquely mapped gene
        if seen_genes.len() == 1 && promote_index.is_some() {
            rescued = true
        } else {
            promote_index = ori_primary;
        }

        // Promote a single alignment
        let promote_index = promote_index.expect("No primary alignment found!");
        let pair = &mut pairs[promote_index];
        pair.set_rescued();
        let (aln1, aln2) = pair.mut_rec();
        aln1.flags_mut().set(Flags::SECONDARY, false);
        aln2.unwrap().flags_mut().set(Flags::SECONDARY, false);
    }

    rescued
}

/// Which end of a transcript reads come from
#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub enum WhichEnd {
    ThreePrime,
    FivePrime,
}

#[derive(Eq, PartialEq, Ord, PartialOrd, Debug, Clone, Copy, Hash)]
pub enum AnnotationRegion {
    Exonic,
    Intronic,
    Intergenic,
}

#[derive(Eq, PartialEq, Debug)]
pub struct Annotation {
    pub aln_sense: Vec<TranscriptAlignment>,
    pub aln_antisense: Vec<TranscriptAlignment>,
    pub genes: Vec<Gene>,
    pub region: AnnotationRegion,
    pub rescued: bool,
    pub chrom: String,
}

impl Annotation {
    fn new(genome: String) -> Annotation {
        Annotation {
            aln_sense: Vec::new(),
            aln_antisense: Vec::new(),
            genes: Vec::new(),
            region: AnnotationRegion::Intergenic,
            rescued: false,
            chrom: genome,
        }
    }
}

#[derive(Eq, PartialEq, Debug)]
pub struct TranscriptAlignment {
    pub gene: Gene,
    pub strand: Strand,
    pub exon_align: Option<TxAlignProperties>,
}

impl TranscriptAlignment {
    pub fn is_exonic(&self) -> bool {
        self.exon_align.is_some()
    }

    pub fn is_intronic(&self) -> bool {
        !self.is_exonic()
    }
}

// These quantities are well defined for a valid transcriptomic alignment
#[derive(Eq, PartialEq, Debug)]
pub struct TxAlignProperties {
    pub id: String,
    pub pos: i64,
    pub cigar: Cigar,
    pub alen: i64,
    pub se_insert_size: u32,
}

/// Fraction of read interval covered by ref interval
fn get_overlap(read_start: i64, read_end: i64, ref_start: i64, ref_end: i64) -> f64 {
    let mut overlap_bases = cmp::min(ref_end, read_end) - cmp::max(ref_start, read_start);
    if overlap_bases < 0 {
        overlap_bases = 0;
    }
    (overlap_bases as f64) / ((read_end - read_start) as f64)
}
