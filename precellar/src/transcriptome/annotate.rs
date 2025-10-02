/// This module provides utilities for annotating bam records by mapping them to a transcriptome.
/// It supports both single-end and paired-end alignments and uses transcript annotations for gene-level
/// and exon-level classification.
use crate::align::MultiMapR;
use crate::transcriptome::{SplicedRecord, Transcript};

use anyhow::{ensure, Result};
use bed_utils::bed::map::GIntervalMap;
use bed_utils::bed::GenomicRange;
use bincode::{Decode, Encode};
use noodles::sam;
use noodles::sam::alignment::record_buf::Cigar;
use noodles::sam::alignment::record_buf::RecordBuf;
use noodles::sam::record::data::field::value::base_modifications::group::Strand;
use std::cmp;
use std::collections::{BTreeMap, HashSet};

#[derive(Eq, PartialEq, Debug)]
pub struct Annotation {
    pub aln_sense: Vec<TranscriptAlignment>,
    pub aln_antisense: Vec<TranscriptAlignment>,
    pub genes: Vec<String>,
    pub region: AnnotationRegion,
    pub rescued: bool,
    pub chrom: String,
}

/// Represents an annotated BAM record.
#[derive(Debug)]
pub enum AnnotatedAlignment {
    /// Represents a single-end mapped read with annotation.
    SeMapped(Annotation),
    /// Represents a paired-end mapped read with annotations for each end.
    PeMapped(
        Annotation,
        Annotation,
        PairAnnotationData,
    ),
}

impl AnnotatedAlignment {
    /// Returns references to the annotations for single-end or paired-end reads.
    pub fn annotation(&self) -> (Option<&Annotation>, Option<&Annotation>) {
        match self {
            AnnotatedAlignment::SeMapped(ref anno) => (Some(anno), None),
            AnnotatedAlignment::PeMapped(ref anno1, ref anno2, _) => {
                (Some(anno1), Some(anno2))
            }
        }
    }

    /// Marks the alignment as rescued.
    fn set_rescued(&mut self) {
        match self {
            AnnotatedAlignment::SeMapped(ref mut anno) => anno.rescued = true,
            AnnotatedAlignment::PeMapped(ref mut anno1, ref mut anno2, _) => {
                anno1.rescued = true;
                anno2.rescued = true;
            }
        }
    }
}

/// Manages the annotation of alignments using transcriptome data.
#[derive(Debug, Clone)]
pub struct AlignmentAnnotator {
    /// Map of genomic intervals to transcripts.
    pub(crate) transcripts: GIntervalMap<Transcript>,
    /// Indicates the strandedness of the chemistry used in the experiment.
    chemistry_strandedness: Strand,
    /// Number of bases to trim from intergenic alignments.
    intergenic_trim_bases: u64,
    /// Number of bases to trim from intronic alignments.
    intronic_trim_bases: u64,
    /// Number of bases to trim from junction alignments.
    junction_trim_bases: u64,
    /// Minimum overlap fraction required for a region to be considered mapped.
    region_min_overlap: f64,
}

impl AlignmentAnnotator {
    /// Creates a new `AlignmentAnnotator` with the provided transcripts.
    pub fn new(transcripts: impl IntoIterator<Item = Transcript>) -> Self {
        let transcripts = transcripts
            .into_iter()
            .map(|x| (GenomicRange::new(&x.chrom, x.start, x.end), x))
            .collect();
        Self {
            transcripts,
            chemistry_strandedness: Strand::Forward,
            intergenic_trim_bases: 0,
            intronic_trim_bases: 0,
            junction_trim_bases: 0,
            region_min_overlap: 0.5,
        }
    }

    /// Annotate the alignments by mapping them to the transcriptome. If multiple
    /// alignments are present, we will try to find the confident ones and promote
    /// them to primary. A read may align to multiple transcripts and genes, but
    /// it is only considered confidently mapped to the transcriptome if it is
    /// mapped to a single gene. The confident alignment will be returned if found.
    ///
    /// # Arguments
    /// * `header` - Reference to the SAM header.
    /// * `rec` - Vector of single-end alignment records.
    pub fn annotate_alignments_se(
        &self,
        header: &sam::Header,
        rec: MultiMapR,
    ) -> Option<AnnotatedAlignment> {
        let results = rec
            .iter()
            .filter_map(|rec| self.annotate_alignment_se(header, rec))
            .collect::<Vec<_>>();
        rescue_alignments_se(results)
    }

    /// Annotates a batch of paired-end alignments.
    ///
    /// # Arguments
    /// * `header` - Reference to the SAM header.
    /// * `rec1` - Vector of first-end alignment records.
    /// * `rec2` - Vector of second-end alignment records.
    pub fn annotate_alignments_pe(
        &self,
        header: &sam::Header,
        rec1: MultiMapR,
        rec2: MultiMapR,
    ) -> Option<AnnotatedAlignment> {
        let _pair_improper = rec1.len() != rec2.len();
        let result: Vec<_> = rec1
            .iter()
            .zip(rec2.iter())
            .filter_map(|(r1, r2)| self.annotate_alignment_pe(header, r1, r2))
            .collect();
        rescue_alignments_pe(result)
    }

    /// Annotates a single-end alignment record.
    fn annotate_alignment_se(&self, header: &sam::Header, rec: &RecordBuf) -> Option<AnnotatedAlignment> {
        if rec.flags().is_unmapped() {
            None
        } else {
            let anno = self.annotate_alignment(header, rec).unwrap();
            Some(AnnotatedAlignment::SeMapped(anno))
        }
    }

    /// Create annotation for a pair of alignments
    fn annotate_alignment_pe(
        &self,
        header: &sam::Header,
        rec1: &RecordBuf,
        rec2: &RecordBuf,
    ) -> Option<AnnotatedAlignment> {
        // STAR _shouldn't_ return pairs where only a single end is mapped,
        //   but if it does, consider the pair unmapped
        if rec1.flags().is_unmapped() || rec2.flags().is_unmapped() {
            None
        } else {
            let anno1 = self.annotate_alignment(header, rec1).unwrap();
            let anno2 = self.annotate_alignment(header, rec2).unwrap();
            let annop = PairAnnotationData::from_pair(&anno1, &anno2);
            Some(AnnotatedAlignment::PeMapped(anno1, anno2, annop))
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
        let annotation_region;
        if alignments.is_empty() {
            annotation_region = AnnotationRegion::Intergenic;
        } else if alignments.iter().any(|x| x.is_exonic()) {
            annotation_region = AnnotationRegion::Exonic;
            // Check if there are transcriptome compatible alignments
            alignments.into_iter().rev().for_each(|aln| {
                if let Some(tx_align) = &aln.exon_align {
                    match aln.strand {
                        Strand::Forward => {
                            // Transcript sense alignment
                            seen_genes.insert(aln.gene_id.clone());
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
            annotation_region = AnnotationRegion::Intronic;
            alignments
                .into_iter()
                .rev()
                .for_each(|aln| match aln.strand {
                    Strand::Forward => {
                        seen_genes.insert(aln.gene_id.clone());
                        transcripts.insert(aln.gene_id.clone(), aln);
                    }
                    Strand::Reverse => {
                        antisense.insert(aln.gene_id.clone(), aln);
                    }
                });
        }

        let mut annotation = Annotation {
            aln_sense: transcripts.into_values().collect::<Vec<_>>(),
            aln_antisense: antisense.into_values().collect::<Vec<_>>(),
            genes: seen_genes.into_iter().collect::<Vec<_>>(),
            region: annotation_region,
            rescued: false,
            chrom,
        };
        // Sorting this makes life easier later.
        annotation.genes.sort_unstable();

        Ok(annotation)
    }

    /// Aligns a read to a transcript and determines the region type (exonic, intronic, or intergenic).
    fn align_to_transcript(
        &self,
        read: &RecordBuf,
        transcript: &Transcript,
    ) -> Option<TranscriptAlignment> {
        // figure out coordinates
        let tx_start = transcript.start;
        let tx_end = transcript.end;
        let genomic_start = read.alignment_start().unwrap().get() as u64;
        let genomic_end = read.alignment_end().unwrap().get() as u64;
        let splice_segments = SplicedRecord::from(read);

        let is_exonic = splice_segments.is_exonic(transcript, self.region_min_overlap);
        if is_exonic || get_overlap(genomic_start, genomic_end, tx_start, tx_end) >= 1.0 {
            // compute strand
            let tx_reverse_strand = transcript.strand == Strand::Reverse;
            let flags = read.flags();
            let mut read_reverse_strand = flags.is_reverse_complemented();
            if flags.is_segmented() && flags.is_last_segment() {
                read_reverse_strand = !read_reverse_strand;
            };
            let is_antisense = match self.chemistry_strandedness {
                Strand::Forward => tx_reverse_strand != read_reverse_strand,
                Strand::Reverse => tx_reverse_strand == read_reverse_strand,
            };
            let tx_strand = if is_antisense {
                Strand::Reverse
            } else {
                Strand::Forward
            };

            let mut alignment = TranscriptAlignment {
                gene_id: transcript.gene_id.clone(),
                strand: tx_strand,
                exon_align: None,
            };

            if is_exonic {
                // compute offsets
                let mut tx_offset = transcript.get_offset(genomic_start).unwrap().max(0) as u64;
                let tx_length = transcript.exon_len();

                // align the read to the exons
                if let Some((mut tx_cigar, tx_aligned_bases)) = splice_segments.align_junctions(
                    transcript,
                    self.junction_trim_bases,
                    self.intergenic_trim_bases,
                    self.intronic_trim_bases,
                ) {
                    // flip reverse strand
                    if tx_reverse_strand {
                        tx_offset = tx_length - (tx_offset + tx_aligned_bases);
                        tx_cigar.as_mut().reverse();
                    };
                    alignment = TranscriptAlignment {
                        gene_id: transcript.gene_id.clone(),
                        strand: tx_strand,
                        exon_align: Some(TxAlignProperties {
                            id: transcript.id.clone(),
                            pos: tx_offset,
                            cigar: tx_cigar,
                            alen: tx_aligned_bases,
                        }),
                    };
                }
            }
            Some(alignment)
        } else {
            None
        }
    }
}

#[derive(Eq, PartialEq, Debug)]
pub struct PairAnnotationData {
    /// Genes associated with the pair of alignments.
    pub genes: Vec<String>,
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
fn rescue_alignments_se(mut recs: Vec<AnnotatedAlignment>) -> Option<AnnotatedAlignment> {
    let n = recs.len();
    if n == 0 {
        None
    } else if n == 1 {
        recs.pop()
    } else {
        let mut promote_index: Option<usize> = None;
        let mut seen_genes = HashSet::new();

        for (i, rec) in recs.iter().enumerate() {
            match rec {
                AnnotatedAlignment::SeMapped(anno) => {
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
            let mut rec = recs.swap_remove(promote_index.unwrap());
            rec.set_rescued();
            Some(rec)
        } else {
            None
        }
    }
}

/// Attempts to rescue paired-end alignments using transcript annotations.
/// Returns true if rescue took place.
fn rescue_alignments_pe(mut pairs: Vec<AnnotatedAlignment>) -> Option<AnnotatedAlignment> {
    let n = pairs.len();
    if n == 0 {
        None
    } else if n == 1 {
        pairs.pop()
    } else {
        // Check if rescue is appropriate and determine which record to promote
        let mut seen_genes = HashSet::new();
        let mut promote_index: Option<usize> = None;

        for (i, pair) in pairs.iter_mut().enumerate() {
            match pair {
                AnnotatedAlignment::PeMapped(_, _, anno) => {
                    let genes = &anno.genes;

                    // Track which record/record-pair we should promote;
                    // Take the first record/pair with 1 gene
                    if genes.len() == 1 {
                        promote_index = promote_index.or(Some(i));
                    }

                    // Track number of distinct genes we're aligned to
                    seen_genes.extend(genes);
                }
                _ => unimplemented!(),
            }
        }

        // The alignment can be rescued if there is only one uniquely mapped gene
        if seen_genes.len() == 1 && promote_index.is_some() {
            let mut pair = pairs.swap_remove(promote_index.unwrap());
            pair.set_rescued();
            Some(pair)
        } else {
            None
        }
    }
}

#[derive(Encode, Decode, Eq, PartialEq, Ord, PartialOrd, Debug, Clone, Copy, Hash)]
pub enum AnnotationRegion {
    Exonic,
    Intronic,
    Intergenic,
}

#[derive(Eq, PartialEq, Debug)]
pub struct TranscriptAlignment {
    pub gene_id: String,
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
    pub pos: u64,
    pub cigar: Cigar,
    pub alen: u64,
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
