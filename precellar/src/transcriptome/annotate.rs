/// This module provides utilities for annotating bam records by mapping them to a transcriptome.
/// It supports both single-end and paired-end alignments and uses transcript annotations for gene-level
/// and exon-level classification.
use crate::align::MultiMapR;
use crate::transcriptome::align::JunctionAlignOptions;
use crate::transcriptome::align::SplicedRecord;
use crate::transcriptome::align::TranscriptAlignment;
use crate::transcriptome::Transcript;

use anyhow::Result;
use bed_utils::bed::map::GIntervalMap;
use bed_utils::bed::GenomicRange;
use bincode::{Decode, Encode};
use noodles::sam;
use noodles::sam::record::data::field::value::base_modifications::group::Strand;
use std::collections::{BTreeMap, HashSet};

#[derive(Encode, Decode, Eq, PartialEq, Ord, PartialOrd, Debug, Copy, Clone, Hash)]
pub enum RegionType {
    Exonic,
    Intronic,
    Intergenic,
}

#[derive(Eq, PartialEq, Debug, Clone)]
pub struct Annotation {
    aln_sense: Vec<TranscriptAlignment>,
    aln_antisense: Vec<TranscriptAlignment>,
    genes: HashSet<String>,
    region_type: RegionType,
}

#[derive(Debug)]
pub struct AnnotatedAlignment {
    read1: SplicedRecord,
    read2: Option<SplicedRecord>,
    mapped_genes: Vec<String>, // genes that this read maps to (transcriptomically)
    region_type: RegionType,
    rescued: bool,
}

impl AnnotatedAlignment {
    pub fn from_se(rec: SplicedRecord, anno: Annotation) -> Self {
        Self {
            read1: rec,
            read2: None,
            mapped_genes: anno.genes.into_iter().collect(),
            region_type: anno.region_type,
            rescued: false,
        }
    }

    pub fn from_pe(
        read1: SplicedRecord,
        anno1: Annotation,
        read2: SplicedRecord,
        anno2: Annotation,
    ) -> Self {
        let region = match (anno1.region_type, anno2.region_type) {
            (RegionType::Intronic, _) => RegionType::Intronic,
            (_, RegionType::Intronic) => RegionType::Intronic,
            (RegionType::Intergenic, RegionType::Intergenic) => RegionType::Intergenic,
            _ => RegionType::Exonic,
        };
        Self {
            read1,
            read2: Some(read2),
            mapped_genes: find_common_genes(&anno1, &anno2).into_iter().collect(),
            region_type: region,
            rescued: false,
        }
    }

    pub fn is_paired(&self) -> bool {
        self.read2.is_some()
    }

    /// Returns the gene if the read is confidently mapped to a single gene.
    pub fn confidently_mapped_gene(&self) -> Option<&str> {
        if self.mapped_genes.len() == 1 {
            Some(&self.mapped_genes[0])
        } else {
            None
        }
    }

    pub fn read1(&self) -> &SplicedRecord {
        &self.read1
    }

    pub fn read2(&self) -> Option<&SplicedRecord> {
        self.read2.as_ref()
    }

    pub fn region_type(&self) -> RegionType {
        self.region_type
    }
}

/// Manages the annotation of alignments using transcriptome data.
#[derive(Debug, Clone)]
pub struct AlignmentAnnotator {
    /// Map of genomic intervals to transcripts.
    pub(crate) transcripts: GIntervalMap<Transcript>,
    options: JunctionAlignOptions,
}

impl AlignmentAnnotator {
    /// Creates a new `AlignmentAnnotator` with the provided transcripts.
    pub fn new(transcripts: impl IntoIterator<Item = Transcript>, options: JunctionAlignOptions) -> Self {
        let transcripts = transcripts
            .into_iter()
            .map(|x| (GenomicRange::new(&x.chrom, x.start, x.end), x))
            .collect();
        Self {
            transcripts,
            options,
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
        multi_map: MultiMapR,
    ) -> Option<AnnotatedAlignment> {
        let results = multi_map.iter().flat_map(|rec| {
            let spliced_rec = SplicedRecord::new(rec, header).unwrap()?;
            let anno = self.annotate_alignment(&spliced_rec).unwrap();
            Some((spliced_rec, anno))
        });
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
        multi_map1: MultiMapR,
        multi_map2: MultiMapR,
    ) -> Option<AnnotatedAlignment> {
        let result = multi_map1
            .iter()
            .zip(multi_map2.iter())
            .flat_map(|(r1, r2)| {
                let rec1 = SplicedRecord::new(r1, header).unwrap()?;
                let rec2 = SplicedRecord::new(r2, header).unwrap()?;
                let anno1 = self.annotate_alignment(&rec1).unwrap();
                let anno2 = self.annotate_alignment(&rec2).unwrap();
                Some(((rec1, anno1), (rec2, anno2)))
            });
        rescue_alignments_pe(result)
    }

    fn annotate_alignment(&self, spliced_rec: &SplicedRecord) -> Result<Annotation> {
        let region = GenomicRange::new(
            &spliced_rec.chrom,
            spliced_rec.start() as u64,
            spliced_rec.end() as u64,
        );
        let alignments = self
            .transcripts
            .find(&region)
            .flat_map(|(_, transcript)| spliced_rec.align_transcript(transcript, &self.options))
            .collect::<Vec<_>>();

        let mut seen_genes = HashSet::new();
        let mut transcripts = BTreeMap::new();
        let mut antisense = BTreeMap::new();
        let annotation_region;
        if alignments.is_empty() {
            annotation_region = RegionType::Intergenic;
        } else if alignments.iter().any(|x| x.is_exonic()) {
            annotation_region = RegionType::Exonic;
            // Check if there are transcriptome compatible alignments
            alignments.into_iter().rev().for_each(|aln| {
                if aln.is_exonic() {
                    match aln.strand {
                        Strand::Forward => {
                            // Transcript sense alignment
                            seen_genes.insert(aln.gene_id.clone());
                            transcripts.insert(aln.transcript_id.clone(), aln);
                        }
                        Strand::Reverse => {
                            // Transcript anti-sense alignment
                            antisense.insert(aln.transcript_id.clone(), aln);
                        }
                    }
                }
            });
        } else {
            annotation_region = RegionType::Intronic;
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

        let annotation = Annotation {
            aln_sense: transcripts.into_values().collect::<Vec<_>>(),
            aln_antisense: antisense.into_values().collect::<Vec<_>>(),
            genes: seen_genes,
            region_type: annotation_region,
        };

        Ok(annotation)
    }
}

/// When multiple genomic alignments exist for a single-end read, attempts to rescue one
/// using transcript annotations. A read can be rescued if it maps to a single gene
/// transcriptomically, even if it maps to multiple genomic loci,
/// e.g., due to pseudogenes or repetitive elements.
fn rescue_alignments_se(
    recs: impl IntoIterator<Item = (SplicedRecord, Annotation)>,
) -> Option<AnnotatedAlignment> {
    let mut recs = recs.into_iter().peekable();
    let first_rec = recs.peek()?.clone();

    let mut n = 0;
    let mut rescued = None;
    let mut seen_genes = HashSet::new();
    recs.for_each(|(rec, anno)| {
        n += 1;
        // Only consider transcriptomic alignments for rescue
        if anno.aln_sense.iter().any(|x| x.is_exonic()) {
            // Track number of distinct genes we're aligned to
            seen_genes.extend(anno.genes.iter().cloned());

            // Track which record/record-pair we should promote;
            // Take the first record/pair with 1 gene
            if anno.genes.len() == 1 && rescued.is_none() {
                rescued = Some((rec, anno));
            }
        }
    });

    if n == 1 {
        Some(AnnotatedAlignment::from_se(first_rec.0, first_rec.1))
    } else if seen_genes.len() == 1 && rescued.is_some() {
        let rescued = rescued.unwrap();
        let mut anno = AnnotatedAlignment::from_se(rescued.0, rescued.1);
        anno.rescued = true;
        Some(anno)
    } else {
        None
    }
}

/// Attempts to rescue paired-end alignments using transcript annotations.
/// Returns true if rescue took place.
fn rescue_alignments_pe(
    recs: impl IntoIterator<Item = ((SplicedRecord, Annotation), (SplicedRecord, Annotation))>,
) -> Option<AnnotatedAlignment> {
    let mut recs = recs.into_iter().peekable();
    let first_rec = recs.peek()?.clone();

    let mut n = 0;
    let mut rescued = None;
    let mut seen_genes = HashSet::new();
    recs.for_each(|((rec1, anno1), (rec2, anno2))| {
        n += 1;
        let genes = find_common_genes(&anno1, &anno2);

        // Track which record/record-pair we should promote;
        // Take the first record/pair with 1 gene
        if genes.len() == 1 && rescued.is_none() {
            rescued = Some(((rec1, anno1), (rec2, anno2)));
        }

        // Track number of distinct genes we're aligned to
        seen_genes.extend(genes);
    });

    if n == 1 {
        Some(AnnotatedAlignment::from_pe(
            first_rec.0 .0,
            first_rec.0 .1,
            first_rec.1 .0,
            first_rec.1 .1,
        ))
    } else if seen_genes.len() == 1 && rescued.is_some() {
        let rescued = rescued.unwrap();
        let mut anno =
            AnnotatedAlignment::from_pe(rescued.0 .0, rescued.0 .1, rescued.1 .0, rescued.1 .1);
        anno.rescued = true;
        Some(anno)
    } else {
        None
    }
}

/// Take the intersection of the non-empty gene sets of the mates
pub fn find_common_genes(anno1: &Annotation, anno2: &Annotation) -> HashSet<String> {
    match (!anno1.genes.is_empty(), !anno2.genes.is_empty()) {
        (true, false) => anno1.genes.clone(),
        (false, true) => anno2.genes.clone(),
        _ => anno1.genes.intersection(&anno2.genes).cloned().collect(),
    }
}
