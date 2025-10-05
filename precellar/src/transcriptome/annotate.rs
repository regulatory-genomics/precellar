/// This module provides utilities for annotating bam records by mapping them to a transcriptome.
/// It supports both single-end and paired-end alignments and uses transcript annotations for gene-level
/// and exon-level classification.
use crate::align::MultiMapR;
use crate::fragment::Fragment;
use crate::transcriptome::align::ChemistryStrandness;
use crate::transcriptome::align::JunctionAlignOptions;
use crate::transcriptome::align::SplicedRecord;
use crate::transcriptome::align::TranscriptAlignment;
use crate::transcriptome::Transcript;

use anyhow::Result;
use bed_utils::bed::map::GIntervalMap;
use bed_utils::bed::GenomicRange;
use bincode::{Decode, Encode};
use log::info;
use noodles::sam;
use noodles::sam::record::data::field::value::base_modifications::group::Strand;
use std::collections::{BTreeMap, HashSet};

#[derive(Encode, Decode, Eq, PartialEq, Ord, PartialOrd, Debug, Copy, Clone, Hash)]
pub enum RegionType {
    /// Aligned to exonic regions of a gene
    Exonic,
    /// Aligned to intronic regions of a gene
    Intronic,
    /// Aligned to regions outside of any gene
    Intergenic,
}

#[derive(Eq, PartialEq, Debug, Clone)]
struct Annotation {
    aln_sense: Vec<TranscriptAlignment>,
    aln_antisense: Vec<TranscriptAlignment>,
    genes: HashSet<String>,
    region_type: RegionType,
}

impl Annotation {
    fn is_antisense(&self) -> bool {
        !self.aln_antisense.is_empty() && self.aln_sense.is_empty()
    }
}

#[derive(Debug, Encode, Decode)]
pub struct AnnotatedAlignment {
    read1: SplicedRecord,
    read2: Option<SplicedRecord>,
    uniquely_mapped_gene: Option<String>, // genes that this read uniquely maps to (transcriptomically)
    region_type: RegionType,
    antisense: bool,
    umi: Option<String>,
}

impl AnnotatedAlignment {
    fn from_se(rec: SplicedRecord, anno: Annotation) -> Self {
        let antisense = anno.is_antisense();
        Self {
            read1: rec,
            read2: None,
            uniquely_mapped_gene: if anno.genes.len() == 1 {
                Some(anno.genes.iter().next().unwrap().to_string())
            } else {
                None
            },
            region_type: anno.region_type,
            antisense,
            umi: None,
        }
    }

    fn from_pe(
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
        let antisense = anno1.is_antisense() && anno2.is_antisense();
        let mapped_genes: Vec<_> = find_common_genes(&anno1, &anno2).into_iter().collect();
        Self {
            read1,
            read2: Some(read2),
            uniquely_mapped_gene: if mapped_genes.len() == 1 {
                Some(mapped_genes[0].clone())
            } else {
                None
            },
            region_type: region,
            antisense,
            umi: None,
        }
    }

    pub fn is_paired(&self) -> bool {
        self.read2.is_some()
    }

    // A read is counted as antisense if it has any alignments that are consistent
    // with an exon of a transcript but antisense to it, and has no sense alignments.
    pub fn is_antisense(&self) -> bool {
        self.antisense
    }

    pub fn is_intergenic(&self) -> bool {
        self.region_type == RegionType::Intergenic
    }

    pub fn is_exonic(&self) -> bool {
        self.region_type == RegionType::Exonic
    }

    pub fn is_intronic(&self) -> bool {
        self.region_type == RegionType::Intronic
    }

    /// Returns true if the read is confidently mapped to a single gene.
    /// It is considered confidently mapped if it maps to exonic or intronic regions of a single gene.
    /// Or if it is intergenic and maps to a single locus.
    pub fn is_confidently_mapped(&self) -> bool {
        self.uniquely_mapped_gene.is_some() || self.is_intergenic() || self.is_antisense()
    }

    /// Returns the gene if the read is confidently mapped to a single gene.
    pub fn uniquely_mapped_gene(&self) -> Option<&str> {
        self.uniquely_mapped_gene.as_deref()
    }

    pub fn read1(&self) -> &SplicedRecord {
        &self.read1
    }

    pub fn read2(&self) -> Option<&SplicedRecord> {
        self.read2.as_ref()
    }

    pub fn umi(&self) -> Option<&str> {
        self.umi.as_deref()
    }

    pub fn region_type(&self) -> RegionType {
        self.region_type
    }

    pub fn to_fragments(&self) -> impl Iterator<Item = Fragment> + '_ {
        self.read1
            .to_fragments()
            .chain(self.read2.iter().flat_map(|r| r.to_fragments()))
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
    pub fn new(
        transcripts: impl IntoIterator<Item = Transcript>,
        options: JunctionAlignOptions,
    ) -> Self {
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
    /// * `multi_map` - Vector of single-end alignment records.
    ///
    /// Returns `Some(AnnotatedAlignment)` if a confident alignment is found, otherwise `None`.
    pub fn annotate_alignments_se(
        &self,
        header: &sam::Header,
        multi_map: &MultiMapR,
    ) -> Option<AnnotatedAlignment> {
        let results = multi_map.iter().flat_map(|rec| {
            let spliced_rec = SplicedRecord::new(rec, header).unwrap()?;
            let anno = self.annotate_alignment(&spliced_rec).unwrap();
            Some((spliced_rec, anno))
        });
        let mut aln = rescue_alignments_se(results)?;
        aln.umi = multi_map.umi().unwrap();
        Some(aln)
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
        multi_map1: &MultiMapR,
        multi_map2: &MultiMapR,
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
        let mut aln = rescue_alignments_pe(result)?;
        aln.umi = multi_map1.umi().unwrap();
        Some(aln)
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

    /// Detects the strandness of the RNA-seq data based on the alignments.
    pub fn set_strandness(
        &mut self,
        header: &sam::Header,
        alignments: &[(Option<MultiMapR>, Option<MultiMapR>)],
    ) {
        if self.options.chemistry_strandedness.is_some() {
            return;
        }

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
                    if let Some(rec) = SplicedRecord::new(&aln.primary, header).unwrap() {
                        let region =
                            GenomicRange::new(&rec.chrom, rec.start() as u64, rec.end() as u64);
                        self.transcripts.find(&region).for_each(|(_, transcript)| {
                            match rec.orientation(transcript, self.options.region_min_overlap) {
                                Some(Strand::Forward) => is_sense = true,
                                Some(Strand::Reverse) => is_antisense = true,
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
            self.options.chemistry_strandedness = Some(ChemistryStrandness::Forward);
        } else if percent_antisense > 75.0 {
            info!("Chemistry strandness is set to: Reverse");
            self.options.chemistry_strandedness = Some(ChemistryStrandness::Reverse);
        } else {
            info!("Chemistry strandness is set to: Unstranded");
            self.options.chemistry_strandedness = Some(ChemistryStrandness::Unstranded);
        }
    }
}

/// When multiple genomic alignments exist for a single-end read, attempts to rescue one
/// using transcript annotations.
/// If a read mapped to exonic loci from a single gene and also to non-exonic loci,
/// it is considered uniquely mapped to one of the exonic loci.
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
        let anno = AnnotatedAlignment::from_se(rescued.0, rescued.1);
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
        let anno =
            AnnotatedAlignment::from_pe(rescued.0 .0, rescued.0 .1, rescued.1 .0, rescued.1 .1);
        Some(anno)
    } else {
        None
    }
}

/// Take the intersection of the non-empty gene sets of the mates
fn find_common_genes(anno1: &Annotation, anno2: &Annotation) -> HashSet<String> {
    match (!anno1.genes.is_empty(), !anno2.genes.is_empty()) {
        (true, false) => anno1.genes.clone(),
        (false, true) => anno2.genes.clone(),
        _ => anno1.genes.intersection(&anno2.genes).cloned().collect(),
    }
}
