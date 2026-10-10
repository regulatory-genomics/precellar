use std::collections::{BTreeMap, HashMap};

use bed_utils::bed::{BEDLike, MergeBed};
use indexmap::IndexMap;
use itertools::Itertools;

use crate::utils::get_directional_umi_mapping;
use crate::{fragment::Fragment, transcriptome::TxAlignment};

type GeneIndex = usize;

pub enum Annotation {
    Spliced,
    Unspliced,
    Ambiguous,
}

/// Per-gene UMI counts of one cell, split by splicing status.
#[derive(Debug, Default)]
pub struct SplicedCounts {
    pub spliced: Vec<(GeneIndex, u32)>,
    pub unspliced: Vec<(GeneIndex, u32)>,
    pub ambiguous: Vec<(GeneIndex, u32)>,
}

/// A group of transcript alignments corresponding to a single UMI.
struct AlignmentGroup<'a>(Vec<&'a TxAlignment>);

impl AlignmentGroup<'_> {
    /// Annotate spliced and unspliced alignments with transcript information.
    ///
    /// 1. A molecule was annotated as spliced if all of the reads in the set supporting a
    ///    given molecule are exonic-only (Caveat: a read contained entirely within an exon
    ///    is likely coming from unspliced transcripts technically)
    /// 2. A molecule was annotated as unspliced if at least one read is spanning.
    /// 3. A molecule was annotated as unspliced if ALL of the compatible transcript models
    ///    had at least one read that maps to introns or exon-intron boundary.
    /// 4. Other molecules are annotated as ambiguous.
    fn to_annotation(&self) -> Annotation {
        let mut exonic = true;
        for aln in self.0.iter() {
            if aln.is_spanning() {
                return Annotation::Unspliced;
            }

            if !aln.is_exonic_only() {
                exonic = false;
            }
        }

        if exonic {
            Annotation::Spliced
        } else {
            let transcript_groups = self
                .0
                .iter()
                .flat_map(|x| x.alignments())
                .sorted_by(|a, b| a.0.cmp(&b.0))
                .chunk_by(|x| x.0.clone());
            if transcript_groups.into_iter().all(|(_, group)| {
                group
                    .into_iter()
                    .any(|x| x.1.is_intronic() || x.1.is_spanning())
            }) {
                Annotation::Unspliced
            } else {
                Annotation::Ambiguous
            }
        }
    }

    fn to_fragments(&self) -> Vec<Fragment> {
        self.0
            .iter()
            .flat_map(|x| x.to_fragments())
            .sorted_by(|a, b| {
                a.strand
                    .cmp(&b.strand)
                    .then(a.chrom.cmp(&b.chrom))
                    .then(a.start.cmp(&b.start))
                    .then(a.end.cmp(&b.end))
            })
            .chunk_by(|x| x.strand)
            .into_iter()
            .flat_map(|(strand, group)| {
                group.merge_sorted_bed().map(move |x| Fragment {
                    chrom: x.chrom().to_string(),
                    start: x.start(),
                    end: x.end(),
                    strand,
                    barcode: None,
                    count: 1,
                    extended: None,
                })
            })
            .collect()
    }
}

impl<'a> From<Vec<&'a TxAlignment>> for AlignmentGroup<'a> {
    fn from(value: Vec<&'a TxAlignment>) -> Self {
        AlignmentGroup(value)
    }
}

/// A list of TxAlignments grouped by gene index and UMI. Each entry contains the gene index
/// and a list of lists of TxAlignments, where each inner list corresponds to a unique
/// UMI within that gene.
pub struct GeneCounter<'a>(Vec<(GeneIndex, Vec<AlignmentGroup<'a>>)>);

impl<'a> GeneCounter<'a> {
    pub fn new(data: &'a Vec<TxAlignment>, genes: &IndexMap<String, String>) -> Self {
        // Group the Tx alignments by gene
        let mut gene_group = BTreeMap::new();
        data.iter().for_each(|x| {
            let idx = genes.get_full(x.uniquely_mapped_gene().unwrap()).unwrap().0;
            gene_group.entry(idx).or_insert_with(Vec::new).push(x);
        });

        let umi_group = gene_group
            .into_iter()
            .map(|(idx, group)| {
                // Within each gene, compute UMI counts
                let mut aln_with_umi = HashMap::new();
                group.iter().for_each(|aln| {
                    if let Some(umi) = aln.umi() {
                        *aln_with_umi.entry(umi.as_bytes().to_vec()).or_insert(0) += 1;
                    }
                });

                let umi_mapping = get_directional_umi_mapping(&aln_with_umi);
                let mut aln_without_umi = Vec::new();

                let mut umi_group = HashMap::new();
                group.into_iter().for_each(|aln| {
                    if let Some(umi) = aln.umi() {
                        let umi = umi.as_bytes();
                        let corrected_umi =
                            umi_mapping.get(umi).map_or(umi.to_vec(), |x| x.clone());
                        umi_group
                            .entry(corrected_umi)
                            .or_insert_with(Vec::new)
                            .push(aln);
                    } else {
                        aln_without_umi.push(vec![aln]);
                    }
                });

                let result = aln_without_umi
                    .into_iter()
                    .chain(umi_group.into_iter().map(|x| x.1))
                    .map(|x| x.into())
                    .collect();
                (idx, result)
            })
            .collect();
        GeneCounter(umi_group)
    }

    pub fn to_counts(&self) -> Vec<(GeneIndex, u32)> {
        self.0
            .iter()
            .map(|(gene_idx, umi_group)| (*gene_idx, umi_group.len() as u32))
            .collect()
    }

    /// Splits the UMI counts of this cell by splicing status.
    ///
    /// The three categories are mutually exclusive and together account for
    /// every UMI, so their sums equal [`GeneCounter::to_counts`].
    pub fn to_spliced_counts(&self) -> SplicedCounts {
        let mut counts = SplicedCounts::default();
        self.0.iter().for_each(|(gene_idx, umi_group)| {
            let mut spliced: u32 = 0;
            let mut unspliced: u32 = 0;
            let mut ambiguous: u32 = 0;
            umi_group
                .iter()
                .for_each(|reads| match reads.to_annotation() {
                    Annotation::Spliced => {
                        spliced = spliced.checked_add(1).unwrap();
                    }
                    Annotation::Unspliced => {
                        unspliced = unspliced.checked_add(1).unwrap();
                    }
                    Annotation::Ambiguous => {
                        ambiguous = ambiguous.checked_add(1).unwrap();
                    }
                });
            if spliced > 0 {
                counts.spliced.push((*gene_idx, spliced));
            }
            if unspliced > 0 {
                counts.unspliced.push((*gene_idx, unspliced));
            }
            if ambiguous > 0 {
                counts.ambiguous.push((*gene_idx, ambiguous));
            }
        });
        counts
    }

    pub fn to_fragments(&self) -> impl Iterator<Item = Fragment> + '_ {
        self.0
            .iter()
            .flat_map(|(_, umi_group)| umi_group.iter().flat_map(|reads| reads.to_fragments()))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::transcriptome::align::{AlignType, SplicedRecord};

    /// One read of a molecule, given as the per-transcript alignment types it
    /// is compatible with. Coordinates are irrelevant here: `to_annotation`
    /// reads only the `(transcript_id, AlignType)` list.
    fn read(umi: &str, alignment: &[(&str, AlignType)]) -> TxAlignment {
        TxAlignment::SeAligned {
            read: SplicedRecord::for_test("chr1", vec![(0, 10)]),
            alignment: alignment
                .iter()
                .map(|(tx, ty)| (tx.to_string(), ty.clone()))
                .collect(),
            gene_id: "gene1".to_string(),
            barcode: "ACGT".to_string(),
            umi: Some(umi.to_string()),
        }
    }

    /// The three categories partition the UMIs of a gene, so the layers sum to
    /// the total count.
    #[test]
    fn splicing_status_partitions_the_umis_of_a_gene() {
        let genes = IndexMap::from([("gene1".to_string(), "GeneName1".to_string())]);

        let alignments = vec![
            // Spliced: every transcript is exonic (rule 2).
            read("AAAAAAAA", &[("tx1", AlignType::Exonic)]),
            // Unspliced: a spanning read decides the molecule (rule 1).
            read("CCCCCCCC", &[("tx1", AlignType::Spanning)]),
            // Ambiguous: tx1 says exonic, tx2 says intronic, and no read spans,
            // so neither rule 2 nor rule 3 holds.
            read(
                "GGGGGGGG",
                &[("tx1", AlignType::Exonic), ("tx2", AlignType::Intronic)],
            ),
        ];

        let counter = GeneCounter::new(&alignments, &genes);
        let counts = counter.to_spliced_counts();

        assert_eq!(counts.spliced, vec![(0, 1)]);
        assert_eq!(counts.unspliced, vec![(0, 1)]);
        assert_eq!(counts.ambiguous, vec![(0, 1)]);

        // The invariant that makes the three layers comparable: nothing is
        // dropped, so `X == spliced + unspliced + ambiguous`.
        let total: u32 = counter.to_counts().iter().map(|(_, c)| c).sum();
        let split: u32 = counts
            .spliced
            .iter()
            .chain(&counts.unspliced)
            .chain(&counts.ambiguous)
            .map(|(_, c)| c)
            .sum();
        assert_eq!(total, 3);
        assert_eq!(total, split);
    }

    /// A molecule with no spanning read is still unspliced when every
    /// compatible transcript has intronic evidence.
    #[test]
    fn intronic_evidence_on_all_transcripts_is_unspliced() {
        let genes = IndexMap::from([("gene1".to_string(), "GeneName1".to_string())]);
        let alignments = vec![read(
            "AAAAAAAA",
            &[("tx1", AlignType::Intronic), ("tx2", AlignType::Intronic)],
        )];

        let counts = GeneCounter::new(&alignments, &genes).to_spliced_counts();

        assert_eq!(counts.unspliced, vec![(0, 1)]);
        assert!(counts.spliced.is_empty());
        assert!(counts.ambiguous.is_empty());
    }
}
