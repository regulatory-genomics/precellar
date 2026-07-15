use std::collections::{BTreeMap, HashMap};

use bed_utils::bed::{BEDLike, MergeBed};
use indexmap::IndexMap;
use itertools::Itertools;

use crate::{fragment::Fragment, transcriptome::TxAlignment};

type GeneIndex = usize;

pub enum Annotation {
    Spliced,
    Unspliced,
    Ambiguous,
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

                let umi_mapping = get_umi_mapping(&aln_with_umi);
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

    pub fn to_spliced_counts(&self) -> (Vec<(GeneIndex, u32)>, Vec<(GeneIndex, u32)>) {
        let mut spliced_counts = Vec::new();
        let mut unspliced_counts = Vec::new();
        self.0.iter().for_each(|(gene_idx, umi_group)| {
            let mut spliced: u32 = 0;
            let mut unspliced: u32 = 0;
            umi_group
                .iter()
                .for_each(|reads| match reads.to_annotation() {
                    Annotation::Spliced => {
                        spliced = spliced.checked_add(1).unwrap();
                    }
                    Annotation::Unspliced => {
                        unspliced = unspliced.checked_add(1).unwrap();
                    }
                    Annotation::Ambiguous => {}
                });
            if spliced > 0 {
                spliced_counts.push((*gene_idx, spliced));
            }
            if unspliced > 0 {
                unspliced_counts.push((*gene_idx, unspliced));
            }
        });
        (spliced_counts, unspliced_counts)
    }

    pub fn to_fragments(&self) -> impl Iterator<Item = Fragment> + '_ {
        self.0
            .iter()
            .flat_map(|(_, umi_group)| umi_group.iter().flat_map(|reads| reads.to_fragments()))
    }
}

/// Returns a map from each UMI to its corrected UMI by correcting Hamming-distance-one UMIs.
fn get_umi_mapping(umi_count: &HashMap<Vec<u8>, u64>) -> HashMap<Vec<u8>, Vec<u8>> {
    let nucs = b"ACGT";

    let mut corrections = HashMap::new();

    for (umi, orig_count) in umi_count {
        let mut test_umi = umi.clone();

        let mut best_dest_count = *orig_count;
        let mut best_dest_umi = umi.to_vec();

        for pos in 0..umi.len() {
            // Try each nucleotide at this position
            for test_char in nucs {
                if *test_char == umi[pos] {
                    // Skip the identitical nucleotide
                    continue;
                }
                test_umi[pos] = *test_char;

                // Test for the existence of this mutated UMI
                let test_count = *umi_count.get(&test_umi).unwrap_or(&0);

                // If there's a 1-HD UMI w/ greater count, move to that UMI.
                // If there's a 1-HD UMI w/ equal count, move to the lexicographically larger UMI.
                if test_count > best_dest_count
                    || (test_count == best_dest_count && test_umi > best_dest_umi)
                {
                    best_dest_umi = test_umi.clone();
                    best_dest_count = test_count;
                }
            }
            // Reset this position to the unmutated sequence
            test_umi[pos] = umi[pos];
        }
        if *umi != best_dest_umi {
            corrections.insert(umi.to_vec(), best_dest_umi);
        }
    }
    corrections
}
