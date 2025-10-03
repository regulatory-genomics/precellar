use super::quantification::GeneAlignment;
use crate::transcriptome::annotate::RegionType;

use std::collections::{BTreeMap, HashMap, HashSet};

type Gene = usize;

#[derive(Debug, Default)]
pub(crate) struct DeDupResult {
    exon_uniq_umi: BTreeMap<Gene, usize>,
    intron_uniq_umi: BTreeMap<Gene, usize>,
    pub total_umi: u64,
    pub unique_umi: u64,
    pub uniq_exon: u64,
    pub uniq_intron: u64,
    pub uniq_mito: u64,
}

impl DeDupResult {
    /// Merge the intron and exon counts into a single count for each gene
    pub fn into_counts(mut self) -> impl Iterator<Item = (Gene, u64)> {
        self.exon_uniq_umi
            .iter()
            .for_each(|(gene, c)| {
                let val = self.intron_uniq_umi.entry(*gene).or_insert(0);
                *val += *c;
            });
        self.intron_uniq_umi.into_iter().map(|(gene, c)| (gene, c as u64))
    }
}

pub fn count_unique_umi<I>(alignments: I, mito_genes: &HashSet<usize>) -> DeDupResult
where
    I: IntoIterator<Item = GeneAlignment>,
{
    let mut result = DeDupResult::default();
    let mut umigene_counts_exon = HashMap::new();
    let mut umigene_counts_intron = HashMap::new();

    alignments.into_iter().for_each(|alignment| {
        let gene = alignment.idx;
        // We use empty vector to represent empty UMI
        let umi = alignment.umi.map(|x| x.into_bytes()).unwrap_or(Vec::new());
        
        match alignment.align_type {
            RegionType::Exonic => {
                *umigene_counts_exon.entry((umi, gene)).or_insert(0) += 1u64;
            }
            RegionType::Intronic => {
                *umigene_counts_intron.entry((umi, gene)).or_insert(0) += 1u64;
            }
            RegionType::Intergenic => {},
        }
    });

    let umi_correction= correct_umis(&umigene_counts_exon);
    let mut uniq_counts: HashMap<Gene, HashSet<&[u8]>> = HashMap::new();
    umigene_counts_exon.iter().for_each(|((umi, gene), c)| {
        result.total_umi += c;

        // Empty UMI
        if umi.is_empty() {
            result.unique_umi += c;
            if mito_genes.contains(&gene) {
                result.uniq_mito += c;
            }
            result.exon_uniq_umi.insert(*gene, *c as usize);
        } else {
            let corrected_umi = umi_correction.get(&(umi, *gene)).unwrap_or(umi);
            uniq_counts.entry(*gene).or_insert(HashSet::new()).insert(corrected_umi);
        }
    });
    uniq_counts.into_iter().for_each(|(gene, umis)| {
        let c = umis.len();
        result.unique_umi += c as u64;
        result.uniq_exon += c as u64;
        if mito_genes.contains(&gene) {
            result.uniq_mito += c as u64;
        }
        result.exon_uniq_umi.entry(gene).and_modify(|x| *x += c).or_insert(c);
    });

    let umi_correction= correct_umis(&umigene_counts_intron);
    let mut uniq_counts: HashMap<Gene, HashSet<&[u8]>> = HashMap::new();
    umigene_counts_intron.iter().for_each(|((umi, gene), c)| {
        result.total_umi += c;

        // Empty UMI
        if umi.is_empty() {
            result.unique_umi += c;
            if mito_genes.contains(&gene) {
                result.uniq_mito += c;
            }
            result.intron_uniq_umi.insert(*gene, *c as usize);
        } else {
            let corrected_umi = umi_correction.get(&(umi, *gene)).unwrap_or(umi);
            uniq_counts.entry(*gene).or_insert(HashSet::new()).insert(corrected_umi);
        }
    });
    uniq_counts.into_iter().for_each(|(gene, umis)| {
        let c = umis.len();
        result.unique_umi += c as u64;
        result.uniq_intron += c as u64;
        if mito_genes.contains(&gene) {
            result.uniq_mito += c as u64;
        }
        result.intron_uniq_umi.entry(gene).and_modify(|x| *x += c).or_insert(c);
    });

    result
}

/// Within each gene, correct Hamming-distance-one UMIs
fn correct_umis<'a>(umigene_counts: &'a HashMap<(Vec<u8>, Gene), u64>) -> HashMap<(&'a [u8], Gene), Vec<u8>> {
    let nucs = b"ACGT";

    let mut corrections = HashMap::new();

    for ((umi, gene), orig_count) in umigene_counts {
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
                let test_count = *umigene_counts.get(&(test_umi.clone(), *gene)).unwrap_or(&0u64);

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
            corrections.insert((umi.as_slice(), *gene), best_dest_umi);
        }
    }
    corrections
}