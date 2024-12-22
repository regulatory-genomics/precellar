use std::collections::{BTreeMap, HashMap, HashSet};

use super::quantification::GeneAlignment;

type Gene = usize;

pub fn count_unique_umi<I>(alignments: I) -> BTreeMap<Gene, usize>
where
    I: IntoIterator<Item = GeneAlignment>,
{
    fn get_uniq_counts(counts: HashMap<(Vec<u8>, Gene), u64>) -> BTreeMap<Gene, usize> {
        let umi_correction= correct_umis(&counts);
        
        let mut uniq_counts: HashMap<Gene, HashSet<&[u8]>> = HashMap::new();
        counts.keys().for_each(|(umi, gene)| {
            let corrected_umi = umi_correction.get(&(umi, *gene)).unwrap_or(umi);
            uniq_counts.entry(*gene).or_insert(HashSet::new()).insert(corrected_umi);
        });

        uniq_counts.into_iter().map(|(gene, umis)| (gene, umis.len())).collect()
    }

    let mut umigene_counts = HashMap::new();
    alignments.into_iter().for_each(|alignment| {
        let gene = alignment.idx;
        let umi = alignment.umi.unwrap().into_bytes();
        *umigene_counts.entry((umi, gene)).or_insert(0) += 1u64;
    });

    get_uniq_counts(umigene_counts)
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