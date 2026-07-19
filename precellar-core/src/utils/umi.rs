use std::collections::HashMap;

/// Returns a map from each UMI to its corrected UMI by correcting Hamming-distance-one UMIs.
pub(crate) fn get_umi_mapping(umi_count: &HashMap<Vec<u8>, u64>) -> HashMap<Vec<u8>, Vec<u8>> {
    let nucs = b"ACGT";

    let mut corrections = HashMap::new();

    for (umi, orig_count) in umi_count {
        let mut test_umi = umi.clone();

        let mut best_dest_count = *orig_count;
        let mut best_dest_umi = umi.to_vec();

        for pos in 0..umi.len() {
            // Try each nucleotide at this position.
            for test_char in nucs {
                if *test_char == umi[pos] {
                    continue;
                }
                test_umi[pos] = *test_char;

                let test_count = *umi_count.get(&test_umi).unwrap_or(&0);

                // Prefer a higher-count UMI, or the lexicographically larger UMI on ties.
                if test_count > best_dest_count
                    || (test_count == best_dest_count && test_umi > best_dest_umi)
                {
                    best_dest_umi = test_umi.clone();
                    best_dest_count = test_count;
                }
            }
            test_umi[pos] = umi[pos];
        }
        if *umi != best_dest_umi {
            corrections.insert(umi.to_vec(), best_dest_umi);
        }
    }
    corrections
}

/// Returns a directional, transitive UMI correction map.
///
/// A lower-count UMI is attached to a one-edit neighbor only when the neighbor has a
/// strictly greater count and satisfies `parent_count >= 2 * child_count - 1`. Mappings
/// are then resolved to their final root.
pub(crate) fn get_directional_umi_mapping(
    umi_count: &HashMap<Vec<u8>, u64>,
) -> HashMap<Vec<u8>, Vec<u8>> {
    let nucs = b"ACGT";
    let mut direct_mapping = HashMap::new();

    for (umi, child_count) in umi_count {
        let mut test_umi = umi.clone();
        let threshold = child_count.saturating_mul(2).saturating_sub(1);
        let mut best_parent: Option<(u64, Vec<u8>)> = None;

        for pos in 0..umi.len() {
            for test_char in nucs {
                if *test_char == umi[pos] {
                    continue;
                }
                test_umi[pos] = *test_char;

                if let Some(&parent_count) = umi_count.get(&test_umi) {
                    if parent_count > *child_count && parent_count >= threshold {
                        let candidate = (parent_count, test_umi.clone());
                        if best_parent.as_ref().is_none_or(|best| candidate > *best) {
                            best_parent = Some(candidate);
                        }
                    }
                }
            }
            test_umi[pos] = umi[pos];
        }

        if let Some((_, parent)) = best_parent {
            direct_mapping.insert(umi.clone(), parent);
        }
    }

    let mut corrections = HashMap::new();
    for umi in umi_count.keys() {
        let mut root = umi;
        while let Some(parent) = direct_mapping.get(root) {
            root = parent;
        }
        if root != umi {
            corrections.insert(umi.clone(), root.clone());
        }
    }
    corrections
}

#[cfg(test)]
mod tests {
    use super::*;

    const TFSEQ_UMI_COUNTS: &str = "AAGTGGATAT:1;AAGTGGATTT:16;AAGTGGGTTT:4;AAGTGGTATT:1;AAGTGGTTTT:40;AAGTGTATTT:47;AAGTGTGTTT:1;AAGTGTTTAT:2;AAGTGTTTTT:24;ATACAATTGG:1;ATCACCCCCG:4;ATCACTCTCG:1;ATCACTGTCG:1;ATCACTGTGG:136;ATCGCTGTGG:1;ATGCCGACAT:1;ATGCCGTCAT:1;ATGCCGTCGT:15;ATGCCGTCTG:2;ATGCCGTCTT:75;ATGCCGTGAT:1;ATGCCGTGGA:1;ATGCCGTGGG:2;ATGCCGTGGT:4;ATGCCGTGTT:8;ATTCGCCCGC:14;GAATGGCGCG:5;GATATTCCTT:1;GGACCTAGGT:103;GTTAAACTAT:20;GTTAAACTTT:2;TACGTCAATG:1;TCTAACGCAT:20;TCTACCGCAT:1;TGTCAGATAA:194;TGTCGGATAA:1;TTATAAAAAT:1;TTATATAAAT:31";

    fn parse_counts(input: &str) -> HashMap<Vec<u8>, u64> {
        input
            .split(';')
            .map(|entry| {
                let (umi, count) = entry.split_once(':').unwrap();
                (umi.as_bytes().to_vec(), count.parse::<u64>().unwrap())
            })
            .collect()
    }

    fn apply_mapping(
        counts: &HashMap<Vec<u8>, u64>,
        mapping: &HashMap<Vec<u8>, Vec<u8>>,
    ) -> HashMap<Vec<u8>, u64> {
        let mut corrected_counts = HashMap::new();
        for (umi, count) in counts {
            let corrected_umi = mapping.get(umi).unwrap_or(umi);
            *corrected_counts.entry(corrected_umi.clone()).or_insert(0) += count;
        }
        corrected_counts
    }

    #[test]
    fn corrects_lower_count_hamming_distance_one_umi() {
        let counts = HashMap::from([(b"AAAA".to_vec(), 3), (b"AAAT".to_vec(), 1)]);

        assert_eq!(
            get_umi_mapping(&counts).get(b"AAAT" as &[u8]),
            Some(&b"AAAA".to_vec())
        );
    }

    #[test]
    fn directional_mapping_resolves_transitive_chains() {
        let counts = HashMap::from([
            (b"AAAA".to_vec(), 1),
            (b"AAAT".to_vec(), 2),
            (b"AAGT".to_vec(), 4),
        ]);

        assert_eq!(
            get_directional_umi_mapping(&counts).get(b"AAAA" as &[u8]),
            Some(&b"AAGT".to_vec())
        );
    }

    #[test]
    fn directional_mapping_does_not_merge_equal_count_umis() {
        let counts = HashMap::from([(b"AAAA".to_vec(), 1), (b"AAAT".to_vec(), 1)]);

        assert!(get_directional_umi_mapping(&counts).is_empty());
    }

    #[test]
    fn uses_lexicographically_larger_umi_on_equal_count() {
        let counts = HashMap::from([
            (b"AAAA".to_vec(), 1),
            (b"AAAC".to_vec(), 1),
            (b"AAAG".to_vec(), 1),
        ]);

        assert_eq!(
            get_umi_mapping(&counts).get(b"AAAA" as &[u8]),
            Some(&b"AAAG".to_vec())
        );
    }

    #[test]
    fn selects_best_neighbor_for_aagtggattt_from_tfseq_example() {
        let counts = parse_counts(TFSEQ_UMI_COUNTS);
        let mapping = get_umi_mapping(&counts);

        assert_eq!(
            mapping.get(b"AAGTGGATTT" as &[u8]),
            Some(&b"AAGTGTATTT".to_vec())
        );

        let expected = parse_counts("AAGTGGATTT:1;AAGTGGTTTT:45;AAGTGTATTT:88;AAGTGTTTTT:2;ATACAATTGG:1;ATCACCCCCG:4;ATCACTGTCG:1;ATCACTGTGG:138;ATGCCGTCAT:1;ATGCCGTCGT:4;ATGCCGTCTT:101;ATGCCGTGGT:3;ATGCCGTGTT:1;ATTCGCCCGC:14;GAATGGCGCG:5;GATATTCCTT:1;GGACCTAGGT:103;GTTAAACTAT:22;TACGTCAATG:1;TCTAACGCAT:21;TGTCAGATAA:195;TTATATAAAT:32");
        assert_eq!(apply_mapping(&counts, &mapping), expected);
    }

    #[test]
    fn directional_mapping_corrects_complete_tfseq_example() {
        let counts = parse_counts(TFSEQ_UMI_COUNTS);
        let mapping = get_directional_umi_mapping(&counts);
        let expected = parse_counts("AAGTGGTTTT:45;AAGTGTATTT:91;ATACAATTGG:1;ATCACCCCCG:4;ATCACTCTCG:1;ATCACTGTGG:138;ATGCCGACAT:1;ATGCCGTCTT:109;ATTCGCCCGC:14;GAATGGCGCG:5;GATATTCCTT:1;GGACCTAGGT:103;GTTAAACTAT:22;TACGTCAATG:1;TCTAACGCAT:21;TGTCAGATAA:195;TTATATAAAT:32");

        assert_eq!(apply_mapping(&counts, &mapping), expected);
    }
}
