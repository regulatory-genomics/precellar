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
/// A higher-count UMI has a directional edge to each one-edit neighbor satisfying
/// `parent_count >= 2 * child_count - 1`. Roots are processed by descending abundance,
/// and each reachable UMI is assigned to the first root that reaches it. Equal-count,
/// empty, non-ACGT, and zero-count UMIs are left uncorrected.
pub(crate) fn get_directional_umi_mapping(
    umi_count: &HashMap<Vec<u8>, u64>,
) -> HashMap<Vec<u8>, Vec<u8>> {
    let nucs = b"ACGT";
    let mut children: HashMap<Vec<u8>, Vec<Vec<u8>>> = HashMap::new();

    for (umi, child_count) in umi_count {
        if *child_count > 0 && is_canonical_umi(umi) {
            let mut test_umi = umi.clone();

            for pos in 0..umi.len() {
                for test_char in nucs {
                    if *test_char != umi[pos] {
                        test_umi[pos] = *test_char;

                        if let Some(&parent_count) = umi_count.get(&test_umi) {
                            if is_directional_parent(parent_count, *child_count) {
                                children
                                    .entry(test_umi.clone())
                                    .or_default()
                                    .push(umi.clone());
                            }
                        }
                    }
                }
                test_umi[pos] = umi[pos];
            }
        }
    }

    let mut roots = umi_count.keys().collect::<Vec<_>>();
    roots.sort_unstable_by(|a, b| umi_count[*b].cmp(&umi_count[*a]).then_with(|| b.cmp(a)));

    let mut assigned = HashMap::new();
    for root in roots {
        if !assigned.contains_key(root) {
            let mut stack = vec![root.as_slice()];
            while let Some(umi) = stack.pop() {
                if !assigned.contains_key(umi) {
                    assigned.insert(umi.to_vec(), root.clone());
                    if let Some(next) = children.get(umi) {
                        stack.extend(next.iter().map(Vec::as_slice));
                    }
                }
            }
        }
    }

    let mut corrections = HashMap::new();
    for (umi, root) in assigned {
        if umi != root {
            corrections.insert(umi, root);
        }
    }
    corrections
}

fn is_canonical_umi(umi: &[u8]) -> bool {
    !umi.is_empty()
        && umi
            .iter()
            .all(|base| matches!(base, b'A' | b'C' | b'G' | b'T'))
}

fn is_directional_parent(parent_count: u64, child_count: u64) -> bool {
    parent_count > child_count && child_count <= parent_count / 2 + parent_count % 2
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
    fn directional_mapping_uses_highest_reachable_root() {
        let entries = [
            (b"AAAA".to_vec(), 3),
            (b"AAAC".to_vec(), 6),
            (b"AAAG".to_vec(), 5),
            (b"AATG".to_vec(), 10),
        ];
        let counts = HashMap::from(entries.clone());
        let mapping = get_directional_umi_mapping(&counts);

        assert_eq!(mapping.get(b"AAAA" as &[u8]), Some(&b"AATG".to_vec()));
        assert_eq!(mapping.get(b"AAAG" as &[u8]), Some(&b"AATG".to_vec()));
        assert!(!mapping.contains_key(b"AAAC" as &[u8]));

        let reversed = entries.into_iter().rev().collect::<HashMap<_, _>>();
        assert_eq!(mapping, get_directional_umi_mapping(&reversed));
    }

    #[test]
    fn directional_mapping_enforces_abundance_threshold() {
        let counts = HashMap::from([
            (b"AAAA".to_vec(), 2),
            (b"AAAT".to_vec(), 3),
            (b"CCCC".to_vec(), 3),
            (b"CCCT".to_vec(), 4),
        ]);
        let mapping = get_directional_umi_mapping(&counts);

        assert_eq!(mapping.get(b"AAAA" as &[u8]), Some(&b"AAAT".to_vec()));
        assert!(!mapping.contains_key(b"CCCC" as &[u8]));
    }

    #[test]
    fn directional_mapping_leaves_invalid_umis_uncorrected() {
        let counts = HashMap::from([
            (Vec::new(), 1),
            (b"NAAA".to_vec(), 1),
            (b"aaaa".to_vec(), 1),
            (b"AAAA".to_vec(), 10),
            (b"CCCC".to_vec(), 0),
            (b"CCCT".to_vec(), 10),
        ]);
        let mapping = get_directional_umi_mapping(&counts);

        assert!(!mapping.contains_key(b"" as &[u8]));
        assert!(!mapping.contains_key(b"NAAA" as &[u8]));
        assert!(!mapping.contains_key(b"aaaa" as &[u8]));
        assert!(!mapping.contains_key(b"CCCC" as &[u8]));
    }

    #[test]
    fn directional_mapping_handles_threshold_without_overflow() {
        let counts = HashMap::from([
            (b"AAAA".to_vec(), u64::MAX / 2 + 2),
            (b"AAAT".to_vec(), u64::MAX),
        ]);

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
        let expected = parse_counts("AAGTGGTTTT:41;AAGTGTATTT:95;ATACAATTGG:1;ATCACCCCCG:4;ATCACTCTCG:1;ATCACTGTGG:138;ATGCCGACAT:1;ATGCCGTCTT:109;ATTCGCCCGC:14;GAATGGCGCG:5;GATATTCCTT:1;GGACCTAGGT:103;GTTAAACTAT:22;TACGTCAATG:1;TCTAACGCAT:21;TGTCAGATAA:195;TTATATAAAT:32");

        assert_eq!(apply_mapping(&counts, &mapping), expected);
        assert!(mapping
            .values()
            .all(|destination| !mapping.contains_key(destination)));
    }
}
