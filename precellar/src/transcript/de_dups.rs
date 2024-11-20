use std::collections::{HashMap, HashSet};

type Gene = usize;

/// Within each gene, correct Hamming-distance-one UMIs
fn correct_umis<'a>(umigene_counts: &'a HashMap<(&'a [u8], Gene), u64>) -> HashMap<(&'a [u8], Gene), Vec<u8>> {
    let nucs = b"ACGT";

    let mut corrections = HashMap::new();

    for ((umi, gene), orig_count) in umigene_counts {
        let mut test_umi = umi.to_vec();

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
                let test_count = *umigene_counts.get(&(&test_umi, *gene)).unwrap_or(&0u64);

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
            corrections.insert((*umi, *gene), best_dest_umi);
        }
    }
    corrections
}

/// Do duplicate marking on a single barcode's worth of data
/// 1) Correct UMI sequences
/// 2) Mark records as low-support (putatively chimeric) UMIs
/// 3) Mark records as PCR duplicates
/// 4) update ReadAnnotation with DupInfo that will be used for metrics and BAM tags
pub struct BarcodeDupMarker<'a> {
    umigene_counts: HashMap<(&'a [u8], Gene), u64>,
    low_support_umigenes: HashSet<(&'a [u8], Gene)>,
    umi_corrections: HashMap<(&'a [u8], Gene), Vec<u8>>,
    umigene_min_key: HashMap<(&'a [u8], Gene), UmiSelectKey>,
    /// If this is Some(r), we filter out (UMI, genes) pairs
    /// with **less than** r reads and not include them in UMI counts.
    targeted_umi_min_read_count: Option<u64>,
}