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

impl<'a> BarcodeDupMarker<'a> {
    pub fn new(
        mut umigene_counts: HashMap<(&'a [u8], Gene), u64>,
        mut umigene_min_key: HashMap<(&'a [u8], Gene), UmiSelectKey>,
        filter_umis: bool,
        targeted_umi_min_read_count: Option<u64>,
    ) -> Self {
        // Determine which UMIs need to be corrected
        let umi_corrections: HashMap<(UmiSeq, Gene), UmiSeq> = match umi_correction {
            UmiCorrection::Enable => correct_umis(&umigene_counts),
            UmiCorrection::Disable => HashMap::new(),
        };

        let umi_correction_counts = umi_corrections
            .iter()
            .map(|(raw_key, corrected_umi)| {
                (
                    raw_key,
                    (*corrected_umi, raw_key.1),
                    umigene_counts[raw_key],
                )
            })
            .collect::<Vec<_>>();

        // Before determining low-support UMI-genes, count one read of each corrected UMI
        // to match the behaviour of Cell Ranger 3.
        for (raw_key, corrected_key, _raw_count) in &umi_correction_counts {
            // One read has been counted before determining low-support UMI-genes.
            *umigene_counts.get_mut(raw_key).unwrap() -= 1;
            *umigene_counts.get_mut(corrected_key).unwrap() += 1;
        }

        // Determine low-support UMI-genes.
        let low_support_umigenes = if filter_umis {
            determine_low_support_umigenes(&umigene_counts)
        } else {
            HashSet::new()
        };

        // After determining low-support UMI-genes, count the remaining reads of each corrected UMI.
        for (raw_key, corrected_key, raw_count) in &umi_correction_counts {
            // One read has already been counted before determining low-support UMI-genes.
            *umigene_counts.get_mut(raw_key).unwrap() -= raw_count - 1;
            *umigene_counts.get_mut(corrected_key).unwrap() += raw_count - 1;
        }

        // Which is the lowest raw UMI lexicographicaly that would be corrected to (UmiSeq, Gene)
        // and would potentially be a "UMI count".
        let mut min_raw_umis: HashMap<(UmiSeq, Gene), UmiSeq> = HashMap::new();
        for ((raw_seq, gene), corr_seq) in &umi_corrections {
            if raw_seq < corr_seq || umi_corrections.contains_key(&(*corr_seq, *gene)) {
                let min_umi = match min_raw_umis.remove(&(*corr_seq, *gene)) {
                    Some(prev_min) => prev_min.min(*raw_seq),
                    None => *raw_seq,
                };
                min_raw_umis.insert((*corr_seq, *gene), min_umi);
            }
        }
        // The min UmiSelectKey after correction
        let mut min_umi_key_corrections = HashMap::new();
        for ((corr_seq, gene), raw_seq) in min_raw_umis {
            min_umi_key_corrections
                .insert((corr_seq, gene), umigene_min_key[&(raw_seq, gene)].clone());
        }
        for (key, umi_key) in min_umi_key_corrections {
            umigene_min_key.insert(key, umi_key.clone());
        }

        Self {
            umigene_counts,
            umi_corrections,
            low_support_umigenes,
            umigene_min_key,
            targeted_umi_min_read_count,
        }
    }
}