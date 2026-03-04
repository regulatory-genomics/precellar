use indexmap::IndexSet;

use super::sequence_aligner::fitting_alignment_distance;

/// K-mer size for indexing (6-mer)
const KMER_SIZE: usize = 6;
/// Number of possible k-mers (4^6 = 4096)
const NUM_KMERS: usize = 4096;
/// Minimum number of top candidates to consider
const MIN_TOP_CANDIDATES: usize = 50;
/// Minimum confidence threshold for barcode matching
const MIN_CONFIDENCE: f64 = 0.7;

/// Pre-built index for fast barcode candidate filtering using k-mer voting.
///
/// Instead of comparing a candidate sequence against all barcodes in the whitelist,
/// this index uses 6-mer voting to quickly identify the most likely matches,
/// then only runs expensive fitting alignment on the top candidates.
#[derive(Debug, Clone)]
pub struct BarcodeIndex {
    /// Original whitelist barcodes
    barcodes: Vec<Vec<u8>>,
    /// Inverted index: kmer_hash -> list of barcode indices containing this k-mer
    kmer_to_barcodes: Vec<Vec<u32>>,
}

impl BarcodeIndex {
    /// Build index from whitelist
    pub fn new(whitelist: &IndexSet<Vec<u8>>) -> Self {
        let barcodes: Vec<Vec<u8>> = whitelist.iter().cloned().collect();
        let mut kmer_to_barcodes: Vec<Vec<u32>> = vec![Vec::new(); NUM_KMERS];

        // Build inverted index
        for (barcode_idx, barcode) in barcodes.iter().enumerate() {
            if barcode.len() < KMER_SIZE {
                continue;
            }

            for i in 0..=(barcode.len() - KMER_SIZE) {
                if let Some(hash) = encode_kmer(&barcode[i..i + KMER_SIZE]) {
                    kmer_to_barcodes[hash as usize].push(barcode_idx as u32);
                }
            }
        }

        Self {
            barcodes,
            kmer_to_barcodes,
        }
    }

    /// Find best barcode match using k-mer voting followed by fitting alignment.
    ///
    /// Returns the matched barcode and confidence score, or None if no good match found.
    pub fn find_best_match(&self, candidate: &[u8]) -> Option<(Vec<u8>, f64)> {
        if self.barcodes.is_empty() {
            return None;
        }

        // For small whitelists, skip k-mer voting and do direct comparison
        if self.barcodes.len() <= MIN_TOP_CANDIDATES || candidate.len() < KMER_SIZE {
            return self.find_best_match_linear(candidate);
        }

        // Vote counting: barcode_idx -> vote count
        let mut votes: Vec<u32> = vec![0; self.barcodes.len()];

        // Extract k-mers from candidate and vote
        for i in 0..=(candidate.len() - KMER_SIZE) {
            if let Some(hash) = encode_kmer(&candidate[i..i + KMER_SIZE]) {
                for &barcode_idx in &self.kmer_to_barcodes[hash as usize] {
                    votes[barcode_idx as usize] += 1;
                }
            }
        }

        // Select top candidates with ties
        let top_candidates = self.select_top_candidates_with_ties(&votes);

        // If no candidates from k-mer voting, fall back to linear search
        if top_candidates.is_empty() {
            return self.find_best_match_linear(candidate);
        }

        // Run fitting alignment only on top candidates
        let mut best_barcode_idx = None;
        let mut min_distance = usize::MAX;

        for barcode_idx in top_candidates {
            let barcode = &self.barcodes[barcode_idx];
            let distance = fitting_alignment_distance(barcode, candidate);
            if distance < min_distance {
                min_distance = distance;
                best_barcode_idx = Some(barcode_idx);
            }
        }

        // Convert to result with confidence
        best_barcode_idx.and_then(|idx| {
            let barcode = &self.barcodes[idx];
            let barcode_length = barcode.len().max(1);
            let confidence = 1.0 - (min_distance as f64 / barcode_length as f64);

            if confidence >= MIN_CONFIDENCE {
                Some((barcode.clone(), confidence))
            } else {
                None
            }
        })
    }

    /// Linear search fallback for small whitelists or when k-mer voting fails
    fn find_best_match_linear(&self, candidate: &[u8]) -> Option<(Vec<u8>, f64)> {
        let mut best_barcode_idx = None;
        let mut min_distance = usize::MAX;

        for (idx, barcode) in self.barcodes.iter().enumerate() {
            let distance = fitting_alignment_distance(barcode, candidate);
            if distance < min_distance {
                min_distance = distance;
                best_barcode_idx = Some(idx);
            }
        }

        best_barcode_idx.and_then(|idx| {
            let barcode = &self.barcodes[idx];
            let barcode_length = barcode.len().max(1);
            let confidence = 1.0 - (min_distance as f64 / barcode_length as f64);

            if confidence >= MIN_CONFIDENCE {
                Some((barcode.clone(), confidence))
            } else {
                None
            }
        })
    }

    /// Select top candidates with ties.
    ///
    /// Adds candidates until count > MIN_TOP_CANDIDATES, including all barcodes
    /// with the same vote count as the cutoff.
    fn select_top_candidates_with_ties(&self, votes: &[u32]) -> Vec<usize> {
        // Collect candidates with non-zero votes
        let mut candidates: Vec<(usize, u32)> = votes
            .iter()
            .enumerate()
            .filter(|(_, &v)| v > 0)
            .map(|(i, &v)| (i, v))
            .collect();

        if candidates.is_empty() {
            return Vec::new();
        }

        // Sort by vote count descending
        candidates.sort_unstable_by(|a, b| b.1.cmp(&a.1));

        // Find cutoff vote count
        let mut result = Vec::new();
        let mut cutoff_votes = 0u32;

        for (idx, vote_count) in candidates {
            if result.len() < MIN_TOP_CANDIDATES {
                // Still under minimum, add this candidate
                result.push(idx);
                cutoff_votes = vote_count;
            } else if vote_count == cutoff_votes {
                // Same vote count as cutoff, include it (tie)
                result.push(idx);
            } else {
                // Lower vote count and we have enough candidates, stop
                break;
            }
        }

        result
    }

    /// Get the number of barcodes in the index
    pub fn len(&self) -> usize {
        self.barcodes.len()
    }

    /// Check if the index is empty
    pub fn is_empty(&self) -> bool {
        self.barcodes.is_empty()
    }
}

/// Encode a 6-mer as a 12-bit integer.
///
/// Each base is encoded as 2 bits: A=0, C=1, G=2, T=3.
/// Returns None if sequence contains non-ACGT characters.
#[inline]
fn encode_kmer(kmer: &[u8]) -> Option<u16> {
    debug_assert_eq!(kmer.len(), KMER_SIZE);

    let mut hash: u16 = 0;
    for &base in kmer {
        let bits = match base {
            b'A' | b'a' => 0,
            b'C' | b'c' => 1,
            b'G' | b'g' => 2,
            b'T' | b't' => 3,
            _ => return None, // Skip k-mers with N or other chars
        };
        hash = (hash << 2) | bits;
    }
    Some(hash)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_encode_kmer() {
        // Test basic encoding
        assert_eq!(encode_kmer(b"AAAAAA"), Some(0b000000000000)); // 0
        assert_eq!(encode_kmer(b"CCCCCC"), Some(0b010101010101)); // 1365
        assert_eq!(encode_kmer(b"GGGGGG"), Some(0b101010101010)); // 2730
        assert_eq!(encode_kmer(b"TTTTTT"), Some(0b111111111111)); // 4095

        // Test mixed sequence
        assert_eq!(encode_kmer(b"ACGTAC"), Some(0b000110110001)); // A=00, C=01, G=10, T=11, A=00, C=01

        // Test with N (should return None)
        assert_eq!(encode_kmer(b"ACNTGA"), None);

        // Test lowercase
        assert_eq!(encode_kmer(b"acgtac"), Some(0b000110110001));
    }

    #[test]
    fn test_barcode_index_basic() {
        let mut whitelist = IndexSet::new();
        whitelist.insert(b"AAAAAAAAAA".to_vec()); // 10bp barcode
        whitelist.insert(b"CCCCCCCCCC".to_vec());
        whitelist.insert(b"GGGGGGGGGG".to_vec());
        whitelist.insert(b"TTTTTTTTTT".to_vec());

        let index = BarcodeIndex::new(&whitelist);
        assert_eq!(index.len(), 4);

        // Test exact match
        let result = index.find_best_match(b"AAAAAAAAAA");
        assert!(result.is_some());
        let (barcode, confidence) = result.unwrap();
        assert_eq!(barcode, b"AAAAAAAAAA".to_vec());
        assert_eq!(confidence, 1.0);

        // Test with 1 mismatch
        let result = index.find_best_match(b"AAAAAAAAAT");
        assert!(result.is_some());
        let (barcode, confidence) = result.unwrap();
        assert_eq!(barcode, b"AAAAAAAAAA".to_vec());
        assert!(confidence >= 0.7);
    }

    #[test]
    fn test_barcode_index_with_indels() {
        let mut whitelist = IndexSet::new();
        whitelist.insert(b"ACGTACGTAC".to_vec()); // 10bp barcode

        let index = BarcodeIndex::new(&whitelist);

        // Test with insertion in candidate
        let result = index.find_best_match(b"ACGTAACGTAC"); // extra A
        assert!(result.is_some());
        let (barcode, _) = result.unwrap();
        assert_eq!(barcode, b"ACGTACGTAC".to_vec());

        // Test with deletion in candidate
        let result = index.find_best_match(b"ACGTCGTAC"); // missing A
        assert!(result.is_some());
        let (barcode, _) = result.unwrap();
        assert_eq!(barcode, b"ACGTACGTAC".to_vec());
    }

    #[test]
    fn test_select_top_candidates_with_ties() {
        let mut whitelist = IndexSet::new();
        for i in 0..100 {
            whitelist.insert(format!("ACGT{:06}", i).into_bytes());
        }

        let index = BarcodeIndex::new(&whitelist);

        // Create votes with ties
        let mut votes = vec![0u32; 100];
        // 10 candidates with 5 votes each
        for i in 0..10 {
            votes[i] = 5;
        }
        // 50 candidates with 3 votes each
        for i in 10..60 {
            votes[i] = 3;
        }
        // Rest with 1 vote
        for i in 60..100 {
            votes[i] = 1;
        }

        let top = index.select_top_candidates_with_ties(&votes);

        // Should include all 10 with 5 votes + all 50 with 3 votes (ties at cutoff)
        assert_eq!(top.len(), 60);
    }

    #[test]
    fn test_empty_whitelist() {
        let whitelist = IndexSet::new();
        let index = BarcodeIndex::new(&whitelist);

        assert!(index.is_empty());
        assert_eq!(index.find_best_match(b"ACGTACGTAC"), None);
    }
}
