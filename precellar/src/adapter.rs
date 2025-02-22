/// Trim polyA tail from the 3' end of the sequence.
/// The aim of the --poly-A trimming algorithm is to find a suffix of the read that contains a high number of A nucleotides. Conceptually, we consider all possible suffixes of the read. For each suffix, we count how many A and non-A nucleotides it contains. We then exclude all suffixes from consideration that have more than 20% non-A because we assume these are not actually poly-A tails. For each remaining suffix, we compute a score: Non-A nucleotides get -2 and A nucleotides get +1, which we add up to get the score for the suffix. Finally, we choose the suffix that maximizes that score and remove it from the read. Shorter suffixes win if there is a tie.
pub fn trim_poly_nucleotide(nucl: u8, seq: impl Iterator<Item = u8>) -> Option<usize> {
    let mut best_index = None;
    let mut best_score: i32 = 0;
    let mut score: i32 = 0;
    let mut errors: usize = 0;
    for (i, base) in seq.enumerate() {
        if base == nucl {
            score += 1;
        } else {
            score -= 2;
            errors += 1;
        }

        if score > best_score && 5 * errors <= i {
            best_index = Some(i);
            best_score = score;
        }
    }

    best_index.map(|i| i + 1)
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_trim() {
        assert_eq!(
            trim_poly_nucleotide(b'A', b"AAAAAA".iter().copied()),
            Some(6)
        );

        assert_eq!(
            trim_poly_nucleotide(b'A', b"ATCGTC".iter().copied()),
            Some(1),
        );

        assert_eq!(
            trim_poly_nucleotide(b'A', b"CTCGAA".iter().copied()),
            None,
        );
    }
}
