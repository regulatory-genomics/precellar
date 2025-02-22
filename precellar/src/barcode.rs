use anyhow::Result;
use core::f64;
use noodles::sam::alignment::{
    record::data::field::{Tag, Value},
    Record,
};
use std::{
    collections::HashMap,
    ops::{Deref, DerefMut},
};

const BC_MAX_QV: u8 = 66; // This is the illumina quality value
pub(crate) const BASE_OPTS: [u8; 4] = [b'A', b'C', b'G', b'T'];

/// A map of oligo species to their frequency in a given library.
#[derive(Debug, Clone)]
pub struct OligoFrequncy(HashMap<Vec<u8>, usize>);

impl Deref for OligoFrequncy {
    type Target = HashMap<Vec<u8>, usize>;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl DerefMut for OligoFrequncy {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

impl FromIterator<(Vec<u8>, usize)> for OligoFrequncy {
    fn from_iter<I: IntoIterator<Item = (Vec<u8>, usize)>>(iter: I) -> Self {
        Self(iter.into_iter().collect())
    }
}

// Implement default for OligoFrequncy
impl Default for OligoFrequncy {
    fn default() -> Self {
        Self::new()
    }
}

impl OligoFrequncy {
    pub fn new() -> Self {
        Self(HashMap::new())
    }

    /// The likelihood of a query oligo being generated by the library.
    /// If the query is present in the library, the likelihood is 1.0.
    /// Otherwise, the likelihood is calculated as
    pub fn likelihood<'a>(
        &'a self,
        query: &'a [u8],
        qual: &[u8],
        n_mismatch: usize,
    ) -> (&'a [u8], f64) {
        if n_mismatch == 0 {
            if self.0.contains_key(query) {
                (query, 1.0)
            } else {
                (query, 0.0)
            }
        } else if n_mismatch == 1 {
            self.likelihood1(query, qual)
        } else if n_mismatch == 2 {
            self.likelihood2(query, qual)
        } else {
            todo!()
        }
    }

    /// The likelihood up to 2 mismatches.
    fn likelihood2<'a>(&'a self, query: &'a [u8], qual: &[u8]) -> (&'a [u8], f64) {
        if self.0.contains_key(query) {
            return (query, 1.0);
        }

        let mut best_option = None;
        let mut total_likelihood = 0.0;
        let mut query_bytes = query.to_vec();

        // Single mismatch loop
        for (pos1, &qv1) in qual.iter().enumerate() {
            let qv1 = qv1.min(BC_MAX_QV);
            let original1 = query_bytes[pos1];

            for base1 in BASE_OPTS {
                if base1 != original1 {
                    query_bytes[pos1] = base1;

                    // Check for 1-mismatch barcode match
                    if let Some((key, raw_count)) = self.0.get_key_value(&query_bytes) {
                        let bc_count = 1 + raw_count;
                        let likelihood = bc_count as f64 * error_probability(qv1);
                        update_best_option(&mut best_option, likelihood, key);
                        total_likelihood += likelihood;
                    }

                    // Loop for the second mismatch
                    for (pos2, &qv2) in qual.iter().enumerate().skip(pos1 + 1) {
                        let qv2 = qv2.min(BC_MAX_QV);
                        let original2 = query_bytes[pos2];

                        for val2 in BASE_OPTS {
                            if val2 != original2 {
                                query_bytes[pos2] = val2;

                                // Check for 2-mismatch barcode match
                                if let Some((key, raw_count)) = self.0.get_key_value(&query_bytes) {
                                    let bc_count = 1 + raw_count;
                                    let likelihood = bc_count as f64
                                        * error_probability(qv1)
                                        * error_probability(qv2);
                                    update_best_option(&mut best_option, likelihood, key);
                                    total_likelihood += likelihood;
                                }
                            }
                        }
                        // Restore original value for second position
                        query_bytes[pos2] = original2;
                    }
                }
            }
            // Restore original value for first position
            query_bytes[pos1] = original1;
        }

        if let Some((best_like, best_bc)) = best_option {
            (best_bc, best_like / total_likelihood)
        } else {
            (query, 0.0)
        }
    }

    /// The likehood up to 1 mismatch.
    fn likelihood1<'a>(&'a self, query: &'a [u8], qual: &[u8]) -> (&'a [u8], f64) {
        if self.0.contains_key(query) {
            return (query, 1.0);
        }

        let mut best_option = None;
        let mut total_likelihood = 0.0;
        let mut query_bytes = query.to_vec();
        for (pos, &qv) in qual.iter().enumerate() {
            let qv = qv.min(BC_MAX_QV);
            let existing = query_bytes[pos];
            for val in BASE_OPTS {
                if val != existing {
                    query_bytes[pos] = val;
                    if let Some((key, raw_count)) = self.0.get_key_value(&query_bytes) {
                        let bc_count = 1 + raw_count;
                        let likelihood = bc_count as f64 * error_probability(qv);
                        update_best_option(&mut best_option, likelihood, key);
                        total_likelihood += likelihood;
                    }
                }
            }
            query_bytes[pos] = existing;
        }

        if let Some((best_like, best_bc)) = best_option {
            (best_bc, best_like / total_likelihood)
        } else {
            (query, 0.0)
        }
    }
}

// Helper function to update the best option
fn update_best_option<'a>(
    best_option: &mut Option<(f64, &'a [u8])>,
    likelihood: f64,
    key: &'a [u8],
) {
    match best_option {
        None => *best_option = Some((likelihood, key)),
        Some(ref old_best) => {
            if old_best.0 < likelihood {
                *best_option = Some((likelihood, key));
            }
        }
    }
}

#[derive(Debug)]
pub struct Whitelist {
    whitelist_exists: bool,
    barcode_counts: OligoFrequncy,
    mismatch_count: usize,
    pub(crate) total_count: usize,
    pub(crate) total_base_count: u64,
    q30_base_count: u64,
    base_qual_sum: i64,
}

impl Whitelist {
    pub fn empty() -> Self {
        Self {
            whitelist_exists: false,
            barcode_counts: OligoFrequncy::new(),
            mismatch_count: 0,
            total_count: 0,
            total_base_count: 0,
            q30_base_count: 0,
            base_qual_sum: 0,
        }
    }

    /// Create a new whitelist from an iterator of strings.
    pub fn new<I: IntoIterator<Item = S>, S: Into<Vec<u8>>>(iter: I) -> Self {
        let mut whitelist = Self::empty();
        whitelist.whitelist_exists = true;
        whitelist.barcode_counts = iter.into_iter().map(|x| (x.into(), 0)).collect();
        whitelist
    }

    /// Update the barcode counter with a barcode and its quality scores.
    pub fn count_barcode(&mut self, barcode: &[u8], barcode_qual: &[u8]) {
        if self.whitelist_exists {
            if let Some(count) = self.barcode_counts.get_mut(barcode) {
                *count += 1;
            } else {
                self.mismatch_count += 1;
            }
        } else if barcode.len() > 1 {
            *self.barcode_counts.entry(barcode.to_vec()).or_insert(0) += 1;
        } else {
            self.mismatch_count += 1;
        }

        self.total_count += 1;

        for &qual in barcode_qual {
            let qual_int = (qual as u32) - 33;
            self.base_qual_sum += qual_int as i64;
            if qual_int >= 30 {
                self.q30_base_count += 1;
            }
            self.total_base_count += 1;
        }
    }

    pub fn num_seen_barcodes(&self) -> usize {
        self.barcode_counts.values().filter(|&&x| x > 0).count()
    }

    pub fn get_barcode_counts(&self) -> &OligoFrequncy {
        &self.barcode_counts
    }

    pub fn mean_base_quality_score(&self) -> f64 {
        if self.total_base_count == 0 {
            // u64 never < 0
            0.0
        } else {
            self.base_qual_sum as f64 / self.total_base_count as f64
        }
    }

    pub fn frac_q30_bases(&self) -> f64 {
        if self.total_base_count == 0 {
            0.0
        } else {
            self.q30_base_count as f64 / self.total_base_count as f64
        }
    }

    pub fn frac_exact_match(&self) -> f64 {
        if self.total_count == 0 {
            0.0
        } else {
            1.0 - (self.mismatch_count as f64 / self.total_count as f64)
        }
    }
}

/// A barcode validator that uses a barcode counter to validate barcodes.
#[derive(Debug, Clone)]
pub struct BarcodeCorrector {
    /// threshold for sum of probability of error on barcode QVs. Barcodes exceeding
    /// this threshold will be marked as not valid.
    max_expected_errors: f64,
    /// if the posterior probability of a correction
    /// exceeds this threshold, the barcode will be corrected.
    bc_confidence_threshold: f64,
    /// The number of mismatches allowed in barcode
    max_mismatch: usize,
}

impl Default for BarcodeCorrector {
    fn default() -> Self {
        Self {
            max_expected_errors: f64::MAX,
            bc_confidence_threshold: 0.975,
            max_mismatch: 1,
        }
    }
}

impl BarcodeCorrector {
    pub fn with_bc_confidence_threshold(mut self, threshold: f64) -> Self {
        self.bc_confidence_threshold = threshold;
        self
    }

    pub fn with_max_missmatch(mut self, max_mismatch: usize) -> Self {
        self.max_mismatch = max_mismatch;
        self
    }
}

#[derive(Copy, Clone, Debug)]
pub enum BarcodeError {
    ExceedExpectedError(f64),
    LowConfidence(f64),
    NoMatch,
}

impl BarcodeCorrector {
    /// Determine if a barcode is valid. A barcode is valid if any of the following conditions are met:
    /// 1) It is in the whitelist and the number of expected errors is less than the max_expected_errors.
    /// 2) It is not in the whitelist, but the number of expected errors is less than the max_expected_errors and the corrected barcode is in the whitelist.
    /// 3) If the whitelist does not exist, the barcode is always valid.
    ///
    /// Return the corrected barcode
    pub fn correct<'a>(
        &'a self,
        barcode_counts: &'a OligoFrequncy,
        barcode: &'a [u8],
        qual: &[u8],
    ) -> Result<&[u8], BarcodeError> {
        let expected_errors: f64 = qual.iter().map(|&q| error_probability(q)).sum();
        if expected_errors >= self.max_expected_errors {
            return Err(BarcodeError::ExceedExpectedError(expected_errors));
        }

        let (bc, prob) = barcode_counts.likelihood(barcode, qual, self.max_mismatch);
        if prob <= 0.0 {
            Err(BarcodeError::NoMatch)
        } else if prob >= self.bc_confidence_threshold {
            Ok(bc)
        } else {
            Err(BarcodeError::LowConfidence(prob))
        }
    }
}

/// Convert Illumina quality scores to base-calling error probabilities, i.e.,
/// the probability of an incorrect base call.
#[inline(always)]
fn error_probability(qual: u8) -> f64 {
    let offset = 33.0; // Illumina quality score offset
    10f64.powf(-((qual as f64 - offset) / 10.0))
}

pub(crate) fn get_barcode<R: Record>(rec: &R) -> Result<Option<String>> {
    Ok(rec
        .data()
        .get(&Tag::CELL_BARCODE_ID)
        .transpose()?
        .and_then(|x| match x {
            Value::String(barcode) => Some(barcode.to_string()),
            _ => None,
        }))
}

pub(crate) fn get_umi<R: Record>(rec: &R) -> Result<Option<String>> {
    Ok(rec
        .data()
        .get(&Tag::UMI_SEQUENCE)
        .transpose()?
        .and_then(|x| match x {
            Value::String(umi) => Some(umi.to_string()),
            _ => None,
        }))
}
