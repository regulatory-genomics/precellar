use bed_utils::bed::BEDLike;
use noodles::sam;
use noodles::sam::alignment::{record::data::field::tag::Tag, Record};
use std::collections::{HashMap, HashSet};
use std::fmt::Display;
use std::ops::{Deref, DerefMut};

use crate::fragment::Fragment;

#[derive(Debug, Default, Clone)]
pub struct Metrics(HashMap<String, f64>);

impl From<HashMap<String, f64>> for Metrics {
    fn from(map: HashMap<String, f64>) -> Self {
        Metrics(map)
    }
}

impl From<Metrics> for HashMap<String, f64> {
    fn from(val: Metrics) -> Self {
        val.0
    }
}

impl Deref for Metrics {
    type Target = HashMap<String, f64>;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl DerefMut for Metrics {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

impl Display for Metrics {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for (key, value) in &self.0 {
            writeln!(f, "{}\t{}", key, value)?;
        }
        Ok(())
    }
}

/// Alignment record statistics.
#[derive(Debug, Default)]
pub struct FlagStat {
    pub read: u64,
    pub primary: u64,
    pub secondary: u64,
    pub supplementary: u64,
    pub duplicate: u64,
    pub primary_duplicate: u64,
    pub mapped: u64,
    pub primary_mapped: u64,
    pub paired: u64,
    pub read_1: u64,
    pub read_2: u64,
    pub proper_pair: u64,
    pub mate_mapped: u64,
    pub singleton: u64,
    pub mate_reference_sequence_id_mismatch: u64,
}

impl FlagStat {
    pub fn add(&mut self, other: &FlagStat) {
        self.read += other.read;
        self.primary += other.primary;
        self.secondary += other.secondary;
        self.supplementary += other.supplementary;
        self.duplicate += other.duplicate;
        self.primary_duplicate += other.primary_duplicate;
        self.mapped += other.mapped;
        self.primary_mapped += other.primary_mapped;
        self.paired += other.paired;
        self.read_1 += other.read_1;
        self.read_2 += other.read_2;
        self.proper_pair += other.proper_pair;
        self.mate_mapped += other.mate_mapped;
        self.singleton += other.singleton;
        self.mate_reference_sequence_id_mismatch += other.mate_reference_sequence_id_mismatch;
    }

    pub fn update<R: Record>(&mut self, header: &sam::Header, record: &R) {
        self.read += 1;
        let flags = record.flags().unwrap();

        if flags.is_duplicate() {
            self.duplicate += 1;
        }

        if !flags.is_unmapped() {
            self.mapped += 1;
        }

        if flags.is_secondary() {
            self.secondary += 1;
        } else if flags.is_supplementary() {
            self.supplementary += 1;
        } else {
            self.primary += 1;

            if !flags.is_unmapped() {
                self.primary_mapped += 1;
            }

            if flags.is_duplicate() {
                self.primary_duplicate += 1;
            }

            if flags.is_segmented() {
                self.paired += 1;

                if flags.is_first_segment() {
                    self.read_1 += 1;
                }

                if flags.is_last_segment() {
                    self.read_2 += 1;
                }

                if !flags.is_unmapped() {
                    if flags.is_properly_segmented() {
                        self.proper_pair += 1;
                    }

                    if flags.is_mate_unmapped() {
                        self.singleton += 1;
                    } else {
                        self.mate_mapped += 1;
                        let rec_id = record.mate_reference_sequence_id(header).unwrap().unwrap();
                        let mat_id = record.reference_sequence_id(header).unwrap().unwrap();

                        if mat_id != rec_id {
                            self.mate_reference_sequence_id_mismatch += 1;
                        }
                    }
                }
            }
        }
    }
}

#[derive(Debug, Default)]
pub struct AlignQC {
    pub(crate) mito_dna: HashSet<usize>, // Mitochondrial DNA reference sequence IDs
    pub(crate) all_reads_flagstat: FlagStat,
    pub(crate) barcoded_reads_flagstat: FlagStat,
    pub(crate) hq_flagstat: FlagStat,
    pub(crate) mito_flagstat: FlagStat,
    pub(crate) num_read1_bases: u64,
    pub(crate) num_read1_q30_bases: u64,
    pub(crate) num_read2_bases: u64,
    pub(crate) num_read2_q30_bases: u64,
}

impl AlignQC {
    pub fn add_mito_dna(&mut self, mito_dna: usize) {
        self.mito_dna.insert(mito_dna);
    }

    pub fn update<R: Record>(&mut self, record: &R, header: &sam::Header) {
        let mut flagstat = FlagStat::default();
        flagstat.update(header, record);
        if flagstat.paired == 1 && flagstat.read_2 == 1 {
            self.num_read2_bases += record.sequence().len() as u64;
            self.num_read2_q30_bases += record
                .quality_scores()
                .as_ref()
                .iter()
                .filter(|x| *x >= 30)
                .count() as u64;
        } else {
            self.num_read1_bases += record.sequence().len() as u64;
            self.num_read1_q30_bases += record
                .quality_scores()
                .as_ref()
                .iter()
                .filter(|x| *x >= 30)
                .count() as u64;
        }

        self.all_reads_flagstat.add(&flagstat);
        let is_hq = record
            .mapping_quality()
            .map_or(true, |x| x.unwrap().get() >= 30);
        if is_hq {
            self.hq_flagstat.add(&flagstat);
        }

        if record
            .data()
            .get(&Tag::CELL_BARCODE_ID)
            .transpose()
            .unwrap()
            .is_some()
        {
            self.barcoded_reads_flagstat.add(&flagstat);
            if let Some(rid) = record.reference_sequence_id(header) {
                if is_hq && self.mito_dna.contains(&rid.unwrap()) {
                    self.mito_flagstat.add(&flagstat);
                }
            }
        }
    }

    pub fn report(&self, metric: &mut Metrics) {
        let flagstat_all = &self.all_reads_flagstat;
        let flagstat_barcoded = &self.barcoded_reads_flagstat;
        let num_reads = flagstat_all.read;
        let num_pairs = flagstat_all.paired / 2;
        let num_barcoded_reads = flagstat_barcoded.read;
        let num_barcoded_pairs = flagstat_barcoded.paired / 2;
        let mapped_pairs = flagstat_barcoded.mate_mapped / 2;
        let is_paired = num_pairs > 0;

        let fraction_unmapped = if is_paired {
            1.0 - mapped_pairs as f64 / num_barcoded_pairs as f64
        } else {
            1.0 - flagstat_barcoded.mapped as f64 / num_barcoded_reads as f64
        };
        let valid_barcode = if is_paired {
            num_barcoded_pairs as f64 / num_pairs as f64
        } else {
            num_barcoded_reads as f64 / num_reads as f64
        };
        let fraction_confidently_mapped = if is_paired {
            (self.hq_flagstat.paired / 2) as f64 / num_pairs as f64
        } else {
            self.hq_flagstat.read as f64 / num_reads as f64
        };
        let fraction_nonnuclear = if is_paired {
            (self.mito_flagstat.paired / 2) as f64 / num_pairs as f64
        } else {
            self.mito_flagstat.read as f64 / num_reads as f64
        };

        metric.insert("sequenced_reads".to_string(), num_reads as f64);
        metric.insert("sequenced_read_pairs".to_string(), num_pairs as f64);
        metric.insert(
            "frac_q30_bases_read1".to_string(),
            self.num_read1_q30_bases as f64 / self.num_read1_bases as f64,
        );
        metric.insert(
            "frac_q30_bases_read2".to_string(),
            self.num_read2_q30_bases as f64 / self.num_read2_bases as f64,
        );
        metric.insert(
            "frac_confidently_mapped".to_string(),
            fraction_confidently_mapped,
        );
        metric.insert("frac_unmapped".to_string(), fraction_unmapped);
        metric.insert("frac_valid_barcode".to_string(), valid_barcode);
        metric.insert("frac_nonnuclear".to_string(), fraction_nonnuclear);
    }
}

#[derive(Debug, Clone, Default)]
pub struct FragmentQC {
    mito_dna: HashSet<String>,
    num_pcr_duplicates: u64,
    num_unique_fragments: u64,
    num_frag_nfr: u64,
    num_frag_single: u64,
}

impl FragmentQC {
    pub fn add_mito_dna<S: Into<String>>(&mut self, mito_dna: S) {
        self.mito_dna.insert(mito_dna.into());
    }

    pub fn update(&mut self, fragment: &Fragment) {
        self.num_pcr_duplicates += fragment.count as u64 - 1;
        self.num_unique_fragments += 1;
        let size = fragment.len();
        if self.mito_dna.contains(fragment.chrom()) {
            if size < 147 {
                self.num_frag_nfr += 1;
            } else if size <= 294 {
                self.num_frag_single += 1;
            }
        }
    }

    pub fn report(&self, metric: &mut Metrics) {
        metric.insert(
            "frac_duplicates".to_string(),
            self.num_pcr_duplicates as f64
                / (self.num_unique_fragments + self.num_pcr_duplicates) as f64,
        );
        metric.insert(
            "frac_fragment_in_nucleosome_free_region".to_string(),
            self.num_frag_nfr as f64 / self.num_unique_fragments as f64,
        );
        metric.insert(
            "frac_fragment_flanking_single_nucleosome".to_string(),
            self.num_frag_single as f64 / self.num_unique_fragments as f64,
        );
    }
}
