use bed_utils::bed::BEDLike;
use noodles::sam;
use noodles::sam::alignment::{record::data::field::tag::Tag, Record};
use std::collections::{HashMap, HashSet};
use std::fmt::Display;
use std::ops::{Deref, DerefMut};
use anyhow::Result;

use crate::align::MultiMap;
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

#[derive(Debug, Default)]
struct AlignStat {
    total: u64, // Total number of reads
    mapped: u64, // Number of mapped reads
    high_quality: u64, // Number of high-quality mapped reads: unique, non-duplicate, and mapping quality >= 30
    multimapped: u64, // Number of reads with multiple alignments
    duplicate: u64, // Number of duplicate reads
}

impl AlignStat {
    pub fn add<R: Record>(&mut self, record: &MultiMap<R>) -> Result<()> {
        self.total += 1;
        let flags = record.primary.flags()?;
        if flags.is_duplicate() {
            self.duplicate += 1;
        }
        if !flags.is_unmapped() {
            self.mapped += 1;
            if record.others.is_some() {
                self.multimapped += 1;
            } else {
                let q = record.primary.mapping_quality().transpose()?.map(|x| x.get()).unwrap_or(60);
                if q >= 30 {
                    self.high_quality += 1;
                }
            }
        }
        Ok(())
    }

    pub fn combine(&mut self, other: &Self) {
        self.total += other.total;
        self.mapped += other.mapped;
        self.high_quality += other.high_quality;
        self.multimapped += other.multimapped;
        self.duplicate += other.duplicate;
    }
}

#[derive(Debug, Default)]
struct PairAlignStat {
    read1: AlignStat,
    read2: AlignStat,
    proper_pairs: u64,
}

impl PairAlignStat {
    fn total_reads(&self) -> u64 {
        self.read1.total + self.read2.total
    }

    fn total_pairs(&self) -> u64 {
        self.read2.total
    }

    fn total_mapped(&self) -> u64 {
        self.read1.mapped + self.read2.mapped
    }

    fn total_high_quality(&self) -> u64 {
        self.read1.high_quality + self.read2.high_quality
    }

    fn total_duplicate(&self) -> u64 {
        self.read1.duplicate + self.read2.duplicate
    }

    fn add_read1<R: Record>(&mut self, record: &MultiMap<R>) -> Result<()> {
        self.read1.add(record)
    }

    fn add_read2<R: Record>(&mut self, record: &MultiMap<R>) -> Result<()> {
        self.read2.add(record)
    }

    fn add_pair<R: Record>(&mut self, record1: &MultiMap<R>, record2: &MultiMap<R>) -> Result<()> {
        self.read1.add(record1)?;
        self.read2.add(record2)?;
        if record1.primary.flags()?.is_properly_segmented() {
            self.proper_pairs += 1;
        }
        Ok(())
    }

    fn combine(&mut self, other: &Self) {
        self.read1.combine(&other.read1);
        self.read2.combine(&other.read2);
        self.proper_pairs += other.proper_pairs;
    }
}

#[derive(Debug, Default)]
pub struct AlignQC {
    pub(crate) mito_dna: HashSet<usize>, // Mitochondrial DNA reference sequence IDs
    stat_all: PairAlignStat,
    stat_barcoded: PairAlignStat,
    stat_mito: PairAlignStat,
    pub(crate) num_read1_bases: u64,
    pub(crate) num_read1_q30_bases: u64,
    pub(crate) num_read2_bases: u64,
    pub(crate) num_read2_q30_bases: u64,
}

impl AlignQC {
    pub fn add_mito_dna(&mut self, mito_dna: usize) {
        self.mito_dna.insert(mito_dna);
    }

    pub fn add_pair<R: Record>(&mut self, header: &sam::Header, record1: &MultiMap<R>, record2: &MultiMap<R>) -> Result<()> {
        let mut stat= PairAlignStat::default();

        self.num_read1_bases += record1.primary.sequence().len() as u64;
        self.num_read2_bases += record2.primary.sequence().len() as u64;

        self.num_read1_q30_bases += record1.primary
            .quality_scores()
            .iter()
            .filter(|s| s.as_ref().map(|x| *x >= 30).unwrap_or(false))
            .count() as u64;
        self.num_read2_q30_bases += record2.primary
            .quality_scores()
            .iter()
            .filter(|s| s.as_ref().map(|x| *x >= 30).unwrap_or(false))
            .count() as u64;

        stat.add_pair(record1, record2)?;

        self.stat_all.combine(&stat);
 
        if record1.primary
            .data()
            .get(&Tag::CELL_BARCODE_ID)
            .transpose()
            .unwrap()
            .is_some()
        {
            self.stat_barcoded.combine(&stat);
            if let Some(rid) = record1.primary.reference_sequence_id(header) {
                if self.mito_dna.contains(&rid.unwrap()) {
                    self.stat_mito.combine(&stat);
                }
            }
        }
        Ok(())
    }

    pub fn add_read1<R: Record>(&mut self, header: &sam::Header, record: &MultiMap<R>) -> Result<()> {
        let mut stat= PairAlignStat::default();

        self.num_read1_bases += record.primary.sequence().len() as u64;
        self.num_read1_q30_bases += record.primary
            .quality_scores()
            .iter()
            .filter(|s| s.as_ref().map(|x| *x >= 30).unwrap_or(false))
            .count() as u64;
        stat.add_read1(record)?;

        self.stat_all.combine(&stat);
 
        if record.primary
            .data()
            .get(&Tag::CELL_BARCODE_ID)
            .transpose()
            .unwrap()
            .is_some()
        {
            self.stat_barcoded.combine(&stat);
            if let Some(rid) = record.primary.reference_sequence_id(header) {
                if self.mito_dna.contains(&rid.unwrap()) {
                    self.stat_mito.combine(&stat);
                }
            }
        }
        Ok(())
    }

    pub fn add_read2<R: Record>(&mut self, header: &sam::Header, record: &MultiMap<R>) -> Result<()> {
        let mut stat= PairAlignStat::default();

        self.num_read2_bases += record.primary.sequence().len() as u64;
        self.num_read2_q30_bases += record.primary
            .quality_scores()
            .iter()
            .filter(|s| s.as_ref().map(|x| *x >= 30).unwrap_or(false))
            .count() as u64;
        stat.add_read2(record)?;

        self.stat_all.combine(&stat);
 
        if record.primary
            .data()
            .get(&Tag::CELL_BARCODE_ID)
            .transpose()
            .unwrap()
            .is_some()
        {
            self.stat_barcoded.combine(&stat);
            if let Some(rid) = record.primary.reference_sequence_id(header) {
                if self.mito_dna.contains(&rid.unwrap()) {
                    self.stat_mito.combine(&stat);
                }
            }
        }
        Ok(())
    }

    pub fn report(&self, metric: &mut Metrics) {
        let stat_all = &self.stat_all;
        let stat_barcoded = &self.stat_barcoded;

        let fraction_unmapped = 1.0 - stat_barcoded.total_mapped() as f64 / stat_barcoded.total_reads() as f64;
        let valid_barcode = stat_barcoded.total_reads() as f64 / stat_all.total_reads() as f64;
        let fraction_confidently_mapped = stat_barcoded.total_high_quality() as f64 / stat_barcoded.total_reads() as f64;
        let fraction_nonnuclear = self.stat_mito.total_reads() as f64 / stat_barcoded.total_reads() as f64;

        metric.insert("sequenced_reads".to_string(), stat_all.total_reads() as f64);
        metric.insert("sequenced_read_pairs".to_string(), stat_all.total_pairs() as f64);
        if self.num_read1_bases > 0 {
            metric.insert(
                "frac_q30_bases_read1".to_string(),
                self.num_read1_q30_bases as f64 / self.num_read1_bases as f64,
            );
        }
        if self.num_read2_bases > 0 {
            metric.insert(
                "frac_q30_bases_read2".to_string(),
                self.num_read2_q30_bases as f64 / self.num_read2_bases as f64,
            );
        }
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
        if !self.mito_dna.contains(fragment.chrom()) {
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
