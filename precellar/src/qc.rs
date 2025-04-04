use crate::align::{AnnotatedFastq, MultiMap};
use crate::fragment::Fragment;

use anyhow::Result;
use bed_utils::bed::BEDLike;
use noodles::sam;
use noodles::sam::alignment::{record::data::field::tag::Tag, Record};
use serde_json::{json, Value};
use std::collections::{HashMap, HashSet};

/// The trait for quality control metrics.
pub trait Metric: Sized + Extend<Self> {
    fn to_json(&self) -> Value;

    fn into_json(self) -> Value {
        self.to_json()
    }
}

#[derive(Debug, Default)]
pub struct QcFastq {
    pub(crate) num_reads: HashMap<String, u64>,
    pub(crate) num_defect: HashMap<String, u64>, // Number of reads with defects, e.g., misformed structure.
    pub(crate) num_q30_bases: HashMap<String, u64>,
    pub(crate) num_total_bases: HashMap<String, u64>,
}

impl QcFastq {
    pub fn update(&mut self, fq: &AnnotatedFastq) {
        if let Some(umi) = &fq.umi {
            *self.num_total_bases.entry("umi".to_string()).or_insert(0) +=
                umi.sequence().len() as u64;
            *self.num_q30_bases.entry("umi".to_string()).or_insert(0) +=
                umi.quality_scores()
                    .iter()
                    .filter(|s| **s - 33 >= 30)
                    .count() as u64;
        }
        if let Some(barcode) = &fq.barcode {
            *self.num_total_bases.entry("barcode".to_string()).or_insert(0) +=
                barcode.raw.sequence().len() as u64;
            *self.num_q30_bases.entry("barcode".to_string()).or_insert(0) +=
                barcode.raw.quality_scores()
                    .iter()
                    .filter(|s| **s - 33 >= 30)
                    .count() as u64;
        }
        if let Some(read1) = &fq.read1 {
            *self.num_reads.entry("read1".to_string()).or_insert(0) += 1;
            *self.num_total_bases.entry("read1".to_string()).or_insert(0) +=
                read1.sequence().len() as u64;
            *self.num_q30_bases.entry("read1".to_string()).or_insert(0) +=
                read1.quality_scores()
                    .iter()
                    .filter(|s| **s - 33 >= 30)
                    .count() as u64;
        }
        if let Some(read2) = &fq.read2 {
            *self.num_reads.entry("read2".to_string()).or_insert(0) += 1;
            *self.num_total_bases.entry("read2".to_string()).or_insert(0) +=
                read2.sequence().len() as u64;
            *self.num_q30_bases.entry("read2".to_string()).or_insert(0) +=
                read2.quality_scores()
                    .iter()
                    .filter(|s| **s - 33 >= 30)
                    .count() as u64;
        }
    }
}

impl Extend<Self> for QcFastq {
    fn extend<T: IntoIterator<Item = Self>>(&mut self, iter: T) {
        for qc in iter {
            for (k, v) in qc.num_reads {
                *self.num_reads.entry(k).or_insert(0) += v;
            }
            for (k, v) in qc.num_defect {
                *self.num_defect.entry(k).or_insert(0) += v;
            }
            for (k, v) in qc.num_q30_bases {
                *self.num_q30_bases.entry(k).or_insert(0) += v;
            }
            for (k, v) in qc.num_total_bases {
                *self.num_total_bases.entry(k).or_insert(0) += v;
            }
        }
    }
}

impl Metric for QcFastq {
    fn to_json(&self) -> Value {
        let mut map = serde_json::Map::new();
        let defect = self
            .num_defect
            .iter()
            .filter(|x| *x.1 > 0)
            .map(|(k, v)| (k.clone(), (*v as f64 / self.num_reads[k] as f64).into()))
            .collect::<serde_json::Map<String, Value>>();
        if !defect.is_empty() {
            map.insert(
                "frac_reads_with_misformed_structure".to_string(),
                defect.into(),
            );
        }
        map.insert(
            "frac_q30_bases".to_string(),
            self.num_q30_bases
                .iter()
                .map(|(k, v)| {
                    (
                        k.clone(),
                        (*v as f64 / self.num_total_bases[k] as f64).into(),
                    )
                })
                .collect::<serde_json::Map<String, Value>>()
                .into(),
        );
        map.into()
    }
}

#[derive(Debug, Default)]
struct AlignStat {
    total: u64,        // Total number of reads
    mapped: u64,       // Number of mapped reads
    high_quality: u64, // Number of high-quality mapped reads: unique, non-duplicate, and mapping quality >= 30
    multimapped: u64,  // Number of reads with multiple alignments
    duplicate: u64,    // Number of duplicate reads
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
                let q = record
                    .primary
                    .mapping_quality()
                    .transpose()?
                    .map(|x| x.get())
                    .unwrap_or(60);
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
        self.read2.total.min(self.read1.total)
    }

    fn total_mapped(&self) -> u64 {
        self.read1.mapped + self.read2.mapped
    }

    fn total_high_quality(&self) -> u64 {
        self.read1.high_quality + self.read2.high_quality
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
pub struct QcAlign {
    pub(crate) mito_dna: HashSet<usize>, // Mitochondrial DNA reference sequence IDs
    stat_all: PairAlignStat,
    stat_barcoded: PairAlignStat,
    stat_mito: PairAlignStat,
}

impl Extend<Self> for QcAlign {
    fn extend<T: IntoIterator<Item = Self>>(&mut self, iter: T) {
        for qc in iter {
            self.stat_all.combine(&qc.stat_all);
            self.stat_barcoded.combine(&qc.stat_barcoded);
            self.stat_mito.combine(&qc.stat_mito);
        }
    }
}

impl Metric for QcAlign {
    fn to_json(&self) -> Value {
        let mut map = serde_json::Map::new();

        let stat_all = &self.stat_all;
        let stat_barcoded = &self.stat_barcoded;

        let fraction_confidently_mapped =
            stat_barcoded.total_high_quality() as f64 / stat_barcoded.total_reads() as f64;

        map.insert("sequenced_reads".to_string(), stat_all.total_reads().into());
        map.insert(
            "sequenced_read_pairs".to_string(),
            stat_all.total_pairs().into(),
        );
        if stat_all.total_pairs() > 0 {
            map.insert(
                "frac_properly_paired".to_string(),
                (stat_all.proper_pairs as f64 / stat_all.total_pairs() as f64).into(),
            );
        }
        map.insert(
            "frac_confidently_mapped".to_string(),
            fraction_confidently_mapped.into(),
        );
        map.insert("frac_unmapped".to_string(), self.frac_unmapped().into());
        map.insert(
            "frac_valid_barcode".to_string(),
            self.frac_valid_barcode().into(),
        );
        map.insert(
            "frac_mitochondrial".to_string(),
            self.frac_mitochondrial().into(),
        );

        map.into()
    }
}

impl QcAlign {
    pub fn add_mito_dna(&mut self, mito_dna: usize) {
        self.mito_dna.insert(mito_dna);
    }

    pub fn add_pair<R: Record>(
        &mut self,
        header: &sam::Header,
        record1: &MultiMap<R>,
        record2: &MultiMap<R>,
    ) -> Result<()> {
        let mut stat = PairAlignStat::default();

        stat.add_pair(record1, record2)?;

        self.stat_all.combine(&stat);

        if record1
            .primary
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

    pub fn add_read1<R: Record>(
        &mut self,
        header: &sam::Header,
        record: &MultiMap<R>,
    ) -> Result<()> {
        let mut stat = PairAlignStat::default();

        stat.add_read1(record)?;

        self.stat_all.combine(&stat);

        if record
            .primary
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

    pub fn add_read2<R: Record>(
        &mut self,
        header: &sam::Header,
        record: &MultiMap<R>,
    ) -> Result<()> {
        let mut stat = PairAlignStat::default();

        stat.add_read2(record)?;

        self.stat_all.combine(&stat);

        if record
            .primary
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

    /// Fraction of read pairs with barcodes that match the whitelist after error correction.
    pub fn frac_valid_barcode(&self) -> f64 {
        self.stat_barcoded.total_reads() as f64 / self.stat_all.total_reads() as f64
    }

    /// Fraction of sequenced read pairs with a valid barcode that could not be
    /// mapped to the genome, defined as the number of unmapped
    /// barcoded reads divided by the total number of barcoded reads.
    pub fn frac_unmapped(&self) -> f64 {
        1.0 - self.stat_barcoded.total_mapped() as f64 / self.stat_barcoded.total_reads() as f64
    }

    /// Estimated fraction of sequenced read pairs with a valid barcode that map to mitochondria,
    /// defined as the number of barcoded reads mapped to mitochondria divided by the total number of mapped barcoded reads.
    pub fn frac_mitochondrial(&self) -> f64 {
        self.stat_mito.total_reads() as f64 / self.stat_barcoded.total_mapped() as f64
    }
}

#[derive(Debug, Clone, Default)]
pub struct QcFragment {
    mito_dna: HashSet<String>,
    num_pcr_duplicates: u64,
    num_unique_fragments: u64,
    num_frag_nfr: u64,
    num_frag_single: u64,
}

impl From<QcFragment> for Value {
    fn from(qc: QcFragment) -> Self {
        let mut map = serde_json::Map::new();
        map.insert(
            "num_unique_fragments".to_string(),
            qc.num_unique_fragments.into(),
        );
        map.insert(
            "frac_duplicates".to_string(),
            (qc.num_pcr_duplicates as f64
                / (qc.num_unique_fragments + qc.num_pcr_duplicates) as f64)
                .into(),
        );
        map.insert(
            "frac_fragment_in_nucleosome_free_region".to_string(),
            (qc.num_frag_nfr as f64 / qc.num_unique_fragments as f64).into(),
        );
        map.insert(
            "frac_fragment_flanking_single_nucleosome".to_string(),
            (qc.num_frag_single as f64 / qc.num_unique_fragments as f64).into(),
        );

        map.into()
    }
}

impl QcFragment {
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
}

#[derive(Debug, Default)]
pub struct QcGeneQuant {
    pub(crate) total_reads: u64,
    pub(crate) total_umi: u64,
    pub(crate) unique_umi: u64,
}

impl From<QcGeneQuant> for Value {
    fn from(qc: QcGeneQuant) -> Self {
        json!({
            "frac_transcriptome": qc.total_umi as f64 / qc.total_reads as f64,
            "frac_duplicates": 1.0 - qc.unique_umi as f64 / qc.total_umi as f64,
        })
    }
}
