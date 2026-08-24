use crate::align::{AnnotatedFastq, MultiMap};
use crate::fragment::Fragment;
use crate::transcriptome::TxAlignment;

use anyhow::Result;
use bed_utils::bed::BEDLike;
use noodles_sam as sam;
use noodles_sam::alignment::{record::data::field::tag::Tag, Record};
use serde_json::{json, Value};
use std::collections::{HashMap, HashSet};

/// The trait for quality control metrics.
pub trait Metric: Sized + Extend<Self> {
    fn to_json(&self) -> Value;
}

#[derive(Debug, Clone, Default)]
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
            *self
                .num_total_bases
                .entry("barcode".to_string())
                .or_insert(0) += barcode.raw.sequence().len() as u64;
            *self.num_q30_bases.entry("barcode".to_string()).or_insert(0) += barcode
                .raw
                .quality_scores()
                .iter()
                .filter(|s| **s - 33 >= 30)
                .count()
                as u64;
        }
        if let Some(read1) = &fq.read1 {
            *self.num_reads.entry("read1".to_string()).or_insert(0) += 1;
            *self.num_total_bases.entry("read1".to_string()).or_insert(0) +=
                read1.sequence().len() as u64;
            *self.num_q30_bases.entry("read1".to_string()).or_insert(0) += read1
                .quality_scores()
                .iter()
                .filter(|s| **s - 33 >= 30)
                .count()
                as u64;
        }
        if let Some(read2) = &fq.read2 {
            *self.num_reads.entry("read2".to_string()).or_insert(0) += 1;
            *self.num_total_bases.entry("read2".to_string()).or_insert(0) +=
                read2.sequence().len() as u64;
            *self.num_q30_bases.entry("read2".to_string()).or_insert(0) += read2
                .quality_scores()
                .iter()
                .filter(|s| **s - 33 >= 30)
                .count()
                as u64;
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

#[derive(Debug, Clone, Default)]
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

#[derive(Debug, Clone, Default)]
struct PairAlignStat {
    read1: AlignStat,
    read2: AlignStat,
    sequenced_pairs: u64,
    proper_pairs: u64,
}

impl PairAlignStat {
    fn total_reads(&self) -> u64 {
        self.read1.total + self.read2.total
    }

    fn total_pairs(&self) -> u64 {
        self.sequenced_pairs
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
        // Opposite slots can hold independent single-end reads, so only this path establishes a pair.
        self.sequenced_pairs += 1;
        if record1.primary.flags()?.is_properly_segmented() {
            self.proper_pairs += 1;
        }
        Ok(())
    }

    fn combine(&mut self, other: &Self) {
        self.read1.combine(&other.read1);
        self.read2.combine(&other.read2);
        self.sequenced_pairs += other.sequenced_pairs;
        self.proper_pairs += other.proper_pairs;
    }
}

#[derive(Debug, Clone, Default)]
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
        let stat_all = &self.stat_all;
        let stat_barcoded = &self.stat_barcoded;
        let fraction_confidently_mapped =
            stat_barcoded.total_high_quality() as f64 / stat_barcoded.total_reads() as f64;
        json!({
            "sequenced_reads": stat_all.total_reads(),
            "sequenced_read_pairs": stat_all.total_pairs(),
            "frac_properly_paired": if stat_all.total_pairs() > 0 { stat_all.proper_pairs as f64 / stat_all.total_pairs() as f64 } else { 0.0 },
            "frac_confidently_mapped": fraction_confidently_mapped,
            "frac_unmapped": self.frac_unmapped(),
            "frac_valid_barcode": self.frac_valid_barcode(),
            "frac_mitochondrial": self.frac_mitochondrial(),
        })
    }
}

impl QcAlign {
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
    num_frag_nfr: u64,    // Nucleosome-free region fragments (<147 bp)
    num_frag_single: u64, // Flanking single nucleosome fragments (147-294 bp)
}

impl From<QcFragment> for Value {
    fn from(qc: QcFragment) -> Self {
        json!({
            "num_unique_fragments": qc.num_unique_fragments,
            "frac_duplicates": qc.num_pcr_duplicates as f64 / (qc.num_unique_fragments + qc.num_pcr_duplicates) as f64,
            "frac_fragment_in_nucleosome_free_region": qc.num_frag_nfr as f64 / qc.num_unique_fragments as f64,
            "frac_fragment_flanking_single_nucleosome": qc.num_frag_single as f64 / qc.num_unique_fragments as f64,
        })
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
    total_raw_count: u64,
    num_antisense: u64,
    num_intergenic: u64,
    num_multimapped: u64,
    num_discordant: u64,
    num_exonic: u64,
    num_spanning: u64,
    num_intronic: u64,
    num_mix: u64,
    pub(crate) num_unique_umi: u64,
    pub(crate) num_spliced: u64,
    pub(crate) num_unspliced: u64,
}

impl QcGeneQuant {
    pub fn update(&mut self, alignment: Option<&TxAlignment>) {
        self.total_raw_count += 1;
        match alignment {
            None => {}
            Some(TxAlignment::Antisense) => self.num_antisense += 1,
            Some(TxAlignment::Discordant) => self.num_discordant += 1,
            Some(TxAlignment::Intergenic) => self.num_intergenic += 1,
            Some(TxAlignment::Multimapped) => self.num_multimapped += 1,
            Some(aln) => {
                if aln.is_spanning() {
                    self.num_spanning += 1;
                } else if aln.is_exonic_only() {
                    self.num_exonic += 1;
                } else if aln.is_intronic_only() {
                    self.num_intronic += 1;
                } else {
                    self.num_mix += 1;
                }
            }
        }
    }
}

impl From<QcGeneQuant> for Value {
    fn from(qc: QcGeneQuant) -> Self {
        let num_transcriptomic = qc.num_exonic + qc.num_intronic + qc.num_mix + qc.num_spanning;
        let mapping = json!({
            "frac_intronic": qc.num_intronic as f64 / qc.total_raw_count as f64,
            "frac_exonic": qc.num_exonic as f64 / qc.total_raw_count as f64,
            "frac_spanning": qc.num_spanning as f64 / qc.total_raw_count as f64,
            "frac_mixed": qc.num_mix as f64 / qc.total_raw_count as f64,
            "frac_intergenic": qc.num_intergenic as f64 / qc.total_raw_count as f64,
            "frac_multimapped": qc.num_multimapped as f64 / qc.total_raw_count as f64,
            "frac_discordant": qc.num_discordant as f64 / qc.total_raw_count as f64,
            "frac_antisense": qc.num_antisense as f64 / qc.total_raw_count as f64,
        });
        json!({
            "frac_transcriptomic": num_transcriptomic as f64 / qc.total_raw_count as f64,
            "mapping": mapping,
            "quantification": {
                "num_unique_umi": qc.num_unique_umi,
                "frac_duplicates": 1.0 - qc.num_unique_umi as f64 / num_transcriptomic as f64,
                "frac_spliced": qc.num_spliced as f64 / qc.num_unique_umi as f64,
                "frac_unspliced": qc.num_unspliced as f64 / qc.num_unique_umi as f64,
            },
        })
    }
}

// Long-Read QC
/// Per-barcode-group consensus resolution statistics.
#[derive(Debug, Default, Clone)]
pub struct ConsensusStats {
    /// Number of reads where this group had entries from multiple ends
    pub attempted: u64,
    /// Number of times the candidate intersection was non-empty
    pub intersection_hit: u64,
    /// Number of times the candidate intersection was empty (fallback to best)
    pub intersection_miss: u64,
}

impl ConsensusStats {
    pub fn combine(&mut self, other: &Self) {
        self.attempted += other.attempted;
        self.intersection_hit += other.intersection_hit;
        self.intersection_miss += other.intersection_miss;
    }
}

/// QC metrics specific to long-read barcode extraction.
#[derive(Debug, Default)]
pub struct QcLongRead {
    // Orientation detection
    pub orientation_forward: u64,
    pub orientation_reverse: u64,
    pub orientation_undetermined: u64,

    // Composite alignment quality gate
    pub composite_alignment_pass: u64,
    pub composite_alignment_fail: u64,

    // Barcode extraction outcome
    pub barcode_extracted: u64,
    pub barcode_failed: u64,

    // UMI extraction outcome (only enabled for supported UMI layouts)
    pub umi_extracted: u64,
    pub umi_failed: u64,
    umi_extraction_enabled: bool,

    // Consensus resolution (keyed by barcode group name)
    pub consensus_stats: HashMap<String, ConsensusStats>,
}

impl QcLongRead {
    pub fn record_orientation(&mut self, is_reverse: bool, is_undetermined: bool) {
        if is_undetermined {
            self.orientation_undetermined += 1;
        } else if is_reverse {
            self.orientation_reverse += 1;
        } else {
            self.orientation_forward += 1;
        }
    }

    pub fn record_composite_alignment(&mut self, pass: bool) {
        if pass {
            self.composite_alignment_pass += 1;
        } else {
            self.composite_alignment_fail += 1;
        }
    }

    pub fn record_barcode_extraction(&mut self, success: bool) {
        if success {
            self.barcode_extracted += 1;
        } else {
            self.barcode_failed += 1;
        }
    }

    pub fn enable_umi_extraction(&mut self) {
        self.umi_extraction_enabled = true;
    }

    pub fn record_umi_extraction(&mut self, success: bool) {
        self.umi_extraction_enabled = true;
        if success {
            self.umi_extracted += 1;
        } else {
            self.umi_failed += 1;
        }
    }

    pub fn record_consensus(&mut self, group_name: &str, intersection_hit: bool) {
        let stats = self
            .consensus_stats
            .entry(group_name.to_string())
            .or_default();
        stats.attempted += 1;
        if intersection_hit {
            stats.intersection_hit += 1;
        } else {
            stats.intersection_miss += 1;
        }
    }
}

impl Extend<Self> for QcLongRead {
    fn extend<T: IntoIterator<Item = Self>>(&mut self, iter: T) {
        for other in iter {
            self.orientation_forward += other.orientation_forward;
            self.orientation_reverse += other.orientation_reverse;
            self.orientation_undetermined += other.orientation_undetermined;
            self.composite_alignment_pass += other.composite_alignment_pass;
            self.composite_alignment_fail += other.composite_alignment_fail;
            self.barcode_extracted += other.barcode_extracted;
            self.barcode_failed += other.barcode_failed;
            self.umi_extracted += other.umi_extracted;
            self.umi_failed += other.umi_failed;
            self.umi_extraction_enabled |= other.umi_extraction_enabled;
            for (k, v) in other.consensus_stats {
                self.consensus_stats.entry(k).or_default().combine(&v);
            }
        }
    }
}

impl Metric for QcLongRead {
    fn to_json(&self) -> Value {
        let total_reads =
            self.orientation_forward + self.orientation_reverse + self.orientation_undetermined;

        let consensus: serde_json::Map<String, Value> = self
            .consensus_stats
            .iter()
            .map(|(name, stats)| {
                (
                    name.clone(),
                    json!({
                        "attempted": stats.attempted,
                        "intersection_hit": stats.intersection_hit,
                        "intersection_miss": stats.intersection_miss,
                    }),
                )
            })
            .collect();

        let mut result = json!({
            "total_reads": total_reads,
            "orientation": {
                "forward": self.orientation_forward,
                "reverse": self.orientation_reverse,
                "undetermined": self.orientation_undetermined,
            },
            "composite_alignment": {
                "pass": self.composite_alignment_pass,
                "fail": self.composite_alignment_fail,
            },
            "barcode_extraction": {
                "success": self.barcode_extracted,
                "fail": self.barcode_failed,
            },
            "consensus": consensus,
        });
        if self.umi_extraction_enabled {
            result.as_object_mut().unwrap().insert(
                "umi_extraction".to_string(),
                json!({
                    "success": self.umi_extracted,
                    "fail": self.umi_failed,
                }),
            );
        }
        result
    }
}

#[cfg(test)]
mod tests {
    use noodles_sam::alignment::record::Flags;
    use noodles_sam::alignment::record_buf::RecordBuf;

    use super::*;

    fn alignment(flags: Flags) -> MultiMap<RecordBuf> {
        let mut record = RecordBuf::default();
        *record.flags_mut() = flags;
        record.into()
    }

    #[test]
    fn independent_single_end_slots_are_not_counted_as_a_pair() {
        let header = sam::Header::default();
        let mut qc = QcAlign::default();

        qc.add_read1(&header, &alignment(Flags::empty())).unwrap();
        qc.add_read2(&header, &alignment(Flags::empty())).unwrap();

        let report = qc.to_json();
        assert_eq!(report["sequenced_reads"], json!(2));
        assert_eq!(report["sequenced_read_pairs"], json!(0));
        assert_eq!(report["frac_properly_paired"], json!(0.0));
    }

    #[test]
    fn paired_end_record_is_counted_as_one_pair() {
        let header = sam::Header::default();
        let mut qc = QcAlign::default();
        let read1 = alignment(Flags::SEGMENTED | Flags::PROPERLY_SEGMENTED);
        let read2 = alignment(Flags::SEGMENTED);

        qc.add_pair(&header, &read1, &read2).unwrap();

        let report = qc.to_json();
        assert_eq!(report["sequenced_reads"], json!(2));
        assert_eq!(report["sequenced_read_pairs"], json!(1));
        assert_eq!(report["frac_properly_paired"], json!(1.0));
    }
}
