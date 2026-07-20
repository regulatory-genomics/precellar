//! TFseq floating-barcode middleware.

use super::{AnnotatedFastqBatch, FastqStage, MiddlewareQcReport};
use crate::align::AnnotatedFastq;
use crate::barcode::{BarcodeCorrectOptions, Whitelist, WhitelistBuilder};
use crate::utils::get_directional_umi_mapping;
use crate::utils::insertion_extractor::InsertionExtractor;
use anyhow::Result;
use bitcode::{Decode, Encode};
use rayon::prelude::*;
use serde_json::json;
use std::collections::{BTreeMap, HashMap};
use std::io::Write;

const CORRECTION_BATCH_SIZE: usize = 4_194_304;

#[derive(Debug, Default)]
pub struct FloatingBarcodeQc {
    pub num_records: u64,
    pub num_extracted: u64,
    pub num_matched: u64,
    pub num_extracted_single: Vec<u64>,
    pub num_extracted_all: u64,
    pub num_extracted_with_valid_cell_barcode: u64,
}

impl FloatingBarcodeQc {
    pub fn frac_valid_cell_barcode(&self) -> f64 {
        if self.num_extracted == 0 {
            0.0
        } else {
            self.num_extracted_with_valid_cell_barcode as f64 / self.num_extracted as f64
        }
    }

    pub fn to_json(&self) -> serde_json::Value {
        let mut metrics = json!({
            "num_records": self.num_records,
            "num_extracted": self.num_extracted,
            "num_matched": self.num_matched,
            "num_extracted_all": self.num_extracted_all,
            "frac_valid_cell_barcode": self.frac_valid_cell_barcode(),
        });
        if let Some(metrics) = metrics.as_object_mut() {
            for (index, count) in self.num_extracted_single.iter().enumerate() {
                metrics.insert(
                    format!("num_extracted_single_{}", index + 1),
                    (*count).into(),
                );
            }
        }
        metrics
    }
}

#[derive(Decode, Encode)]
struct PendingInsertion {
    insertion_index: usize,
    sequence: Vec<u8>,
    quality_scores: Vec<u8>,
}

#[derive(Decode, Encode)]
struct PendingFloatingBarcode {
    barcode: Vec<u8>,
    umi: Vec<u8>,
    insertions: Vec<PendingInsertion>,
}

enum ProcessedFloatingBarcode {
    Forward(AnnotatedFastq),
    Extracted {
        barcode: Option<Vec<u8>>,
        umi: Vec<u8>,
        insertions: Vec<PendingInsertion>,
    },
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct FloatingBarcodeEntry {
    pub name: Vec<u8>,
    pub barcode_1: Vec<u8>,
    pub barcode_2: Vec<u8>,
}

impl FloatingBarcodeEntry {
    pub fn new<N, B1, B2>(name: N, barcode_1: B1, barcode_2: B2) -> Self
    where
        N: Into<Vec<u8>>,
        B1: Into<Vec<u8>>,
        B2: Into<Vec<u8>>,
    {
        Self {
            name: name.into(),
            barcode_1: barcode_1.into(),
            barcode_2: barcode_2.into(),
        }
    }
}

/// Floating-barcode entries and their inferred insertion-length profiles.
#[derive(Clone)]
pub struct FloatingBarcodeTable {
    entries: Vec<FloatingBarcodeEntry>,
    length_profiles: Vec<Vec<usize>>,
    barcode_1_rows: HashMap<Vec<u8>, Option<usize>>,
    barcode_2_rows: HashMap<Vec<u8>, Option<usize>>,
    pair_rows: HashMap<Vec<u8>, HashMap<Vec<u8>, usize>>,
}

impl FloatingBarcodeTable {
    /// Builds a table and records each unique barcode-length pair in first-seen order.
    pub fn new(mut entries: Vec<FloatingBarcodeEntry>) -> Result<Self> {
        if entries.is_empty() {
            anyhow::bail!("barcode_table must contain at least one row");
        }

        let mut length_profiles = Vec::new();
        let mut name_rows = HashMap::new();
        let mut barcode_1_rows = HashMap::new();
        let mut barcode_2_rows = HashMap::new();
        let mut pair_rows: HashMap<Vec<u8>, HashMap<Vec<u8>, usize>> = HashMap::new();
        for (index, entry) in entries.iter_mut().enumerate() {
            if entry.name.is_empty() {
                anyhow::bail!("barcode table names must not be empty");
            }
            if entry
                .name
                .iter()
                .any(|byte| matches!(byte, b'\t' | b'\n' | b'\r' | 0..=31))
            {
                anyhow::bail!("barcode table names must not contain control characters");
            }
            if let Some(first_index) = name_rows.insert(entry.name.clone(), index) {
                anyhow::bail!(
                    "barcode table row {} has duplicate name '{}' first seen on row {}",
                    index + 1,
                    String::from_utf8_lossy(&entry.name),
                    first_index + 1
                );
            }
            entry.barcode_1.make_ascii_uppercase();
            entry.barcode_2.make_ascii_uppercase();
            if !is_dna(&entry.barcode_1) || !is_dna(&entry.barcode_2) {
                anyhow::bail!(
                    "barcode table barcodes contain invalid DNA sequences: {}",
                    String::from_utf8_lossy(&entry.barcode_1)
                );
            }
            let profile = vec![entry.barcode_1.len(), entry.barcode_2.len()];
            if !length_profiles.contains(&profile) {
                length_profiles.push(profile);
            }

            insert_unique_index(&mut barcode_1_rows, &entry.barcode_1, index);
            insert_unique_index(&mut barcode_2_rows, &entry.barcode_2, index);
            if pair_rows
                .entry(entry.barcode_1.clone())
                .or_default()
                .insert(entry.barcode_2.clone(), index)
                .is_some()
            {
                anyhow::bail!("barcode table contains a duplicate barcode combination");
            }
        }

        Ok(Self {
            entries,
            length_profiles,
            barcode_1_rows,
            barcode_2_rows,
            pair_rows,
        })
    }

    pub fn length_profiles(&self) -> &[Vec<usize>] {
        &self.length_profiles
    }

    fn resolve(
        &self,
        barcode_1: Option<&[u8]>,
        barcode_2: Option<&[u8]>,
    ) -> Option<&FloatingBarcodeEntry> {
        let index = match (barcode_1, barcode_2) {
            (Some(barcode_1), Some(barcode_2)) => *self.pair_rows.get(barcode_1)?.get(barcode_2)?,
            (Some(barcode_1), None) => self.barcode_1_rows.get(barcode_1).copied().flatten()?,
            (None, Some(barcode_2)) => self.barcode_2_rows.get(barcode_2).copied().flatten()?,
            (None, None) => return None,
        };
        self.entries.get(index)
    }
}

fn is_dna(sequence: &[u8]) -> bool {
    !sequence.is_empty()
        && sequence
            .iter()
            .all(|base| matches!(base, b'A' | b'C' | b'G' | b'T'))
}

fn insert_unique_index(rows: &mut HashMap<Vec<u8>, Option<usize>>, barcode: &[u8], index: usize) {
    rows.entry(barcode.to_vec())
        .and_modify(|existing| *existing = None)
        .or_insert(Some(index));
}

/// Extracts floating barcodes and removes extracted records from downstream.
///
/// Floating-barcode correction uses an explicit, non-empty whitelist. Every
/// extracted floating sequence contributes to whitelist frequencies, but output
/// is emitted only when the cell barcode and at least one floating barcode have
/// corrected values. Unavailable or uncorrectable floating barcodes are `-`.
pub struct FloatingBarcodeStage<W> {
    use_read1: bool,
    extractor: InsertionExtractor,
    barcode_table: FloatingBarcodeTable,
    whitelist_builders: Option<Vec<WhitelistBuilder>>,
    correction_options: BarcodeCorrectOptions,
    pending: Option<bed_utils::extsort::ExternalChunkBuilder<PendingFloatingBarcode>>,
    output: W,
    thread_pool: rayon::ThreadPool,
    qc: FloatingBarcodeQc,
    finished: bool,
}

impl<W> FloatingBarcodeStage<W>
where
    W: Write + Send,
{
    pub fn new(
        use_read1: bool,
        extractor: InsertionExtractor,
        barcode_table: FloatingBarcodeTable,
        correction_options: BarcodeCorrectOptions,
        output: W,
    ) -> Result<Self> {
        if extractor.length_profiles() != barcode_table.length_profiles() {
            anyhow::bail!("barcode table length profiles do not match extractor length profiles");
        }
        let num_insertions = extractor.num_gaps();
        let valid_barcodes = [
            barcode_table
                .entries
                .iter()
                .map(|entry| entry.barcode_1.clone())
                .collect::<Vec<_>>(),
            barcode_table
                .entries
                .iter()
                .map(|entry| entry.barcode_2.clone())
                .collect::<Vec<_>>(),
        ];

        Ok(Self {
            use_read1,
            extractor,
            barcode_table,
            whitelist_builders: Some(valid_barcodes.into_iter().map(Whitelist::builder).collect()),
            correction_options,
            pending: Some(bed_utils::extsort::ExternalChunkBuilder::new(
                tempfile::tempfile()?,
                2,
            )?),
            output,
            thread_pool: rayon::ThreadPoolBuilder::new().num_threads(1).build()?,
            qc: FloatingBarcodeQc {
                num_extracted_single: vec![0; num_insertions],
                ..FloatingBarcodeQc::default()
            },
            finished: false,
        })
    }

    pub fn qc(&self) -> &FloatingBarcodeQc {
        &self.qc
    }

    pub fn into_output(self) -> W {
        self.output
    }

    pub fn with_num_threads(mut self, num_threads: usize) -> Result<Self> {
        if num_threads == 0 {
            anyhow::bail!("num_threads must be greater than zero");
        }
        self.thread_pool = rayon::ThreadPoolBuilder::new()
            .num_threads(num_threads)
            .build()?;
        Ok(self)
    }

    fn selected_sequence<'a>(
        use_read1: bool,
        record: &'a AnnotatedFastq,
    ) -> Option<(&'a [u8], &'a [u8])> {
        if use_read1 {
            record
                .read1
                .as_ref()
                .map(|read| (read.sequence(), read.quality_scores()))
        } else {
            record
                .read2
                .as_ref()
                .map(|read| (read.sequence(), read.quality_scores()))
        }
    }

    fn write_result(
        output: &mut W,
        (cell_barcode, name): (Vec<u8>, Vec<u8>),
        umis: BTreeMap<Vec<u8>, u64>,
    ) -> std::io::Result<()> {
        output.write_all(&cell_barcode)?;
        output.write_all(b"\t")?;
        output.write_all(&name)?;
        output.write_all(b"\t")?;
        for (index, (umi, count)) in umis.into_iter().enumerate() {
            if index > 0 {
                output.write_all(b";")?;
            }
            output.write_all(&umi)?;
            write!(output, ":{count}")?;
        }
        output.write_all(b"\n")?;
        Ok(())
    }

    fn process_record(
        extractor: &InsertionExtractor,
        use_read1: bool,
        record: AnnotatedFastq,
    ) -> Result<ProcessedFloatingBarcode> {
        let Some((sequence, quality_scores)) = Self::selected_sequence(use_read1, &record) else {
            return Ok(ProcessedFloatingBarcode::Forward(record));
        };
        if sequence.len() != quality_scores.len() {
            anyhow::bail!("selected FASTQ sequence and quality lengths differ");
        }
        let Some(ranges) = extractor.extract_ranges(sequence) else {
            return Ok(ProcessedFloatingBarcode::Forward(record));
        };
        if ranges.is_empty() {
            return Ok(ProcessedFloatingBarcode::Forward(record));
        }

        let barcode = record
            .barcode
            .as_ref()
            .and_then(|barcode| barcode.corrected.clone());
        let umi = if barcode.is_some() {
            record
                .umi
                .as_ref()
                .map(|umi| umi.sequence().to_vec())
                .unwrap_or_default()
        } else {
            Vec::new()
        };
        let insertions = ranges
            .into_iter()
            .map(|(insertion_index, range)| PendingInsertion {
                insertion_index,
                sequence: sequence[range.clone()].to_vec(),
                quality_scores: quality_scores[range].to_vec(),
            })
            .collect();
        Ok(ProcessedFloatingBarcode::Extracted {
            barcode,
            umi,
            insertions,
        })
    }

    fn resolve_pending(
        record: PendingFloatingBarcode,
        whitelists: &[Whitelist],
        barcode_table: &FloatingBarcodeTable,
        correction_options: &BarcodeCorrectOptions,
    ) -> Option<(Vec<u8>, Vec<u8>, Vec<u8>)> {
        let PendingFloatingBarcode {
            barcode,
            umi,
            insertions,
        } = record;
        let mut floating_barcodes: [Option<Vec<u8>>; 2] = [None, None];
        for insertion in insertions {
            let PendingInsertion {
                insertion_index,
                sequence,
                quality_scores,
            } = insertion;
            let whitelist = whitelists.get(insertion_index)?;
            if let Ok(corrected) =
                whitelist.correct_barcode(&sequence, &quality_scores, correction_options)
            {
                floating_barcodes[insertion_index] = Some(corrected.to_vec());
            }
        }
        let entry = barcode_table.resolve(
            floating_barcodes[0].as_deref(),
            floating_barcodes[1].as_deref(),
        )?;
        Some((barcode, entry.name.clone(), umi))
    }
}

impl<W> FastqStage for FloatingBarcodeStage<W>
where
    W: Write + Send,
{
    fn process(&mut self, batch: AnnotatedFastqBatch) -> Result<AnnotatedFastqBatch> {
        let num_records = batch.len() as u64;
        let extractor = &self.extractor;
        let use_read1 = self.use_read1;
        let results: Vec<_> = self.thread_pool.install(|| {
            batch
                .into_par_iter()
                .map(|record| Self::process_record(extractor, use_read1, record))
                .collect()
        });
        let results: Vec<_> = results.into_iter().collect::<Result<_>>()?;
        self.qc.num_records += num_records;

        let mut forwarded = Vec::with_capacity(results.len());
        for result in results {
            match result {
                ProcessedFloatingBarcode::Forward(record) => forwarded.push(record),
                ProcessedFloatingBarcode::Extracted {
                    barcode,
                    umi,
                    insertions,
                } => {
                    self.qc.num_extracted += 1;
                    if barcode.is_some() {
                        self.qc.num_extracted_with_valid_cell_barcode += 1;
                    }
                    if insertions.len() == 1 {
                        let insertion_index = insertions[0].insertion_index;
                        *self
                            .qc
                            .num_extracted_single
                            .get_mut(insertion_index)
                            .expect("extractor returned an invalid insertion index") += 1;
                    } else if insertions.len() == self.qc.num_extracted_single.len() {
                        self.qc.num_extracted_all += 1;
                    }
                    for insertion in &insertions {
                        self.whitelist_builders
                            .as_mut()
                            .expect("floating barcode stage already finalized")
                            .get_mut(insertion.insertion_index)
                            .expect("extractor returned an invalid insertion index")
                            .add(&insertion.sequence);
                    }
                    if let Some(barcode) = barcode {
                        self.pending
                            .as_mut()
                            .expect("floating barcode stage already finalized")
                            .add(PendingFloatingBarcode {
                                barcode,
                                umi,
                                insertions,
                            })?;
                    }
                }
            }
        }
        Ok(forwarded)
    }

    fn finish(&mut self) -> Result<()> {
        if self.finished {
            return Ok(());
        }

        let whitelists: Vec<Whitelist> = self
            .whitelist_builders
            .take()
            .expect("floating barcode stage already finalized")
            .into_iter()
            .map(WhitelistBuilder::finish)
            .collect();
        let pending = self
            .pending
            .take()
            .expect("floating barcode stage already finalized")
            .finish()?;
        let whitelists = &whitelists;
        let barcode_table = &self.barcode_table;
        let correction_options = &self.correction_options;
        let pool = &self.thread_pool;
        let mut pending = pending;
        let mut matches: BTreeMap<_, BTreeMap<Vec<u8>, u64>> = BTreeMap::new();
        loop {
            let mut batch = Vec::with_capacity(CORRECTION_BATCH_SIZE);
            while batch.len() < CORRECTION_BATCH_SIZE {
                match pending.next() {
                    Some(Ok(record)) => batch.push(record),
                    Some(Err(error)) => return Err(error.into()),
                    None => break,
                }
            }
            if batch.is_empty() {
                break;
            }
            let resolved = pool.install(|| {
                batch
                    .into_par_iter()
                    .filter_map(|record| {
                        Self::resolve_pending(record, whitelists, barcode_table, correction_options)
                    })
                    .collect::<Vec<_>>()
            });
            self.qc.num_matched += resolved.len() as u64;
            for (barcode, name, umi) in resolved {
                *matches
                    .entry((barcode, name))
                    .or_default()
                    .entry(umi)
                    .or_default() += 1;
            }
        }

        for umis in matches.values_mut() {
            let umi_counts = std::mem::take(umis).into_iter().collect::<HashMap<_, _>>();
            let umi_mapping = get_directional_umi_mapping(&umi_counts);
            for (umi, count) in umi_counts {
                let corrected_umi = umi_mapping.get(&umi).unwrap_or(&umi).clone();
                *umis.entry(corrected_umi).or_default() += count;
            }
        }

        self.output.write_all(b"cell_barcode\tname\tumi_counts\n")?;
        for (key, umis) in matches {
            Self::write_result(&mut self.output, key, umis)?;
        }
        self.output.flush()?;
        self.finished = true;
        Ok(())
    }

    fn report(&self) -> Option<MiddlewareQcReport> {
        Some(MiddlewareQcReport {
            name: "floating_barcode".to_string(),
            metrics: self.qc.to_json(),
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use noodles_fastq as fastq;

    fn record_with_metadata(
        sequence: &[u8],
        cell_barcode: Option<&[u8]>,
        umi: Option<&[u8]>,
    ) -> AnnotatedFastq {
        let quality_scores = vec![b'I'; sequence.len()];
        AnnotatedFastq {
            barcode: cell_barcode.map(|barcode| crate::align::Barcode {
                raw: fastq::Record::new(
                    fastq::record::Definition::default(),
                    barcode,
                    vec![b'I'; barcode.len()],
                ),
                corrected: Some(barcode.to_vec()),
            }),
            umi: umi.map(|umi| {
                fastq::Record::new(
                    fastq::record::Definition::default(),
                    umi,
                    vec![b'I'; umi.len()],
                )
            }),
            read1: Some(fastq::Record::new(
                fastq::record::Definition::default(),
                sequence,
                quality_scores,
            )),
            read2: None,
        }
    }

    fn record(sequence: &[u8]) -> AnnotatedFastq {
        record_with_metadata(sequence, Some(b"CELL"), None)
    }

    fn extractor() -> InsertionExtractor {
        InsertionExtractor::new(
            vec!["AAAA".into(), "CCCC".into(), "GGGG".into()],
            vec![vec![2, 3]],
            3,
            0.0,
        )
        .unwrap()
    }

    fn table() -> FloatingBarcodeTable {
        FloatingBarcodeTable::new(vec![FloatingBarcodeEntry::new("TF1", "TT", "AAA")]).unwrap()
    }

    fn full_sequence() -> &'static [u8] {
        b"AAAATTCCCCAAAGGGG"
    }

    #[test]
    fn infers_unique_length_profiles_in_table_order() {
        let table = FloatingBarcodeTable::new(vec![
            FloatingBarcodeEntry::new("TF1", "TT", "AAA"),
            FloatingBarcodeEntry::new("TF2", "GGGG", "CCCCC"),
            FloatingBarcodeEntry::new("TF3", "AA", "GGG"),
        ])
        .unwrap();

        assert_eq!(table.length_profiles(), &[vec![2, 3], vec![4, 5]]);
    }

    #[test]
    fn processes_barcode_rows_with_different_length_profiles() {
        let table = FloatingBarcodeTable::new(vec![
            FloatingBarcodeEntry::new("TF1", "TT", "AAA"),
            FloatingBarcodeEntry::new("TF2", "TTTT", "AAAAA"),
        ])
        .unwrap();
        let extractor = InsertionExtractor::new(
            vec!["AAAA".into(), "CCCC".into(), "GGGG".into()],
            table.length_profiles().to_vec(),
            3,
            0.0,
        )
        .unwrap();
        let mut stage = FloatingBarcodeStage::new(
            true,
            extractor,
            table,
            BarcodeCorrectOptions::default(),
            Vec::new(),
        )
        .unwrap();

        stage
            .process(vec![record(b"AAAATTTTCCCCAAAAAGGGG")])
            .unwrap();
        stage.finish().unwrap();
        assert_eq!(
            stage.into_output(),
            b"cell_barcode\tname\tumi_counts\nCELL\tTF2\t:1\n"
        );
    }

    #[test]
    fn serializes_position_specific_extraction_qc() {
        let qc = FloatingBarcodeQc {
            num_records: 10,
            num_extracted: 9,
            num_matched: 8,
            num_extracted_single: vec![1, 2, 3],
            num_extracted_all: 4,
            num_extracted_with_valid_cell_barcode: 6,
        };

        assert_eq!(
            qc.to_json(),
            json!({
                "num_records": 10,
                "num_extracted": 9,
                "num_matched": 8,
                "num_extracted_single_1": 1,
                "num_extracted_single_2": 2,
                "num_extracted_single_3": 3,
                "num_extracted_all": 4,
                "frac_valid_cell_barcode": 6.0 / 9.0,
            })
        );
    }

    #[test]
    fn valid_cell_barcode_fraction_is_zero_without_extracted_reads() {
        assert_eq!(FloatingBarcodeQc::default().frac_valid_cell_barcode(), 0.0);
    }

    #[test]
    fn removes_extracted_reads_and_writes_corrected_barcode() {
        let mut stage = FloatingBarcodeStage::new(
            true,
            extractor(),
            table(),
            BarcodeCorrectOptions::default(),
            Vec::new(),
        )
        .unwrap();

        let forwarded = stage.process(vec![record(full_sequence())]).unwrap();
        assert!(forwarded.is_empty());
        assert_eq!(stage.qc().frac_valid_cell_barcode(), 1.0);
        stage.finish().unwrap();
        assert_eq!(
            stage.into_output(),
            b"cell_barcode\tname\tumi_counts\nCELL\tTF1\t:1\n"
        );
    }

    #[test]
    fn removes_extracted_reads_when_barcode_is_not_correctable() {
        let mut stage = FloatingBarcodeStage::new(
            true,
            extractor(),
            FloatingBarcodeTable::new(vec![FloatingBarcodeEntry::new("TF2", "GG", "CCC")]).unwrap(),
            BarcodeCorrectOptions::default(),
            Vec::new(),
        )
        .unwrap();

        let forwarded = stage.process(vec![record(full_sequence())]).unwrap();
        assert!(forwarded.is_empty());
        stage.finish().unwrap();
        assert_eq!(stage.into_output(), b"cell_barcode\tname\tumi_counts\n");
    }

    #[test]
    fn skips_matches_without_a_corrected_cell_barcode() {
        let mut stage = FloatingBarcodeStage::new(
            true,
            extractor(),
            table(),
            BarcodeCorrectOptions::default(),
            Vec::new(),
        )
        .unwrap();
        let mut record = record(full_sequence());
        record.barcode.as_mut().unwrap().corrected = None;

        stage.process(vec![record]).unwrap();
        assert_eq!(stage.qc().frac_valid_cell_barcode(), 0.0);
        stage.finish().unwrap();
        assert_eq!(stage.into_output(), b"cell_barcode\tname\tumi_counts\n");
    }

    #[test]
    fn groups_matches_by_cell_barcode_and_umi() {
        let mut stage = FloatingBarcodeStage::new(
            true,
            extractor(),
            table(),
            BarcodeCorrectOptions::default(),
            Vec::new(),
        )
        .unwrap();
        let sequence = full_sequence();
        stage
            .process(vec![
                record_with_metadata(sequence, Some(b"CCCC"), Some(b"UMI2")),
                record_with_metadata(sequence, Some(b"AAAA"), None),
                record_with_metadata(sequence, Some(b"CCCC"), Some(b"UMI1")),
                record_with_metadata(sequence, Some(b"CCCC"), Some(b"UMI2")),
            ])
            .unwrap();
        stage.finish().unwrap();
        assert_eq!(
            stage.into_output(),
            b"cell_barcode\tname\tumi_counts\nAAAA\tTF1\t:1\nCCCC\tTF1\tUMI1:1;UMI2:2\n"
        );
    }

    #[test]
    fn corrects_umis_and_preserves_read_counts() {
        let mut stage = FloatingBarcodeStage::new(
            true,
            extractor(),
            table(),
            BarcodeCorrectOptions::default(),
            Vec::new(),
        )
        .unwrap();
        let sequence = full_sequence();
        stage
            .process(vec![
                record_with_metadata(sequence, Some(b"CELL"), Some(b"AAAA")),
                record_with_metadata(sequence, Some(b"CELL"), Some(b"AAAA")),
                record_with_metadata(sequence, Some(b"CELL"), Some(b"AAAT")),
            ])
            .unwrap();
        stage.finish().unwrap();
        assert_eq!(
            stage.into_output(),
            b"cell_barcode\tname\tumi_counts\nCELL\tTF1\tAAAA:3\n"
        );
    }

    #[test]
    fn resolves_transitive_umi_corrections_to_the_root() {
        let mut stage = FloatingBarcodeStage::new(
            true,
            extractor(),
            table(),
            BarcodeCorrectOptions::default(),
            Vec::new(),
        )
        .unwrap();
        let sequence = full_sequence();
        let mut records = Vec::new();
        for (umi, count) in [
            (b"AAAA".as_slice(), 1),
            (b"AAAT".as_slice(), 2),
            (b"AAGT".as_slice(), 4),
        ] {
            for _ in 0..count {
                records.push(record_with_metadata(sequence, Some(b"CELL"), Some(umi)));
            }
        }

        stage.process(records).unwrap();
        stage.finish().unwrap();
        assert_eq!(
            stage.into_output(),
            b"cell_barcode\tname\tumi_counts\nCELL\tTF1\tAAGT:7\n"
        );
    }

    #[test]
    fn corrects_umis_within_each_cell_and_barcode_group() {
        let mut stage = FloatingBarcodeStage::new(
            true,
            extractor(),
            table(),
            BarcodeCorrectOptions::default(),
            Vec::new(),
        )
        .unwrap();
        let sequence = full_sequence();
        stage
            .process(vec![
                record_with_metadata(sequence, Some(b"CELL1"), Some(b"AAAT")),
                record_with_metadata(sequence, Some(b"CELL2"), Some(b"AAAA")),
                record_with_metadata(sequence, Some(b"CELL2"), Some(b"AAAA")),
            ])
            .unwrap();
        stage.finish().unwrap();
        assert_eq!(
            stage.into_output(),
            b"cell_barcode\tname\tumi_counts\nCELL1\tTF1\tAAAT:1\nCELL2\tTF1\tAAAA:2\n"
        );
    }

    fn run_with_threads(num_threads: usize) -> (Vec<u8>, Vec<Vec<u8>>, serde_json::Value) {
        let mut stage = FloatingBarcodeStage::new(
            true,
            extractor(),
            table(),
            BarcodeCorrectOptions::default(),
            Vec::new(),
        )
        .unwrap()
        .with_num_threads(num_threads)
        .unwrap();
        let forwarded = stage
            .process(vec![
                record(b"GATTGATT"),
                record(full_sequence()),
                record(b"TCTCTCTC"),
                record(b"AAAATTCCCC"),
            ])
            .unwrap()
            .into_iter()
            .filter_map(|record| record.read1.map(|read| read.sequence().to_vec()))
            .collect();
        stage.finish().unwrap();
        let metrics = stage.qc().to_json();
        (stage.into_output(), forwarded, metrics)
    }

    #[test]
    fn parallel_processing_matches_serial_processing() {
        assert_eq!(run_with_threads(1), run_with_threads(2));
    }

    #[test]
    fn writes_multiple_insertions_as_one_combined_key() {
        let extractor = InsertionExtractor::new(
            vec!["AAAA".into(), "CCCC".into(), "GGGG".into()],
            vec![vec![2, 3]],
            3,
            0.0,
        )
        .unwrap();
        let mut stage = FloatingBarcodeStage::new(
            true,
            extractor,
            table(),
            BarcodeCorrectOptions::default(),
            Vec::new(),
        )
        .unwrap();

        let forwarded = stage.process(vec![record(b"AAAATTCCCCAAAGGGG")]).unwrap();
        assert!(forwarded.is_empty());
        assert_eq!(stage.qc().num_extracted, 1);
        assert_eq!(stage.qc().num_extracted_single, vec![0, 0]);
        assert_eq!(stage.qc().num_extracted_all, 1);
        stage.finish().unwrap();
        assert_eq!(
            stage.into_output(),
            b"cell_barcode\tname\tumi_counts\nCELL\tTF1\t:1\n"
        );
    }

    #[test]
    fn writes_dash_for_an_unavailable_insertion() {
        let extractor = InsertionExtractor::new(
            vec!["AAAA".into(), "CCCC".into(), "GGGG".into()],
            vec![vec![2, 3]],
            3,
            0.0,
        )
        .unwrap();
        let mut stage = FloatingBarcodeStage::new(
            true,
            extractor,
            table(),
            BarcodeCorrectOptions::default(),
            Vec::new(),
        )
        .unwrap();

        stage
            .process(vec![record(b"AAAATTCCCCA"), record(b"TCCCCAAAGGGG")])
            .unwrap();
        assert_eq!(stage.qc().num_extracted, 2);
        assert_eq!(stage.qc().num_extracted_single, vec![1, 1]);
        assert_eq!(stage.qc().num_extracted_all, 0);
        stage.finish().unwrap();
        assert_eq!(
            stage.into_output(),
            b"cell_barcode\tname\tumi_counts\nCELL\tTF1\t:2\n"
        );
    }

    #[test]
    fn writes_a_read_when_at_least_one_insertion_corrects() {
        let extractor = InsertionExtractor::new(
            vec!["AAAA".into(), "CCCC".into(), "GGGG".into()],
            vec![vec![2, 3]],
            3,
            0.0,
        )
        .unwrap();
        let mut stage = FloatingBarcodeStage::new(
            true,
            extractor,
            table(),
            BarcodeCorrectOptions::default(),
            Vec::new(),
        )
        .unwrap();

        stage.process(vec![record(b"AAAATTCCCCTTTGGGG")]).unwrap();
        stage.finish().unwrap();
        assert_eq!(stage.qc().num_extracted, 1);
        assert_eq!(stage.qc().num_matched, 1);
        assert_eq!(
            stage.into_output(),
            b"cell_barcode\tname\tumi_counts\nCELL\tTF1\t:1\n"
        );
    }

    #[test]
    fn rejects_an_empty_barcode_table() {
        assert!(FloatingBarcodeTable::new(Vec::new()).is_err());
    }

    #[test]
    fn rejects_duplicate_barcode_table_names() {
        let error = FloatingBarcodeTable::new(vec![
            FloatingBarcodeEntry::new("TF1", "TT", "AAA"),
            FloatingBarcodeEntry::new("TF1", "GG", "CCC"),
        ])
        .err()
        .expect("duplicate names must be rejected");

        assert!(error.to_string().contains("duplicate name 'TF1'"));
        assert!(error.to_string().contains("row 2"));
    }

    #[test]
    fn rejects_an_ambiguous_single_barcode_lookup() {
        let mut stage = FloatingBarcodeStage::new(
            true,
            extractor(),
            FloatingBarcodeTable::new(vec![
                FloatingBarcodeEntry::new("TF1", "TT", "AAA"),
                FloatingBarcodeEntry::new("TF2", "TT", "CCC"),
            ])
            .unwrap(),
            BarcodeCorrectOptions::default(),
            Vec::new(),
        )
        .unwrap();

        stage.process(vec![record(b"AAAATTCCCC")]).unwrap();
        stage.finish().unwrap();
        assert_eq!(stage.qc().num_matched, 0);
        assert_eq!(stage.into_output(), b"cell_barcode\tname\tumi_counts\n");
    }

    #[test]
    fn forwards_reads_without_an_extracted_barcode() {
        let mut stage = FloatingBarcodeStage::new(
            true,
            extractor(),
            table(),
            BarcodeCorrectOptions::default(),
            Vec::new(),
        )
        .unwrap();

        let forwarded = stage.process(vec![record(b"GATTGATTGATT")]).unwrap();
        assert_eq!(forwarded.len(), 1);
        stage.finish().unwrap();
        assert_eq!(stage.into_output(), b"cell_barcode\tname\tumi_counts\n");
    }
}
