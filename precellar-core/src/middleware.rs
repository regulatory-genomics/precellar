//! FASTQ processing middleware.

use crate::align::AnnotatedFastq;
use crate::barcode::{BarcodeCorrectOptions, Whitelist, WhitelistBuilder};
pub use crate::pipeline::{AnnotatedFastqBatch, FastqStage, MiddlewareQcReport};
use crate::utils::insertion_extractor::InsertionExtractor;
use anyhow::Result;
use rayon::prelude::*;
use serde_json::json;
use std::collections::{BTreeMap, HashMap};
use std::io::Write;

#[derive(Debug, Default)]
pub struct FloatingBarcodeQc {
    pub num_records: u64,
    pub num_extracted: u64,
    pub num_matched: u64,
}

impl FloatingBarcodeQc {
    pub fn to_json(&self) -> serde_json::Value {
        json!({
            "num_records": self.num_records,
            "num_extracted": self.num_extracted,
            "num_matched": self.num_matched,
        })
    }
}

struct PendingInsertion {
    insertion_index: usize,
    sequence: Vec<u8>,
    quality_scores: Vec<u8>,
}

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
    pub name_1: Vec<u8>,
    pub name_2: Vec<u8>,
    pub barcode_1: Vec<u8>,
    pub barcode_2: Vec<u8>,
}

impl FloatingBarcodeEntry {
    pub fn new<N1, N2, B1, B2>(name_1: N1, name_2: N2, barcode_1: B1, barcode_2: B2) -> Self
    where
        N1: Into<Vec<u8>>,
        N2: Into<Vec<u8>>,
        B1: Into<Vec<u8>>,
        B2: Into<Vec<u8>>,
    {
        Self {
            name_1: name_1.into(),
            name_2: name_2.into(),
            barcode_1: barcode_1.into(),
            barcode_2: barcode_2.into(),
        }
    }
}

#[derive(Clone)]
pub struct FloatingBarcodeTable {
    entries: Vec<FloatingBarcodeEntry>,
    expected_lens: [usize; 2],
    barcode_1_rows: HashMap<Vec<u8>, Option<usize>>,
    barcode_2_rows: HashMap<Vec<u8>, Option<usize>>,
    pair_rows: HashMap<Vec<u8>, HashMap<Vec<u8>, usize>>,
}

impl FloatingBarcodeTable {
    pub fn new(mut entries: Vec<FloatingBarcodeEntry>, expected_lens: &[usize]) -> Result<Self> {
        if expected_lens.len() != 2 {
            anyhow::bail!("barcode tables require exactly two insertion gaps");
        }
        if entries.is_empty() {
            anyhow::bail!("barcode_table must contain at least one row");
        }

        let expected_lens = [expected_lens[0], expected_lens[1]];
        let mut barcode_1_rows = HashMap::new();
        let mut barcode_2_rows = HashMap::new();
        let mut pair_rows: HashMap<Vec<u8>, HashMap<Vec<u8>, usize>> = HashMap::new();
        for (index, entry) in entries.iter_mut().enumerate() {
            if entry.name_1.is_empty() || entry.name_2.is_empty() {
                anyhow::bail!("barcode table names must not be empty");
            }
            if entry
                .name_1
                .iter()
                .chain(&entry.name_2)
                .any(|byte| matches!(byte, b'\t' | b'\n' | b'\r' | 0..=31))
            {
                anyhow::bail!("barcode table names must not contain control characters");
            }
            entry.barcode_1.make_ascii_uppercase();
            entry.barcode_2.make_ascii_uppercase();
            if !is_dna(&entry.barcode_1) || !is_dna(&entry.barcode_2) {
                anyhow::bail!("barcode table barcodes must be non-empty DNA sequences");
            }
            if entry.barcode_1.len() != expected_lens[0]
                || entry.barcode_2.len() != expected_lens[1]
            {
                anyhow::bail!(
                    "barcode table row {} has lengths ({}, {}), expected ({}, {})",
                    index + 1,
                    entry.barcode_1.len(),
                    entry.barcode_2.len(),
                    expected_lens[0],
                    expected_lens[1]
                );
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
            expected_lens,
            barcode_1_rows,
            barcode_2_rows,
            pair_rows,
        })
    }

    pub fn expected_lens(&self) -> &[usize; 2] {
        &self.expected_lens
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
    pending: Vec<PendingFloatingBarcode>,
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
        if extractor.expected_lens() != &barcode_table.expected_lens()[..] {
            anyhow::bail!("barcode table lengths do not match extractor insertion lengths");
        }
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
            pending: Vec::new(),
            output,
            thread_pool: rayon::ThreadPoolBuilder::new().num_threads(1).build()?,
            qc: FloatingBarcodeQc::default(),
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

    fn write_results(
        output: &mut W,
        matches: BTreeMap<(Vec<u8>, Vec<u8>, Vec<u8>), BTreeMap<Vec<u8>, u64>>,
    ) -> std::io::Result<()> {
        output.write_all(b"cell_barcode\tname_1\tname_2\tumi_counts\n")?;
        for ((cell_barcode, name_1, name_2), umis) in matches {
            output.write_all(&cell_barcode)?;
            output.write_all(b"\t")?;
            output.write_all(&name_1)?;
            output.write_all(b"\t")?;
            output.write_all(&name_2)?;
            output.write_all(b"\t")?;
            for (index, (umi, count)) in umis.into_iter().enumerate() {
                if index > 0 {
                    output.write_all(b";")?;
                }
                output.write_all(&umi)?;
                write!(output, ":{count}")?;
            }
            output.write_all(b"\n")?;
        }
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
    ) -> Option<(Vec<u8>, Vec<u8>, Vec<u8>, Vec<u8>)> {
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
        Some((barcode, entry.name_1.clone(), entry.name_2.clone(), umi))
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
                    for insertion in &insertions {
                        self.whitelist_builders
                            .as_mut()
                            .expect("floating barcode stage already finalized")
                            .get_mut(insertion.insertion_index)
                            .expect("extractor returned an invalid insertion index")
                            .add(&insertion.sequence);
                    }
                    if let Some(barcode) = barcode {
                        self.pending.push(PendingFloatingBarcode {
                            barcode,
                            umi,
                            insertions,
                        });
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
        let pending = std::mem::take(&mut self.pending);
        let whitelists = &whitelists;
        let barcode_table = &self.barcode_table;
        let correction_options = &self.correction_options;
        let resolved: Vec<_> = self.thread_pool.install(|| {
            pending
                .into_par_iter()
                .filter_map(|record| {
                    Self::resolve_pending(record, whitelists, barcode_table, correction_options)
                })
                .collect()
        });
        self.qc.num_matched += resolved.len() as u64;
        let mut matches: BTreeMap<_, BTreeMap<Vec<u8>, u64>> = BTreeMap::new();
        for (barcode, name_1, name_2, umi) in resolved {
            *matches
                .entry((barcode, name_1, name_2))
                .or_default()
                .entry(umi)
                .or_default() += 1;
        }
        Self::write_results(&mut self.output, matches)?;
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
    use noodles::fastq;

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
            vec![2, 3],
            3,
            0.0,
        )
        .unwrap()
    }

    fn table() -> FloatingBarcodeTable {
        FloatingBarcodeTable::new(
            vec![FloatingBarcodeEntry::new("TF1", "SITE1", "TT", "AAA")],
            &[2, 3],
        )
        .unwrap()
    }

    fn full_sequence() -> &'static [u8] {
        b"AAAATTCCCCAAAGGGG"
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
        stage.finish().unwrap();
        assert_eq!(
            stage.into_output(),
            b"cell_barcode\tname_1\tname_2\tumi_counts\nCELL\tTF1\tSITE1\t:1\n"
        );
    }

    #[test]
    fn removes_extracted_reads_when_barcode_is_not_correctable() {
        let mut stage = FloatingBarcodeStage::new(
            true,
            extractor(),
            FloatingBarcodeTable::new(
                vec![FloatingBarcodeEntry::new("TF2", "SITE2", "GG", "CCC")],
                &[2, 3],
            )
            .unwrap(),
            BarcodeCorrectOptions::default(),
            Vec::new(),
        )
        .unwrap();

        let forwarded = stage.process(vec![record(full_sequence())]).unwrap();
        assert!(forwarded.is_empty());
        stage.finish().unwrap();
        assert_eq!(
            stage.into_output(),
            b"cell_barcode\tname_1\tname_2\tumi_counts\n"
        );
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
        stage.finish().unwrap();
        assert_eq!(
            stage.into_output(),
            b"cell_barcode\tname_1\tname_2\tumi_counts\n"
        );
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
            b"cell_barcode\tname_1\tname_2\tumi_counts\nAAAA\tTF1\tSITE1\t:1\nCCCC\tTF1\tSITE1\tUMI1:1;UMI2:2\n"
        );
    }

    fn run_with_threads(num_threads: usize) -> (Vec<u8>, Vec<Vec<u8>>, (u64, u64, u64)) {
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
        let qc = stage.qc();
        let metrics = (qc.num_records, qc.num_extracted, qc.num_matched);
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
            vec![2, 3],
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
        stage.finish().unwrap();
        assert_eq!(
            stage.into_output(),
            b"cell_barcode\tname_1\tname_2\tumi_counts\nCELL\tTF1\tSITE1\t:1\n"
        );
    }

    #[test]
    fn writes_dash_for_an_unavailable_insertion() {
        let extractor = InsertionExtractor::new(
            vec!["AAAA".into(), "CCCC".into(), "GGGG".into()],
            vec![2, 3],
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
        stage.finish().unwrap();
        assert_eq!(
            stage.into_output(),
            b"cell_barcode\tname_1\tname_2\tumi_counts\nCELL\tTF1\tSITE1\t:2\n"
        );
    }

    #[test]
    fn writes_a_read_when_at_least_one_insertion_corrects() {
        let extractor = InsertionExtractor::new(
            vec!["AAAA".into(), "CCCC".into(), "GGGG".into()],
            vec![2, 3],
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
            b"cell_barcode\tname_1\tname_2\tumi_counts\nCELL\tTF1\tSITE1\t:1\n"
        );
    }

    #[test]
    fn rejects_an_empty_barcode_table() {
        assert!(FloatingBarcodeTable::new(Vec::new(), &[2, 3]).is_err());
    }

    #[test]
    fn rejects_an_ambiguous_single_barcode_lookup() {
        let mut stage = FloatingBarcodeStage::new(
            true,
            extractor(),
            FloatingBarcodeTable::new(
                vec![
                    FloatingBarcodeEntry::new("TF1", "SITE1", "TT", "AAA"),
                    FloatingBarcodeEntry::new("TF2", "SITE2", "TT", "CCC"),
                ],
                &[2, 3],
            )
            .unwrap(),
            BarcodeCorrectOptions::default(),
            Vec::new(),
        )
        .unwrap();

        stage.process(vec![record(b"AAAATTCCCC")]).unwrap();
        stage.finish().unwrap();
        assert_eq!(stage.qc().num_matched, 0);
        assert_eq!(
            stage.into_output(),
            b"cell_barcode\tname_1\tname_2\tumi_counts\n"
        );
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
        assert_eq!(
            stage.into_output(),
            b"cell_barcode\tname_1\tname_2\tumi_counts\n"
        );
    }
}
