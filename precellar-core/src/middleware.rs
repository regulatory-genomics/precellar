//! FASTQ processing middleware.

use crate::align::AnnotatedFastq;
use crate::barcode::{BarcodeCorrectOptions, Whitelist, WhitelistBuilder};
pub use crate::pipeline::{AnnotatedFastqBatch, FastqStage, MiddlewareQcReport};
use crate::utils::insertion_extractor::InsertionExtractor;
use anyhow::Result;
use serde_json::json;
use std::collections::BTreeMap;
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

/// Extracts floating barcodes and removes extracted records from downstream.
///
/// Floating-barcode correction uses an explicit, non-empty whitelist. Every
/// extracted floating sequence contributes to whitelist frequencies, but output
/// is emitted only when the cell barcode and at least one floating barcode have
/// corrected values. Unavailable or uncorrectable floating barcodes are `-`.
pub struct FloatingBarcodeStage<W> {
    use_read1: bool,
    extractor: InsertionExtractor,
    whitelist_builders: Option<Vec<WhitelistBuilder>>,
    correction_options: BarcodeCorrectOptions,
    pending: Vec<PendingFloatingBarcode>,
    output: W,
    qc: FloatingBarcodeQc,
    finished: bool,
}

impl<W> FloatingBarcodeStage<W>
where
    W: Write + Send,
{
    pub fn new<B, V, S>(
        use_read1: bool,
        extractor: InsertionExtractor,
        valid_barcodes: B,
        correction_options: BarcodeCorrectOptions,
        output: W,
    ) -> Result<Self>
    where
        B: IntoIterator<Item = S>,
        S: IntoIterator<Item = V>,
        V: Into<Vec<u8>>,
    {
        let valid_barcodes: Vec<Vec<Vec<u8>>> = valid_barcodes
            .into_iter()
            .map(|barcodes| barcodes.into_iter().map(Into::into).collect())
            .collect();
        if valid_barcodes.len() != extractor.num_gaps() {
            anyhow::bail!(
                "valid_barcodes must contain one whitelist per insertion gap (expected {}, got {})",
                extractor.num_gaps(),
                valid_barcodes.len()
            );
        }
        if valid_barcodes.iter().any(Vec::is_empty) {
            anyhow::bail!("each insertion gap must have a non-empty barcode whitelist");
        }

        Ok(Self {
            use_read1,
            extractor,
            whitelist_builders: Some(valid_barcodes.into_iter().map(Whitelist::builder).collect()),
            correction_options,
            pending: Vec::new(),
            output,
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

    fn pending_record(
        record: &AnnotatedFastq,
        insertions: Vec<PendingInsertion>,
    ) -> Option<PendingFloatingBarcode> {
        Some(PendingFloatingBarcode {
            barcode: record.barcode.as_ref()?.corrected.clone()?,
            umi: record
                .umi
                .as_ref()
                .map(|umi| umi.sequence().to_vec())
                .unwrap_or_default(),
            insertions,
        })
    }

    fn write_results(
        output: &mut W,
        num_gaps: usize,
        matches: BTreeMap<(Vec<u8>, Vec<Vec<u8>>), BTreeMap<Vec<u8>, u64>>,
    ) -> std::io::Result<()> {
        output.write_all(b"cell_barcode")?;
        for index in 0..num_gaps {
            write!(output, "\tfloating_barcode{}", index + 1)?;
        }
        output.write_all(b"\tumi_counts\n")?;
        for ((cell_barcode, floating_barcodes), umis) in matches {
            output.write_all(&cell_barcode)?;
            for floating_barcode in floating_barcodes {
                output.write_all(b"\t")?;
                output.write_all(&floating_barcode)?;
            }
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
}

impl<W> FastqStage for FloatingBarcodeStage<W>
where
    W: Write + Send,
{
    fn process(&mut self, batch: AnnotatedFastqBatch) -> Result<AnnotatedFastqBatch> {
        let mut forwarded = Vec::with_capacity(batch.len());
        for record in batch {
            self.qc.num_records += 1;
            let Some((sequence, quality_scores)) = Self::selected_sequence(self.use_read1, &record)
            else {
                forwarded.push(record);
                continue;
            };
            if sequence.len() != quality_scores.len() {
                anyhow::bail!("selected FASTQ sequence and quality lengths differ");
            }
            let Some(ranges) = self.extractor.extract_ranges(sequence) else {
                forwarded.push(record);
                continue;
            };
            if ranges.is_empty() {
                forwarded.push(record);
                continue;
            }

            let retain_pending = record
                .barcode
                .as_ref()
                .and_then(|barcode| barcode.corrected.as_ref())
                .is_some();
            let mut insertions = retain_pending.then(|| Vec::with_capacity(ranges.len()));
            for (insertion_index, range) in ranges {
                self.qc.num_extracted += 1;
                let extracted = &sequence[range.clone()];
                self.whitelist_builders
                    .as_mut()
                    .expect("floating barcode stage already finalized")
                    .get_mut(insertion_index)
                    .expect("extractor returned an invalid insertion index")
                    .add(extracted);
                if let Some(insertions) = insertions.as_mut() {
                    insertions.push(PendingInsertion {
                        insertion_index,
                        sequence: extracted.to_vec(),
                        quality_scores: quality_scores[range].to_vec(),
                    });
                }
            }
            if let Some(insertions) = insertions {
                if let Some(pending) = Self::pending_record(&record, insertions) {
                    self.pending.push(pending);
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
        let num_gaps = whitelists.len();
        let mut matches: BTreeMap<_, BTreeMap<Vec<u8>, u64>> = BTreeMap::new();
        for record in std::mem::take(&mut self.pending) {
            let PendingFloatingBarcode {
                barcode,
                umi,
                insertions,
            } = record;
            let mut floating_barcodes = vec![b"-".to_vec(); num_gaps];
            let mut matched = false;
            for insertion in insertions {
                let PendingInsertion {
                    insertion_index,
                    sequence,
                    quality_scores,
                } = insertion;
                let Some(whitelist) = whitelists.get(insertion_index) else {
                    continue;
                };
                match whitelist.correct_barcode(
                    &sequence,
                    &quality_scores,
                    &self.correction_options,
                ) {
                    Ok(corrected) => {
                        floating_barcodes[insertion_index] = corrected.to_vec();
                        matched = true;
                    }
                    Err(_) => {}
                }
            }
            if matched {
                *matches
                    .entry((barcode, floating_barcodes))
                    .or_default()
                    .entry(umi)
                    .or_default() += 1;
                self.qc.num_matched += 1;
            }
        }
        Self::write_results(&mut self.output, num_gaps, matches)?;
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
            vec!["ACGTCAGTGGCA".into(), "TTGGAACCTTGG".into()],
            vec![4],
            12,
            0,
        )
        .unwrap()
    }

    #[test]
    fn removes_extracted_reads_and_writes_corrected_barcode() {
        let mut stage = FloatingBarcodeStage::new(
            true,
            extractor(),
            [[b"ACGT".to_vec()]],
            BarcodeCorrectOptions::default(),
            Vec::new(),
        )
        .unwrap();

        let forwarded = stage
            .process(vec![record(b"ACGTCAGTGGCAACGTTTGGAACCTTGG")])
            .unwrap();
        assert!(forwarded.is_empty());
        stage.finish().unwrap();
        assert_eq!(
            stage.into_output(),
            b"cell_barcode\tfloating_barcode1\tumi_counts\nCELL\tACGT\t:1\n"
        );
    }

    #[test]
    fn removes_extracted_reads_when_barcode_is_not_correctable() {
        let mut stage = FloatingBarcodeStage::new(
            true,
            extractor(),
            [[b"GGGG".to_vec()]],
            BarcodeCorrectOptions::default(),
            Vec::new(),
        )
        .unwrap();

        let forwarded = stage
            .process(vec![record(b"ACGTCAGTGGCAACGTTTGGAACCTTGG")])
            .unwrap();
        assert!(forwarded.is_empty());
        stage.finish().unwrap();
        assert_eq!(
            stage.into_output(),
            b"cell_barcode\tfloating_barcode1\tumi_counts\n"
        );
    }

    #[test]
    fn skips_matches_without_a_corrected_cell_barcode() {
        let mut stage = FloatingBarcodeStage::new(
            true,
            extractor(),
            [[b"ACGT".to_vec()]],
            BarcodeCorrectOptions::default(),
            Vec::new(),
        )
        .unwrap();
        let mut record = record(b"ACGTCAGTGGCAACGTTTGGAACCTTGG");
        record.barcode.as_mut().unwrap().corrected = None;

        stage.process(vec![record]).unwrap();
        stage.finish().unwrap();
        assert_eq!(
            stage.into_output(),
            b"cell_barcode\tfloating_barcode1\tumi_counts\n"
        );
    }

    #[test]
    fn groups_matches_by_cell_barcode_and_umi() {
        let mut stage = FloatingBarcodeStage::new(
            true,
            extractor(),
            [[b"ACGT".to_vec()]],
            BarcodeCorrectOptions::default(),
            Vec::new(),
        )
        .unwrap();
        let sequence = b"ACGTCAGTGGCAACGTTTGGAACCTTGG";
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
            b"cell_barcode\tfloating_barcode1\tumi_counts\nAAAA\tACGT\t:1\nCCCC\tACGT\tUMI1:1;UMI2:2\n"
        );
    }

    #[test]
    fn rejects_an_empty_whitelist() {
        assert!(FloatingBarcodeStage::new(
            true,
            extractor(),
            Vec::<Vec<Vec<u8>>>::new(),
            BarcodeCorrectOptions::default(),
            Vec::new(),
        )
        .is_err());
    }

    #[test]
    fn writes_multiple_insertions_as_one_combined_key() {
        let extractor = InsertionExtractor::new(
            vec!["AAAA".into(), "CCCC".into(), "GGGG".into()],
            vec![2, 3],
            3,
            0,
        )
        .unwrap();
        let mut stage = FloatingBarcodeStage::new(
            true,
            extractor,
            [[b"TT".to_vec()], [b"AAA".to_vec()]],
            BarcodeCorrectOptions::default(),
            Vec::new(),
        )
        .unwrap();

        let forwarded = stage.process(vec![record(b"AAAATTCCCCAAAGGGG")]).unwrap();
        assert!(forwarded.is_empty());
        assert_eq!(stage.qc().num_extracted, 2);
        stage.finish().unwrap();
        assert_eq!(
            stage.into_output(),
            b"cell_barcode\tfloating_barcode1\tfloating_barcode2\tumi_counts\nCELL\tTT\tAAA\t:1\n"
        );
    }

    #[test]
    fn writes_dash_for_an_unavailable_insertion() {
        let extractor = InsertionExtractor::new(
            vec!["AAAA".into(), "CCCC".into(), "GGGG".into()],
            vec![2, 3],
            3,
            0,
        )
        .unwrap();
        let mut stage = FloatingBarcodeStage::new(
            true,
            extractor,
            [[b"TT".to_vec()], [b"AAA".to_vec()]],
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
            b"cell_barcode\tfloating_barcode1\tfloating_barcode2\tumi_counts\nCELL\t-\tAAA\t:1\nCELL\tTT\t-\t:1\n"
        );
    }

    #[test]
    fn writes_a_read_when_at_least_one_insertion_corrects() {
        let extractor = InsertionExtractor::new(
            vec!["AAAA".into(), "CCCC".into(), "GGGG".into()],
            vec![2, 3],
            3,
            0,
        )
        .unwrap();
        let mut stage = FloatingBarcodeStage::new(
            true,
            extractor,
            [[b"TT".to_vec()], [b"GGG".to_vec()]],
            BarcodeCorrectOptions::default(),
            Vec::new(),
        )
        .unwrap();

        stage.process(vec![record(b"AAAATTCCCCTTTGGGG")]).unwrap();
        stage.finish().unwrap();
        assert_eq!(stage.qc().num_extracted, 2);
        assert_eq!(stage.qc().num_matched, 1);
        assert_eq!(
            stage.into_output(),
            b"cell_barcode\tfloating_barcode1\tfloating_barcode2\tumi_counts\nCELL\tTT\t-\t:1\n"
        );
    }

    #[test]
    fn rejects_an_empty_per_gap_whitelist() {
        assert!(FloatingBarcodeStage::new(
            true,
            extractor(),
            [Vec::<Vec<u8>>::new()],
            BarcodeCorrectOptions::default(),
            Vec::new(),
        )
        .is_err());
    }

    #[test]
    fn forwards_reads_without_an_extracted_barcode() {
        let mut stage = FloatingBarcodeStage::new(
            true,
            extractor(),
            [[b"ACGT".to_vec()]],
            BarcodeCorrectOptions::default(),
            Vec::new(),
        )
        .unwrap();

        let forwarded = stage.process(vec![record(b"GATTGATTGATT")]).unwrap();
        assert_eq!(forwarded.len(), 1);
        stage.finish().unwrap();
        assert_eq!(
            stage.into_output(),
            b"cell_barcode\tfloating_barcode1\tumi_counts\n"
        );
    }
}
