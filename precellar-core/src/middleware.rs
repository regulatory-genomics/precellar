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

struct PendingFloatingBarcode {
    barcode: Option<Vec<u8>>,
    umi: Vec<u8>,
    sequence: Vec<u8>,
    quality_scores: Vec<u8>,
}

/// Extracts floating barcodes and removes extracted records from downstream.
///
/// Floating-barcode correction uses an explicit, non-empty whitelist. Every
/// extracted floating sequence contributes to whitelist frequencies, but output
/// is emitted only when the cell barcode has a corrected value.
pub struct FloatingBarcodeStage<W> {
    use_read1: bool,
    extractor: InsertionExtractor,
    whitelist_builder: Option<WhitelistBuilder>,
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
    pub fn new<B, S>(
        use_read1: bool,
        extractor: InsertionExtractor,
        valid_barcodes: B,
        correction_options: BarcodeCorrectOptions,
        output: W,
    ) -> Result<Self>
    where
        B: IntoIterator<Item = S>,
        S: Into<Vec<u8>>,
    {
        let valid_barcodes: Vec<Vec<u8>> = valid_barcodes.into_iter().map(Into::into).collect();
        if valid_barcodes.is_empty() {
            anyhow::bail!("valid_barcodes must contain at least one barcode");
        }

        Ok(Self {
            use_read1,
            extractor,
            whitelist_builder: Some(Whitelist::builder(valid_barcodes)),
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

    fn pending_record(record: &AnnotatedFastq) -> PendingFloatingBarcode {
        PendingFloatingBarcode {
            barcode: record
                .barcode
                .as_ref()
                .and_then(|barcode| barcode.corrected.clone()),
            umi: record
                .umi
                .as_ref()
                .map(|umi| umi.sequence().to_vec())
                .unwrap_or_default(),
            sequence: Vec::new(),
            quality_scores: Vec::new(),
        }
    }

    fn write_results(
        output: &mut W,
        matches: BTreeMap<(Vec<u8>, Vec<u8>), BTreeMap<Vec<u8>, u64>>,
    ) -> std::io::Result<()> {
        output.write_all(b"cell_barcode\tfloating_barcode\tumi_counts\n")?;
        for ((cell_barcode, floating_barcode), umis) in matches {
            output.write_all(&cell_barcode)?;
            output.write_all(b"\t")?;
            output.write_all(&floating_barcode)?;
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
            let Some(range) = self.extractor.extract_range(sequence) else {
                forwarded.push(record);
                continue;
            };

            self.qc.num_extracted += 1;
            let mut pending = Self::pending_record(&record);
            pending.sequence = sequence[range.clone()].to_vec();
            pending.quality_scores = quality_scores[range].to_vec();
            self.whitelist_builder
                .as_mut()
                .expect("floating barcode stage already finalized")
                .add(&pending.sequence);
            self.pending.push(pending);
        }
        Ok(forwarded)
    }

    fn finish(&mut self) -> Result<()> {
        if self.finished {
            return Ok(());
        }

        let whitelist: Whitelist = self
            .whitelist_builder
            .take()
            .expect("floating barcode stage already finalized")
            .finish();
        let mut matches: BTreeMap<_, BTreeMap<Vec<u8>, u64>> = BTreeMap::new();
        for record in std::mem::take(&mut self.pending) {
            let PendingFloatingBarcode {
                barcode,
                umi,
                sequence,
                quality_scores,
            } = record;
            let Some(cell_barcode) = barcode else {
                continue;
            };
            match whitelist.correct_barcode(&sequence, &quality_scores, &self.correction_options) {
                Ok(corrected) => {
                    let floating_barcode = corrected.to_vec();
                    *matches
                        .entry((cell_barcode, floating_barcode))
                        .or_default()
                        .entry(umi)
                        .or_default() += 1;
                    self.qc.num_matched += 1;
                }
                Err(_) => {}
            }
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
        InsertionExtractor::new("ACGTCAGTGGCA", "TTGGAACCTTGG", 12, 4, 0)
    }

    #[test]
    fn removes_extracted_reads_and_writes_corrected_barcode() {
        let mut stage = FloatingBarcodeStage::new(
            true,
            extractor(),
            [b"ACGT".to_vec()],
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
            b"cell_barcode\tfloating_barcode\tumi_counts\nCELL\tACGT\t:1\n"
        );
    }

    #[test]
    fn removes_extracted_reads_when_barcode_is_not_correctable() {
        let mut stage = FloatingBarcodeStage::new(
            true,
            extractor(),
            [b"GGGG".to_vec()],
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
            b"cell_barcode\tfloating_barcode\tumi_counts\n"
        );
    }

    #[test]
    fn skips_matches_without_a_corrected_cell_barcode() {
        let mut stage = FloatingBarcodeStage::new(
            true,
            extractor(),
            [b"ACGT".to_vec()],
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
            b"cell_barcode\tfloating_barcode\tumi_counts\n"
        );
    }

    #[test]
    fn groups_matches_by_cell_barcode_and_umi() {
        let mut stage = FloatingBarcodeStage::new(
            true,
            extractor(),
            [b"ACGT".to_vec()],
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
            b"cell_barcode\tfloating_barcode\tumi_counts\nAAAA\tACGT\t:1\nCCCC\tACGT\tUMI1:1;UMI2:2\n"
        );
    }

    #[test]
    fn rejects_an_empty_whitelist() {
        assert!(FloatingBarcodeStage::new(
            true,
            extractor(),
            Vec::<Vec<u8>>::new(),
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
            [b"ACGT".to_vec()],
            BarcodeCorrectOptions::default(),
            Vec::new(),
        )
        .unwrap();

        let forwarded = stage.process(vec![record(b"GATTGATTGATT")]).unwrap();
        assert_eq!(forwarded.len(), 1);
        stage.finish().unwrap();
        assert_eq!(
            stage.into_output(),
            b"cell_barcode\tfloating_barcode\tumi_counts\n"
        );
    }
}
