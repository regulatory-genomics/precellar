//! FASTQ processing middleware.

use crate::align::AnnotatedFastq;
use crate::barcode::{BarcodeCorrectOptions, Whitelist, WhitelistBuilder};
pub use crate::pipeline::{AnnotatedFastqBatch, FastqStage, MiddlewareQcReport};
use crate::utils::insertion_extractor::InsertionExtractor;
use anyhow::Result;
use serde_json::json;
use std::io::Write;

#[derive(Debug, Default)]
pub struct FloatingBarcodeQc {
    pub num_records: u64,
    pub num_extracted: u64,
    pub num_forwarded: u64,
    pub num_corrected: u64,
    pub num_skipped: u64,
}

impl FloatingBarcodeQc {
    pub fn to_json(&self) -> serde_json::Value {
        json!({
            "num_records": self.num_records,
            "num_extracted": self.num_extracted,
            "num_forwarded": self.num_forwarded,
            "num_corrected": self.num_corrected,
            "num_skipped": self.num_skipped,
        })
    }
}

struct PendingFloatingBarcode {
    barcode: Vec<u8>,
    umi: Vec<u8>,
    sequence: Vec<u8>,
    quality_scores: Vec<u8>,
}

/// Extracts floating barcodes and removes extracted records from downstream.
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
    ) -> Self
    where
        B: IntoIterator<Item = S>,
        S: Into<Vec<u8>>,
    {
        Self {
            use_read1,
            extractor,
            whitelist_builder: Some(Whitelist::builder(valid_barcodes)),
            correction_options,
            pending: Vec::new(),
            output,
            qc: FloatingBarcodeQc::default(),
            finished: false,
        }
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
                .and_then(|barcode| {
                    barcode
                        .corrected
                        .as_deref()
                        .or(Some(barcode.raw.sequence()))
                })
                .unwrap_or_default()
                .to_vec(),
            umi: record
                .umi
                .as_ref()
                .map(|umi| umi.sequence().to_vec())
                .unwrap_or_default(),
            sequence: Vec::new(),
            quality_scores: Vec::new(),
        }
    }

    fn write_result(
        output: &mut W,
        record: &PendingFloatingBarcode,
        corrected: &[u8],
    ) -> std::io::Result<()> {
        output.write_all(&record.barcode)?;
        output.write_all(b"\t")?;
        output.write_all(&record.umi)?;
        output.write_all(b"\t")?;
        output.write_all(corrected)?;
        output.write_all(b"\n")
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
                self.qc.num_forwarded += 1;
                forwarded.push(record);
                continue;
            };
            let Some(range) = self.extractor.extract_range(sequence) else {
                self.qc.num_forwarded += 1;
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
        for record in std::mem::take(&mut self.pending) {
            match whitelist.correct_barcode(
                &record.sequence,
                &record.quality_scores,
                &self.correction_options,
            ) {
                Ok(corrected) => {
                    Self::write_result(&mut self.output, &record, corrected)?;
                    self.qc.num_corrected += 1;
                }
                Err(_) => self.qc.num_skipped += 1,
            }
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
    use noodles::fastq;

    fn record(sequence: &[u8]) -> AnnotatedFastq {
        let quality_scores = vec![b'I'; sequence.len()];
        AnnotatedFastq {
            barcode: None,
            umi: None,
            read1: Some(fastq::Record::new(
                fastq::record::Definition::default(),
                sequence,
                quality_scores,
            )),
            read2: None,
        }
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
        );

        let forwarded = stage
            .process(vec![record(b"ACGTCAGTGGCAACGTTTGGAACCTTGG")])
            .unwrap();
        assert!(forwarded.is_empty());
        stage.finish().unwrap();
        assert_eq!(stage.into_output(), b"\t\tACGT\n");
    }

    #[test]
    fn removes_extracted_reads_when_barcode_is_not_correctable() {
        let mut stage = FloatingBarcodeStage::new(
            true,
            extractor(),
            [b"GGGG".to_vec()],
            BarcodeCorrectOptions::default(),
            Vec::new(),
        );

        let forwarded = stage
            .process(vec![record(b"ACGTCAGTGGCAACGTTTGGAACCTTGG")])
            .unwrap();
        assert!(forwarded.is_empty());
        stage.finish().unwrap();
        assert!(stage.into_output().is_empty());
    }

    #[test]
    fn forwards_reads_without_an_extracted_barcode() {
        let mut stage = FloatingBarcodeStage::new(
            true,
            extractor(),
            [b"ACGT".to_vec()],
            BarcodeCorrectOptions::default(),
            Vec::new(),
        );

        let forwarded = stage.process(vec![record(b"GATTGATTGATT")]).unwrap();
        assert_eq!(forwarded.len(), 1);
        stage.finish().unwrap();
        assert!(stage.into_output().is_empty());
    }
}
