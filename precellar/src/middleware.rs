//! FASTQ processing middleware.
//!
//! Middleware in this module consumes annotated FASTQ records, performs an
//! operation on them, and either forwards or removes each record.

use crate::align::AnnotatedFastq;
use crate::barcode::{BarcodeCorrectOptions, Whitelist, WhitelistBuilder};
use crate::utils::insertion_extractor::InsertionExtractor;
use std::io::{self, Write};

pub type AnnotatedFastqStream = Box<dyn Iterator<Item = Vec<AnnotatedFastq>> + Send>;
pub type MiddlewareFactory = Box<dyn FnOnce(AnnotatedFastqStream) -> AnnotatedFastqStream + Send>;

struct PendingFloatingBarcode {
    barcode: Vec<u8>,
    umi: Vec<u8>,
    sequence: Vec<u8>,
    quality_scores: Vec<u8>,
}

/// Extracts floating barcodes and removes extracted records from downstream.
///
/// Records without an extractable floating barcode are forwarded immediately.
/// Extracted records are removed immediately but retained in memory until the
/// input is exhausted, allowing their frequencies to be included in the
/// whitelist before correction. Only successfully corrected barcodes are
/// written as headerless, tab-separated `barcode`, `umi`, and floating-barcode
/// rows.
pub struct FloatingBarcodeExtracter<I, W> {
    input: I,
    use_read1: bool,
    extractor: InsertionExtractor,
    whitelist_builder: Option<WhitelistBuilder>,
    whitelist: Option<Whitelist>,
    correction_options: BarcodeCorrectOptions,
    pending: Vec<PendingFloatingBarcode>,
    output: W,
    finalized: bool,
    error: Option<io::Error>,
}

impl<I, W> FloatingBarcodeExtracter<I, W>
where
    I: Iterator<Item = Vec<AnnotatedFastq>>,
    W: Write,
{
    pub fn new<B, S>(
        input: I,
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
            input,
            use_read1,
            extractor,
            whitelist_builder: Some(Whitelist::builder(valid_barcodes)),
            whitelist: None,
            correction_options,
            pending: Vec::new(),
            output,
            finalized: false,
            error: None,
        }
    }

    /// Returns a write error encountered while processing the stream.
    pub fn error(&self) -> Option<&io::Error> {
        self.error.as_ref()
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
                .map(|x| (x.sequence(), x.quality_scores()))
        } else {
            record
                .read2
                .as_ref()
                .map(|x| (x.sequence(), x.quality_scores()))
        }
    }

    fn pending_record(record: &AnnotatedFastq) -> PendingFloatingBarcode {
        PendingFloatingBarcode {
            barcode: record
                .barcode
                .as_ref()
                .and_then(|x| x.corrected.as_deref().or(Some(x.raw.sequence())))
                .unwrap_or_default()
                .to_vec(),
            umi: record
                .umi
                .as_ref()
                .map(|x| x.sequence().to_vec())
                .unwrap_or_default(),
            sequence: Vec::new(),
            quality_scores: Vec::new(),
        }
    }

    fn write_result(
        output: &mut W,
        record: &PendingFloatingBarcode,
        corrected: &[u8],
    ) -> io::Result<()> {
        output.write_all(&record.barcode)?;
        output.write_all(b"\t")?;
        output.write_all(&record.umi)?;
        output.write_all(b"\t")?;
        output.write_all(corrected)?;
        output.write_all(b"\n")
    }

    fn finalize(&mut self) {
        if self.finalized {
            return;
        }
        self.finalized = true;

        let whitelist = self
            .whitelist_builder
            .take()
            .expect("floating barcode whitelist builder already finalized")
            .finish();
        self.whitelist = Some(whitelist);

        let pending = std::mem::take(&mut self.pending);
        let whitelist = self.whitelist.as_ref().unwrap();
        for record in pending {
            let corrected = whitelist.correct_barcode(
                &record.sequence,
                &record.quality_scores,
                &self.correction_options,
            );
            if let Ok(corrected) = corrected {
                if let Err(error) = Self::write_result(&mut self.output, &record, corrected) {
                    self.error = Some(error);
                    break;
                }
            }
        }
    }
}

impl<I, W> Iterator for FloatingBarcodeExtracter<I, W>
where
    I: Iterator<Item = Vec<AnnotatedFastq>>,
    W: Write,
{
    type Item = Vec<AnnotatedFastq>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.finalized {
            return None;
        }

        loop {
            let Some(chunk) = self.input.next() else {
                self.finalize();
                return None;
            };

            let mut forwarded = Vec::with_capacity(chunk.len());
            for record in chunk {
                let Some((sequence, quality_scores)) =
                    Self::selected_sequence(self.use_read1, &record)
                else {
                    forwarded.push(record);
                    continue;
                };
                let Some(range) = self.extractor.extract_range(sequence) else {
                    forwarded.push(record);
                    continue;
                };

                let mut pending = Self::pending_record(&record);
                pending.sequence = sequence[range.clone()].to_vec();
                pending.quality_scores = quality_scores[range].to_vec();
                self.whitelist_builder
                    .as_mut()
                    .expect("floating barcode whitelist builder already finalized")
                    .add(&pending.sequence);
                self.pending.push(pending);
            }

            if !forwarded.is_empty() {
                return Some(forwarded);
            }
        }
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
        let input = vec![vec![record(b"ACGTCAGTGGCAACGTTTGGAACCTTGG")]];
        let mut middleware = FloatingBarcodeExtracter::new(
            input.into_iter(),
            true,
            extractor(),
            [b"ACGT".to_vec()],
            BarcodeCorrectOptions::default(),
            Vec::new(),
        );

        assert!(middleware.next().is_none());
        assert_eq!(middleware.into_output(), b"\t\tACGT\n");
    }

    #[test]
    fn removes_extracted_reads_when_barcode_is_not_correctable() {
        let input = vec![vec![record(b"ACGTCAGTGGCAACGTTTGGAACCTTGG")]];
        let mut middleware = FloatingBarcodeExtracter::new(
            input.into_iter(),
            true,
            extractor(),
            [b"GGGG".to_vec()],
            BarcodeCorrectOptions::default(),
            Vec::new(),
        );

        assert!(middleware.next().is_none());
        assert!(middleware.into_output().is_empty());
    }

    #[test]
    fn forwards_reads_without_an_extracted_barcode() {
        let input = vec![vec![record(b"GATTGATTGATT")]];
        let mut middleware = FloatingBarcodeExtracter::new(
            input.into_iter(),
            true,
            extractor(),
            [b"ACGT".to_vec()],
            BarcodeCorrectOptions::default(),
            Vec::new(),
        );

        assert_eq!(middleware.next().unwrap().len(), 1);
        assert!(middleware.next().is_none());
        assert!(middleware.into_output().is_empty());
    }
}
