use std::path::PathBuf;

use anyhow::{Context, Result};
use bitcode::{Decode, Encode};
use noodles_sam::alignment::record::cigar::op::Kind;
use noodles_sam::alignment::record::data::field::Tag;
use noodles_sam::alignment::record::Cigar;
use noodles_sam::alignment::record_buf::data::field::Value;
use noodles_sam::alignment::RecordBuf;

use crate::barcode::get_barcode;

/// Which alignments of a record enter the UCF payload.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AlignmentPolicy {
    /// Store the primary alignment only; count the rest as excluded.
    Primary,
    /// Store every mapped alignment, including secondary and supplementary.
    All,
}

/// What to do with a mapped record that has no cell barcode.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MissingBarcodePolicy {
    /// Skip the record and count it.
    Drop,
    /// Fail the conversion.
    Error,
}

/// One cell-aware RNA alignment, already reduced to what UCF stores.
#[derive(Encode, Decode, Debug, Clone, PartialEq, Eq)]
pub struct ScRnaRow {
    pub chrom_id: u32,
    pub start: u32,
    pub mate_chrom_id: Option<u32>,
    pub mate_start: Option<u32>,
    pub template_length: i32,
    pub flags: u16,
    pub nh: u16,
    pub mapq: u8,
    pub cell_id: u32,
    pub umi_id: u32,
    pub reference_span: u32,
    /// Number of aligned blocks derived from `cigar`.
    pub block_count: u16,
    /// BAM-compatible `(operation_length, operation_code)` pairs.
    pub cigar: Vec<(u32, u8)>,
}

/// What the converter did with one alignment record.
#[derive(Debug)]
pub enum ScRnaDisposition {
    Stored(ScRnaRow),
    Unmapped,
    Secondary,
    Supplementary,
    MissingBarcode,
    MissingReference,
    EmptyCigar,
}

/// Counts of stored and excluded alignment records.
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub struct ScRnaRetention {
    pub stored: u64,
    pub unmapped: u64,
    pub secondary: u64,
    pub supplementary: u64,
    pub missing_barcode: u64,
    pub missing_reference: u64,
    pub empty_cigar: u64,
}

impl ScRnaRetention {
    pub fn record(&mut self, disposition: &ScRnaDisposition) {
        match disposition {
            ScRnaDisposition::Stored(_) => self.stored += 1,
            ScRnaDisposition::Unmapped => self.unmapped += 1,
            ScRnaDisposition::Secondary => self.secondary += 1,
            ScRnaDisposition::Supplementary => self.supplementary += 1,
            ScRnaDisposition::MissingBarcode => self.missing_barcode += 1,
            ScRnaDisposition::MissingReference => self.missing_reference += 1,
            ScRnaDisposition::EmptyCigar => self.empty_cigar += 1,
        }
    }

    pub fn excluded(&self) -> u64 {
        self.unmapped
            + self.secondary
            + self.supplementary
            + self.missing_barcode
            + self.missing_reference
            + self.empty_cigar
    }
}

/// Reduce one alignment record to a UCF row, or classify why it is excluded.
pub fn classify_record(
    record: &RecordBuf,
    mut cell_id: impl FnMut(&str) -> u32,
    mut umi_id: impl FnMut(&str) -> u32,
    alignments: AlignmentPolicy,
    missing_barcode: MissingBarcodePolicy,
) -> Result<ScRnaDisposition> {
    let flags = record.flags();
    if flags.is_unmapped() {
        return Ok(ScRnaDisposition::Unmapped);
    }
    if alignments == AlignmentPolicy::Primary {
        if flags.is_secondary() {
            return Ok(ScRnaDisposition::Secondary);
        }
        if flags.is_supplementary() {
            return Ok(ScRnaDisposition::Supplementary);
        }
    }
    let Some(barcode) = get_barcode(record)? else {
        if missing_barcode == MissingBarcodePolicy::Error {
            let name = record
                .name()
                .map(|name| name.to_string())
                .unwrap_or_else(|| "<unnamed>".to_owned());
            anyhow::bail!("mapped read {name} has no cell barcode tag");
        }
        return Ok(ScRnaDisposition::MissingBarcode);
    };
    let Some(reference_id) = record.reference_sequence_id() else {
        return Ok(ScRnaDisposition::MissingReference);
    };
    let Some(start) = record.alignment_start() else {
        return Ok(ScRnaDisposition::MissingReference);
    };
    let mut cigar = Vec::new();
    let mut reference_span = 0u32;
    for op in record.cigar().iter() {
        let op = op?;
        let len = u32::try_from(op.len())?;
        cigar.push((len, cigar_code(op.kind())));
        if matches!(
            op.kind(),
            Kind::Match | Kind::Deletion | Kind::Skip | Kind::SequenceMatch | Kind::SequenceMismatch
        ) {
            reference_span = reference_span.saturating_add(len);
        }
    }
    if cigar.is_empty() {
        return Ok(ScRnaDisposition::EmptyCigar);
    }
    let mate = match (
        record.mate_reference_sequence_id(),
        record.mate_alignment_start(),
    ) {
        (Some(chrom), Some(pos)) => (Some(chrom), Some(pos)),
        _ => (None, None),
    };
    let umi = record
        .data()
        .get(&Tag::UMI_SEQUENCE)
        .and_then(|value| match value {
            Value::String(umi) => Some(umi.to_string()),
            _ => None,
        });
    let umi_id = match umi {
        Some(umi) => umi_id(&umi),
        None => u32::MAX,
    };
    Ok(ScRnaDisposition::Stored(ScRnaRow {
        chrom_id: u32::try_from(reference_id)?,
        start: u32::try_from(start.get() - 1)?,
        mate_chrom_id: mate.0.map(u32::try_from).transpose()?,
        mate_start: mate.1.map(|pos| u32::try_from(pos.get() - 1)).transpose()?,
        template_length: record.template_length(),
        flags: flags.bits(),
        nh: read_nh(record)?,
        mapq: record.mapping_quality().map_or(255, |value| value.get()),
        cell_id: cell_id(&barcode),
        umi_id,
        reference_span,
        block_count: count_blocks(&cigar)?,
        cigar,
    }))
}

/// Count aligned blocks the same way `ucf-rna` derives them from a CIGAR:
/// `M/=/X` extend coverage; `D/N` split blocks.
fn count_blocks(cigar: &[(u32, u8)]) -> Result<u16> {
    let mut blocks = 0u16;
    let mut inside = false;
    for &(_, code) in cigar {
        match code {
            0 | 7 | 8 => {
                if !inside {
                    blocks = blocks
                        .checked_add(1)
                        .context("scRNA alignment has too many blocks")?;
                    inside = true;
                }
            }
            2 | 3 => inside = false,
            _ => {}
        }
    }
    Ok(blocks)
}

/// Sort rows into the genomic order the UCF writer validates.
pub fn sort_rows(
    rows: impl Iterator<Item = ScRnaRow>,
    temp_dir: Option<PathBuf>,
    chunk_size: usize,
) -> anyhow::Result<impl ExactSizeIterator<Item = anyhow::Result<ScRnaRow>>> {
    super::external_sort(rows, temp_dir, chunk_size, |a, b| {
        a.chrom_id
            .cmp(&b.chrom_id)
            .then_with(|| a.start.cmp(&b.start))
            .then_with(|| a.cell_id.cmp(&b.cell_id))
            .then_with(|| a.umi_id.cmp(&b.umi_id))
    })
}

fn read_nh(record: &RecordBuf) -> Result<u16> {
    let Some(value) = record.data().get(&Tag::ALIGNMENT_HIT_COUNT) else {
        return Ok(0);
    };
    let value = match value {
        Value::Int8(n) => i64::from(*n),
        Value::UInt8(n) => i64::from(*n),
        Value::Int16(n) => i64::from(*n),
        Value::UInt16(n) => i64::from(*n),
        Value::Int32(n) => i64::from(*n),
        Value::UInt32(n) => i64::from(*n),
        _ => return Ok(0),
    };
    Ok(u16::try_from(value.max(0)).unwrap_or(0))
}

const fn cigar_code(kind: Kind) -> u8 {
    match kind {
        Kind::Match => 0,
        Kind::Insertion => 1,
        Kind::Deletion => 2,
        Kind::Skip => 3,
        Kind::SoftClip => 4,
        Kind::HardClip => 5,
        Kind::Pad => 6,
        Kind::SequenceMatch => 7,
        Kind::SequenceMismatch => 8,
    }
}
