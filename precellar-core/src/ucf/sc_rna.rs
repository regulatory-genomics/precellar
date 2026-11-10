use std::path::PathBuf;

use anyhow::{Context, Result};
use bitcode::{Decode, Encode};
use noodles_sam::alignment::record::cigar::op::Kind;
use noodles_sam::alignment::record::data::field::Tag;
use noodles_sam::alignment::record::Cigar;
use noodles_sam::alignment::record_buf::data::field::Value;
use noodles_sam::alignment::RecordBuf;

use crate::barcode::get_barcode;

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

/// Reduce a primary alignment to a UCF row. Returns `None` for unmapped or
/// secondary/supplementary records and for records without a cell barcode.
pub fn row_from_record(
    record: &RecordBuf,
    mut cell_id: impl FnMut(&str) -> u32,
    mut umi_id: impl FnMut(&str) -> u32,
) -> Result<Option<ScRnaRow>> {
    let flags = record.flags();
    if flags.is_unmapped() {
        return Ok(None);
    }
    let Some(barcode) = get_barcode(record)? else {
        return Ok(None);
    };
    let Some(reference_id) = record.reference_sequence_id() else {
        return Ok(None);
    };
    let Some(start) = record.alignment_start() else {
        return Ok(None);
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
        return Ok(None);
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
        Some(umi) => {
            let umi = umi.to_string();
            umi_id(&umi)
        }
        None => u32::MAX,
    };
    Ok(Some(ScRnaRow {
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
) -> impl ExactSizeIterator<Item = ScRnaRow> {
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
