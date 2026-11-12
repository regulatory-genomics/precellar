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
    /// Base-level mismatches and indels, empty unless evidence is requested.
    pub events: Vec<RnaEvent>,
    /// True when evidence was requested but the record had no usable `MD:Z`.
    pub missing_md: bool,
}

/// One base-level difference between a read and the reference.
#[derive(Encode, Decode, Debug, Clone, Copy, PartialEq, Eq)]
pub struct RnaEvent {
    /// 0-based position in the original read sequence, soft clips included.
    pub read_pos: u32,
    /// 0-based offset from the alignment start.
    pub ref_offset: u32,
    /// `b'X'` substitution, `b'I'` insertion, `b'D'` deletion.
    pub op: u8,
    /// Reference base, `b'.'` when not applicable.
    pub ref_base: u8,
    /// Read base, `b'.'` when not applicable.
    pub read_base: u8,
    /// Phred score of the read base; 0 for deletions.
    pub qual: u8,
}

/// Evidence extracted from one alignment record.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct RnaEvidence {
    pub events: Vec<RnaEvent>,
    /// True when the record had no usable `MD:Z` tag.
    pub missing_md: bool,
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
    collect_evidence: bool,
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
            Kind::Match
                | Kind::Deletion
                | Kind::Skip
                | Kind::SequenceMatch
                | Kind::SequenceMismatch
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
    let (events, missing_md) = if collect_evidence {
        let evidence = rna_events(record)?;
        (evidence.events, evidence.missing_md)
    } else {
        (Vec::new(), false)
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
        events,
        missing_md,
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

/// Extract substitutions, insertions, and deletions from one alignment.
///
/// The reference path comes from the SAM `MD:Z` tag; the read path and base
/// qualities come from the record itself. Rows without `MD` still report
/// insertions, deletions, and explicit `X` CIGAR mismatches.
pub fn rna_events(record: &RecordBuf) -> Result<RnaEvidence> {
    use crate::align::snv::{parse_md_tag_str, MDKind};
    use std::collections::VecDeque;

    let mut evidence = RnaEvidence::default();
    if record.flags().is_unmapped() {
        return Ok(evidence);
    }
    let sequence = record.sequence().as_ref();
    let qualities = record.quality_scores().as_ref();
    let md = match record.data().get(&Tag::MISMATCHED_POSITIONS) {
        Some(Value::String(value)) => Some(VecDeque::from(
            parse_md_tag_str(value.as_ref()).collect::<Result<Vec<_>>>()?,
        )),
        _ => {
            evidence.missing_md = true;
            None
        }
    };
    let mut md = md;

    let mut read_pos = 0u32; // position in the original read
    let mut ref_offset = 0u32; // offset from the alignment start
    let mut aligned = false; // leading clips are consumed before this flips

    for op in record.cigar().iter() {
        let op = op?;
        let len = op.len();
        match op.kind() {
            Kind::SoftClip | Kind::HardClip | Kind::Pad => {
                if !aligned && op.kind() == Kind::SoftClip {
                    read_pos += len as u32;
                }
            }
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                aligned = true;
                let mut remaining = len;
                while remaining > 0 {
                    match md.as_mut().and_then(|md| md.front_mut()) {
                        Some(MDKind::Match(matches)) => {
                            let take = (*matches as usize).min(remaining);
                            *matches -= take as u16;
                            read_pos += take as u32;
                            ref_offset += take as u32;
                            remaining -= take;
                            if *matches == 0 {
                                md.as_mut().unwrap().pop_front();
                            }
                        }
                        Some(MDKind::Substitution(base)) => {
                            let base = *base;
                            if let Some(entry) = md.as_mut().unwrap().pop_front() {
                                debug_assert!(matches!(entry, MDKind::Substitution(_)));
                            }
                            push_event(
                                &mut evidence.events,
                                sequence,
                                qualities,
                                read_pos,
                                ref_offset,
                                b'X',
                                base,
                            );
                            read_pos += 1;
                            ref_offset += 1;
                            remaining -= 1;
                        }
                        Some(MDKind::Deletion(_)) => {
                            // A deletion entry belongs to the following `D` op.
                            break;
                        }
                        None => {
                            if op.kind() == Kind::SequenceMismatch {
                                push_event(
                                    &mut evidence.events,
                                    sequence,
                                    qualities,
                                    read_pos,
                                    ref_offset,
                                    b'X',
                                    b'N',
                                );
                            }
                            read_pos += 1;
                            ref_offset += 1;
                            remaining -= 1;
                        }
                    }
                }
            }
            Kind::Insertion => {
                aligned = true;
                for _ in 0..len {
                    push_event(
                        &mut evidence.events,
                        sequence,
                        qualities,
                        read_pos,
                        ref_offset,
                        b'I',
                        b'.',
                    );
                    read_pos += 1;
                }
            }
            Kind::Deletion => {
                aligned = true;
                let bases = match md.as_mut() {
                    Some(md) => match md.front() {
                        Some(MDKind::Deletion(bases)) => {
                            let bases = bases.clone();
                            md.pop_front();
                            bases
                        }
                        _ => vec![b'N'; len],
                    },
                    None => vec![b'N'; len],
                };
                for index in 0..len {
                    let base = bases.get(index).copied().unwrap_or(b'N');
                    evidence.events.push(RnaEvent {
                        read_pos,
                        ref_offset,
                        op: b'D',
                        ref_base: base.to_ascii_uppercase(),
                        read_base: b'.',
                        qual: 0,
                    });
                    ref_offset += 1;
                }
            }
            Kind::Skip => {
                aligned = true;
                ref_offset += len as u32;
            }
            _ => {
                aligned = true;
            }
        }
    }
    Ok(evidence)
}

fn push_event(
    events: &mut Vec<RnaEvent>,
    sequence: &[u8],
    qualities: &[u8],
    read_pos: u32,
    ref_offset: u32,
    op: u8,
    ref_base: u8,
) {
    let index = read_pos as usize;
    events.push(RnaEvent {
        read_pos,
        ref_offset,
        op,
        ref_base: ref_base.to_ascii_uppercase(),
        read_base: sequence
            .get(index)
            .copied()
            .unwrap_or(b'N')
            .to_ascii_uppercase(),
        qual: qualities.get(index).copied().unwrap_or(0),
    });
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

/// Sort rows into the `sc_rna_alignment_v2` layout:
/// `(chromosome, genomic block, cell tile, start, cell rank, UMI)`.
pub fn sort_rows_v2(
    rows: impl Iterator<Item = ScRnaRow>,
    temp_dir: Option<PathBuf>,
    chunk_size: usize,
    block_size_bp: u32,
    cell_tile_size: u32,
) -> anyhow::Result<impl ExactSizeIterator<Item = anyhow::Result<ScRnaRow>>> {
    if block_size_bp == 0 || cell_tile_size == 0 {
        anyhow::bail!("scRNA v2 block and cell tile sizes must be greater than zero");
    }
    let block = u64::from(block_size_bp);
    let tile = cell_tile_size;
    super::external_sort(rows, temp_dir, chunk_size, move |a, b| {
        a.chrom_id
            .cmp(&b.chrom_id)
            .then_with(|| (a.start as u64 / block).cmp(&(b.start as u64 / block)))
            .then_with(|| (a.cell_id / tile).cmp(&(b.cell_id / tile)))
            .then_with(|| a.start.cmp(&b.start))
            .then_with(|| (a.cell_id % tile).cmp(&(b.cell_id % tile)))
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

#[cfg(test)]
mod tests {
    use super::*;
    use noodles_core::Position;
    use noodles_sam::alignment::record::cigar::op::Op;
    use noodles_sam::alignment::record::data::field::Tag;
    use noodles_sam::alignment::record_buf::data::field::Value;

    fn record(seq: &str, cigar: Vec<(Kind, usize)>, md: Option<&str>) -> RecordBuf {
        let mut rec = RecordBuf::default();
        *rec.flags_mut() = noodles_sam::alignment::record::Flags::empty();
        *rec.reference_sequence_id_mut() = Some(0);
        *rec.alignment_start_mut() = Position::new(101);
        *rec.sequence_mut() = seq.as_bytes().to_vec().into();
        *rec.quality_scores_mut() = (0..seq.len())
            .map(|i| (i + 30) as u8)
            .collect::<Vec<_>>()
            .into();
        *rec.cigar_mut() = cigar
            .into_iter()
            .map(|(kind, len)| Op::new(kind, len))
            .collect();
        if let Some(md) = md {
            rec.data_mut()
                .insert(Tag::MISMATCHED_POSITIONS, Value::String(md.into()));
        }
        rec
    }

    #[test]
    fn substitution_uses_md_reference_base() {
        let rec = record("ACGT", vec![(Kind::Match, 4)], Some("1T2"));
        let evidence = rna_events(&rec).unwrap();
        assert!(!evidence.missing_md);
        assert_eq!(evidence.events.len(), 1);
        let event = evidence.events[0];
        assert_eq!(event.read_pos, 1);
        assert_eq!(event.ref_offset, 1);
        assert_eq!(event.op, b'X' as u8);
        assert_eq!(event.ref_base, b'T');
        assert_eq!(event.read_base, b'C');
        assert_eq!(event.qual, 31);
    }

    #[test]
    fn soft_clip_offsets_read_position() {
        let rec = record(
            "GGACGT",
            vec![(Kind::SoftClip, 2), (Kind::Match, 4)],
            Some("1T2"),
        );
        let evidence = rna_events(&rec).unwrap();
        assert_eq!(evidence.events.len(), 1);
        assert_eq!(evidence.events[0].read_pos, 3);
        assert_eq!(evidence.events[0].ref_offset, 1);
    }

    #[test]
    fn insertion_emits_every_inserted_base() {
        let rec = record(
            "ACGTAC",
            vec![(Kind::Match, 2), (Kind::Insertion, 2), (Kind::Match, 2)],
            Some("4"),
        );
        let evidence = rna_events(&rec).unwrap();
        assert_eq!(evidence.events.len(), 2);
        assert_eq!(
            (
                evidence.events[0].read_pos,
                evidence.events[0].ref_offset,
                evidence.events[0].op
            ),
            (2, 2, b'I' as u8)
        );
        assert_eq!(
            (evidence.events[0].read_base, evidence.events[0].qual),
            (b'G', 32)
        );
        assert_eq!(
            (evidence.events[1].read_pos, evidence.events[1].read_base),
            (3, b'T')
        );
    }

    #[test]
    fn deletion_uses_md_bases() {
        let rec = record(
            "ACGT",
            vec![(Kind::Match, 2), (Kind::Deletion, 2), (Kind::Match, 2)],
            Some("2^AG2"),
        );
        let evidence = rna_events(&rec).unwrap();
        assert_eq!(evidence.events.len(), 2);
        assert_eq!(
            (
                evidence.events[0].ref_offset,
                evidence.events[0].ref_base,
                evidence.events[0].op
            ),
            (2, b'A', b'D' as u8)
        );
        assert_eq!(
            (evidence.events[1].ref_offset, evidence.events[1].ref_base),
            (3, b'G')
        );
        assert_eq!(evidence.events[0].read_base, b'.');
        assert_eq!(evidence.events[0].qual, 0);
    }

    #[test]
    fn explicit_mismatch_operation_is_reported() {
        let rec = record(
            "ACGT",
            vec![
                (Kind::Match, 2),
                (Kind::SequenceMismatch, 1),
                (Kind::Match, 1),
            ],
            Some("2A1"),
        );
        let evidence = rna_events(&rec).unwrap();
        assert_eq!(evidence.events.len(), 1);
        assert_eq!(
            (evidence.events[0].read_pos, evidence.events[0].ref_offset),
            (2, 2)
        );
        assert_eq!(evidence.events[0].op, b'X' as u8);
    }

    #[test]
    fn missing_md_still_reports_indels() {
        let rec = record(
            "ACGTAC",
            vec![(Kind::Match, 2), (Kind::Insertion, 2), (Kind::Match, 2)],
            None,
        );
        let evidence = rna_events(&rec).unwrap();
        assert!(evidence.missing_md);
        assert_eq!(evidence.events.len(), 2);
        assert!(evidence.events.iter().all(|event| event.op == b'I' as u8));
    }
}
