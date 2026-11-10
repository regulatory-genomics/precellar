use std::collections::HashMap;
use std::path::PathBuf;

use bitcode::{Decode, Encode};

/// One compact scATAC row in UCF `genome_start` order.
#[derive(Encode, Decode, Debug, Clone, Copy, PartialEq, Eq)]
pub struct ScAtacRow {
    pub chrom_id: u32,
    pub start: u32,
    pub end: u32,
    pub cell_id: u32,
    pub multiplicity: u32,
    pub flags: u8,
}

/// Strand bits already used by `ucf-sc-atac` for `+` and `-`.
pub const STRAND_FLAG_PLUS: u8 = 0b0000_0001;
pub const STRAND_FLAG_MINUS: u8 = 0b0000_0010;

pub fn strand_flags(strand: Option<bed_utils::bed::Strand>) -> u8 {
    match strand {
        Some(bed_utils::bed::Strand::Forward) => STRAND_FLAG_PLUS,
        Some(bed_utils::bed::Strand::Reverse) => STRAND_FLAG_MINUS,
        None => 0,
    }
}

/// Assigns stable contiguous IDs to strings in first-seen order.
#[derive(Default)]
pub struct Dictionary {
    ids: HashMap<String, u32>,
    entries: Vec<String>,
}

impl Dictionary {
    pub fn id(&mut self, value: &str) -> u32 {
        if let Some(id) = self.ids.get(value) {
            return *id;
        }
        let id = u32::try_from(self.entries.len()).expect("dictionary exceeds u32 entries");
        self.ids.insert(value.to_owned(), id);
        self.entries.push(value.to_owned());
        id
    }

    pub fn entries(&self) -> &[String] {
        &self.entries
    }

    pub fn len(&self) -> usize {
        self.entries.len()
    }

    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }
}

/// Sort fragment rows for `ScAtacWriter`, which validates
/// `(chrom_id, block, start, cell_id, end, flags)` order.
pub fn sort_rows(
    rows: impl Iterator<Item = ScAtacRow>,
    temp_dir: Option<PathBuf>,
    chunk_size: usize,
) -> impl ExactSizeIterator<Item = ScAtacRow> {
    super::external_sort(rows, temp_dir, chunk_size, compare_rows)
}

/// A scATAC row plus its optional TSV SNV line, sorted together so the SNV
/// sidecar and the UCF payload share one order.
#[derive(Encode, Decode, Debug, Clone)]
pub struct ScAtacSnvRow {
    pub row: ScAtacRow,
    pub line: String,
}

pub fn sort_snv_rows(
    rows: impl Iterator<Item = ScAtacSnvRow>,
    temp_dir: Option<PathBuf>,
    chunk_size: usize,
) -> impl ExactSizeIterator<Item = ScAtacSnvRow> {
    super::external_sort(rows, temp_dir, chunk_size, |a, b| compare_rows(&a.row, &b.row))
}

fn compare_rows(a: &ScAtacRow, b: &ScAtacRow) -> std::cmp::Ordering {
    a.chrom_id
        .cmp(&b.chrom_id)
        .then_with(|| a.start.cmp(&b.start))
        .then_with(|| a.cell_id.cmp(&b.cell_id))
        .then_with(|| a.end.cmp(&b.end))
        .then_with(|| a.flags.cmp(&b.flags))
}
