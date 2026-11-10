//! Conversion helpers for UCF sinks.
//!
//! These keep the format-agnostic alignment/fragment stream out of the PyO3
//! layer: one external sort per output, no fragment TSV round trip.

use std::cmp::Ordering;
use std::path::PathBuf;

use bed_utils::extsort::ExternalSorterBuilder;
use bitcode::{DecodeOwned, Encode};

pub mod sc_atac;
pub mod sc_rna;

/// Sort records on disk so a UCF writer can stream them in genome order.
///
/// Returns the merged iterator. The sorter builds its runs eagerly, so any
/// dictionary filled while producing `items` is complete once this returns.
pub(crate) fn external_sort<T, F>(
    items: impl Iterator<Item = T>,
    temp_dir: Option<PathBuf>,
    chunk_size: usize,
    cmp: F,
) -> impl ExactSizeIterator<Item = T>
where
    T: Encode + DecodeOwned + Send + 'static,
    F: Fn(&T, &T) -> Ordering + Send + Sync + Copy + 'static,
{
    let mut sorter = ExternalSorterBuilder::new()
        .with_chunk_size(chunk_size)
        .with_compression(2);
    if let Some(tmp) = temp_dir {
        sorter = sorter.with_tmp_dir(tmp);
    }
    sorter
        .build()
        .expect("build external sorter")
        .sort_by_async(items, cmp)
        .expect("external sort failed")
        .map(|item| item.expect("external sort record"))
}
