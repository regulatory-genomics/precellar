//! Conversion helpers for UCF sinks.
//!
//! These keep the format-agnostic alignment/fragment stream out of the PyO3
//! layer: one external sort per output, no fragment TSV round trip.

use std::cmp::Ordering;
use std::path::{Path, PathBuf};

use anyhow::{Context, Result};
use bed_utils::extsort::ExternalSorterBuilder;
use bitcode::{DecodeOwned, Encode};

pub mod sc_atac;
pub mod sc_rna;

/// Create the external-sort directory when the caller supplied one.
///
/// The sorter treats a missing directory as a hard failure and does not create
/// it, so validate once here and return a normal error instead of a panic.
pub fn prepare_temp_dir(dir: Option<&Path>) -> Result<()> {
    let Some(dir) = dir else {
        return Ok(());
    };
    std::fs::create_dir_all(dir)
        .with_context(|| format!("cannot create temporary directory {}", dir.display()))?;
    if !dir.is_dir() {
        anyhow::bail!("temporary path is not a directory: {}", dir.display());
    }
    Ok(())
}

/// Sort records on disk so a UCF writer can stream them in genome order.
///
/// The sorter builds its runs eagerly, so any dictionary filled while producing
/// `items` is complete once this returns. Chunk read failures surface as
/// per-item errors instead of panics.
pub(crate) fn external_sort<T, F>(
    items: impl Iterator<Item = T>,
    temp_dir: Option<PathBuf>,
    chunk_size: usize,
    cmp: F,
) -> Result<impl ExactSizeIterator<Item = Result<T>>>
where
    T: Encode + DecodeOwned + Send + 'static,
    F: Fn(&T, &T) -> Ordering + Send + Sync + Copy + 'static,
{
    if chunk_size == 0 {
        anyhow::bail!("external sort chunk size must be greater than zero");
    }
    let mut sorter = ExternalSorterBuilder::new()
        .with_chunk_size(chunk_size)
        .with_compression(2);
    if let Some(tmp) = temp_dir {
        prepare_temp_dir(Some(&tmp))?;
        sorter = sorter.with_tmp_dir(tmp);
    }
    let sorted = sorter
        .build()
        .context("cannot build external sorter")?
        .sort_by_async(items, cmp)
        .context("external sort failed")?;
    Ok(sorted.map(|item| item.map_err(anyhow::Error::from)))
}
