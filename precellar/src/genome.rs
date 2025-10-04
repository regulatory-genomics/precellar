use std::ops::Range;

use bed_utils::bed::{GenomicRange, Strand};
use indexmap::IndexMap;
use noodles::sam;
use polars::{frame::DataFrame, prelude::Column, series::Series};

use crate::fragment::Fragment;

/// 0-based index that maps genomic loci to integers.
#[derive(Debug, Clone)]
pub struct GenomeBaseIndex {
    pub(crate) chrom_sizes: IndexMap<String, u64>,
    pub(crate) base_accum_len: Vec<u64>,
    pub(crate) binned_accum_len: Vec<u64>,
}

impl GenomeBaseIndex {
    pub fn new(header: &sam::Header) -> Self {
        let chrom_sizes: IndexMap<String, u64> = header
            .reference_sequences()
            .iter()
            .map(|(name, seq)| {
                let len = seq.length().get() as u64;
                (name.to_string(), len)
            })
            .collect();

        let mut acc = 0;
        let base_accum_len = chrom_sizes
            .iter()
            .map(|(_, length)| {
                acc += length;
                acc
            })
            .collect::<Vec<_>>();
        Self {
            chrom_sizes,
            binned_accum_len: base_accum_len.clone(),
            base_accum_len,
        }
    }

    pub fn get_chrom_sizes(&self) -> DataFrame {
        DataFrame::new(vec![
            Column::new(
                "reference_seq_name".into(),
                self.chrom_sizes
                    .iter()
                    .map(|x| x.0.clone())
                    .collect::<Series>(),
            ),
            Column::new(
                "reference_seq_length".into(),
                self.chrom_sizes.iter().map(|x| x.1).collect::<Series>(),
            ),
        ])
        .unwrap()
    }

    /// Retreive the range of a chromosome.
    pub fn get_range(&self, chr: &str) -> Option<Range<usize>> {
        let i = self.chrom_sizes.get_index_of(chr)?;
        let end = self.binned_accum_len[i];
        let start = if i == 0 {
            0
        } else {
            self.binned_accum_len[i - 1]
        };
        Some(start as usize..end as usize)
    }

    pub fn to_index(&self) -> anndata::data::index::Index {
        self.chrom_sizes()
            .map(|(chrom, length)| {
                let i = anndata::data::index::Interval {
                    start: 0,
                    end: length as usize,
                    size: 1,
                    step: 1,
                };
                (chrom.to_owned(), i)
            })
            .collect()
    }

    /// Number of indices.
    pub fn len(&self) -> usize {
        self.binned_accum_len
            .last()
            .map(|x| *x as usize)
            .unwrap_or(0)
    }

    pub fn chrom_sizes(&self) -> impl Iterator<Item = (&String, u64)> + '_ {
        let mut prev = 0;
        self.chrom_sizes
            .keys()
            .zip(self.base_accum_len.iter())
            .map(move |(chrom, acc)| {
                let length = acc - prev;
                prev = *acc;
                (chrom, length)
            })
    }

    /// Check if the index contains the given chromosome.
    pub fn contain_chrom(&self, chrom: &str) -> bool {
        self.chrom_sizes.contains_key(chrom)
    }

    /// Given a genomic position, return the corresponding index.
    pub fn get_position_rev(&self, chrom: &str, pos: u64) -> usize {
        let i = self
            .chrom_sizes
            .get_index_of(chrom)
            .expect(format!("Chromosome {} not found", chrom).as_str());
        let size = if i == 0 {
            self.base_accum_len[i]
        } else {
            self.base_accum_len[i] - self.base_accum_len[i - 1]
        };
        if pos as u64 >= size {
            panic!("Position {} is out of range for chromosome {}", pos, chrom);
        }
        if i == 0 {
            pos as usize
        } else {
            self.binned_accum_len[i - 1] as usize + pos as usize
        }
    }

    /// O(log(N)). Given a index, find the corresponding chromosome.
    pub fn get_chrom(&self, pos: usize) -> &String {
        let i = pos as u64;
        let j = match self.binned_accum_len.binary_search(&i) {
            Ok(j) => j + 1,
            Err(j) => j,
        };
        self.chrom_sizes.get_index(j).unwrap().0
    }

    /// O(log(N)). Given a index, find the corresponding chromosome and position.
    pub fn get_position(&self, pos: usize) -> (&String, u64) {
        let i = pos as u64;
        match self.binned_accum_len.binary_search(&i) {
            Ok(j) => (self.chrom_sizes.get_index(j + 1).unwrap().0, 0),
            Err(j) => {
                let chr = self.chrom_sizes.get_index(j).unwrap().0;
                let prev = if j == 0 {
                    0
                } else {
                    self.binned_accum_len[j - 1]
                };
                let start = i - prev;
                (chr, start)
            }
        }
    }

    /// O(log(N)). Given a index, find the corresponding chromosome and position.
    pub fn get_region(&self, pos: usize) -> GenomicRange {
        let i = pos as u64;
        match self.binned_accum_len.binary_search(&i) {
            Ok(j) => {
                let chr = self.chrom_sizes.get_index(j + 1).unwrap().0;
                let acc = self.base_accum_len[j + 1];
                let size = acc - self.base_accum_len[j];
                let start = 0;
                let end = (start + 1).min(size);
                GenomicRange::new(chr, start, end)
            }
            Err(j) => {
                let chr = self.chrom_sizes.get_index(j).unwrap().0;
                let acc = self.base_accum_len[j];
                let size = if j == 0 {
                    acc
                } else {
                    acc - self.base_accum_len[j - 1]
                };
                let prev = if j == 0 {
                    0
                } else {
                    self.binned_accum_len[j - 1]
                };
                let start = i - prev;
                let end = (start + 1).min(size);
                GenomicRange::new(chr, start, end)
            }
        }
    }

    pub fn count_fragments<V>(
        &self,
        fragments: impl IntoIterator<Item = Fragment>,
    ) -> Vec<(usize, V)>
    where
        V: TryFrom<i64> + Ord,
        <V as TryFrom<i64>>::Error: std::fmt::Debug,
    {
        let mut values = Vec::new();
        fragments.into_iter().for_each(|f| {
            let chrom = &f.chrom;
            if self.contain_chrom(chrom) {
                let start = f.start as i64;
                let end = f.end as i64;
                let size = end - start;
                let pos;
                let shift: V;
                match f.strand {
                    Some(Strand::Reverse) => {
                        pos = self.get_position_rev(chrom, (end - 1) as u64);
                        shift = (-size).try_into().expect(
                            format!(
                                "cannot convert size {} to {}",
                                -size,
                                std::any::type_name::<V>()
                            )
                            .as_str(),
                        );
                    }
                    _ => {
                        pos = self.get_position_rev(chrom, start as u64);
                        shift = size.try_into().expect(
                            format!(
                                "cannot convert size {} to {}",
                                size,
                                std::any::type_name::<V>()
                            )
                            .as_str(),
                        );
                    }
                }
                values.push((pos, shift));
            }
        });
        values.sort();
        values
    }
}
