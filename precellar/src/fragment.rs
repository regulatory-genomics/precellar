mod deduplicate;

use anyhow::Result;
use bed_utils::{
    bed::{BEDLike, ParseError, Strand},
    extsort::ExternalSorterBuilder,
};
use deduplicate::{remove_duplicates, AlignmentInfo};
use either::Either;
use itertools::Itertools;
use noodles::sam::{
    alignment::{record::Flags, Record},
    Header,
};
use rayon::prelude::ParallelSliceMut;
use serde::{Deserialize, Serialize};
use std::path::PathBuf;

pub type CellBarcode = String;

/// Fragments from single-cell ATAC-seq experiment. Each fragment is represented
/// by a genomic coordinate, cell barcode and a integer value.
#[derive(Serialize, Deserialize, Debug)]
pub struct Fragment {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    pub barcode: Option<CellBarcode>,
    pub count: u32,
    pub strand: Option<Strand>,
}

impl BEDLike for Fragment {
    fn chrom(&self) -> &str {
        &self.chrom
    }
    fn set_chrom(&mut self, chrom: &str) -> &mut Self {
        self.chrom = chrom.to_string();
        self
    }
    fn start(&self) -> u64 {
        self.start
    }
    fn set_start(&mut self, start: u64) -> &mut Self {
        self.start = start;
        self
    }
    fn end(&self) -> u64 {
        self.end
    }
    fn set_end(&mut self, end: u64) -> &mut Self {
        self.end = end;
        self
    }
    fn name(&self) -> Option<&str> {
        self.barcode.as_deref()
    }
    fn score(&self) -> Option<bed_utils::bed::Score> {
        Some(self.count.try_into().unwrap())
    }
    fn strand(&self) -> Option<Strand> {
        self.strand
    }
}

impl core::fmt::Display for Fragment {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(
            f,
            "{}\t{}\t{}\t{}\t{}",
            self.chrom(),
            self.start(),
            self.end(),
            self.barcode.as_deref().unwrap_or("."),
            self.count,
        )?;
        if let Some(strand) = self.strand() {
            write!(f, "\t{}", strand)?;
        }
        Ok(())
    }
}

impl std::str::FromStr for Fragment {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut fields = s.split('\t');
        let chrom = fields
            .next()
            .ok_or(ParseError::MissingReferenceSequenceName)?
            .to_string();
        let start = fields
            .next()
            .ok_or(ParseError::MissingStartPosition)
            .and_then(|s| lexical::parse(s).map_err(ParseError::InvalidStartPosition))?;
        let end = fields
            .next()
            .ok_or(ParseError::MissingEndPosition)
            .and_then(|s| lexical::parse(s).map_err(ParseError::InvalidEndPosition))?;
        let barcode = fields
            .next()
            .ok_or(ParseError::MissingName)
            .map(|s| match s {
                "." => None,
                _ => Some(s.into()),
            })?;
        let count = fields.next().map_or(Ok(1), |s| {
            if s == "." {
                Ok(1)
            } else {
                lexical::parse(s).map_err(ParseError::InvalidStartPosition)
            }
        })?;
        let strand = fields.next().map_or(Ok(None), |s| {
            if s == "." {
                Ok(None)
            } else {
                s.parse().map(Some).map_err(ParseError::InvalidStrand)
            }
        })?;
        Ok(Fragment {
            chrom,
            start,
            end,
            barcode,
            count,
            strand,
        })
    }
}

#[derive(Debug)]
pub struct FragmentGenerator {
    shift_left: i64,  // `shift_left` - Insertion site correction for the left end.
    shift_right: i64, // `shift_right` - Insertion site correction for the right end.
    mapq: u8,
    chunk_size: usize,
    temp_dir: Option<PathBuf>,
}

impl Default for FragmentGenerator {
    fn default() -> Self {
        Self {
            shift_left: 4,
            shift_right: -5,
            mapq: 30,
            chunk_size: 30000000,
            temp_dir: None,
        }
    }
}

impl FragmentGenerator {
    pub fn set_temp_dir<P: Into<PathBuf>>(&mut self, temp_dir: P) {
        self.temp_dir = Some(temp_dir.into());
    }

    pub fn set_shift_left(&mut self, shift_left: i64) {
        self.shift_left = shift_left;
    }

    pub fn set_shift_right(&mut self, shift_right: i64) {
        self.shift_right = shift_right;
    }

    pub fn gen_unique_fragments<'a, I, R>(
        &'a self,
        header: &'a Header,
        records: I,
    ) -> UniqueFragments<
        impl Iterator<Item = AlignmentInfo> + 'a,
        impl FnMut(&AlignmentInfo) -> String + 'a,
    >
    where
        I: Iterator<Item = Either<Vec<R>, Vec<(R, R)>>> + 'a,
        R: Record + 'a,
    {
        let data = records.flat_map(|x| match x {
            Either::Left(chunk) => Box::new(chunk.into_iter().flat_map(|r| {
                if filter_read(&r, self.mapq) {
                    AlignmentInfo::from_read(&r, header).unwrap()
                } else {
                    None
                }
            })) as Box<dyn Iterator<Item = AlignmentInfo>>,
            Either::Right(chunk) => Box::new(chunk.into_iter().flat_map(|(r1, r2)| {
                if filter_read_pair((&r1, &r2), self.mapq) {
                    AlignmentInfo::from_read_pair((&r1, &r2), header).unwrap()
                } else {
                    None
                }
            })) as Box<dyn Iterator<Item = AlignmentInfo>>,
        });

        let sorted = sort_by_barcode(data, self.temp_dir.clone(), self.chunk_size);
        UniqueFragments {
            shift_left: self.shift_left,
            shift_right: self.shift_right,
            chunks: sorted.chunk_by(|x| x.barcode.clone()),
        }
    }
}

pub struct UniqueFragments<I: Iterator, F> {
    shift_left: i64,
    shift_right: i64,
    chunks: itertools::ChunkBy<String, I, F>,
}

impl<'a, I: Iterator<Item = AlignmentInfo>, F: FnMut(&AlignmentInfo) -> String> IntoIterator
    for &'a UniqueFragments<I, F>
{
    type Item = Vec<Fragment>;
    type IntoIter = UniqueFragmentsIter<'a, I, F>;

    fn into_iter(self) -> Self::IntoIter {
        UniqueFragmentsIter {
            shift_left: self.shift_left,
            shift_right: self.shift_right,
            iter: self.chunks.into_iter(),
        }
    }
}

pub struct UniqueFragmentsIter<'a, I: Iterator, F> {
    shift_left: i64,
    shift_right: i64,
    iter: itertools::structs::Groups<'a, String, I, F>,
}

impl<'a, I: Iterator<Item = AlignmentInfo>, F: FnMut(&AlignmentInfo) -> String> Iterator
    for UniqueFragmentsIter<'a, I, F>
{
    type Item = Vec<Fragment>;

    fn next(&mut self) -> Option<Self::Item> {
        let (_, chunk) = self.iter.next()?;
        let mut fragments: Vec<_> = remove_duplicates(chunk)
            .drain()
            .flat_map(|(_, mut frag)| {
                if frag.strand().is_none() {
                    // perform fragment length correction for paired-end reads
                    frag.set_start(frag.start().saturating_add_signed(self.shift_left));
                    frag.set_end(frag.end().saturating_add_signed(self.shift_right));
                }
                if frag.len() > 0 {
                    Some(frag)
                } else {
                    None
                }
            })
            .collect();
        fragments.par_sort_unstable_by(|a, b| {
            a.chrom()
                .cmp(b.chrom())
                .then_with(|| a.start().cmp(&b.start()))
                .then_with(|| a.end().cmp(&b.end()))
                .then_with(|| a.score().as_deref().cmp(&b.score().as_deref()))
        });
        Some(fragments)
    }
}

fn sort_by_barcode<I, P>(
    reads: I,
    temp_dir: Option<P>,
    chunk_size: usize,
) -> impl ExactSizeIterator<Item = AlignmentInfo>
where
    I: Iterator<Item = AlignmentInfo>,
    P: AsRef<std::path::Path>,
{
    let mut sorter = ExternalSorterBuilder::new()
        .with_chunk_size(chunk_size)
        .with_compression(2);
    if let Some(tmp) = temp_dir {
        sorter = sorter.with_tmp_dir(tmp);
    }
    sorter
        .build()
        .unwrap()
        .sort_by(reads, |a, b| a.barcode.cmp(&b.barcode))
        .unwrap()
        .map(|x| x.unwrap())
}

fn filter_read<R: Record>(record: &R, min_q: u8) -> bool {
    // flag (3852) meaning:
    //   - read unmapped
    //   - mate unmapped
    //   - not primary alignment
    //   - read fails platform/vendor quality checks
    //   - read is PCR or optical duplicate
    //   - supplementary alignment
    !record
        .flags()
        .unwrap()
        .intersects(Flags::from_bits_retain(3852))
        && record.mapping_quality().map_or(255, |x| x.unwrap().get()) >= min_q
}

fn filter_read_pair<R: Record>(pair: (&R, &R), min_q: u8) -> bool {
    filter_read(pair.0, min_q)
        && filter_read(pair.1, min_q)
        && pair.0.flags().unwrap().is_properly_segmented()
}
