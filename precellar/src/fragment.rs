mod de_dups;

use anyhow::Result;
use bed_utils::{
    bed::{BEDLike, ParseError, Strand},
    extsort::ExternalSorterBuilder,
};
use bincode::{Decode, Encode};
use de_dups::AlignmentInfo;
use itertools::Itertools;
use noodles::sam::{
    alignment::{record::Flags, Record},
    Header,
};
use rayon::prelude::ParallelSliceMut;
use std::path::PathBuf;

use crate::{
    align::{MultiMap, MultiMapR},
    fragment::de_dups::{MutationCount, RemoveDuplicates},
};

pub type CellBarcode = String;

/// Fragments from single-cell ATAC-seq experiment. Each fragment is represented
/// by a genomic coordinate, cell barcode and a integer value.
#[derive(Encode, Decode, Debug)]
pub struct Fragment {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    pub barcode: Option<CellBarcode>,
    pub count: u32,
    pub strand: Option<Strand>,
    pub snps: Option<MutationCount>,
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
        } else {
            write!(f, "\t.")?;
        }
        if let Some(snp) = &self.snps {
            if snp.len() == 0 {
                write!(f, "\t.")?;
            } else {
                write!(f, "\t{}", snp)?;
            }
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
            snps: None, // FIXME: handle SNPs
        })
    }
}

#[derive(Debug, Clone)]
pub struct IntoFragOpts {
    pub shift_left: i64,  // Insertion site correction for the left end.
    pub shift_right: i64, // Insertion site correction for the right end.
    pub min_mapq: u8,
    pub chunk_size: usize,
    pub temp_dir: Option<PathBuf>,
    pub compute_snv: bool, // Whether to compute SNVs for fragments.
}

impl Default for IntoFragOpts {
    fn default() -> Self {
        Self {
            shift_left: 0,
            shift_right: 0,
            min_mapq: 30,
            chunk_size: 30000000,
            temp_dir: None,
            compute_snv: false,
        }
    }
}

pub trait IntoFragments: Iterator {
    fn into_fragments<'a>(
        self,
        header: &'a Header,
        opts: IntoFragOpts,
    ) -> UniqueFragments<
        impl Iterator<Item = AlignmentInfo> + 'a,
        impl FnMut(&AlignmentInfo) -> String + 'a,
    >
    where
        Self: Iterator<Item = Vec<(Option<MultiMapR>, Option<MultiMapR>)>> + 'a + Sized,
    {
        let data = self.flat_map(move |chunk| {
            chunk.into_iter().flat_map(move |(r1, r2)| {
                if r1.is_some() && r2.is_some() {
                    let r1 = r1.unwrap();
                    let r2 = r2.unwrap();
                    if filter_read_pair((&r1.primary, &r2.primary), opts.min_mapq) {
                        AlignmentInfo::from_read_pair((&r1.primary, &r2.primary), header, opts.compute_snv).unwrap()
                    } else {
                        None
                    }
                } else {
                    let r = r1.or(r2).unwrap();
                    if filter_read(&r.primary, opts.min_mapq) {
                        AlignmentInfo::from_read(&r.primary, header, opts.compute_snv).unwrap()
                    } else {
                        None
                    }
                }
            })
        });

        let sorted = sort_by_barcode(data, opts.temp_dir.clone(), opts.chunk_size);
        UniqueFragments {
            shift_left: opts.shift_left,
            shift_right: opts.shift_right,
            chunks: sorted.chunk_by(|x| x.barcode.clone()),
        }
    }
}

impl<T, R> IntoFragments for T
where
    T: Iterator<Item = Vec<(Option<MultiMap<R>>, Option<MultiMap<R>>)>> + Sized,
    R: Record,
{
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
        let mut fragments: Vec<_> = chunk
            .into_iter()
            .remove_duplicates()
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

/*
pub fn txxx(cigar: Vec<u8>) -> Result<(u64, Vec<(u64, u64)>)>{
    let cigar = noodles::sam::record::Cigar::new(&cigar);
    let mut left_clip = 0;
    let mut segments = Vec::new();
    let mut seen_nonclips = false; // whether we've seen non-clip bases yet

    let mut cur_start = 0u64;
    let mut cur_end = 0u64;
    for c in cigar.iter() {
        let c = c?;
        match c.kind() {
            Kind::HardClip | Kind::SoftClip => {
                if !seen_nonclips {
                    left_clip += c.len() as u64;
                }
            }
            Kind::Skip => {
                seen_nonclips = true;
                let next_start = cur_end + c.len() as u64;
                segments.push((cur_start, cur_end));
                cur_start = next_start;
                cur_end = next_start;
            }
            Kind::Insertion => {
                seen_nonclips = true;
            }
            Kind::Match | Kind::Deletion | Kind::SequenceMatch | Kind::SequenceMismatch => {
                seen_nonclips = true;
                cur_end += c.len() as u64;
            }
            Kind::Pad => unreachable!(),
        }
    }
    if cur_start < cur_end {
        segments.push((cur_start, cur_end));
    }
    Ok((left_clip, segments))
}
*/
