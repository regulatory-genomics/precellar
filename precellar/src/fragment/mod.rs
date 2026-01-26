mod dedups;

use anyhow::Result;
use bed_utils::{
    bed::{BEDLike, ParseError, Strand},
    extsort::ExternalSorterBuilder,
};
use bitcode::{Decode, Encode};
use dedups::AlignmentInfo;
use itertools::Itertools;
use noodles::sam::{
    alignment::{record::Flags, Record},
    Header,
};
use rayon::prelude::ParallelSliceMut;
use std::{collections::BTreeMap, path::PathBuf};

use crate::{
    align::{MultiMap, MultiMapR, SNVs},
    fragment::dedups::{AlignmentMini, RemoveDuplicates},
};

pub type CellBarcode = String;

#[derive(Encode, Decode, Debug)]
pub struct ExtendedFields {
    pub first_segment: ReadInfo,
    pub second_segment: Option<ReadInfo>,
}

impl ExtendedFields {
    pub fn new(ali: &AlignmentInfo) -> Option<Self> {
        let first_segment = ReadInfo::new(&ali.read1)?;
        let second_segment = if let Some(ali2) = &ali.read2 {
            Some(ReadInfo::new(ali2)?)
        } else {
            None
        };
        Some(Self {
            first_segment,
            second_segment,
        })
    }

    pub fn add(&mut self, ali: &AlignmentInfo) {
        self.first_segment.update(&ali.read1);

        if let Some(read2) = ali.read2.as_ref() {
            self.second_segment.as_mut().unwrap().update(read2);
        }
    }
}

#[derive(Encode, Decode, Debug)]
pub struct ReadInfo {
    pub read_start: u32,
    pub read_end: u32,
    pub snps: Option<Vec<SNVs>>,  // We use Option because most reads will not have SNPs
}

impl ReadInfo {
    pub(crate) fn new(ali: &AlignmentMini) -> Option<Self> {
        let snv = ali.snps.as_ref()?;
        let snps = if snv.is_empty() {
            None
        } else {
            Some(vec![snv.clone()])
        };
        Some(Self {
            read_start: ali.alignment_start - 1,
            read_end: ali.alignment_end,
            snps,
        })
    }

    fn len(&self) -> u32 {
        self.read_end - self.read_start
    }

    fn update(&mut self, ali: &AlignmentMini) {
        if ali.alignment_start - 1 < self.read_start {
            self.read_start = ali.alignment_start - 1;
        }
        if ali.alignment_end > self.read_end {
            self.read_end = ali.alignment_end;
        }

        if let Some(snp) = &ali.snps {
            if !snp.is_empty() {
                self.add_snp(snp.clone());
            }
        }
    }

    fn add_snp(&mut self, snp: SNVs) {
        if let Some(snps) = &mut self.snps {
            snps.push(snp);
        } else {
            self.snps = Some(vec![snp]);
        }
    }

    fn count_snp(&self) -> BTreeMap<String, u32> {
        let mut count = BTreeMap::new();
        if let Some(snps) = self.snps.as_ref() {
            snps.iter().for_each(|snp| {
                let key = snp.to_string(self.read_start, self.read_start + self.len()).unwrap();
                count.entry(key)
                    .and_modify(|e| *e += 1u32)
                    .or_insert(1);
            });
        }
        count
    }

    fn stringify(&self, n: u32) -> String {
        let mut count = self.count_snp();

        let n_snps = self.snps.as_ref().map(|x| x.len()).unwrap_or(0);
        let n_match = n - n_snps as u32;
        if n_match > 0 {
            count.insert(self.len().to_string(), n_match);
        }

        let snps = count.into_iter()
            .map(|(snv, c)| {
                if c == n {
                    snv
                } else {
                    format!("{}:{}", snv, c)
                }
            })
            .join(";");
        format!("{}\t{}", self.read_start, snps)
    }
}

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
    pub extended: Option<ExtendedFields>,
}

impl Fragment {
    fn from_alignment(value: AlignmentInfo) -> Fragment {
        let extended = ExtendedFields::new(&value);
        if value.read2.is_none() {
            Fragment {
                chrom: value.reference_sequence.clone(),
                start: value.read1.alignment_start as u64 - 1,
                end: value.read1.alignment_end as u64,
                barcode: Some(value.barcode.clone()),
                count: 1,
                strand: Some(if value.read1.is_reverse_complemented() {
                    Strand::Reverse
                } else {
                    Strand::Forward
                }),
                extended,
            }
        } else {
            let rec1_5p = value.read1.alignment_5p();
            let rec2_5p = value.read2.as_ref().unwrap().alignment_5p();
            let (start, end) = if rec1_5p < rec2_5p {
                (rec1_5p, rec2_5p)
            } else {
                (rec2_5p, rec1_5p)
            };

            Fragment {
                chrom: value.reference_sequence.clone(),
                start: start as u64 - 1,
                end: end as u64,
                barcode: Some(value.barcode.clone()),
                count: 1,
                strand: None,
                extended,
            }
        }
    }
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
        if let Some(ext) = &self.extended {
            write!(f, "\t{}", ext.first_segment.stringify(self.count))?;
            if let Some(second_segment) = ext.second_segment.as_ref() {
                write!(f, "\t{}", second_segment.stringify(self.count))?;
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
            extended: None, // FIXME: handle SNPs
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
            chunk_size: 10000000,
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
        .sort_by_async(reads, |a, b| a.barcode.cmp(&b.barcode))
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