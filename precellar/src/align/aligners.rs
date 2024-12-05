use super::fastq::AnnotatedFastq;
use anyhow::{bail, ensure};
pub use bwa_mem2::BurrowsWheelerAligner;
use noodles::sam::alignment::Record;
pub use star_aligner::StarAligner;

use noodles::sam::alignment::record_buf::{data::field::value::Value, RecordBuf};
use noodles::sam::alignment::record::data::field::tag::Tag;
use noodles::sam;
use rayon::iter::ParallelIterator;
use rayon::slice::ParallelSlice;
use bwa_mem2::{AlignerOpts, FMIndex, PairedEndStats};
use star_aligner::StarOpts;

pub type MultiMapR = MultiMap<RecordBuf>;

/// `MultiMap` represents a set of alignments for a single read.
#[derive(Debug, Clone)]
pub struct MultiMap<R> {
    pub primary: R,
    pub others: Option<Vec<R>>,
}

impl<R> MultiMap<R> {
    pub fn new(primary: R, others: Option<Vec<R>>) -> Self {
        Self { primary, others }
    }

    pub fn into_primary(self) -> R {
        self.primary
    }

    pub fn iter(&self) -> impl Iterator<Item = &R> {
        std::iter::once(&self.primary).chain(self.others.iter().flatten())
    }
}

impl<R> From<R> for MultiMap<R> {
    fn from(record: R) -> Self {
        Self {
            primary: record,
            others: None,
        }
    }
}

impl<R: Record> TryFrom<Vec<R>> for MultiMap<R> {
    type Error = anyhow::Error;

    fn try_from(mut vec: Vec<R>) -> Result<Self, Self::Error> {
        let n = vec.len();
        if n == 0 {
            Err(anyhow::anyhow!("No alignments"))
        } else if n == 1 {
            Ok(MultiMap::from(vec.into_iter().next().unwrap()))
        } else {
            let mut primary = None;
            vec.iter().enumerate().try_for_each(|(i, rec)| {
                if !rec.flags()?.is_secondary() {
                    if primary.is_some() {
                        bail!("Multiple primary alignments");
                    } else {
                        primary = Some(i);
                    }
                }
                Ok(())
            })?;
            ensure!(primary.is_some(), "No primary alignment");

            Ok(MultiMap::new(
                vec.swap_remove(primary.unwrap()),
                Some(vec),
            ))
        }
    }
}

pub trait Aligner {
    fn header(&self) -> sam::Header;

    fn align_reads(
        &mut self,
        num_threads: u16,
        records: Vec<AnnotatedFastq>,
    ) -> Vec<MultiMapR>;

    fn align_read_pairs(
        &mut self,
        num_threads: u16,
        records: Vec<AnnotatedFastq>,
    ) -> Vec<(MultiMapR, MultiMapR)>;
}

pub struct DummyAligner;

impl Aligner for DummyAligner {
    fn header(&self) -> sam::Header {
        sam::Header::default()
    }

    fn align_reads(&mut self, _: u16, _: Vec<AnnotatedFastq>) -> Vec<MultiMapR> {
        Vec::new()
    }

    fn align_read_pairs(
        &mut self,
        _: u16,
        _: Vec<AnnotatedFastq>,
    ) -> Vec<(MultiMapR, MultiMapR)> {
        Vec::new()
    }
}

impl Aligner for BurrowsWheelerAligner {
    fn header(&self) -> sam::Header {
        self.get_sam_header()
    }

    fn align_reads(
        &mut self,
        num_threads: u16,
        records: Vec<AnnotatedFastq>,
    ) -> Vec<MultiMapR> {
        let (info, mut reads): (Vec<_>, Vec<_>) = records
            .into_iter()
            .map(|rec| ((rec.barcode.unwrap(), rec.umi), rec.read1.unwrap()))
            .unzip();

        self.align_reads(num_threads, reads.as_mut_slice())
            .enumerate()
            .map(|(i, mut alignment)| {
                let (bc, umi) = info.get(i).unwrap();
                add_cell_barcode(
                    &mut alignment,
                    bc.raw.sequence(),
                    bc.raw.quality_scores(),
                    bc.corrected.as_deref(),
                );
                if let Some(umi) = umi {
                    add_umi(&mut alignment, umi.sequence(), umi.quality_scores());
                }
                alignment.into()
            })
            .collect()
    }

    fn align_read_pairs(
        &mut self,
        num_threads: u16,
        records: Vec<AnnotatedFastq>,
    ) -> Vec<(MultiMapR, MultiMapR)> {
        let (info, mut reads): (Vec<_>, Vec<_>) = records
            .into_iter()
            .map(|rec| {
                (
                    (rec.barcode.unwrap(), rec.umi),
                    (rec.read1.unwrap(), rec.read2.unwrap()),
                )
            })
            .unzip();
        self.align_read_pairs(num_threads, &mut reads)
            .enumerate()
            .map(|(i, (mut ali1, mut ali2))| {
                let (bc, umi) = info.get(i).unwrap();
                add_cell_barcode(
                    &mut ali1,
                    bc.raw.sequence(),
                    bc.raw.quality_scores(),
                    bc.corrected.as_deref(),
                );
                add_cell_barcode(
                    &mut ali2,
                    bc.raw.sequence(),
                    bc.raw.quality_scores(),
                    bc.corrected.as_deref(),
                );
                if let Some(umi) = umi {
                    add_umi(&mut ali1, umi.sequence(), umi.quality_scores());
                    add_umi(&mut ali2, umi.sequence(), umi.quality_scores());
                }
                (ali1.into(), ali2.into())
            })
            .collect()
    }
}

impl Aligner for StarAligner {
    fn header(&self) -> sam::Header {
        self.get_header().clone()
    }

    fn align_reads(
        &mut self,
        num_threads: u16,
        records: Vec<AnnotatedFastq>,
    ) -> Vec<MultiMapR> {
        let chunk_size = get_chunk_size(records.len(), num_threads as usize);

        records.par_chunks(chunk_size).flat_map_iter(|chunk| {
            let mut aligner = self.clone();
            chunk.iter().map(move |rec| {
                let bc = rec.barcode.as_ref().unwrap();
                let mut ali = aligner.align_read(&rec.read1.as_ref().unwrap()).unwrap();
                ali.iter_mut().for_each(|alignment| {
                    add_cell_barcode(
                        alignment,
                        bc.raw.sequence(),
                        bc.raw.quality_scores(),
                        bc.corrected.as_deref(),
                    );
                    if let Some(umi) = &rec.umi {
                        add_umi(alignment, umi.sequence(), umi.quality_scores());
                    };
                });
                ali.try_into().unwrap()
            })
        }).collect()
    }

    fn align_read_pairs(
        &mut self,
        num_threads: u16,
        records: Vec<AnnotatedFastq>,
    ) -> Vec<(MultiMapR, MultiMapR)> {
        let chunk_size = get_chunk_size(records.len(), num_threads as usize);

        records.par_chunks(chunk_size).flat_map_iter(|chunk| {
            let mut aligner = self.clone();
            chunk.iter().map(move |rec| {
                let bc = rec.barcode.as_ref().unwrap();
                let (mut ali1, mut ali2) = aligner.align_read_pair(
                    &rec.read1.as_ref().unwrap(),
                    &rec.read2.as_ref().unwrap()
                ).unwrap();
                ali1.iter_mut().chain(ali2.iter_mut()).for_each(|alignment| {
                    add_cell_barcode(
                        alignment,
                        bc.raw.sequence(),
                        bc.raw.quality_scores(),
                        bc.corrected.as_deref(),
                    );
                    if let Some(umi) = &rec.umi {
                        add_umi(alignment, umi.sequence(), umi.quality_scores());
                    };
                });
                (ali1.try_into().unwrap(), ali2.try_into().unwrap())
            }).collect::<Vec<_>>()
        }).collect()
    }
}

fn get_chunk_size(total_length: usize, num_threads: usize) -> usize {
    let chunk_size = total_length / num_threads;
    if chunk_size == 0 {
        1
    } else {
        chunk_size
    }
}

fn add_cell_barcode(
    record_buf: &mut RecordBuf,
    ori_barcode: &[u8],
    ori_qual: &[u8],
    correct_barcode: Option<&[u8]>,
) {
    let data = record_buf.data_mut();
    data.insert(
        Tag::CELL_BARCODE_SEQUENCE,
        Value::String(ori_barcode.into()),
    );
    data.insert(
        Tag::CELL_BARCODE_QUALITY_SCORES,
        Value::String(ori_qual.into()),
    );
    if let Some(barcode) = correct_barcode {
        data.insert(Tag::CELL_BARCODE_ID, Value::String(barcode.into()));
    }
}

fn add_umi(
    record_buf: &mut RecordBuf,
    umi: &[u8],
    qual: &[u8],
) {
    let data = record_buf.data_mut();
    data.insert(
        Tag::UMI_SEQUENCE,
        Value::String(umi.into()),
    );
    data.insert(
        Tag::UMI_QUALITY_SCORES,
        Value::String(qual.into()),
    );
}

pub trait AlignerBuilder: Aligner {
    fn from_path<P: AsRef<std::path::Path>>(path: P) -> Self;
}

impl AlignerBuilder for BurrowsWheelerAligner {
    fn from_path<P: AsRef<std::path::Path>>(path: P) -> Self {
        BurrowsWheelerAligner::new(
            FMIndex::read(path).unwrap(),
            AlignerOpts::default(),
            PairedEndStats::default(),
        )
    }
}

impl AlignerBuilder for StarAligner {
    fn from_path<P: AsRef<std::path::Path>>(path: P) -> Self {
        let opts = StarOpts::new(path);
        StarAligner::new(opts).unwrap()
    }
}