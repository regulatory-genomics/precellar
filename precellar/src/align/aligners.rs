use super::AnnotatedFastq;
pub use bwa_mem2::BurrowsWheelerAligner;
pub use star_aligner::StarAligner;

use noodles::sam::alignment::record_buf::{data::field::value::Value, RecordBuf};
use noodles::sam::alignment::record::data::field::tag::Tag;
use noodles::sam;
use rayon::iter::ParallelIterator;
use rayon::slice::ParallelSlice;
use bwa_mem2::{AlignerOpts, FMIndex, PairedEndStats};
use star_aligner::StarOpts;

pub trait AsIterator {
    type Item;
    type AsIter<'a>: Iterator<Item = &'a Self::Item>
    where
        Self: 'a;

    fn as_iter(&self) -> Self::AsIter<'_>;
}

impl AsIterator for RecordBuf {
    type Item = RecordBuf;
    type AsIter<'a> = std::iter::Once<&'a RecordBuf>;

    fn as_iter(&self) -> Self::AsIter<'_> {
        std::iter::once(&self)
    }
}

impl AsIterator for Vec<RecordBuf> {
    type Item = RecordBuf;
    type AsIter<'a> = std::slice::Iter<'a, RecordBuf>;

    fn as_iter(&self) -> Self::AsIter<'_> {
        self.iter()
    }
}

pub trait Aligner {
    type AlignOutput: AsIterator<Item = RecordBuf>;

    fn from_path<P: AsRef<std::path::Path>>(path: P) -> Self;

    fn header(&self) -> sam::Header;

    fn align_reads(
        &mut self,
        num_threads: u16,
        records: Vec<AnnotatedFastq>,
    ) -> Vec<Self::AlignOutput>;

    fn align_read_pairs(
        &mut self,
        num_threads: u16,
        records: Vec<AnnotatedFastq>,
    ) -> Vec<(Self::AlignOutput, Self::AlignOutput)>;
}

pub struct DummyAligner;

impl Aligner for DummyAligner {
    type AlignOutput = RecordBuf;

    fn header(&self) -> sam::Header {
        sam::Header::default()
    }

    fn from_path<P: AsRef<std::path::Path>>(_path: P) -> Self {
        Self
    }

    fn align_reads(&mut self, _: u16, _: Vec<AnnotatedFastq>) -> Vec<Self::AlignOutput> {
        Vec::new()
    }

    fn align_read_pairs(
        &mut self,
        _: u16,
        _: Vec<AnnotatedFastq>,
    ) -> Vec<(Self::AlignOutput, Self::AlignOutput)> {
        Vec::new()
    }
}

impl Aligner for BurrowsWheelerAligner {
    type AlignOutput = RecordBuf;

    fn header(&self) -> sam::Header {
        self.get_sam_header()
    }

    fn from_path<P: AsRef<std::path::Path>>(path: P) -> Self {
        BurrowsWheelerAligner::new(
            FMIndex::read(path).unwrap(),
            AlignerOpts::default(),
            PairedEndStats::default(),
        )
    }

    fn align_reads(
        &mut self,
        num_threads: u16,
        records: Vec<AnnotatedFastq>,
    ) -> Vec<Self::AlignOutput> {
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
                alignment
            })
            .collect()
    }

    fn align_read_pairs(
        &mut self,
        num_threads: u16,
        records: Vec<AnnotatedFastq>,
    ) -> Vec<(Self::AlignOutput, Self::AlignOutput)> {
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
                (ali1, ali2)
            })
            .collect()
    }
}

impl Aligner for StarAligner {
    type AlignOutput = Vec<RecordBuf>;

    fn header(&self) -> sam::Header {
        self.get_header().clone()
    }

    fn from_path<P: AsRef<std::path::Path>>(path: P) -> Self {
        let opts = StarOpts::new(path);
        StarAligner::new(opts).unwrap()
    }

    fn align_reads(
        &mut self,
        num_threads: u16,
        records: Vec<AnnotatedFastq>,
    ) -> Vec<Self::AlignOutput> {
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
                ali
            })
        }).collect()
    }

    fn align_read_pairs(
        &mut self,
        num_threads: u16,
        records: Vec<AnnotatedFastq>,
    ) -> Vec<(Self::AlignOutput, Self::AlignOutput)> {
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
                (ali1, ali2)
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