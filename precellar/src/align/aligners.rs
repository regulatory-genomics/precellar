/// This module provides an abstraction for aligning sequencing reads using different alignment tools like BWA and STAR.
mod minimap2;
pub use minimap2::{Minimap2Aligner, Minimap2Opts};

use super::fastq::AnnotatedFastq;
use crate::barcode::{get_barcode, get_umi};

use anyhow::{bail, ensure, Result};
pub use bwa_mem2::BurrowsWheelerAligner;
use noodles::sam::alignment::Record;
pub use star_aligner::StarAligner;

use log;
use noodles::sam;
use noodles::sam::alignment::record::data::field::tag::Tag;
use noodles::sam::alignment::record_buf::{data::field::value::Value, RecordBuf};
use rayon::iter::ParallelIterator;
use rayon::slice::ParallelSlice;

pub type MultiMapR = MultiMap<RecordBuf>;

/// Represents a set of alignments (primary and optional secondary alignments) for a single sequencing read.
#[derive(Debug, Clone)]
pub struct MultiMap<R> {
    /// The primary alignment for the read.
    pub primary: R,
    /// Optional secondary alignments for the read.
    pub others: Option<Vec<R>>,
}

impl<R: Record> MultiMap<R> {
    /// Constructs a new `MultiMap`.
    ///
    /// # Arguments
    /// * `primary` - The primary alignment for the read.
    /// * `others` - Optional secondary alignments for the read.
    pub fn new(primary: R, others: Option<Vec<R>>) -> Self {
        Self { primary, others }
    }

    /// Return the number of records.
    pub fn len(&self) -> usize {
        self.others.as_ref().map_or(0, |x| x.len()) + 1
    }

    /// Return the cell barcode if it exists.
    pub fn barcode(&self) -> Result<Option<String>> {
        get_barcode(&self.primary)
    }

    /// Return the UMI if it exists.
    pub fn umi(&self) -> Result<Option<String>> {
        get_umi(&self.primary)
    }

    /// Consumes the `MultiMap` and returns the primary alignment.
    pub fn into_primary(self) -> R {
        self.primary
    }

    /// Whether the read is confidently mapped. A read is confidently mapped if it
    /// is mapped to a single location.
    pub fn is_confidently_mapped(&self) -> bool {
        self.others.is_none() && !self.primary.flags().unwrap().is_unmapped()
    }

    /// Returns an iterator over all alignments (primary and secondary).
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

            Ok(MultiMap::new(vec.swap_remove(primary.unwrap()), Some(vec)))
        }
    }
}

/// Trait defining the behavior of aligners like BWA and STAR.
pub trait Aligner {
    /// Retrieves the SAM header associated with the aligner.
    fn header(&self) -> sam::Header;

    /// Aligns a batch of sequencing reads.
    ///
    /// # Arguments
    /// * `num_threads` - Number of threads to use for alignment.
    /// * `records` - Vector of annotated FASTQ records to align.
    ///
    /// # Returns
    /// A vector of tuples where each tuple contains the alignments.
    fn align_reads(
        &mut self,
        num_threads: u16,
        records: Vec<AnnotatedFastq>,
    ) -> Vec<(Option<MultiMapR>, Option<MultiMapR>)>;
}

impl Aligner for BurrowsWheelerAligner {
    fn header(&self) -> sam::Header {
        self.get_sam_header()
    }

    fn align_reads(
        &mut self,
        num_threads: u16,
        records: Vec<AnnotatedFastq>,
    ) -> Vec<(Option<MultiMapR>, Option<MultiMapR>)> {
        if records[0].read2.is_some() {
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
                    (Some(ali1.into()), Some(ali2.into()))
                })
                .collect()
        } else {
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
                    (Some(alignment.into()), None)
                })
                .collect()
        }
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
    ) -> Vec<(Option<MultiMapR>, Option<MultiMapR>)> {
        let chunk_size = get_chunk_size(records.len(), num_threads as usize);

        records
            .par_chunks(chunk_size)
            .flat_map_iter(|chunk| {
                let mut aligner = self.clone();
                chunk.iter().map(move |rec| {
                    let bc = rec.barcode.as_ref().unwrap();
                    let read1 = rec.read1.as_ref();
                    let read2 = rec.read2.as_ref();

                    if read1.is_some() && read2.is_some() {
                        let (mut ali1, mut ali2) =
                            aligner.align_read_pair(&read1.unwrap(), &read2.unwrap()).unwrap();
                        ali1.iter_mut()
                            .chain(ali2.iter_mut())
                            .for_each(|alignment| {
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
                        (Some(ali1.try_into().unwrap()), Some(ali2.try_into().unwrap()))
                    } else if let Some(read) = read1.or(read2) {
                        let mut ali = aligner.align_read(read).unwrap();
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
                        if read1.is_some() {
                            (Some(ali.try_into().unwrap()), None)
                        } else {
                            (None, Some(ali.try_into().unwrap()))
                        }
                    } else {
                        log::warn!("Found record with no reads (read1 and read2 are both None). Barcode: {:?}", 
                                  String::from_utf8_lossy(bc.raw.sequence()));
                        (None, None)
                    }
                })
            })
            .collect()
    }
}

impl Aligner for Minimap2Aligner {
    fn header(&self) -> sam::Header {
        self.get_header().clone()
    }

    fn align_reads(
        &mut self,
        num_threads: u16,
        records: Vec<AnnotatedFastq>,
    ) -> Vec<(Option<MultiMapR>, Option<MultiMapR>)> {
        let chunk_size = get_chunk_size(records.len(), num_threads as usize);

        // Use Rayon for parallel processing with chunks
        records
            .par_chunks(chunk_size)
            .flat_map_iter(|chunk| {
                // Clone aligner for this thread (efficient: only clones Arc pointers to shared index)
                let mut thread_aligner = self.clone();

                chunk.iter().map(move |rec| {
                    let bc = rec.barcode.as_ref().unwrap();
                    let read1 = rec.read1.as_ref();
                    let read2 = rec.read2.as_ref();

                    if read1.is_some() && read2.is_some() {
                        let (mut ali1, mut ali2) =
                            thread_aligner.align_read_pair(&read1.unwrap(), &read2.unwrap()).unwrap();
                        ali1.iter_mut()
                            .chain(ali2.iter_mut())
                            .for_each(|alignment| {
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
                        (Some(ali1.try_into().unwrap()), Some(ali2.try_into().unwrap()))
                    } else if let Some(read) = read1.or(read2) {
                        let mut ali = thread_aligner.align_read(read).unwrap();
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
                        if read1.is_some() {
                            (Some(ali.try_into().unwrap()), None)
                        } else {
                            (None, Some(ali.try_into().unwrap()))
                        }
                    } else {
                        log::warn!("Found record with no reads (read1 and read2 are both None). Barcode: {:?}",
                                  String::from_utf8_lossy(bc.raw.sequence()));
                        (None, None)
                    }
                })
            })
            .collect()
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

// Additional helper functions for adding metadata like cell barcodes and UMIs to alignments.
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

fn add_umi(record_buf: &mut RecordBuf, umi: &[u8], qual: &[u8]) {
    let data = record_buf.data_mut();
    data.insert(Tag::UMI_SEQUENCE, Value::String(umi.into()));
    data.insert(Tag::UMI_QUALITY_SCORES, Value::String(qual.into()));
}
