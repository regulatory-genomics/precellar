use crate::seqspec::SeqSpec;

use noodles::{sam, fastq};
use noodles::sam::alignment::{
    Record, record_buf::RecordBuf, record::data::field::tag::Tag,
};
use noodles::sam::alignment::record_buf::data::field::value::Value;

pub trait Alinger {
    fn align_reads<'a, I>(&'a self, records: I) -> impl Iterator<Item = sam::Record> + '_
    where I: Iterator<Item = fastq::Record> + 'a;

    fn align_read_pairs<'a, I>(&'a self, records: I) -> impl Iterator<Item = (sam::Record, sam::Record)> + '_
    where I: Iterator<Item = (fastq::Record, fastq::Record)> + 'a;
}

pub struct FastqProcessor<A> {
    seqspec: SeqSpec,
    aligner: A,
    read1: Option<String>,
    read2: Option<String>,
}

impl<A: Alinger> FastqProcessor<A> {
    pub fn new(seqspec: SeqSpec, aligner: A) -> Self {
        Self { seqspec, aligner, read1: None, read2: None }
    }
}

pub fn add_cell_barcode<R: Record>(
    header: &sam::Header,
    record: &R,
    ori_barcode: &str,
    ori_qual: &[u8],
    correct_barcode: Option<&str>,
) -> std::io::Result<RecordBuf> {
    let mut record_buf = RecordBuf::try_from_alignment_record(header, record)?;
    let data = record_buf.data_mut();
    data.insert(Tag::CELL_BARCODE_SEQUENCE, Value::String(ori_barcode.into()));
    data.insert(Tag::CELL_BARCODE_QUALITY_SCORES, Value::String(ori_qual.into()));
    if let Some(barcode) = correct_barcode {
        data.insert(Tag::CELL_BARCODE_ID, Value::String(barcode.into()));
    }
    Ok(record_buf)
}

fn slice_fastq_record(mut record: fastq::Record, start: usize, end: usize) -> fastq::Record {
    *record.sequence_mut() = record.sequence()[start..end].to_vec();
    *record.quality_scores_mut() = record.quality_scores()[start..end].to_vec();
    record
}