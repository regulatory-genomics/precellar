mod fastq;
mod aligners;

pub use fastq::{extend_fastq_record, Barcode, FastqProcessor, NameCollatedRecords};
pub use aligners::{
    Aligner, AlignerBuilder, BurrowsWheelerAligner, DummyAligner, MultiMap, MultiMapR, StarAligner,
};