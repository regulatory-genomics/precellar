mod aligners;
mod fastq;

pub use aligners::{Aligner, BurrowsWheelerAligner, MultiMap, MultiMapR, StarAligner};
pub use fastq::{extend_fastq_record, Barcode, FastqProcessor, NameCollatedRecords};
