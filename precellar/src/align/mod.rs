mod aligners;
mod fastq;

pub use aligners::{Aligner, BurrowsWheelerAligner, MultiMap, MultiMapR, StarAligner, read_transcriptome_star};
pub use fastq::{extend_fastq_record, Barcode, FastqProcessor, NameCollatedRecords};
