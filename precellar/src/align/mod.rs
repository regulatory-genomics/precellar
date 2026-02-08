//! Raw fastq processing and alignment module.

mod aligners;
mod fastq;
mod snv;

pub use aligners::{Aligner, BurrowsWheelerAligner, MultiMap, MultiMapR, StarAligner, Minimap2Aligner, Minimap2Opts};
pub use fastq::{extend_fastq_record, AnnotatedFastq, Barcode, FastqProcessor, AlignmentResult, NameCollatedRecords};
pub use snv::{SNV, SNVs};