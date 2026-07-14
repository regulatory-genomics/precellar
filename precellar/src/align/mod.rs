//! Raw fastq processing and alignment module.

mod aligners;
mod fastq;
mod snv;

pub use aligners::{
    Aligner, BurrowsWheelerAligner, MiniBwaIndex, MiniBwaOptions, MiniBwaSR, Minimap2Aligner,
    Minimap2Opts, MultiMap, MultiMapR, StarAligner,
};
pub use fastq::{
    extend_fastq_record, AlignmentBatch, AlignmentInput, AlignmentResult, AlignmentRunner,
    AlignmentStream, AnnotatedFastq, Barcode, BarcodeCorrectionConfig, FastqExecution, FastqPlan,
    FastqReport, MultiAnnotatedFqReader, NameCollatedRecords, ReadMetadata, RunReport,
};
pub use snv::{SNVs, SNV};
