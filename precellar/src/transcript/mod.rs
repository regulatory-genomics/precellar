mod quantification;
mod annotate;
mod transcriptome;
pub(crate) mod de_dups;

pub use quantification::Quantifier;
pub use transcriptome::{Transcript, Gene, SpliceSegments, Exon, Exons};
pub use annotate::{AlignmentAnnotator, AnnotatedAlignment};