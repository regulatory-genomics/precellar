mod quantification;
mod annotate;
mod transcriptome;
pub(crate) mod de_dups;
#[cfg(test)]
mod validation_test;

pub use quantification::Quantifier;
pub use transcriptome::{Transcript, Gene, SpliceSegments, Exon, Exons,};
pub use annotate::{AlignmentAnnotator, AnnotatedAlignment};

#[cfg(test)]
pub use validation_test::{IntronValidationTest, ValidatedIntron};