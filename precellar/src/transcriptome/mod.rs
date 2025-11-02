//! Transcriptome module: handling transcripts, spliced alignments, and quantification.
//! 
//! Key structs:
//! * SplicedRecord: represents a spliced BAM/SAM record, split at refskips.
//! * TranscriptAlignment: represents the alignment of a SplicedRecord to a Transcript.
//! * AnnotatedAlignment: final annotated alignment, with gene information.

mod quantification;
mod align;
mod annotate;

pub use quantification::Quantifier;
pub use align::{TxAlignResult, TxAligner, TxAlignment};
pub use annotate::{GeneCounter, CorrectUMI};

use anyhow::{bail, ensure, Result};
use bed_utils::bed::Strand;

/// Represents a transcript. 
#[derive(Debug, Clone)]
pub struct Transcript {
    pub id: String,
    pub chrom: String,
    pub start: u64, // 0-based, half-open
    pub end: u64, // exclusive
    pub strand: Strand,
    pub gene_id: String,
    pub gene_name: String,
    exons: Exons,
}

impl TryFrom<star_aligner::transcript::Transcript> for Transcript {
    type Error = anyhow::Error;

    fn try_from(transcript: star_aligner::transcript::Transcript) -> Result<Self> {
        let start = transcript.start;
        let end = transcript.end;
        let strand = match transcript.strand {
            star_aligner::transcript::Strand::Forward => Strand::Forward,
            star_aligner::transcript::Strand::Reverse => Strand::Reverse,
            _ => bail!("Strand must be Forward or Reverse"),
        };
        let exons = Exons::new(transcript.exons.iter().map(|exon| {
            assert!(exon.start < exon.end, "Exon start must be less than exon end");
            assert!(exon.start >= start, "Exon start must be greater than transcript start");
            assert!(exon.end <= end, "Exon end must be less than transcript end");
            (exon.start, exon.end)
        }))?;
        Ok(Self {
            id: transcript.id,
            chrom: transcript.chrom,
            start,
            end,
            strand,
            gene_id: transcript.gene_id,
            gene_name: transcript.gene_name,
            exons,
        })
    }
}

impl Transcript {
    /// The sum of the lengths of the exons.
    pub fn exon_len(&self) -> u64 {
        self.exons().iter().map(|x| x.len()).sum()
    }

    pub fn exons(&self) -> &[Exon] {
        self.exons.as_ref()
    }

    /// Convert a coordinate in the genome to a coordinate in the exons/transcript.
    /// The coordinate starts at the beginning of the transcript.
    pub fn get_offset(&self, coord: u64) -> Option<i64> {
        let mut cum_len = 0;
        for exon in &self.exons.0 {
            if coord < exon.end {
                return Some((cum_len + coord) as i64 - exon.start as i64);
            }
            cum_len += exon.len();
        }
        None
    }
}

/// Represents a list of exons. Exons must be non-overlapping and in order.
#[derive(Debug, Clone)]
pub struct Exons(Vec<Exon>);

impl Exons {
    pub fn new(iter: impl IntoIterator<Item = (u64, u64)>) -> Result<Self> {
        let mut prev_end = None;
        let exon = iter
            .into_iter()
            .map(|(start, end)| {
                ensure!(
                    prev_end.is_none() || start >= prev_end.unwrap(),
                    "Exons must be non-overlapping and in order"
                );
                ensure!(
                    end > start,
                    "End coordinate must be greater than start coordinate"
                );
                prev_end = Some(end);
                Ok(Exon { start, end })
            })
            .collect::<Result<Vec<_>>>()?;
        Ok(Self(exon))
    }
}

impl AsRef<[Exon]> for Exons {
    fn as_ref(&self) -> &[Exon] {
        &self.0
    }
}

/// Exon coordinates are 0-based, half-open.
#[derive(Eq, PartialEq, Debug, Clone, Ord, PartialOrd)]
pub struct Exon {
    start: u64,
    end: u64,
}

impl Exon {
    pub fn len(&self) -> u64 {
        self.end - self.start
    }
}