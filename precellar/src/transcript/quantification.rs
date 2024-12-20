use anyhow::Result;
use bed_utils::{
    bed::{BEDLike, ParseError, Strand},
    extsort::ExternalSorterBuilder,
};
use itertools::Itertools;
use noodles::sam::{
    alignment::{record::Flags, Record, RecordBuf},
    Header,
};
use rayon::prelude::ParallelSliceMut;
use serde::{Deserialize, Serialize};
use std::path::PathBuf;

use crate::align::{MultiMap, MultiMapR};

use super::AlignmentAnnotator;

#[derive(Debug)]
pub struct Quantifier {
    annotator: AlignmentAnnotator,
}

impl Quantifier {
    pub fn new(annotator: AlignmentAnnotator) -> Self {
        Self { annotator }
    }

    pub fn quantify<'a, I>(
        &'a self,
        header: &'a Header,
        records: I,
    )
    where
        I: Iterator<Item = Vec<(MultiMapR, Option<MultiMapR>)>> + 'a,
    {
        records.for_each(|r| self.quantify_chunk(header, r));
    }

    fn quantify_chunk(
        &self,
        header: &Header,
        records: Vec<(MultiMapR, Option<MultiMapR>)>
    )
    where
    {
        records.into_iter().for_each(|(r1, r2)| {
            let read1: Vec<_> = r1.iter().cloned().collect();
            let alignments = if let Some(r2) = r2 {
                let read2: Vec<_> = r2.iter().cloned().collect();
                self.annotator.annotate_alignments_pe(header, read1, read2)
            } else {
                self.annotator.annotate_alignments_se(header, read1)
            };
            println!("{:?}", alignments);
        })
    }
 
}