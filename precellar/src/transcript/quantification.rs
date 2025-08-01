use anndata::{
    data::utils::{from_csr_data, to_csr_data},
    AnnData, AnnDataOp,
};
use anndata_hdf5::H5;
use anyhow::Result;
use bed_utils::extsort::ExternalSorterBuilder;
use bincode::{Decode, Encode};
use indexmap::IndexMap;
use itertools::Itertools;
use noodles::sam::Header;
use polars::df;
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use std::{collections::HashSet, path::PathBuf};

use crate::{align::MultiMapR, qc::QcGeneQuant, transcript::Gene};
use bed_utils::bed::map::GIntervalMap;
use crate::transcript::annotate::GIntervalMapExt;
use super::{
    annotate::{AnnotationRegion, TranscriptAlignment, SpliceState}, de_dups::count_unique_umi, AlignmentAnnotator, AnnotatedAlignment
};



/// I use this to collect intron validation information outside
#[derive(Debug, Clone, Default)]
pub struct IntronValidationCollector {
    /// Map of transcript_id -> set of intron indices that need validation
    pub validation_requests: std::collections::HashMap<String, HashSet<usize>>,
}

impl IntronValidationCollector {
    pub fn new() -> Self {
        Self {
            validation_requests: std::collections::HashMap::new(),
        }
    }

    /// Add a validation request for a specific transcript and intron
    pub fn add_validation_request(&mut self, transcript_id: String, intron_index: usize) {
        self.validation_requests
            .entry(transcript_id)
            .or_insert_with(HashSet::new)
            .insert(intron_index);
    }

    /// Add multiple validation requests from transcript alignments
    pub fn collect_from_alignments(&mut self, alignments: &[TranscriptAlignment]) {
        for alignment in alignments {
            if let Some(transcript_id) = alignment.transcript_id() {
                for &intron_index in &alignment.validation_requests() {
                    self.add_validation_request(transcript_id.to_string(), intron_index);
                }
            }
        }
    }

    /// Get the total number of validation requests
    pub fn total_requests(&self) -> usize {
        self.validation_requests.values().map(|set| set.len()).sum()
    }

    /// Check if there are any validation requests
    pub fn is_empty(&self) -> bool {
        self.validation_requests.is_empty()
    }

    /// Get all validation requests as a vector of (transcript_id, intron_index) pairs
    pub fn get_all_requests(&self) -> Vec<(String, usize)> {
        let mut requests = Vec::new();
        for (transcript_id, intron_indices) in &self.validation_requests {
            for &intron_index in intron_indices {
                requests.push((transcript_id.clone(), intron_index));
            }
        }
        requests
    }


    /// Extract validation requests from ALL transcript alignments (both intronic and exonic)
    /// This aggregates validation information regardless of final alignment classification
    pub fn extract_validation_requests(alignments: &[TranscriptAlignment]) -> Vec<(String, usize)> {
        let mut validation_requests = Vec::new();
        for alignment in alignments {
            // For exonic alignments, use the transcript ID from exon_align
            if let Some(transcript_id) = alignment.transcript_id() {
                for &intron_index in &alignment.validation_requests() {
                    validation_requests.push((transcript_id.to_string(), intron_index));
                }
            } else {
                // For intronic alignments, we still want to collect validation requests
                // but we don't have a transcript ID. We'll use the gene ID as a fallback.
                // This ensures we don't lose validation information from intronic alignments.
                for &intron_index in &alignment.validation_requests() {
                    validation_requests.push((alignment.gene.id.clone(), intron_index));
                }
            }
        }
        validation_requests
    }
}


#[derive(Debug, Encode, Decode)]
pub struct GeneAlignment {
    pub idx: usize,
    pub umi: Option<String>,
    pub align_type: AnnotationRegion,
}

#[derive(Debug)]
pub struct Quantifier {
    pub annotator: AlignmentAnnotator, // make it public for testing purposes
    genes: IndexMap<String, Gene>,
    temp_dir: Option<PathBuf>,
    chunk_size: usize,
    mito_genes: HashSet<usize>,
}

impl Quantifier {
    pub fn new(annotator: AlignmentAnnotator) -> Self {
        let genes = annotator
            .transcripts
            .iter()
            .map(|(_, t)| {
                let g = t.gene.clone();
                (g.id.clone(), g)
            })
            .collect();
        Self {
            annotator,
            genes,
            temp_dir: None,
            chunk_size: 50000000,
            mito_genes: HashSet::new(),
        }
    }

    pub fn add_mito_dna(&mut self, mito_chr: &str) {
        let iter = self.annotator.transcripts.iter().flat_map(|(_, t)| {
            if t.chrom == mito_chr {
                Some(self.genes.get_full(&t.gene.id).unwrap().0)
            } else {
                None
            }
        });
        self.mito_genes.extend(iter);
    }

    pub fn quantify<'a, I, P>(
        &'a mut self,
        header: &'a Header,
        records: I,
        output: P,
        splice_aware: bool,
    ) -> Result<QcGeneQuant>
    where
        I: Iterator<Item = Vec<(Option<MultiMapR>, Option<MultiMapR>)>> + 'a,
        P: AsRef<std::path::Path>,
    {
        let mut qc = QcGeneQuant::default();

        let adata: AnnData<H5> = AnnData::new(output)?;
        let num_cols = self.genes.len();

        let mut barcodes = Vec::new();
        let mut exon_count: Vec<u64> = Vec::new();
        let mut intron_count: Vec<u64> = Vec::new();
        let mut mito_count: Vec<u64> = Vec::new();
        let mut validation_collector = IntronValidationCollector::new();

        {
            if splice_aware {
                // Extract transcripts, modify them, and rebuild the map
                let transcripts = std::mem::replace(&mut self.annotator.transcripts, GIntervalMap::new());
                let updated_transcripts = transcripts.map_values_mut(|transcript| {
                    transcript.make_intron_by_exons();
                });

                // Update the annotator's transcripts
                self.annotator.transcripts = updated_transcripts;
            }
            let temp_dir = self.temp_dir.clone();
            let chunk_size = self.chunk_size;
            let mito_genes = self.mito_genes.clone();
            let (gene_alignments, all_validation_requests) = if splice_aware {
                // When splice_aware is true, collect both alignments and validation requests
                let tx_alignments = records.flat_map(|recs| {
                    qc.total_reads += recs.len() as u64;
                    recs.into_par_iter()
                        .filter_map(|(r1, r2)| {
                            self.make_gene_alignment(header, r1, r2, splice_aware).map(|(barcode, gene_alignment, validation_requests, splice_state, intron_mapping)| {
                                (barcode, gene_alignment, validation_requests)
                            })
                        })
                        .collect::<Vec<_>>()
                });

                let mut all_validation_requests = Vec::new();
                let gene_alignments: Vec<_> = tx_alignments.into_iter().map(|(barcode, gene_alignment, validation_requests)| {
                    all_validation_requests.extend(validation_requests);
                    (barcode, gene_alignment)
                }).collect();

                (gene_alignments, all_validation_requests)
            } else {
                // When splice_aware is false, only collect alignments (no validation overhead)
                // There might be a minimal time cost compared to the original function, but I think it's acceptable
                let gene_alignments: Vec<_> = records.flat_map(|recs| {
                    qc.total_reads += recs.len() as u64;
                    recs.into_par_iter()
                        .filter_map(|(r1, r2)| {
                            self.make_gene_alignment(header, r1, r2, splice_aware).map(|(barcode, gene_alignment, _, _, _)| {
                                (barcode, gene_alignment)
                            })
                        })
                        .collect::<Vec<_>>()
                }).collect();
                (gene_alignments, Vec::new())
            };

            // Process validation requests only if splice_aware is enabled
            if splice_aware {
                for (transcript_id, intron_index) in all_validation_requests {
                    validation_collector.add_validation_request(transcript_id, intron_index);
                }
            }
            let tx_alignments_chunks =
                sort_alignments(gene_alignments.into_iter(), temp_dir.as_ref(), chunk_size)
                    .chunk_by(|x| x.0.clone());
            let tx_alignments_chunks = tx_alignments_chunks
                .into_iter()
                .map(|(barcode, group)| {
                    barcodes.push(barcode);
                    group.map(|x| x.1).collect::<Vec<_>>()
                })
                .chunks(500);
            let counts = tx_alignments_chunks.into_iter().map(|chunk| {
                let results: Vec<_> = chunk
                    .collect::<Vec<_>>()
                    .into_par_iter()
                    .map(|alignments| count_unique_umi(alignments, &mito_genes))
                    .collect();
                results.iter().for_each(|r| {
                    qc.total_umi += r.total_umi;
                    qc.unique_umi += r.unique_umi;
                    exon_count.push(r.uniq_exon);
                    intron_count.push(r.uniq_intron);
                    mito_count.push(r.uniq_mito);
                });
                let (nrows, ncols, indptr, indices, data) = to_csr_data(
                    results
                        .into_iter()
                        .map(|r| r.into_counts().collect()),
                    num_cols,
                );
                from_csr_data(nrows, ncols, indptr, indices, data).unwrap()
            });

            adata.set_x_from_iter(counts)?;
        }

        adata.set_obs_names(barcodes.into())?;
        adata.set_obs(df!(
            "exon_count" => exon_count,
            "intron_count" => intron_count,
            "mitochondrial_count" => mito_count
        )?)?;

        adata.set_var_names(self.genes.values().map(|g| g.id.clone()).collect())?;
        adata.set_var(df!(
            "gene_name" => self.genes.values().map(|g| g.name.clone()).collect::<Vec<_>>()
        )?)?;

        
        // Log validation summary only if splice_aware is enabled
        if splice_aware && !validation_collector.is_empty() {
            log::info!(
                "Collected {} validation requests for {} transcripts",
                validation_collector.total_requests(),
                validation_collector.validation_requests.len()
            );
        }

        adata.close()?;
        Ok(qc)
    }

    /// I make it public for testing purposes. Better to be private later
    pub fn make_gene_alignment(
        &self,
        header: &Header,
        rec1: Option<MultiMapR>,
        rec2: Option<MultiMapR>,
        splice_aware: bool,
    ) -> Option<(String, GeneAlignment, Vec<(String, usize)>,  SpliceState, std::collections::HashMap<String, std::collections::HashSet<usize>>)> {
        let barcode;
        let umi;
        let anno = if rec1.is_some() && rec2.is_some() {
            let rec1 = rec1.unwrap();
            let rec2 = rec2.unwrap();
            barcode = rec1.barcode().unwrap()?;
            umi = rec1.umi().unwrap();
            self.annotator.annotate_alignments_pe(header, rec1, rec2,splice_aware)
        } else {
            let rec = rec1.or(rec2).unwrap();
            barcode = rec.barcode().unwrap()?;
            umi = rec.umi().unwrap();
            self.annotator.annotate_alignments_se(header, rec,splice_aware)
        }?;

        let gene_id;
        let align_type;
        let mut validation_requests = Vec::new();
        let mut splice_state = SpliceState::Unspliced;
        let mut intron_mapping = std::collections::HashMap::new();

        match anno {
            AnnotatedAlignment::PeMapped(a1, a2, anno) => {
                let genes = anno.genes;
                if genes.len() != 1 {
                    return None;
                }
                let gene = genes.iter().next().unwrap();
                align_type = match (a1.region, a2.region) {
                    (AnnotationRegion::Intronic, _) => AnnotationRegion::Intronic,
                    (_, AnnotationRegion::Intronic) => AnnotationRegion::Intronic,
                    _ => AnnotationRegion::Exonic,
                };
                gene_id = self.genes.get_full(&gene.id).unwrap().0;
                
               
                // Extract validation information only if splice_aware is enabled
                if splice_aware {
                    // Extract validation requests from aggregated annotation data
                    validation_requests.extend(a1.validation_requests.iter().cloned());
                    validation_requests.extend(a2.validation_requests.iter().cloned());

                    // Aggregate splice states: if any is Spliced, result is Spliced
                    let states = vec![a1.splice_state, a2.splice_state];
                    splice_state = if states.iter().any(|&s| s == SpliceState::Spliced) {
                        SpliceState::Spliced
                    } else if states.iter().all(|&s| s == SpliceState::Unspliced) {
                        SpliceState::Unspliced
                    } else {
                        SpliceState::Ambiguous
                    };

                    // Merge intron mappings from both annotations
                    for (transcript_id, intron_set) in &a1.intron_mapping {
                        intron_mapping.entry(transcript_id.clone())
                            .or_insert_with(std::collections::HashSet::new)
                            .extend(intron_set);
                    }
                    for (transcript_id, intron_set) in &a2.intron_mapping {
                        intron_mapping.entry(transcript_id.clone())
                            .or_insert_with(std::collections::HashSet::new)
                            .extend(intron_set);
                    }
                }
            }
            AnnotatedAlignment::SeMapped(anno) => {
                let genes = anno.genes;
                if genes.len() != 1 {
                    return None;
                }
                let gene = genes.iter().next().unwrap();
                align_type = anno.region;
                gene_id = self.genes.get_full(&gene.id).unwrap().0;
                // Extract validation information only if splice_aware is enabled
                // Extract validation information only if splice_aware is enabled
                if splice_aware {
                    // Extract validation requests from aggregated annotation data
                    validation_requests.extend(anno.validation_requests.iter().cloned());

                    // Extract splice state
                    splice_state = anno.splice_state;

                    // Extract intron mapping
                    intron_mapping = anno.intron_mapping.clone();
                }
            }
        }

        let alignment = GeneAlignment {
            idx: gene_id,
            umi,
            align_type,
        };
        Some((barcode, alignment, validation_requests, splice_state, intron_mapping))
    }
}

fn sort_alignments<I, P>(
    alignments: I,
    temp_dir: Option<P>,
    chunk_size: usize,
) -> impl ExactSizeIterator<Item = (String, GeneAlignment)>
where
    I: Iterator<Item = (String, GeneAlignment)>,
    P: AsRef<std::path::Path>,
{
    let mut sorter = ExternalSorterBuilder::new()
        .with_chunk_size(chunk_size)
        .with_compression(2);
    if let Some(tmp) = temp_dir {
        sorter = sorter.with_tmp_dir(tmp);
    }
    sorter
        .build()
        .unwrap()
        .sort_by(alignments, |a, b| a.0.cmp(&b.0))
        .unwrap()
        .map(|x| x.unwrap())
}

