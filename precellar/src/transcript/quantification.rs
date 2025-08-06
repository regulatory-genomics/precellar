use anndata::{
    data::utils::{from_csr_data, to_csr_data},
    AnnData, AnnDataOp,AxisArraysOp
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
use super::{
    annotate::{AnnotationRegion, SpliceState}, de_dups::count_unique_umi, AlignmentAnnotator, AnnotatedAlignment
};

#[derive(Debug, Encode, Decode)]
pub struct GeneAlignment {
    pub idx: usize,
    pub umi: Option<String>,
    pub align_type: AnnotationRegion,
    pub splice_state: SpliceState,
}

#[derive(Debug)]
pub struct Quantifier {
    annotator: AlignmentAnnotator, 
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

    pub fn quantify<I, P>(
        self,
        header: &Header,
        records: I,
        output: P,
        splice_aware: bool,
    ) -> Result<QcGeneQuant>
    where
        I: Iterator<Item = Vec<(Option<MultiMapR>, Option<MultiMapR>)>>,
        P: AsRef<std::path::Path>,
    {
        let mut qc = QcGeneQuant::default();

        let adata: AnnData<H5> = AnnData::new(output)?;
        let num_cols = self.genes.len();

        let mut barcodes = Vec::new();
        let mut exon_count: Vec<u64> = Vec::new();
        let mut intron_count: Vec<u64> = Vec::new();
        let mut mito_count: Vec<u64> = Vec::new();
        let mut spliced_count: Vec<u64> = Vec::new();
        let mut unspliced_count: Vec<u64> = Vec::new();
        let mut ambiguous_count: Vec<u64> = Vec::new();

        {
            let temp_dir = self.temp_dir.clone();
            let chunk_size = self.chunk_size;
            let mito_genes = self.mito_genes.clone();
            let (gene_alignments) = if splice_aware {
                let tx_alignments = records.flat_map(|recs| {
                    qc.total_reads += recs.len() as u64;
                    recs.into_par_iter()
                        .filter_map(|(r1, r2)| {
                            self.make_gene_alignment(header, r1, r2, splice_aware).map(|(barcode, gene_alignment, intron_mapping)| {
                                (barcode, gene_alignment, intron_mapping)
                            })
                        })
                        .collect::<Vec<_>>()
                });
                let gene_alignments: Vec<_> = tx_alignments.into_iter().map(|(barcode, gene_alignment, _)| {
                    (barcode, gene_alignment)
                }).collect();
                (gene_alignments)
            } else {
                // When splice_aware is false, only collect alignments (no validation overhead)
                // There might be a minimal time cost compared to the original function, but I think it's acceptable
                let gene_alignments: Vec<_> = records.flat_map(|recs| {
                    qc.total_reads += recs.len() as u64;
                    recs.into_par_iter()
                        .filter_map(|(r1, r2)| {
                            self.make_gene_alignment(header, r1, r2, splice_aware).map(|(barcode, gene_alignment, _)| {
                                (barcode, gene_alignment)
                            })
                        })
                        .collect::<Vec<_>>()
                }).collect();
                (gene_alignments)
            };


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
            // Collect all results first
            let all_results: Vec<Vec<_>> = tx_alignments_chunks.into_iter().map(|chunk| {
                chunk
                    .collect::<Vec<_>>()
                    .into_par_iter()
                    .map(|alignments| count_unique_umi(alignments, &mito_genes, splice_aware))
                    .collect()
            }).collect(); 

            // Process results for QC metrics
            all_results.iter().for_each(|chunk_results| {
                chunk_results.iter().for_each(|r| {
                    qc.total_umi += r.total_umi;
                    qc.unique_umi += r.unique_umi;
                    exon_count.push(r.uniq_exon);
                    intron_count.push(r.uniq_intron);
                    mito_count.push(r.uniq_mito);
                    if splice_aware {
                        spliced_count.push(r.uniq_spliced);
                        unspliced_count.push(r.uniq_unspliced);
                        ambiguous_count.push(r.uniq_ambiguous);
                    }
                });
            });

            // Create main count matrix (combined exon + intron)
            let counts = all_results.iter().map(|chunk_results| {
                let (nrows, ncols, indptr, indices, data) = to_csr_data(
                    chunk_results
                        .iter()
                        .map(|r| r.get_combined_counts()),
                    num_cols,
                );
                from_csr_data(nrows, ncols, indptr, indices, data).unwrap()
            });
            adata.set_x_from_iter(counts)?;

            // Create splice-state-specific layers if splice_aware is enabled
            if splice_aware {
                // Create spliced layer
                let spliced_counts: Vec<_> = all_results.iter().map(|chunk_results| {
                    let (nrows, ncols, indptr, indices, data) = to_csr_data(
                        chunk_results.iter().map(|r| r.get_spliced_counts()),
                        num_cols,
                    );
                    from_csr_data(nrows, ncols, indptr, indices, data).unwrap()
                }).collect();
                adata.layers().add_iter("spliced", spliced_counts.into_iter())?; // It seems adata-rs don't have layer function now. THis is a placeholder for future implementation

                // Create unspliced layer
                let unspliced_counts: Vec<_> = all_results.iter().map(|chunk_results| {
                    let (nrows, ncols, indptr, indices, data) = to_csr_data(
                        chunk_results
                            .iter()
                            .map(|r| r.get_unspliced_counts()),
                        num_cols,
                    );
                    from_csr_data(nrows, ncols, indptr, indices, data).unwrap()
                }).collect();
                adata.layers().add_iter("unspliced", unspliced_counts.into_iter())?; // It seems adata-rs don't have layer function now. THis is a placeholder for future implementation

                // Create ambiguous layer
                let ambiguous_counts = all_results.iter().map(|chunk_results| {
                    let (nrows, ncols, indptr, indices, data) = to_csr_data(
                        chunk_results
                            .iter()
                            .map(|r| r.get_ambiguous_counts()),
                        num_cols,
                    );
                    from_csr_data(nrows, ncols, indptr, indices, data).unwrap()
                });
                adata.layers().add_iter("ambiguous", ambiguous_counts.into_iter())?; // It seems adata-rs don't have layer function now. THis is a placeholder for future implementation
            }
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

        adata.close()?;
        Ok(qc)
    }

    fn make_gene_alignment(
        &self,
        header: &Header,
        rec1: Option<MultiMapR>,
        rec2: Option<MultiMapR>,
        splice_aware: bool,
    ) -> Option<(String, GeneAlignment,  std::collections::HashMap<String, std::collections::HashSet<usize>>)> {
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
        let mut splice_state = SpliceState::Undetermined;
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
                if splice_aware {
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
            splice_state,
        };
        Some((barcode, alignment, intron_mapping))
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


#[cfg(test)]
mod tests {
    use super::*;
    use crate::transcript::{Gene};
    use crate::transcript::annotate::{AlignmentAnnotator, SpliceState, AnnotationRegion};
    use noodles::sam::{self as sam, header::record::value::{map::ReferenceSequence, Map}};
    use bstr::BString;
    use std::num::NonZeroUsize;
    use crate::transcript::transcriptome::{Transcript, Exons};
    use noodles::sam::record::data::field::value::base_modifications::group::Strand as NoodlesStrand;

    fn create_test_gene() -> Gene {
        Gene {
            id: "gene1".to_string(),
            name: "GENE1".to_string(),
        }
    }

    fn create_test_transcript() -> Transcript {
        // Create a transcript with 3 exons that have gaps between them (will create introns)
        let exons = Exons::new(vec![
            (100, 200),  // Exon 1: 100-200
            (300, 400),  // Exon 2: 300-400 (gap: 200-300)
            (500, 600),  // Exon 3: 500-600 (gap: 400-500)
        ]).unwrap();

        Transcript {
            id: "transcript1".to_string(),
            chrom: "chr1".to_string(),
            start: 100,
            end: 600,
            strand: NoodlesStrand::Forward,
            gene: create_test_gene(),
            exons,
        }
    }

    fn create_test_header() -> sam::Header {
        let mut builder = sam::Header::builder();
        let reference_sequences = [(
            BString::from("chr1"),
            Map::<ReferenceSequence>::new(NonZeroUsize::try_from(1000).unwrap()),
        )]
        .into_iter()
        .collect();
        builder = builder.set_reference_sequences(reference_sequences);
        builder.build()
    }

    #[test]
    fn test_splice_state_specific_counting() {
        // Create test data
        let transcript = create_test_transcript();
        let annotator = AlignmentAnnotator::new(vec![transcript]);

        // Create test alignments with different splice states
        let alignments = vec![
            GeneAlignment {
                idx: 0, // gene index
                umi: Some("AAAA".to_string()),
                align_type: AnnotationRegion::Exonic,
                splice_state: SpliceState::Spliced,
            },
            GeneAlignment {
                idx: 0,
                umi: Some("TTTT".to_string()),
                align_type: AnnotationRegion::Exonic,
                splice_state: SpliceState::Unspliced,
            },
            GeneAlignment {
                idx: 0,
                umi: Some("CCCC".to_string()),
                align_type: AnnotationRegion::Intronic,
                splice_state: SpliceState::Ambiguous,
            },
            GeneAlignment {
                idx: 0,
                umi: Some("GGGG".to_string()),
                align_type: AnnotationRegion::Exonic,
                splice_state: SpliceState::Spliced,
            },
        ];

        // Test the count_unique_umi function with splice_aware = true
        let mito_genes = std::collections::HashSet::new();
        let result = count_unique_umi(alignments, &mito_genes, true);

        // Verify splice-state-specific counts
        assert_eq!(result.uniq_spliced, 2, "Should have 2 spliced UMIs");
        assert_eq!(result.uniq_unspliced, 1, "Should have 1 unspliced UMI");
        assert_eq!(result.uniq_ambiguous, 1, "Should have 1 ambiguous UMI");

        // Verify total counts
        assert_eq!(result.unique_umi, 4, "Should have 4 total unique UMIs");

        // Test the accessor methods
        let spliced_counts: Vec<_> = result.get_spliced_counts();
        let unspliced_counts: Vec<_> = result.get_unspliced_counts();
        let ambiguous_counts: Vec<_> = result.get_ambiguous_counts();

        assert_eq!(spliced_counts.len(), 1, "Should have counts for 1 gene");
        assert_eq!(spliced_counts[0], (0, 2), "Gene 0 should have 2 spliced counts");

        assert_eq!(unspliced_counts.len(), 1, "Should have counts for 1 gene");
        assert_eq!(unspliced_counts[0], (0, 1), "Gene 0 should have 1 unspliced count");

        assert_eq!(ambiguous_counts.len(), 1, "Should have counts for 1 gene");
        assert_eq!(ambiguous_counts[0], (0, 1), "Gene 0 should have 1 ambiguous count");
    }

    #[test]
    fn test_splice_aware_disabled() {
        // Test that when splice_aware is false, splice-state-specific counts are zero
        let alignments = vec![
            GeneAlignment {
                idx: 0,
                umi: Some("AAAA".to_string()),
                align_type: AnnotationRegion::Exonic,
                splice_state: SpliceState::Spliced,
            },
            GeneAlignment {
                idx: 0,
                umi: Some("TTTT".to_string()),
                align_type: AnnotationRegion::Intronic,
                splice_state: SpliceState::Unspliced,
            },
        ];

        let mito_genes = std::collections::HashSet::new();
        let result = count_unique_umi(alignments, &mito_genes, false);

        // When splice_aware is false, splice-state-specific counts should be zero
        assert_eq!(result.uniq_spliced, 0, "Spliced count should be 0 when splice_aware is false");
        assert_eq!(result.uniq_unspliced, 0, "Unspliced count should be 0 when splice_aware is false");
        assert_eq!(result.uniq_ambiguous, 0, "Ambiguous count should be 0 when splice_aware is false");

        // But traditional exon/intron counts should work
        assert_eq!(result.uniq_exon, 1, "Should have 1 exonic UMI");
        assert_eq!(result.uniq_intron, 1, "Should have 1 intronic UMI");
    }

}
