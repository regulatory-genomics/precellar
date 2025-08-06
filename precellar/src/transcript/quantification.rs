use anndata::{
    data::utils::{from_csr_data, to_csr_data},
    AnnData, AnnDataOp,ArrayData,AxisArraysOp
};
use ndarray::Array2;
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

        /// Validate introns in the provided transcripts based on collected requests
    /// Returns a set of (transcript_id, intron_index) pairs that were validated
    pub fn validate_introns(&self, transcripts: &mut [crate::transcript::Transcript]) -> HashSet<(String, usize)> {
        let mut validated_introns = HashSet::new();

        for transcript in transcripts.iter_mut() {
            if let Some(intron_indices) = self.validation_requests.get(&transcript.id) {
                for &intron_index in intron_indices {
                    if transcript.validate_intron(intron_index) {
                        validated_introns.insert((transcript.id.clone(), intron_index));
                    } 
                }
            }
        }
        validated_introns
    }

    /// Extract validation requests from ALL transcript alignments (both intronic and exonic)
    /// This aggregates validation information regardless of final alignment classification
    pub fn extract_validation_requests(alignments: &[TranscriptAlignment]) -> Vec<(String, usize)> {
        let mut validation_requests = Vec::new();
        for alignment in alignments {
            // For exonic alignments, use the transcript ID from exon_align
            if let Some(transcript_id) = alignment.transcript_id.clone() {
                for &intron_index in &alignment.validation_requests() {
                    validation_requests.push((transcript_id.to_string(), intron_index));
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
    pub splice_state: SpliceState,
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
        let mut spliced_count: Vec<u64> = Vec::new();
        let mut unspliced_count: Vec<u64> = Vec::new();
        let mut ambiguous_count: Vec<u64> = Vec::new();
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
                            self.make_gene_alignment(header, r1, r2, splice_aware).map(|(barcode, gene_alignment, validation_requests, intron_mapping)| {
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
                            self.make_gene_alignment(header, r1, r2, splice_aware).map(|(barcode, gene_alignment, _, _)| {
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

            
            // Step: Update splice states based on validated introns
            let mut updated_gene_alignments = gene_alignments;
            if splice_aware && !validation_collector.is_empty() {
                // Validate introns based on collected requests
                let mut transcripts = self.annotator.transcripts().to_vec();
                let validated_introns = validation_collector.validate_introns(&mut transcripts);

                // Update splice states for alignments that have validated introns
                updated_gene_alignments = updated_gene_alignments.into_iter().map(|(barcode, mut gene_alignment)| {
                    // For intronic alignments, check if any introns are validated for this gene
                    if gene_alignment.splice_state == SpliceState::Ambiguous {
                        // Check if there are any validated introns for this gene
                        let gene = &self.genes.get_index(gene_alignment.idx).unwrap().1;
                        // Check if there are any validated introns for any transcript of this gene
                        let has_validated_introns = validated_introns.iter()
                            .any(|(transcript_id, _)| {
                                // Check if this transcript belongs to the current gene
                                // We need to find the transcript and check its gene
                                transcripts.iter().any(|t| &t.id == transcript_id && t.gene.id == gene.id)
                            });

                        if has_validated_introns {
                            gene_alignment.splice_state = SpliceState::Unspliced;
                        }
                    }
                    (barcode, gene_alignment)
                }).collect();
            }
            
            let tx_alignments_chunks =
                sort_alignments(updated_gene_alignments.into_iter(), temp_dir.as_ref(), chunk_size)
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
    ) -> Option<(String, GeneAlignment, Vec<(String, usize)>,  std::collections::HashMap<String, std::collections::HashSet<usize>>)> {
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
            splice_state,
        };
        Some((barcode, alignment, validation_requests, intron_mapping))
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
    use crate::transcript::transcriptome::{Transcript, Exons, Introns};
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

        let introns = Introns::new(std::iter::empty()).unwrap(); // Start with empty introns

        Transcript {
            id: "transcript1".to_string(),
            chrom: "chr1".to_string(),
            start: 100,
            end: 600,
            strand: NoodlesStrand::Forward,
            gene: create_test_gene(),
            exons,
            introns,
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
        let mut quantifier = Quantifier::new(annotator);

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
    #[test]
    fn test_quantify_with_splice_layers() {
        use noodles::sam::{self as sam, header::record::value::{map::ReferenceSequence, Map}};
        use bstr::BString;
        use std::num::NonZeroUsize;
        use crate::align::MultiMapR;
        use noodles::sam::alignment::record::data::field::Tag;
        use noodles::sam::alignment::record_buf::data::field::Value;
        use tempfile::TempDir;

        // Create a test transcript with multiple exons that will generate introns
        let exons = Exons::new(vec![
            (100, 200),  // Exon 1: 100-200 (101bp)
            (300, 400),  // Exon 2: 300-400 (101bp) 
            (500, 600),  // Exon 3: 500-600 (101bp)
        ]).unwrap();

        let introns = Introns::new(std::iter::empty()).unwrap();

        let mut transcript = Transcript {
            id: "TRANSCRIPT001".to_string(),
            chrom: "chr1".to_string(),
            start: 100,
            end: 600,
            strand: NoodlesStrand::Forward,
            gene: Gene {
                id: "GENE001".to_string(),
                name: "TestGene".to_string(),
            },
            exons,
            introns,
        };

        // Generate introns from exons
        transcript.make_intron_by_exons();

        let annotator = AlignmentAnnotator::new(vec![transcript]);
        let mut quantifier = Quantifier::new(annotator);

        // Create a SAM header with proper reference sequence
        let reference_sequences = [(
            BString::from("chr1"),
            Map::<ReferenceSequence>::new(NonZeroUsize::try_from(248956422).unwrap()),
        )]
        .into_iter()
        .collect();
        let header = sam::Header::builder()
            .set_reference_sequences(reference_sequences)
            .build();

        // Create a temporary directory for output
        let temp_dir = TempDir::new().unwrap();
        let output_path = temp_dir.path().join("test_output.h5ad");

        // Create mock records with realistic SAM alignments
        let mut mock_records = Vec::new();
        
        // Record 1: Spliced alignment spanning all exons
        let sequence1 = "A".repeat(303); // 101 + 101 + 101 bases
        let quality1 = "I".repeat(303);
        let sam_line1 = format!(
            "read1\t0\tchr1\t100\t60\t101M99N101M99N101M\t*\t0\t0\t{}\t{}",
            sequence1, quality1
        );
        
        let sam_record1 = noodles::sam::Record::try_from(sam_line1.as_bytes()).unwrap();
        let mut record1 = noodles::sam::alignment::RecordBuf::try_from_alignment_record(&header, &sam_record1).unwrap();
        record1.data_mut().insert(Tag::from([b'C', b'B']), Value::String(BString::from("AAACCTGAGAAACCAT")));
        record1.data_mut().insert(Tag::from([b'U', b'B']), Value::String(BString::from("AAAAAAAAAA")));
        let multi_map1 = Some(MultiMapR::new(record1, None));

        // Record 2: Unspliced alignment within first exon
        let sequence2 = "T".repeat(50);
        let quality2 = "I".repeat(50);
        let sam_line2 = format!(
            "read2\t0\tchr1\t120\t60\t50M\t*\t0\t0\t{}\t{}",
            sequence2, quality2
        );
        
        let sam_record2 = noodles::sam::Record::try_from(sam_line2.as_bytes()).unwrap();
        let mut record2 = noodles::sam::alignment::RecordBuf::try_from_alignment_record(&header, &sam_record2).unwrap();
        record2.data_mut().insert(Tag::from([b'C', b'B']), Value::String(BString::from("AAACCTGAGAAACCAT")));
        record2.data_mut().insert(Tag::from([b'U', b'B']), Value::String(BString::from("TTTTTTTTTT")));
        let multi_map2 = Some(MultiMapR::new(record2, None));

        // Record 3: Intronic alignment
        let sequence3 = "G".repeat(30);
        let quality3 = "I".repeat(30);
        let sam_line3 = format!(
            "read3\t0\tchr1\t250\t60\t30M\t*\t0\t0\t{}\t{}",
            sequence3, quality3
        );
        
        let sam_record3 = noodles::sam::Record::try_from(sam_line3.as_bytes()).unwrap();
        let mut record3 = noodles::sam::alignment::RecordBuf::try_from_alignment_record(&header, &sam_record3).unwrap();
        record3.data_mut().insert(Tag::from([b'C', b'B']), Value::String(BString::from("AAACCTGAGAAACCAT")));
        record3.data_mut().insert(Tag::from([b'U', b'B']), Value::String(BString::from("CCCCCCCCCC")));
        let multi_map3 = Some(MultiMapR::new(record3, None));

        mock_records.push((multi_map1, None));
        mock_records.push((multi_map2, None));
        mock_records.push((multi_map3, None));

        let records = std::iter::once(mock_records).into_iter();

        // Test quantification with splice_aware = true
        let result = quantifier.quantify(&header, records, &output_path, true);
        
        // The test should complete without errors
        match result {
            Ok(_) => {
                // Verify the output file exists
                assert!(output_path.exists(), "Output file should be created");
                
                // Additional verification that splice-aware quantification worked
                println!("Quantification completed successfully with splice-aware mode enabled");
                println!("Output file created at: {:?}", output_path);
            },
            Err(e) => {
                panic!("Quantification failed with error: {:?}", e);
            }
        }

        // Test quantification with splice_aware = false for comparison
        let output_path_no_splice = temp_dir.path().join("test_output_no_splice.h5ad");
        
        // Recreate records for second test
        let mut mock_records_2 = Vec::new();
        
        let sam_record1_copy = noodles::sam::Record::try_from(sam_line1.as_bytes()).unwrap();
        let mut record1_copy = noodles::sam::alignment::RecordBuf::try_from_alignment_record(&header, &sam_record1_copy).unwrap();
        record1_copy.data_mut().insert(Tag::from([b'C', b'B']), Value::String(BString::from("AAACCTGAGAAACCAT")));
        record1_copy.data_mut().insert(Tag::from([b'U', b'B']), Value::String(BString::from("AAAAAAAAAA")));
        let multi_map1_copy = Some(MultiMapR::new(record1_copy, None));
        
        mock_records_2.push((multi_map1_copy, None));
        let records_no_splice = std::iter::once(mock_records_2).into_iter();

        let result_no_splice = quantifier.quantify(&header, records_no_splice, &output_path_no_splice, false);
        
        match result_no_splice {
            Ok(_) => {
                assert!(output_path_no_splice.exists(), "Output file without splice-aware should be created");
                println!("Quantification completed successfully with splice-aware mode disabled");
            },
            Err(e) => {
                panic!("Quantification without splice-aware failed with error: {:?}", e);
            }
        }
    }
}
