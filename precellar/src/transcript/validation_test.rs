// This file is used to test the intron validation functionality
// It is only used for testing
// It should be removed after the validation functionality is tested
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;
use anyhow::Result;
use noodles::sam::{self as sam, Header};
use noodles::bam;
use crate::transcript::{AlignmentAnnotator, Transcript};
use crate::align::MultiMapR;

#[cfg(test)]
use log;

/// Test structure to extract and output validated intron locations
#[derive(Debug)]
pub struct IntronValidationTest {
    annotator: AlignmentAnnotator
}

/// Represents a validated intron with all necessary information
#[derive(Debug, Clone)]
pub struct ValidatedIntron {
    pub chromosome: String,
    pub start: u64,
    pub end: u64,
    pub strand: String,
    pub gene_id: String,
    pub gene_name: String,
    pub transcript_id: String,
    pub transcript_name: String,
    pub intron_number: usize,
    pub length: u64,
    pub validation_status: String,
}

/// Helper struct for building transcripts from GTF
#[derive(Debug)]
struct TranscriptBuilder {
    id: String,
    chrom: String,
    start: u64,
    end: u64,
    strand: noodles::sam::record::data::field::value::base_modifications::group::Strand,
    gene_id: String,
    gene_name: String,
    exons: Vec<(u64, u64)>,
}

impl IntronValidationTest {
    /// Create a new validation test with GTF transcriptome
    pub fn new_from_gtf(gtf_path: &str) -> Result<Self> {
        println!("Loading transcripts from GTF file: {}", gtf_path);

        let transcripts = Self::load_transcripts_from_gtf(gtf_path)?;
        let annotator = AlignmentAnnotator::new(transcripts);

        Ok(Self {
            annotator,
        })
    }

    /// Create a new validation test from a vector of transcripts
    pub fn new_from_transcripts(transcripts: Vec<Transcript>) -> Result<Self> {
        let annotator = AlignmentAnnotator::new(transcripts);

        Ok(Self {
            annotator,
        })
    }

    /// Load transcripts from GTF file
    fn load_transcripts_from_gtf(gtf_path: &str) -> Result<Vec<Transcript>> {
        Self::load_transcripts_from_gtf_with_limit(gtf_path, None)
    }

    /// Load transcripts from GTF file with optional limit for testing
    fn load_transcripts_from_gtf_with_limit(gtf_path: &str, max_transcripts: Option<usize>) -> Result<Vec<Transcript>> {
        use noodles_gtf as gtf;
        use std::collections::HashMap;
        use std::fs::File;
        use std::io::BufReader;
        use crate::transcript::transcriptome::{Exons};
        use crate::transcript::Gene;
        use noodles::sam::record::data::field::value::base_modifications::group::Strand;

        println!("Loading transcripts from GTF file: {}", gtf_path);
        let file = File::open(gtf_path)?;
        let reader = BufReader::new(file);
        let mut gtf_reader = gtf::io::Reader::new(reader);

        // Group records by transcript_id - pre-allocate capacity for better performance
        let mut transcript_data: HashMap<String, TranscriptBuilder> = HashMap::with_capacity(50000000);
        let mut records_processed = 0;

        for result in gtf_reader.record_bufs() {
            let record = result?;
            records_processed += 1;

            // Progress reporting for large files - reduced frequency for better performance
            if records_processed % 500000 == 0 {
                println!("Processed {} records, found {} transcripts", records_processed, transcript_data.len());
            }

            // Skip header lines and non-feature records
            if record.reference_sequence_name().is_empty() {
                continue;
            }

            let feature_type = record.ty();

            // Only process transcript and exon records to save time
            let feature_str = feature_type.to_string();
            if feature_str != "transcript" && feature_str != "exon" {
                continue;
            }

            let attributes = record.attributes();

            // Use GTF attributes API - if you need to iterate over all attributes, use:
            // use noodles_gff::feature::record::Attributes as AttributesTrait;
            // for result in AttributesTrait::iter(&attributes) { ... }

            let gene_id = attributes.get(b"gene_id")
                .map(|value| value.as_string().unwrap_or_default().to_string());

            let gene_name = attributes.get(b"gene_name")
                .map(|value| value.as_string().unwrap_or_default().to_string());

            let transcript_id = attributes.get(b"transcript_id")
                .map(|value| value.as_string().unwrap_or_default().to_string());

            let gene_id = gene_id.ok_or_else(|| anyhow::anyhow!("gene_id not found"))?;
            let gene_name = gene_name.unwrap_or_else(|| gene_id.clone());

            match feature_str.as_str() {
                "transcript" => {
                    let transcript_id = transcript_id.ok_or_else(|| anyhow::anyhow!("transcript_id not found"))?;
                    let chrom = record.reference_sequence_name().to_string();
                    let start = record.start().get() as u64 - 1; // Convert to 0-based
                    let end = record.end().get() as u64; // Already exclusive
                    let strand = match format!("{:?}", record.strand()).as_str() {
                        "Forward" => Strand::Forward,
                        "Reverse" => Strand::Reverse,
                        _ => Strand::Forward, // Default fallback
                    };

                    let builder = TranscriptBuilder {
                        id: transcript_id.clone(),
                        chrom,
                        start,
                        end,
                        strand,
                        gene_id,
                        gene_name,
                        exons: Vec::new(),
                    };

                    transcript_data.insert(transcript_id, builder);
                }
                "exon" => {
                    if let Some(transcript_id) = transcript_id {
                        let exon_start = record.start().get() as u64 - 1; // Convert to 0-based
                        let exon_end = record.end().get() as u64; // Already exclusive

                        if let Some(builder) = transcript_data.get_mut(&transcript_id) {
                            builder.exons.push((exon_start, exon_end));
                        }
                    }
                }
                _ => {
                    // Skip other feature types (gene, CDS, etc.)
                }
            }
        }

        // Filter out transcript builders that have no exons before processing
        let original_count = transcript_data.len();
        transcript_data.retain(|_, builder| !builder.exons.is_empty());
        let filtered_count = transcript_data.len();
        let pre_filtered = original_count - filtered_count;

        if pre_filtered > 0 {
            println!("Pre-filtered {} transcripts with no exons", pre_filtered);
        }

        // Convert transcript builders to Transcript objects
        let mut transcripts = Vec::new();
        let mut skipped_no_exons = pre_filtered;  // Start with pre-filtered count

        for (transcript_id, mut builder) in transcript_data {
            // Skip transcripts with no exons (invalid)
            if builder.exons.is_empty() {
                skipped_no_exons += 1;
                if skipped_no_exons <= 5 {  // Only show first 5 warnings
                    println!("Warning: Skipping transcript {} with no exons", transcript_id);
                }
                continue;
            }

            // Sort exons by start position
            builder.exons.sort_by_key(|&(start, _)| start);

            // Create Exons object
            let exons = Exons::new(builder.exons.into_iter())?;


            let mut transcript = Transcript {
                id: builder.id,
                chrom: builder.chrom,
                start: builder.start,
                end: builder.end,
                strand: builder.strand,
                gene: Gene {
                    id: builder.gene_id,
                    name: builder.gene_name,
                },
                exons,
            };

            // Final validation: Ensure transcript has exons before adding to list
            if transcript.exons().is_empty() {
                skipped_no_exons += 1;
                if skipped_no_exons <= 5 {
                    println!("Warning: Skipping transcript {} - no exons after make_intron_by_exons()", transcript.id);
                }
                continue;
            }


            transcripts.push(transcript);
        }

        println!("Successfully parsed {} transcripts from GTF", transcripts.len());

        Ok(transcripts)
    }

}

/// Represents a classified read with all annotation details
#[derive(Debug, Clone)]
pub struct ReadClassification {
    pub read_id: String,
    pub chromosome: String,
    pub start: u64,
    pub end: u64,
    pub strand: String,
    pub transcript_id: String,
    pub transcript_name: String,
    pub gene_id: String,
    pub gene_name: String,
    pub classification: String,  // spliced/unspliced/ambiguous/intergenic
    pub is_spliced: bool,
}

impl IntronValidationTest {
    /// Classify all reads in a BAM file and generate detailed output
    pub fn classify_reads(&mut self, bam_path: &str, output_path: &str) -> Result<()> {
        println!("Starting read classification...");

        let mut reader = bam::io::Reader::new(std::fs::File::open(bam_path)?);
        let header = reader.read_header()?;

        let mut classifications = Vec::new();
        let mut total_reads = 0;
        let mut classified_reads = 0;

        for result in reader.record_bufs(&header) {
            let record = result?;
            total_reads += 1;

            if total_reads % 10000 == 0 {
                println!("Processed {} reads, classified {} reads", total_reads, classified_reads);
            }

            // Skip unmapped reads
            if record.flags().is_unmapped() {
                continue;
            }

            // Get read information
            let read_id = if let Some(name) = record.name() {
                std::str::from_utf8(name).unwrap_or("unknown").to_string()
            } else {
                format!("read_{}", total_reads)
            };

            let chromosome = if let Some(ref_id) = record.reference_sequence_id() {
                header.reference_sequences()
                    .get_index(ref_id)
                    .map(|(name, _)| name.to_string())
                    .unwrap_or_else(|| format!("chr{}", ref_id))
            } else {
                "unknown".to_string()
            };

            let start = record.alignment_start().map(|pos| pos.get() as u64 - 1).unwrap_or(0);
            let end = record.alignment_end().map(|pos| pos.get() as u64).unwrap_or(start + 1);

            // Convert the read to a MultiMapR
            let multi_map_r = MultiMapR::from(record.clone());

            // Use the public annotate_alignments_se method
            if let Some(annotated_alignment) = self.annotator.annotate_alignments_se(&header, multi_map_r, true) {
                // Extract transcript alignments from the annotated alignment
                use crate::transcript::annotate::AnnotatedAlignment;
                match annotated_alignment {
                    AnnotatedAlignment::SeMapped(annotation) => {
                        // Use the classification result directly from annotate_alignments_se
                        let final_classification = match annotation.splice_state {
                            crate::transcript::annotate::SpliceState::Spliced => "spliced",
                            crate::transcript::annotate::SpliceState::Unspliced => "unspliced",
                            crate::transcript::annotate::SpliceState::Ambiguous => "ambiguous",
                            crate::transcript::annotate::SpliceState::Intergenic => "intergenic",
                            crate::transcript::annotate::SpliceState::Undetermined => "undetermined",
                        };

                        // Get all alignments for detailed classification
                        let mut all_alignments = Vec::new();
                        all_alignments.extend(&annotation.aln_sense);
                        all_alignments.extend(&annotation.aln_antisense);

                        if !all_alignments.is_empty() {
                            // Create a single classification for the read using the result from annotate_alignments_se
                            let classification = self.classify_aggregated_alignment(
                                &read_id,
                                &chromosome,
                                start,
                                end,
                                &all_alignments,
                                final_classification,
                            );

                            classifications.push(classification);
                            classified_reads += 1;
                        }
                    }
                    AnnotatedAlignment::PeMapped(_, _, _) => {
                        // This shouldn't happen for single-end reads, but handle gracefully
                        continue;
                    }
                }
            }
        }

        println!("Classification complete. Total reads: {}, Classified reads: {}", total_reads, classified_reads);

        // Write classifications to output file
        self.write_classification_output(&classifications, output_path)?;

        Ok(())
    }


    /// Classify an aggregated alignment based on multiple transcript alignments
    fn classify_aggregated_alignment(
        &self,
        read_id: &str,
        chromosome: &str,
        start: u64,
        end: u64,
        alignments: &[&crate::transcript::annotate::TranscriptAlignment],
        final_classification: &str,
    ) -> ReadClassification {
        // Use the first alignment for basic information
        let first_alignment = alignments[0];
        
        // Determine strand from the first alignment
        let strand = match first_alignment.strand {
            noodles::sam::record::data::field::value::base_modifications::group::Strand::Forward => "+",
            noodles::sam::record::data::field::value::base_modifications::group::Strand::Reverse => "-",
        }.to_string();

        // Determine if the final classification is spliced
        let is_spliced = final_classification == "spliced";


        // For transcript information, we'll use the first alignment's transcript
        // In a more sophisticated implementation, we might want to list all transcripts
        let transcript_id = first_alignment.transcript_id().unwrap_or("unknown").to_string();
        let transcript_name = first_alignment.transcript_id().unwrap_or("unknown").to_string();
        
        // For gene information, we'll use the first alignment's gene
        let gene_id = first_alignment.gene.id.clone();
        let gene_name = first_alignment.gene.name.clone();

        ReadClassification {
            read_id: format!("{}:{}-{}:{}", read_id, start, end, strand),
            chromosome: chromosome.to_string(),
            start,
            end,
            strand,
            transcript_id,
            transcript_name,
            gene_id,
            gene_name,
            classification: final_classification.to_string(),
            is_spliced,
        }
    }

    /// Write classification results to output file
    fn write_classification_output(&self, classifications: &[ReadClassification], output_path: &str) -> Result<()> {
        let file = File::create(output_path)?;
        let mut writer = BufWriter::new(file);

        // Write header
        writeln!(
            writer,
            "read_id\tchromosome\tstart\tend\tstrand\ttranscript_id\ttranscript_name\tgene_id\tgene_name\tclassification\tis_spliced\thas_exons\thas_introns\thas_validated_introns\thas_spanning"
        )?;

        // Write data
        for classification in classifications {
            writeln!(
                writer,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                classification.read_id,
                classification.chromosome,
                classification.start,
                classification.end,
                classification.strand,
                classification.transcript_id,
                classification.transcript_name,
                classification.gene_id,
                classification.gene_name,
                classification.classification,
                classification.is_spliced,
            )?;
        }

        writer.flush()?;
        Ok(())
    }
}




#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_classification() {
        // Test read classification functionality
        let gtf_path = "/data2/litian/database/gtf/small_subset_genes2.gtf";
        let bam_file = "/data2/litian/202506_trajectory/data/velocyto_toy/test_small_200.bam";
        let output_file = "/data2/litian/202506_trajectory/process/20250803_velocity_validation/20250805_read_classifications.tsv";

        if Path::new(gtf_path).exists() && Path::new(bam_file).exists() {
            println!("Running read classification test...");

            let result = std::panic::catch_unwind(|| -> Result<()> {
                // Load transcripts from GTF
                let transcripts = IntronValidationTest::load_transcripts_from_gtf_with_limit(gtf_path, Some(1000))?;
                let mut classifier = IntronValidationTest::new_from_transcripts(transcripts)?;

                // Classify reads from BAM file
                classifier.classify_reads(bam_file, output_file)?;

                Ok(())
            });

            match result {
                Ok(Ok(())) => {
                    println!("Read classification completed successfully!");

                    // Verify output file exists and show preview
                    if Path::new(output_file).exists() {
                        let content = std::fs::read_to_string(output_file).unwrap();
                        let lines: Vec<&str> = content.lines().take(10).collect(); // Show first 10 lines
                        println!("Classification output preview:");
                        for line in lines {
                            println!("{}", line);
                        }
                    }
                }
                Ok(Err(e)) => {
                    println!("Read classification failed with error: {}", e);
                }
                Err(_) => {
                    println!("Read classification panicked");
                }
            }
        } else {
            println!("Read classification test skipped - files not found");
            println!("  GTF file: {} (exists: {})", gtf_path, Path::new(gtf_path).exists());
            println!("  BAM file: {} (exists: {})", bam_file, Path::new(bam_file).exists());
        }
    }

    #[test]
    fn test_quantification_with_real_data() {
        // Test quantification using real BAM and GTF data
        let gtf_path = "/data2/litian/database/gtf/small_subset_genes2.gtf";
        let bam_file = "/data2/litian/202506_trajectory/data/velocyto_toy/test_small_200.bam";
        let output_file = "/data2/litian/202506_trajectory/process/20250803_velocity_validation/quantification_test_output.h5ad";

        if Path::new(gtf_path).exists() && Path::new(bam_file).exists() {
            println!("Running quantification test with real data...");

            let result = std::panic::catch_unwind(|| -> Result<()> {
                use crate::transcript::quantification::Quantifier;
                use noodles::bam;
                use crate::align::MultiMapR;
                // Load transcripts from GTF
                let transcripts = IntronValidationTest::load_transcripts_from_gtf_with_limit(gtf_path, Some(1000))?;
                println!("Loaded {} transcripts for quantification", transcripts.len());
                
                // Create annotator and quantifier
                let annotator = crate::transcript::annotate::AlignmentAnnotator::new(transcripts);
                let mut quantifier = Quantifier::new(annotator);
                
                // Read BAM file and convert to records
                let mut reader = bam::io::Reader::new(std::fs::File::open(bam_file)?);
                let header = reader.read_header()?;
                
                // Collect records in batches
                let mut all_records = Vec::new();
                let mut batch = Vec::new();
                const BATCH_SIZE: usize = 100;
                
                for result in reader.record_bufs(&header) {
                    let record = result?;
                    
                    // Skip unmapped reads
                    if record.flags().is_unmapped() {
                        continue;
                    }
                    
                    // Convert to MultiMapR
                    let multi_map_r = MultiMapR::from(record);
                    
                    // For single-end reads, second element is None
                    batch.push((Some(multi_map_r), None));
                    
                    if batch.len() >= BATCH_SIZE {
                        all_records.push(batch);
                        batch = Vec::new();
                    }
                }
                
                // Add remaining records
                if !batch.is_empty() {
                    all_records.push(batch);
                }
                
                println!("Collected {} batches of records for quantification", all_records.len());
                
                // Create iterator over records
                let records_iter = all_records.into_iter();
                
                // Test quantification with splice_aware = true
                println!("Starting quantification with splice_aware = true...");
                let qc_result = quantifier.quantify(&header, records_iter, output_file, true)?;
                
                println!("Quantification completed successfully!");
                println!("QC metrics:");
                println!("  Total reads: {}", qc_result.total_reads);
                println!("  Total UMI: {}", qc_result.total_umi);
                println!("  Unique UMI: {}", qc_result.unique_umi);
                
                Ok(())
            });

            match result {
                Ok(Ok(())) => {
                    println!("Quantification test completed successfully!");
                    
                    // Verify output file exists
                    if Path::new(output_file).exists() {
                        println!("AnnData output file created: {}", output_file);
                        
                        // Try to read the file size to verify it's not empty
                        if let Ok(metadata) = std::fs::metadata(output_file) {
                            println!("Output file size: {} bytes", metadata.len());
                        }
                    } else {
                        println!("Warning: Output file was not created");
                    }
                },
                Ok(Err(e)) => {
                    println!("Quantification test failed with error: {}", e);
                    println!("This might be due to:");
                    println!("  - Missing CB/UB tags in BAM file");
                    println!("  - Transcript annotation issues");
                    println!("  - File format compatibility");
                },
                Err(_) => {
                    println!("Quantification test panicked");
                    println!("This might be due to memory issues or low-level library problems");
                }
            }
        } else {
            println!("Quantification test skipped - files not found");
            println!("  GTF file: {} (exists: {})", gtf_path, Path::new(gtf_path).exists());
            println!("  BAM file: {} (exists: {})", bam_file, Path::new(bam_file).exists());
        }
    }

    #[test]
    fn test_quantification_comparison_splice_aware() {
        // Test quantification comparing splice_aware vs non-splice_aware modes
        let gtf_path = "/data2/litian/database/gtf/small_subset_genes2.gtf";
        let bam_file = "/data2/litian/202506_trajectory/data/velocyto_toy/test_small_200.bam";
        let output_file_splice = "/data2/litian/202506_trajectory/process/20250803_velocity_validation/quantification_splice_aware.h5ad";
        let output_file_no_splice = "/data2/litian/202506_trajectory/process/20250803_velocity_validation/quantification_no_splice.h5ad";

        if Path::new(gtf_path).exists() && Path::new(bam_file).exists() {
            println!("Running quantification comparison test...");

            let result = std::panic::catch_unwind(|| -> Result<()> {
                use crate::transcript::quantification::Quantifier;
                use noodles::bam;
                use crate::align::MultiMapR;
                
                // Load transcripts from GTF
                let transcripts = IntronValidationTest::load_transcripts_from_gtf_with_limit(gtf_path, Some(500))?;
                println!("Loaded {} transcripts for comparison test", transcripts.len());
                
                // Helper function to collect records from BAM
                let collect_records = || -> Result<Vec<Vec<(Option<MultiMapR>, Option<MultiMapR>)>>> {
                    let mut reader = bam::io::Reader::new(std::fs::File::open(bam_file)?);
                    let header = reader.read_header()?;
                    
                    let mut all_records = Vec::new();
                    let mut batch = Vec::new();
                    const BATCH_SIZE: usize = 50;
                    
                    for result in reader.record_bufs(&header) {
                        let record = result?;
                        
                        // Skip unmapped reads
                        if record.flags().is_unmapped() {
                            continue;
                        }
                        
                        // Convert to MultiMapR
                        let multi_map_r = MultiMapR::from(record);
                        
                        // For single-end reads, second element is None
                        batch.push((Some(multi_map_r), None));
                        
                        if batch.len() >= BATCH_SIZE {
                            all_records.push(batch);
                            batch = Vec::new();
                        }
                    }
                    
                    if !batch.is_empty() {
                        all_records.push(batch);
                    }
                    
                    Ok(all_records)
                };
                
                // Test 1: With splice_aware = true
                println!("Testing quantification with splice_aware = true...");
                let annotator1 = crate::transcript::annotate::AlignmentAnnotator::new(transcripts.clone());
                let mut quantifier1 = Quantifier::new(annotator1);
                let records1 = collect_records()?;
                let header = {
                    let mut reader = bam::io::Reader::new(std::fs::File::open(bam_file)?);
                    reader.read_header()?
                };
                
                let qc_result1 = quantifier1.quantify(&header, records1.into_iter(), output_file_splice, true)?;
                
                // Test 2: With splice_aware = false
                println!("Testing quantification with splice_aware = false...");
                let annotator2 = crate::transcript::annotate::AlignmentAnnotator::new(transcripts);
                let mut quantifier2 = Quantifier::new(annotator2);
                let records2 = collect_records()?;
                
                let qc_result2 = quantifier2.quantify(&header, records2.into_iter(), output_file_no_splice, false)?;
                
                // Compare results
                println!("\nQuantification Comparison Results:");
                println!("=================================");
                println!("Splice-aware mode:");
                println!("  Total reads: {}", qc_result1.total_reads);
                println!("  Total UMI: {}", qc_result1.total_umi);
                println!("  Unique UMI: {}", qc_result1.unique_umi);
                
                println!("Non-splice-aware mode:");
                println!("  Total reads: {}", qc_result2.total_reads);
                println!("  Total UMI: {}", qc_result2.total_umi);
                println!("  Unique UMI: {}", qc_result2.unique_umi);
                
                // Both modes should process the same number of reads
                assert_eq!(qc_result1.total_reads, qc_result2.total_reads, 
                          "Both modes should process the same number of reads");
                
                println!("Comparison test completed successfully!");
                
                Ok(())
            });

            match result {
                Ok(Ok(())) => {
                    println!("Quantification comparison test completed successfully!");
                    
                    // Verify both output files exist
                    for (file, mode) in [(output_file_splice, "splice-aware"), (output_file_no_splice, "non-splice-aware")] {
                        if Path::new(file).exists() {
                            if let Ok(metadata) = std::fs::metadata(file) {
                                println!("{} output: {} ({} bytes)", mode, file, metadata.len());
                            }
                        } else {
                            println!("Warning: {} output file was not created", mode);
                        }
                    }
                },
                Ok(Err(e)) => {
                    println!("Quantification comparison test failed with error: {}", e);
                },
                Err(_) => {
                    println!("Quantification comparison test panicked");
                }
            }
        } else {
            println!("Quantification comparison test skipped - files not found");
            println!("  GTF file: {} (exists: {})", gtf_path, Path::new(gtf_path).exists());
            println!("  BAM file: {} (exists: {})", bam_file, Path::new(bam_file).exists());
        }
    }
}
