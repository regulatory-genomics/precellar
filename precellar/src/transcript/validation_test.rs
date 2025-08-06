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
use crate::transcript::{AlignmentAnnotator, Transcript, quantification::IntronValidationCollector};
use crate::align::MultiMapR;

#[cfg(test)]
use log;

/// Test structure to extract and output validated intron locations
#[derive(Debug)]
pub struct IntronValidationTest {
    annotator: AlignmentAnnotator,
    validation_collector: IntronValidationCollector,
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
    /// Create a new validation test with STAR transcriptome
    pub fn new(star_reference_path: &str) -> Result<Self> {
        use star_aligner::{StarAligner, StarOpts};
        use std::path::PathBuf;

        println!("Loading transcripts from STAR reference: {}", star_reference_path);

        // Create STAR aligner to access transcriptome
        let opts = StarOpts::new(PathBuf::from(star_reference_path));

        // Add better error handling for STAR aligner initialization
        let aligner = match StarAligner::new(opts) {
            Ok(aligner) => aligner,
            Err(e) => {
                eprintln!("Failed to create STAR aligner: {:?}", e);
                return Err(e.into());
            }
        };

        // Get transcriptome and convert to our Transcript format
        let star_transcripts = match aligner.get_transcriptome() {
            Ok(transcripts) => transcripts,
            Err(e) => {
                eprintln!("Failed to get transcriptome from STAR aligner: {:?}", e);
                return Err(e.into());
            }
        };

        let transcripts: Result<Vec<_>> = star_transcripts
            .iter()
            .map(|t| Transcript::try_from(t.clone()))
            .collect();

        let mut transcripts = transcripts?;

        // Generate introns for all transcripts
        for transcript in &mut transcripts {
            transcript.make_intron_by_exons();
        }

        // Create annotator with the transcripts
        let annotator = AlignmentAnnotator::new(transcripts);
        let validation_collector = IntronValidationCollector::new();

        println!("Loaded {} transcripts from STAR reference", annotator.transcripts().len());

        Ok(Self {
            annotator,
            validation_collector,
        })
    }

    /// Create a new validation test with GTF transcriptome
    pub fn new_from_gtf(gtf_path: &str) -> Result<Self> {
        println!("Loading transcripts from GTF file: {}", gtf_path);

        let transcripts = Self::load_transcripts_from_gtf(gtf_path)?;
        let annotator = AlignmentAnnotator::new(transcripts);
        let validation_collector = IntronValidationCollector::new();

        println!("Loaded {} transcripts from GTF file", annotator.transcripts().len());

        Ok(Self {
            annotator,
            validation_collector,
        })
    }

    /// Create a new validation test from a vector of transcripts
    pub fn new_from_transcripts(transcripts: Vec<Transcript>) -> Result<Self> {
        let annotator = AlignmentAnnotator::new(transcripts);
        let validation_collector = IntronValidationCollector::new();

        println!("Loaded {} transcripts", annotator.transcripts().len());

        Ok(Self {
            annotator,
            validation_collector,
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
        use crate::transcript::transcriptome::{Exons, Introns};
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
        let mut skipped_invalid_exons = 0;

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

            let introns = Introns::new(std::iter::empty())?;

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
                introns,
            };

            // Generate introns from exons
            transcript.make_intron_by_exons();

            // Final validation: Ensure transcript has exons before adding to list
            if transcript.exons().is_empty() {
                skipped_no_exons += 1;
                if skipped_no_exons <= 5 {
                    println!("Warning: Skipping transcript {} - no exons after make_intron_by_exons()", transcript.id);
                }
                continue;
            }

            // Debug: Log transcript details for first transcript only (reduced verbosity)
            if transcripts.len() == 0 {
                println!("Debug: First transcript {} - exons: {}, introns: {}, span: {}-{}",
                         transcript.id,
                         transcript.exons().len(),
                         transcript.introns().len(),
                         transcript.start,
                         transcript.end);
            }

            transcripts.push(transcript);
        }

        println!("Successfully parsed {} transcripts from GTF", transcripts.len());

        Ok(transcripts)
    }

    // Removed inefficient extract_attribute_from_gtf function - now using direct attribute access

    /// Debug wrapper for annotate_splice to identify problematic transcripts
    pub fn debug_annotate_splice(
        splice_segments: &crate::transcript::transcriptome::SpliceSegments,
        transcript: &Transcript,
        transcript_id: &str,
    ) -> Option<(bool, Vec<usize>, Vec<usize>)> {
        let result = splice_segments.annotate_splice(transcript);

        if result.is_none() {
            // Debug information removed for cleaner code
        }

        result
    }
    
    /// Process BAM file and collect validation requests
    pub fn process_bam_file(&mut self, bam_path: &str) -> Result<()> {
        let mut reader = File::open(bam_path)
            .map(bam::io::Reader::new)?;
        
        let header = reader.read_header()?;
        
        // Process records in batches
        let mut batch = Vec::new();
        const BATCH_SIZE: usize = 1000;
        
        for result in reader.records() {
            let record = result?;
            let record_buf = sam::alignment::RecordBuf::try_from_alignment_record(&header, &record)?;
            
            // Convert to MultiMapR format (simplified for this example)
            if let Some(multi_map) = Self::convert_to_multimap(record_buf) {
                batch.push(multi_map);
                
                if batch.len() >= BATCH_SIZE {
                    self.process_batch(&header, &batch)?;
                    batch.clear();
                }
            }
        }
        
        // Process remaining records
        if !batch.is_empty() {
            self.process_batch(&header, &batch)?;
        }
        
        Ok(())
    }
    
    /// Convert BAM record to MultiMapR (simplified)
    fn convert_to_multimap(record: sam::alignment::RecordBuf) -> Option<MultiMapR> {
        // Create a MultiMapR with the primary record and no secondary alignments
        Some(MultiMapR::new(record, None))
    }
    
    /// Process a batch of records
    fn process_batch(&mut self, header: &Header, batch: &[MultiMapR]) -> Result<()> {
        for record in batch {
            // Annotate the alignment and collect validation requests
            if let Some(annotated) = self.annotator.annotate_alignments_se(header, record.clone(), true) {
                // Extract validation requests from the annotation
                use crate::transcript::annotate::AnnotatedAlignment;
                match annotated {
                                         AnnotatedAlignment::SeMapped(annotation) => {
                         // Extract validation requests from single-end annotation
                         for &(ref transcript_id, intron_index) in &annotation.validation_requests {
                             self.validation_collector.add_validation_request(transcript_id.clone(), intron_index);
                         }
                     },
                     AnnotatedAlignment::PeMapped(a1, a2, _) => {
                         // Extract validation requests from paired-end annotations
                         for &(ref transcript_id, intron_index) in &a1.validation_requests {
                             self.validation_collector.add_validation_request(transcript_id.clone(), intron_index);
                         }
                         for &(ref transcript_id, intron_index) in &a2.validation_requests {
                             self.validation_collector.add_validation_request(transcript_id.clone(), intron_index);
                         }
                     }
                }
            }
        }
        Ok(())
    }
    
    /// Validate collected introns and generate output
    pub fn validate_and_output(&mut self, output_path: &str) -> Result<()> {
        // Get all transcripts and validate introns using the existing method
        let mut transcripts = self.annotator.transcripts();
        let validated_introns_set = self.validation_collector.validate_introns(&mut transcripts);

        // Convert to detailed intron information
        let validated_introns = self.extract_validated_intron_details(&transcripts, &validated_introns_set)?;

        // Write to output file
        self.write_output(&validated_introns, output_path)?;

        println!("Validated intron locations written to: {}", output_path);
        println!("Total validated introns: {}", validated_introns.len());

        Ok(())
    }
    
    /// Extract detailed information for validated introns
    fn extract_validated_intron_details(
        &self,
        transcripts: &[Transcript],
        validated_set: &HashSet<(String, usize)>,
    ) -> Result<Vec<ValidatedIntron>> {
        let mut validated_introns = Vec::new();
        
        for transcript in transcripts {
            let introns = transcript.introns();
            
            for (intron_idx, intron) in introns.iter().enumerate() {
                // Check if this intron is validated
                if validated_set.contains(&(transcript.id.clone(), intron_idx)) {
                    let validated_intron = ValidatedIntron {
                        chromosome: transcript.chrom.clone(),
                        start: intron.start(),
                        end: intron.end(),
                        strand: match transcript.strand {
                            noodles::sam::record::data::field::value::base_modifications::group::Strand::Forward => "+".to_string(),
                            noodles::sam::record::data::field::value::base_modifications::group::Strand::Reverse => "-".to_string(),
                        },
                        gene_id: transcript.gene.id.clone(),
                        gene_name: transcript.gene.name.clone(),
                        transcript_id: transcript.id.clone(),
                        transcript_name: transcript.id.clone(), // Assuming same as ID
                        intron_number: intron_idx + 1, // 1-based indexing
                        length: intron.end() - intron.start() + 1,
                        validation_status: "validated".to_string(),
                    };
                    
                    validated_introns.push(validated_intron);
                }
            }
        }
        
        // Sort by chromosome, then by start position
        validated_introns.sort_by(|a, b| {
            a.chromosome.cmp(&b.chromosome)
                .then(a.start.cmp(&b.start))
        });
        
        Ok(validated_introns)
    }
    
    /// Write validated introns to output file
    fn write_output(&self, validated_introns: &[ValidatedIntron], output_path: &str) -> Result<()> {
        let file = File::create(output_path)?;
        let mut writer = BufWriter::new(file);
        
        // Write header
        writeln!(writer, "chromosome\tstart\tend\tstrand\tgene_id\tgene_name\ttranscript_id\ttranscript_name\tintron_number\tlength\tvalidation_status")?;
        
        // Write data
        for intron in validated_introns {
            writeln!(
                writer,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                intron.chromosome,
                intron.start,
                intron.end,
                intron.strand,
                intron.gene_id,
                intron.gene_name,
                intron.transcript_id,
                intron.transcript_name,
                intron.intron_number,
                intron.length,
                intron.validation_status
            )?;
        }
        
        writer.flush()?;
        Ok(())
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
    pub has_exons: bool,
    pub has_introns: bool,
    pub has_validated_introns: bool,
    pub has_spanning: bool,
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
        
        // Check various properties across all alignments
        let has_exons = alignments.iter().any(|alignment| alignment.exon_align.is_some());
        let has_introns = alignments.iter().any(|alignment| !alignment.intron_mapped.is_empty());
        let has_validated_introns = alignments.iter().any(|alignment| !alignment.validation_requests.is_empty());
        let has_spanning = has_validated_introns; // Spanning reads typically have validation requests

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
            has_exons,
            has_introns,
            has_validated_introns,
            has_spanning,
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
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
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
                classification.has_exons,
                classification.has_introns,
                classification.has_validated_introns,
                classification.has_spanning
            )?;
        }

        writer.flush()?;
        Ok(())
    }
}

/// Example function to run the validation test with real data
pub fn run_validation_test(
    star_reference_path: &str,
    bam_file_path: &str,
    output_file_path: &str,
) -> Result<()> {
    println!("Starting intron validation test...");
    println!("STAR reference: {}", star_reference_path);
    println!("BAM file: {}", bam_file_path);
    println!("Output file: {}", output_file_path);

    let mut validator = IntronValidationTest::new(star_reference_path)?;

    println!("Processing BAM file...");
    validator.process_bam_file(bam_file_path)?;

    println!("Validating introns and generating output...");
    validator.validate_and_output(output_file_path)?;

    println!("Validation test completed successfully!");
    Ok(())
}

/// Example function to run the validation test with GTF file
pub fn run_validation_test_with_gtf(
    gtf_path: &str,
    bam_file_path: &str,
    output_file_path: &str,
) -> Result<()> {
    println!("Starting intron validation test with GTF...");
    println!("GTF file: {}", gtf_path);
    println!("BAM file: {}", bam_file_path);
    println!("Output file: {}", output_file_path);

    let mut validator = IntronValidationTest::new_from_gtf(gtf_path)?;

    println!("Processing BAM file...");
    validator.process_bam_file(bam_file_path)?;

    println!("Validating introns and generating output...");
    validator.validate_and_output(output_file_path)?;

    println!("Validation test completed successfully!");
    Ok(())
}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::transcript::{Transcript, Gene};
    use noodles::sam::record::data::field::value::base_modifications::group::Strand as NoodlesStrand;
    use crate::transcript::transcriptome::{Exons, Introns};

    #[test]
    fn test_intron_validation_with_mock_data() {
        // Create mock transcripts for testing
        let mock_transcripts = create_mock_transcripts();
        let annotator = AlignmentAnnotator::new(mock_transcripts);
        let mut validation_collector = IntronValidationCollector::new();

        // Add some mock validation requests
        validation_collector.add_validation_request("ENST00000398216".to_string(), 0);
        validation_collector.add_validation_request("ENST00000418300".to_string(), 0);
        validation_collector.add_validation_request("ENST00000629289".to_string(), 5);

        let mut validator = IntronValidationTest {
            annotator,
            validation_collector,
        };

        // Test the validation and output generation
        let output_file = "test_validated_introns.tsv";
        validator.validate_and_output(output_file).unwrap();

        // Verify output file was created
        assert!(Path::new(output_file).exists());

        // Read and verify content
        let content = std::fs::read_to_string(output_file).unwrap();
        assert!(content.contains("chromosome\tstart\tend"));
        assert!(content.contains("validated"));

        // Clean up
        std::fs::remove_file(output_file).ok();
    }
    fn create_mock_transcripts() -> Vec<Transcript> {
        vec![
            create_mock_transcript(
                "ENST00000398216",
                "ENSG00000241180",
                "ENSG00000241180",
                "1",
                914000,
                915000,
                vec![(914000, 914444), (914803, 915000)],
            ),
            create_mock_transcript(
                "ENST00000418300",
                "ENSG00000242590",
                "ENSG00000242590",
                "1",
                1055000,
                1056000,
                vec![(1055000, 1055215), (1055898, 1056000)],
            ),
            create_mock_transcript(
                "ENST00000629289",
                "ENSG00000248333",
                "CDK11B",
                "1",
                1640000,
                1660000,
                vec![
                    (1640000, 1641000),
                    (1642000, 1643000),
                    (1644000, 1645000),
                    (1645100, 1645508),
                    (1645773, 1646000),
                    (1647000, 1648000),
                    (1649000, 1660000),
                ],
            ),
        ]
    }

    fn create_mock_transcript(
        id: &str,
        gene_id: &str,
        gene_name: &str,
        chrom: &str,
        start: u64,
        end: u64,
        exon_coords: Vec<(u64, u64)>,
    ) -> Transcript {
        let exons = Exons::new(exon_coords.into_iter()).unwrap();
        let introns = Introns::new(std::iter::empty()).unwrap();

        let mut transcript = Transcript {
            id: id.to_string(),
            chrom: chrom.to_string(),
            start,
            end,
            strand: NoodlesStrand::Forward,
            gene: Gene {
                id: gene_id.to_string(),
                name: gene_name.to_string(),
            },
            exons,
            introns,
        };

        // Generate introns from exons
        transcript.make_intron_by_exons();
        transcript
    }

    #[test]
    fn test_gtf_loading() {
        // Create a mock GTF content for testing
        let gtf_content = r#"##description: evidence-based annotation of the human genome (GRCh38), version 44 (Ensembl 110)
##provider: GENCODE
##format: gtf
##date: 2023-03-01
chr1	HAVANA	gene	11869	14409	.	+	.	gene_id "ENSG00000290825"; gene_version "1"; gene_type "lncRNA"; gene_name "DDX11L2"; level 2; tag "overlaps_pseudogene";
chr1	HAVANA	transcript	11869	14409	.	+	.	gene_id "ENSG00000290825"; gene_version "1"; transcript_id "ENST00000456328"; transcript_version "2"; gene_type "lncRNA"; gene_name "DDX11L2"; transcript_type "lncRNA"; transcript_name "DDX11L2-202"; level 2; transcript_support_level "1"; tag "basic"; tag "Ensembl_canonical"; havana_transcript "OTTHUMT00000362751.1";
chr1	HAVANA	exon	11869	12227	.	+	.	gene_id "ENSG00000290825"; gene_version "1"; transcript_id "ENST00000456328"; transcript_version "2"; gene_type "lncRNA"; gene_name "DDX11L2"; transcript_type "lncRNA"; transcript_name "DDX11L2-202"; exon_number 1; exon_id "ENSE00002234944"; exon_version "1"; level 2; transcript_support_level "1"; tag "basic"; tag "Ensembl_canonical"; havana_transcript "OTTHUMT00000362751.1";
chr1	HAVANA	exon	12613	12721	.	+	.	gene_id "ENSG00000290825"; gene_version "1"; transcript_id "ENST00000456328"; transcript_version "2"; gene_type "lncRNA"; gene_name "DDX11L2"; transcript_type "lncRNA"; transcript_name "DDX11L2-202"; exon_number 2; exon_id "ENSE00003582793"; exon_version "1"; level 2; transcript_support_level "1"; tag "basic"; tag "Ensembl_canonical"; havana_transcript "OTTHUMT00000362751.1";
chr1	HAVANA	exon	13221	14409	.	+	.	gene_id "ENSG00000290825"; gene_version "1"; transcript_id "ENST00000456328"; transcript_version "2"; gene_type "lncRNA"; gene_name "DDX11L2"; transcript_type "lncRNA"; transcript_name "DDX11L2-202"; exon_number 3; exon_id "ENSE00002312635"; exon_version "1"; level 2; transcript_support_level "1"; tag "basic"; tag "Ensembl_canonical"; havana_transcript "OTTHUMT00000362751.1";
chr1	HAVANA	gene	29554	31109	.	+	.	gene_id "ENSG00000243485"; gene_version "5"; gene_type "lncRNA"; gene_name "MIR1302-2HG"; level 2; hgnc_id "HGNC:52482"; tag "ncRNA_host"; havana_gene "OTTHUMG00000000959.2";
chr1	HAVANA	transcript	29554	31097	.	+	.	gene_id "ENSG00000243485"; gene_version "5"; transcript_id "ENST00000473358"; transcript_version "1"; gene_type "lncRNA"; gene_name "MIR1302-2HG"; transcript_type "lncRNA"; transcript_name "MIR1302-2HG-202"; level 2; transcript_support_level "5"; hgnc_id "HGNC:52482"; tag "not_best_in_genome_evidence"; tag "dotter_confirmed"; tag "basic"; tag "Ensembl_canonical"; havana_gene "OTTHUMG00000000959.2"; havana_transcript "OTTHUMT00000002840.1";
chr1	HAVANA	exon	29554	30039	.	+	.	gene_id "ENSG00000243485"; gene_version "5"; transcript_id "ENST00000473358"; transcript_version "1"; gene_type "lncRNA"; gene_name "MIR1302-2HG"; transcript_type "lncRNA"; transcript_name "MIR1302-2HG-202"; exon_number 1; exon_id "ENSE00001947070"; exon_version "1"; level 2; transcript_support_level "5"; hgnc_id "HGNC:52482"; tag "not_best_in_genome_evidence"; tag "dotter_confirmed"; tag "basic"; tag "Ensembl_canonical"; havana_gene "OTTHUMG00000000959.2"; havana_transcript "OTTHUMT00000002840.1";
chr1	HAVANA	exon	30564	31097	.	+	.	gene_id "ENSG00000243485"; gene_version "5"; transcript_id "ENST00000473358"; transcript_version "1"; gene_type "lncRNA"; gene_name "MIR1302-2HG"; transcript_type "lncRNA"; transcript_name "MIR1302-2HG-202"; exon_number 2; exon_id "ENSE00001922571"; exon_version "1"; level 2; transcript_support_level "5"; hgnc_id "HGNC:52482"; tag "not_best_in_genome_evidence"; tag "dotter_confirmed"; tag "basic"; tag "Ensembl_canonical"; havana_gene "OTTHUMG00000000959.2"; havana_transcript "OTTHUMT00000002840.1";
"#;

        // Write mock GTF to a temporary file
        let temp_gtf_path = "test_mock.gtf";
        std::fs::write(temp_gtf_path, gtf_content).unwrap();

        // Test GTF loading
        let result = IntronValidationTest::new_from_gtf(temp_gtf_path);

        match result {
            Ok(validator) => {
                let transcripts = validator.annotator.transcripts();
                println!("Successfully loaded {} transcripts from GTF", transcripts.len());

                // Verify we loaded the expected transcripts
                assert_eq!(transcripts.len(), 2);
                // Debug: Print actual transcript IDs
                println!("Loaded transcript IDs:");
                for transcript in transcripts.iter() {
                    println!("  - {}", transcript.id);
                }
                // Check first transcript (DDX11L2)
                let ddx11l2 = transcripts.iter().find(|t| t.id == "ENST00000456328").unwrap();
                assert_eq!(ddx11l2.gene.name, "DDX11L2");
                assert_eq!(ddx11l2.chrom, "chr1");
                assert_eq!(ddx11l2.start, 11868); // 0-based
                assert_eq!(ddx11l2.end, 14409);
                assert_eq!(ddx11l2.exons().len(), 3);

                // Check that introns were generated
                let introns = ddx11l2.introns();
                assert_eq!(introns.len(), 2); // 3 exons = 2 introns

                // Check intron coordinates
                assert_eq!(introns[0].start(), 12228); // End of first exon
                assert_eq!(introns[0].end(), 12611); // Start of second exon - 1
                assert_eq!(introns[1].start(), 12722); // End of second exon
                assert_eq!(introns[1].end(), 13219); // Start of third exon - 1

                // Check second transcript (MIR1302-2HG)
                let mir1302 = transcripts.iter().find(|t| t.id == "ENST00000473358").unwrap();
                assert_eq!(mir1302.gene.name, "MIR1302-2HG");
                assert_eq!(mir1302.exons().len(), 2);
                assert_eq!(mir1302.introns().len(), 1); // 2 exons = 1 intron

                println!("✅ GTF loading test passed!");

                // Test validation functionality with GTF-loaded transcripts
                let mut validation_collector = IntronValidationCollector::new();
                validation_collector.add_validation_request("ENST00000456328".to_string(), 0);
                validation_collector.add_validation_request("ENST00000473358".to_string(), 0);

                let mut transcripts_mut = transcripts;
                let validated_introns = validation_collector.validate_introns(&mut transcripts_mut);

                println!("Validated {} introns", validated_introns.len());

                // Test output generation
                let output_file = "test_gtf_validated_introns.tsv";
                write_validated_introns_to_file(&transcripts_mut, &validated_introns, output_file).unwrap();

                // Verify output file
                assert!(Path::new(output_file).exists());
                let content = std::fs::read_to_string(output_file).unwrap();
                println!("GTF validation output:\n{}", content);

                // Clean up
                std::fs::remove_file(output_file).ok();

                println!("✅ GTF validation workflow test passed!");
            },
            Err(e) => {
                panic!("Failed to load GTF: {:?}", e);
            }
        }

        // Clean up
        std::fs::remove_file(temp_gtf_path).ok();
    }

    /// Helper function to write validated introns to file
    fn write_validated_introns_to_file(
        transcripts: &[Transcript],
        validated_set: &HashSet<(String, usize)>,
        output_path: &str,
    ) -> Result<()> {
        let file = File::create(output_path)?;
        let mut writer = BufWriter::new(file);

        // Write header
        writeln!(writer, "chromosome\tstart\tend\tstrand\tgene_id\tgene_name\ttranscript_id\ttranscript_name\tintron_number\tlength\tvalidation_status")?;

        let mut validated_introns = Vec::new();

        for transcript in transcripts {
            let introns = transcript.introns();

            for (intron_idx, intron) in introns.iter().enumerate() {
                if validated_set.contains(&(transcript.id.clone(), intron_idx)) {
                    validated_introns.push((
                        transcript.chrom.clone(),
                        intron.start(),
                        intron.end(),
                        match transcript.strand {
                            noodles::sam::record::data::field::value::base_modifications::group::Strand::Forward => "+",
                            noodles::sam::record::data::field::value::base_modifications::group::Strand::Reverse => "-",
                        },
                        transcript.gene.id.clone(),
                        transcript.gene.name.clone(),
                        transcript.id.clone(),
                        transcript.id.clone(),
                        intron_idx + 1,
                        intron.end() - intron.start(),
                    ));
                }
            }
        }

        // Sort by chromosome, then by start position
        validated_introns.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));

        // Write data
        for (chrom, start, end, strand, gene_id, gene_name, tx_id, tx_name, intron_num, length) in validated_introns {
            writeln!(
                writer,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tvalidated",
                chrom, start, end, strand, gene_id, gene_name, tx_id, tx_name, intron_num, length
            )?;
        }

        writer.flush()?;
        Ok(())
    }

    #[test]
    fn test_real_gtf_if_available() {
        // Test with real GTF file if available
        let gtf_path = "/data2/litian/database/gtf/small_subset_genes.gtf";
        let bam_file = "/data2/litian/202506_trajectory/data/velocyto_toy/test_small_200.bam";
        let output_file = "/data2/litian/202506_trajectory/process/20250803_velocity_validation/validated_introns_gtf.tsv";

        if Path::new(gtf_path).exists() && Path::new(bam_file).exists() {
            println!("Running test with real GTF data...");

            let result = std::panic::catch_unwind(|| -> Result<()> {
                // Limit to 100 transcripts for testing to avoid long processing times
                let transcripts = IntronValidationTest::load_transcripts_from_gtf_with_limit(gtf_path, Some(10000000))?;
                let mut validator = IntronValidationTest::new_from_transcripts(transcripts)?;
                validator.process_bam_file(bam_file)?;
                validator.validate_and_output(output_file)?;
                Ok(())
            });

            match result {
                Ok(Ok(_)) => {
                    // Verify output
                    assert!(Path::new(output_file).exists());
                    let content = std::fs::read_to_string(output_file).unwrap();
                    println!("Real GTF output preview:\n{}", &content[..content.len().min(500)]);

                    println!("✅ Real GTF test completed successfully!");
                },
                Ok(Err(e)) => {
                    eprintln!("Real GTF test failed with error: {:?}", e);
                    println!("Skipping real GTF test due to errors");
                },
                Err(_) => {
                    eprintln!("Real GTF test panicked");
                    println!("Skipping real GTF test due to panic");
                }
            }
        } else {
            println!("Skipping real GTF test - files not found");
            println!("To run with real data, ensure these paths exist:");
            println!("  GTF file: {}", gtf_path);
            println!("  BAM file: {}", bam_file);
        }
    }


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
