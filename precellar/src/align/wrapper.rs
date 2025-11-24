//! Wrapper for community minimap2 crate to provide Minimap2Aligner interface

use anyhow::Result;
use noodles::sam::{self, alignment::record_buf::RecordBuf};
use noodles::sam::alignment::record::cigar::op::Kind;
use noodles::sam::alignment::record::cigar::Op;
use noodles::fastq;
use std::path::PathBuf;

/// Parse a CIGAR string (e.g., "10M2I5D") into a vector of CIGAR operations
fn parse_cigar_string(cigar_str: &str) -> Result<Vec<Op>> {
    let mut ops = Vec::new();
    let mut current_len = 0u32;

    for ch in cigar_str.chars() {
        if ch.is_ascii_digit() {
            current_len = current_len * 10 + ch.to_digit(10).unwrap();
        } else {
            let kind = match ch {
                'M' => Kind::Match,
                'I' => Kind::Insertion,
                'D' => Kind::Deletion,
                'N' => Kind::Skip,
                'S' => Kind::SoftClip,
                'H' => Kind::HardClip,
                'P' => Kind::Pad,
                '=' => Kind::SequenceMatch,
                'X' => Kind::SequenceMismatch,
                _ => return Err(anyhow::anyhow!("Invalid CIGAR operation: {}", ch)),
            };

            ops.push(Op::new(kind, current_len as usize));
            current_len = 0;
        }
    }

    Ok(ops)
}

/// Options for Minimap2 aligner
#[derive(Clone)]
pub struct Minimap2Opts {
    index_path: PathBuf,
    preset: Option<minimap2::Preset>,
}

impl Minimap2Opts {
    pub fn new(index_path: PathBuf) -> Self {
        Self {
            index_path,
            preset: None,
        }
    }

    pub fn with_preset(mut self, preset: minimap2::Preset) -> Self {
        self.preset = Some(preset);
        self
    }
}

/// Wrapper around community minimap2::Aligner to match the expected interface
pub struct Minimap2Aligner {
    aligner: minimap2::Aligner<minimap2::Built>,
    header: sam::Header,
    opts: Minimap2Opts,
}

impl Minimap2Aligner {
    pub fn new(opts: Minimap2Opts) -> Result<Self> {

        let builder = minimap2::Aligner::builder();
        // Set preset if provided, otherwise use map-ont for long reads
        let builder = match opts.preset.as_ref() {
            Some(preset) => builder.preset(preset.clone()),
            None => builder.map_ont(),
        };

        // Load the index to reach Aligner<Built> type
        let aligner = builder
            .with_cigar()
            .with_index(&opts.index_path, None)
            .map_err(|e| anyhow::anyhow!("Failed to build minimap2 index: {}", e))?;

        // Build SAM header from minimap2 index
        let header = Self::build_header(&aligner);

        Ok(Self {
            aligner, // ready to align reads
            header,
            opts,
        })
    }

    pub fn get_header(&self) -> &sam::Header {
        &self.header
    }

    pub fn get_opts(&self) -> &Minimap2Opts {
        &self.opts
    }

    /// Align a single read and return alignments
    pub fn align_read(&mut self, record: &fastq::Record) -> Result<Vec<RecordBuf>> {
        let seq = record.sequence();
        let name = record.name();

        // Map the sequence using the minimap2 API
        // Parameters: seq, cs (long cs tag), md (MD tag), max_frag_len, extra_flags, query_name
        let mappings = self.aligner.map(
            seq,
            false,  // cs - long cs tag
            false,  // md - MD tag
            None,   // max_frag_len
            None,   // extra_flags
            Some(name),
        ).map_err(|e| anyhow::anyhow!("Minimap2 mapping failed: {}", e))?;

        // Convert minimap2 mappings to RecordBuf
        let mut results = Vec::new();

        for mapping in mappings {
            let mut record_buf = RecordBuf::default();

            // Set name (Option<BString>)
            *record_buf.name_mut() = Some(name.to_vec().into());
            // Set sequence (using Vec<u8> then into())
            *record_buf.sequence_mut() = seq.to_vec().into();
            // Set quality scores (using Vec<u8> then into())
            *record_buf.quality_scores_mut() = record.quality_scores().to_vec().into();

            // Set mapping information
            if let Some(target_name) = &mapping.target_name {
                let ref_id = self.header
                    .reference_sequences()
                    .get_index_of(target_name.as_bytes())
                    .map(|i| i);

                if let Some(id) = ref_id {
                    *record_buf.reference_sequence_id_mut() = Some(id);
                }
            }

            // Set alignment start position (1-based)
            if let Ok(pos) = noodles::core::Position::try_from(mapping.target_start as usize + 1) {
                *record_buf.alignment_start_mut() = Some(pos);
            }

            // Set mapping quality (mapq is u32, but needs to be u8 for MappingQuality)
            let mapq_u8 = mapping.mapq.min(255) as u8;
            if let Ok(mq) = noodles::sam::alignment::record::MappingQuality::try_from(mapq_u8) {
                *record_buf.mapping_quality_mut() = Some(mq);
            }

            // Set flags
            let mut flags = noodles::sam::alignment::record::Flags::empty();
            if mapping.strand == minimap2::Strand::Reverse {
                flags |= noodles::sam::alignment::record::Flags::REVERSE_COMPLEMENTED;
            }
            if !mapping.is_primary {
                flags |= noodles::sam::alignment::record::Flags::SECONDARY;
            }
            if mapping.is_supplementary {
                flags |= noodles::sam::alignment::record::Flags::SUPPLEMENTARY;
            }
            *record_buf.flags_mut() = flags;

            // Add alignment information if available
            if let Some(ref alignment) = mapping.alignment {
                // Parse and set CIGAR string
                if let Some(cigar_string) = &alignment.cigar_str {
                    match parse_cigar_string(cigar_string) {
                        Ok(ops) => {
                            *record_buf.cigar_mut() = ops.into_iter().collect();
                        }
                        Err(e) => {
                            eprintln!("Failed to parse CIGAR string '{}': {}", cigar_string, e);
                        }
                    }
                }

                // Add alignment score if available
                if let Some(score) = alignment.alignment_score {
                    record_buf.data_mut().insert(
                        noodles::sam::alignment::record::data::field::tag::Tag::ALIGNMENT_SCORE,
                        noodles::sam::alignment::record_buf::data::field::value::Value::Int32(score),
                    );
                }
            }

            results.push(record_buf);
        }

        Ok(results)
    }

    fn build_header(aligner: &minimap2::Aligner<minimap2::Built>) -> sam::Header {
        use std::ffi::CStr;
        use std::num::NonZeroUsize;
        use noodles::sam::header::record::value::{map::ReferenceSequence, Map};

        let mut header = sam::Header::default();

        // Add reference sequences from index using n_seq() and get_seq()
        let n_sequences = aligner.n_seq();
        for i in 0..n_sequences {
            if let Some(seq_info) = aligner.get_seq(i as usize) {
                // Get the sequence name from the C string
                let c_str = unsafe { CStr::from_ptr(seq_info.name) };
                let name = c_str.to_str().unwrap_or("unknown");

                // Get the sequence length
                let length = seq_info.len as usize;

                // Create the reference sequence entry
                if let Ok(length_nz) = NonZeroUsize::try_from(length) {
                    let ref_seq = Map::<ReferenceSequence>::new(length_nz);
                    header.reference_sequences_mut()
                        .insert(name.as_bytes().to_vec().into(), ref_seq);
                }
            }
        }

        header
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs::File;
    use std::io::BufReader;
    use noodles::{fastq, sam};
    use noodles::sam::alignment::io::Write;
    use noodles::sam::alignment::record::Cigar as CigarTrait;

    #[test]
    #[ignore] 
    fn test_alignment() -> Result<()> {
        let ref_path = "/data/Public/genome/GRCh38/GRCh38.primary_assembly.genome.fa.gz";
        let reads_path = "/data/xurui/project/minimap2-rs/data/h38_dna/SRR17666426_20.fastq";
        let output_sam_path = "/data/xurui/project/minimap2-rs/data/h38_dna/SRR17666426_20_rs.sam";

        // Set expected results from minimap2 (RNAME, POS, FLAG)
        let expected_results = vec![
            ("chr13", 48234884, 0),
            ("chr3", 5862304, 16),
            ("chr12", 21126545, 0),
            ("chr7", 24828037, 0),
            ("chr1", 119770793, 16),
        ];

        // Initialize aligner
        let opts = Minimap2Opts::new(PathBuf::from(ref_path))
            .with_preset(minimap2::Preset::MapOnt);
        let mut aligner = Minimap2Aligner::new(opts)?;
        let header = aligner.get_header().clone(); // Clone header for later validation

        // Read FASTQ file
        let mut reader = fastq::io::Reader::new(BufReader::new(File::open(&reads_path)?));
        let mut records = Vec::new();
        for result in reader.records() {
            match result {
                Ok(record) => records.push(record),
                Err(e) => {
                    eprintln!("Warning: Skipping invalid FASTQ record: {}", e);
                    continue;
                }
            }
        }

        println!("Read {} FASTQ records", records.len());

        // Align reads and collect alignments
        let mut all_alignments = Vec::new();
        for record in &records {
            let alignments = aligner.align_read(record)?;
            all_alignments.extend(alignments);
        }

        println!("Generated {} alignments", all_alignments.len());

        // Write alignments to output SAM file
        let output_file = File::create(&output_sam_path)?;
        let mut writer = sam::io::Writer::new(output_file);
        writer.write_header(&header)?;

        for alignment in &all_alignments {
            writer.write_alignment_record(&header, alignment)?;
        }

        println!("Wrote alignments to {}", output_sam_path);

        // Validate first 5 primary alignments against expected results
        let primary_alignments: Vec<_> = all_alignments.iter()
            .filter(|a| !a.flags().is_secondary() && !a.flags().is_supplementary())
            .take(5)
            .collect();

        assert_eq!(primary_alignments.len(), 5, "Expected at least 5 primary alignments");

        for (i, (expected_rname, expected_pos, expected_flag)) in expected_results.iter().enumerate() {
            let alignment = primary_alignments[i];

            // Validate RNAME (reference sequence name)
            if let Some(ref_id) = alignment.reference_sequence_id() {
                let ref_name = header.reference_sequences()
                    .get_index(ref_id)
                    .map(|(name, _)| std::str::from_utf8(name.as_ref()).unwrap());
                assert_eq!(ref_name, Some(*expected_rname),
                    "Alignment {} RNAME mismatch", i);
            } else {
                panic!("Alignment {} missing reference sequence", i);
            }

            // Validate POS (alignment position)
            if let Some(pos) = alignment.alignment_start() {
                assert_eq!(usize::from(pos), *expected_pos,
                    "Alignment {} POS mismatch", i);
            } else {
                panic!("Alignment {} missing alignment start position", i);
            }

            // Validate FLAG
            let flag_value = alignment.flags().bits();
            assert_eq!(flag_value, *expected_flag,
                "Alignment {} FLAG mismatch", i);

            // Validate CIGAR is present and non-empty
            assert!(!alignment.cigar().is_empty(),
                "Alignment {} missing CIGAR string", i);

            println!("Alignment {}: {} @ {} FLAG={} CIGAR={:?} - PASS",
                i, expected_rname, expected_pos, expected_flag, alignment.cigar());
        }

        Ok(())
    }
}
