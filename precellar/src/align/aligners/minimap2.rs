//! Wrapper for community minimap2 crate to provide Minimap2Aligner interface

use anyhow::Result;
use noodles::sam::{self, alignment::record_buf::RecordBuf};
use noodles::sam::alignment::record::cigar::op::Kind;
use noodles::sam::alignment::record::cigar::Op;
use noodles::fastq;
use std::path::PathBuf;

/// Static lookup table for nucleotide complement
/// This is computed at compile time and provides O(1) branch-free complement lookups
static COMPLEMENT_TABLE: [u8; 256] = {
    let mut table = [0u8; 256];
    
    // Initialize all entries to themselves (for non-ACGT characters)
    let mut i = 0;
    while i < 256 {
        table[i] = i as u8;
        i += 1;
    }
    
    // Set complementary bases (uppercase)
    table[b'A' as usize] = b'T';
    table[b'T' as usize] = b'A';
    table[b'C' as usize] = b'G';
    table[b'G' as usize] = b'C';
    
    // Set complementary bases (lowercase)
    table[b'a' as usize] = b't';
    table[b't' as usize] = b'a';
    table[b'c' as usize] = b'g';
    table[b'g' as usize] = b'c';
    
    // N and other characters stay as themselves
    
    table
};

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

    pub fn preset(&self) -> Option<&minimap2::Preset> {
        self.preset.as_ref()
    }
}

/// Wrapper around community minimap2::Aligner to match the expected interface
#[derive(Clone)]
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
        let qual = record.quality_scores();

        // Map the sequence using the minimap2 API
        // Parameters: seq, cs (long cs tag), md (MD tag), max_frag_len, extra_flags, query_name (default values)
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

            // Set flags
            let mut flags = noodles::sam::alignment::record::Flags::empty();
            let is_reverse = mapping.strand == minimap2::Strand::Reverse;

            if is_reverse {
                flags |= noodles::sam::alignment::record::Flags::REVERSE_COMPLEMENTED;
            }
            if !mapping.is_primary {
                flags |= noodles::sam::alignment::record::Flags::SECONDARY;
            }
            if mapping.is_supplementary {
                flags |= noodles::sam::alignment::record::Flags::SUPPLEMENTARY;
            }
            *record_buf.flags_mut() = flags;

            // Set name (Option<BString>)
            *record_buf.name_mut() = Some(name.to_vec().into());

            // For reverse-complemented alignments, sequence and quality scores of SAM record should be reversed as well
            if is_reverse {
                // Use lookup table for fast, branch-free complement
                let rc_seq: Vec<u8> = seq.iter()
                    .rev()
                    .map(|&b| COMPLEMENT_TABLE[b as usize])
                    .collect();
                *record_buf.sequence_mut() = rc_seq.into();

                let rev_qual: Vec<u8> = qual.iter()
                    .rev()
                    .map(|&b| b.saturating_sub(33)) // NOTE: Decode FASTQ ASCII (Phred+33) to raw Phred score required by RecordBuf
                    .collect();
                *record_buf.quality_scores_mut() = rev_qual.into();
            } else {
                *record_buf.sequence_mut() = seq.to_vec().into();
                
                let decoded_qual: Vec<u8> = qual.iter()
                    .map(|&b| b.saturating_sub(33)) // NOTE: Decode FASTQ ASCII (Phred+33) to raw Phred score required by RecordBuf
                    .collect();
                *record_buf.quality_scores_mut() = decoded_qual.into();
            }   

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

            // Add alignment information if available
            if let Some(ref alignment) = mapping.alignment {
                // Set CIGAR from pre-parsed minimap2 cigar operations
                // Note: minimap2's CIGAR doesn't include soft clips for unaligned query regions
                // We need to add them based on query_start and query_end (both 0-based)
                if let Some(cigar_ops) = &alignment.cigar {
                    // Pre-allocate capacity: soft clips (0-2) + cigar operations
                    let estimated_capacity = cigar_ops.len() + 2;
                    let mut ops: Vec<Op> = Vec::with_capacity(estimated_capacity);

                    // Calculate leading and trailing soft clip sizes
                    let query_len = mapping.query_len.map(|l| l.get() as i32).unwrap_or(seq.len() as i32);

                    let leading_clip = if mapping.query_start > 0 {
                        mapping.query_start as usize
                    } else {
                        0
                    };
                    let trailing_clip = if mapping.query_end < query_len {
                        (query_len - mapping.query_end) as usize
                    } else {
                        0
                    };

                    // For reverse-complemented alignments, swap the soft clip positions
                    // Add leading soft clip (or trailing if reverse)
                    if is_reverse && trailing_clip > 0 {
                        ops.push(Op::new(Kind::SoftClip, trailing_clip));
                    } else if !is_reverse && leading_clip > 0 {
                        ops.push(Op::new(Kind::SoftClip, leading_clip));
                    }

                    // Add minimap2 CIGAR operations
                    for (len, op_type) in cigar_ops {
                        if *len == 0 {
                            continue
                        }

                        let kind = match op_type {
                            0 => Kind::Match,           // M - Match/Mismatch
                            1 => Kind::Insertion,       // I - Insertion
                            2 => Kind::Deletion,        // D - Deletion
                            3 => Kind::Skip,            // N - Skip (intron)
                            4 => Kind::SoftClip,        // S - Soft clip
                            5 => Kind::HardClip,        // H - Hard clip
                            6 => Kind::Pad,             // P - Pad
                            7 => Kind::SequenceMatch,   // = - Sequence match
                            8 => Kind::SequenceMismatch, // X - Sequence mismatch
                            _ => Kind::Match,           // Default to match for unknown types
                        };
                        ops.push(Op::new(kind, *len as usize));
                    }

                    // Add trailing soft clip (or leading if reverse)
                    if is_reverse && leading_clip > 0 {
                        ops.push(Op::new(Kind::SoftClip, leading_clip));
                    } else if !is_reverse && trailing_clip > 0 {
                        ops.push(Op::new(Kind::SoftClip, trailing_clip));
                    }

                    *record_buf.cigar_mut() = ops.into_iter().collect();
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
        let reads_path = "/data/xurui/project/archive/minimap2-rs/data/h38_dna/SRR17666426_20.fastq";
        let output_sam_path = "/data/xurui/project/archive/minimap2-rs/data/h38_dna/SRR17666426_20_rs.sam";

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

        for alignment in all_alignments.iter() {
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
