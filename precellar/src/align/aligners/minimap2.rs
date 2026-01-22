//! Wrapper for community minimap2 crate to provide Minimap2Aligner interface

use anyhow::Result;
use noodles::sam::{self, alignment::record_buf::RecordBuf};
use noodles::sam::alignment::record::cigar::op::Kind;
use noodles::sam::alignment::record::cigar::Op;
use noodles::sam::alignment::record::Flags;
use noodles::fastq;
use std::path::PathBuf;

/// Default maximum insert size for paired-end reads (used for proper pair detection)
/// TODO: get this from .yaml config
const DEFAULT_MAX_INSERT_SIZE: i64 = 1000;

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
        let header = build_header(&aligner);

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
            true,  // cs - long cs tag
            true,  // md - MD tag
            None,   // max_frag_len
            None,   // extra_flags
            Some(name),
        ).map_err(|e| anyhow::anyhow!("Minimap2 mapping failed: {}", e))?;

        // Return unmapped record when no alignments found
        if mappings.is_empty() {
            return Ok(vec![create_unmapped_record(record, false, false)]);
        }

        // Convert minimap2 mappings to RecordBuf
        let results = mappings
            .into_iter()
            .map(|mapping| mapping_to_record_buf(&self.header, &mapping, seq, name, qual))
            .collect();

        Ok(results)
    }

    /// Align a read pair and return alignments for both reads.
    pub fn align_read_pair(
        &mut self,
        read1: &fastq::Record,
        read2: &fastq::Record,
    ) -> Result<(Vec<RecordBuf>, Vec<RecordBuf>)> {
        let seq1 = read1.sequence();
        let seq2 = read2.sequence();
        let name = read1.name();
        let qual1 = read1.quality_scores();
        let qual2 = read2.quality_scores();

        let (mappings1, mappings2) = self.aligner.map_pair(
            seq1,
            seq2,
            true,  // cs - long cs tag
            true,  // md - MD tag
            Some(DEFAULT_MAX_INSERT_SIZE as usize),  // max_frag_len
            None,   // extra_flags
            Some(name),
        ).map_err(|e| anyhow::anyhow!("Minimap2 paired mapping failed: {}", e))?;

        // Extract primary mapping info for mate fields
        let r1_info = extract_mapping_info(&mappings1, &self.header);
        let r2_info = extract_mapping_info(&mappings2, &self.header);

        // Create records for both reads
        let ali1 = create_paired_records(
            &self.header, read1, &mappings1, seq1, name, qual1,
            true, &r1_info, &r2_info,
        );
        let ali2 = create_paired_records(
            &self.header, read2, &mappings2, seq2, name, qual2,
            false, &r2_info, &r1_info,
        );

        Ok((ali1, ali2))
    }
}

/// Build SAM header from minimap2 index
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

/// Primary mapping info extracted for mate field
struct MappingInfo {
    ref_id: Option<usize>,
    pos: Option<noodles::core::Position>,
    end_pos: i64,  // 1-based end position (inclusive)
    is_reverse: bool,
    is_mapped: bool,
}

/// Extract primary mapping info from a list of mappings
fn extract_mapping_info(mappings: &[minimap2::Mapping], header: &sam::Header) -> MappingInfo {
    if let Some(primary) = mappings.iter().find(|m| m.is_primary && !m.is_supplementary) {
        let ref_id = primary.target_name.as_ref().and_then(|name| {
            header.reference_sequences().get_index_of(name.as_bytes())
        });
        let pos = noodles::core::Position::try_from(primary.target_start as usize + 1).ok();
        // target_end is 0-based exclusive, so it equals 1-based inclusive end position
        let end_pos = primary.target_end as i64;
        let is_reverse = primary.strand == minimap2::Strand::Reverse;
        MappingInfo { ref_id, pos, end_pos, is_reverse, is_mapped: true }
    } else {
        MappingInfo { ref_id: None, pos: None, end_pos: 0, is_reverse: false, is_mapped: false }
    }
}

/// Create paired-end records from mappings, setting all flags and mate info in one pass
fn create_paired_records(
    header: &sam::Header,
    record: &fastq::Record,
    mappings: &[minimap2::Mapping],
    seq: &[u8],
    name: &[u8],
    qual: &[u8],
    is_first: bool,
    self_info: &MappingInfo,
    mate_info: &MappingInfo,
) -> Vec<RecordBuf> {
    if mappings.is_empty() {
        // Unmapped read
        let mut rec = create_unmapped_record(record, true, is_first);

        // Set mate unmapped flag if mate is also unmapped
        if !mate_info.is_mapped {
            *rec.flags_mut() |= Flags::MATE_UNMAPPED;
        } else {
            // Mate is mapped - set the RNAME/POS and mate info to mate's position
            *rec.reference_sequence_id_mut() = mate_info.ref_id;
            *rec.alignment_start_mut() = mate_info.pos;
            *rec.mate_reference_sequence_id_mut() = mate_info.ref_id;
            *rec.mate_alignment_start_mut() = mate_info.pos;
            if mate_info.is_reverse {
                *rec.flags_mut() |= Flags::MATE_REVERSE_COMPLEMENTED;
            }
        }
        vec![rec]
    } else {
        // Mapped read
        mappings
            .iter()
            .map(|mapping| {
                let mut rec = mapping_to_record_buf(header, mapping, seq, name, qual);

                // Set paired-end segment flags
                *rec.flags_mut() |= Flags::SEGMENTED;
                if is_first {
                    *rec.flags_mut() |= Flags::FIRST_SEGMENT;
                } else {
                    *rec.flags_mut() |= Flags::LAST_SEGMENT;
                }

                // Set mate info
                if !mate_info.is_mapped {
                    // Mate unmapped - mate position points to self's primary
                    *rec.flags_mut() |= Flags::MATE_UNMAPPED;
                    *rec.mate_reference_sequence_id_mut() = self_info.ref_id;
                    *rec.mate_alignment_start_mut() = self_info.pos;
                } else {
                    // Both mapped - set mate info and check proper pair
                    *rec.mate_reference_sequence_id_mut() = mate_info.ref_id;
                    *rec.mate_alignment_start_mut() = mate_info.pos;
                    if mate_info.is_reverse {
                        *rec.flags_mut() |= Flags::MATE_REVERSE_COMPLEMENTED;
                    }

                    // Calculate TLEN and proper pair status for primary alignments
                    if mapping.is_primary && !mapping.is_supplementary {
                        if let (Some(self_pos), Some(mate_pos)) = (self_info.pos, mate_info.pos) {
                            let self_start = usize::from(self_pos) as i64;
                            let mate_start = usize::from(mate_pos) as i64;

                            // TLEN = rightmost base - leftmost base + 1
                            let leftmost = self_start.min(mate_start);
                            let rightmost = self_info.end_pos.max(mate_info.end_pos);
                            let tlen_abs = (rightmost - leftmost + 1) as i32;

                            // Positive for leftmost read, negative for rightmost
                            // When positions are equal, first segment gets positive TLEN
                            let tlen = if self_start < mate_start || (self_start == mate_start && is_first) {
                                tlen_abs
                            } else {
                                -tlen_abs
                            };
                            *rec.template_length_mut() = tlen;

                            // Check proper pair: same ref, FR orientation, within insert size
                            // TODO: if not FR pair, need specific processing
                            let same_ref = self_info.ref_id == mate_info.ref_id;
                            let fr_orientation = self_info.is_reverse != mate_info.is_reverse;
                            let within_insert = (tlen_abs as i64) <= DEFAULT_MAX_INSERT_SIZE;

                            if same_ref && fr_orientation && within_insert {
                                *rec.flags_mut() |= Flags::PROPERLY_SEGMENTED;
                            }
                        }
                    }
                }

                rec
            })
            .collect()
    }
}

/// Convert a minimap2 Mapping to a SAM RecordBuf
fn mapping_to_record_buf(
    header: &sam::Header,
    mapping: &minimap2::Mapping,
    seq: &[u8],
    name: &[u8],
    qual: &[u8],
) -> RecordBuf {
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

    // For reverse-complemented alignments, sequence and
    // quality scores of SAM record should be reversed as well
    set_sequence_and_quality(&mut record_buf, seq, qual, is_reverse);

    // Set reference sequence ID
    if let Some(target_name) = &mapping.target_name {
        if let Some(id) = header
            .reference_sequences()
            .get_index_of(target_name.as_bytes())
        {
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

    // Add alignment information if available (when .with_cigar() is used)
    if let Some(ref alignment) = mapping.alignment {
        set_cigar(&mut record_buf, mapping, alignment, seq.len(), is_reverse);

        // Add alignment score if available
        if let Some(score) = alignment.alignment_score {
            record_buf.data_mut().insert(
                noodles::sam::alignment::record::data::field::tag::Tag::ALIGNMENT_SCORE,
                noodles::sam::alignment::record_buf::data::field::value::Value::Int32(score),
            );
        }

        // Add MD tag if available
        if let Some(ref md) = alignment.md {
            record_buf.data_mut().insert(
                noodles::sam::alignment::record::data::field::tag::Tag::MISMATCHED_POSITIONS,
                noodles::sam::alignment::record_buf::data::field::value::Value::String(md.clone().into()),
            );
        }

        // Add CS tag if available (minimap2-specific difference string)
        if let Some(ref cs) = alignment.cs {
            // CS tag uses lowercase 'cs' - need to create custom tag
            let cs_tag = noodles::sam::alignment::record::data::field::tag::Tag::new(b'c', b's');
            record_buf.data_mut().insert(
                cs_tag,
                noodles::sam::alignment::record_buf::data::field::value::Value::String(cs.clone().into()),
            );
        }
    }

    record_buf
}

/// Set sequence and quality scores on a RecordBuf, handling reverse complement
fn set_sequence_and_quality(
    record_buf: &mut RecordBuf,
    seq: &[u8],
    qual: &[u8],
    is_reverse: bool,
) {
    let len = seq.len();
    if is_reverse {
        // Pre-allocate with exact capacity to avoid reallocations
        let mut rc_seq = Vec::with_capacity(len);
        let mut rev_qual = Vec::with_capacity(len);

        // Use lookup table for fast, branch-free complement
        for i in (0..len).rev() {
            rc_seq.push(COMPLEMENT_TABLE[seq[i] as usize]);
            // Decode FASTQ ASCII (Phred+33) to raw Phred score
            rev_qual.push(qual[i].saturating_sub(33));
        }

        *record_buf.sequence_mut() = rc_seq.into();
        *record_buf.quality_scores_mut() = rev_qual.into();
    } else {
        *record_buf.sequence_mut() = seq.to_vec().into();

        // Decode FASTQ ASCII (Phred+33) to raw Phred score
        let mut decoded_qual = Vec::with_capacity(len);
        for &b in qual {
            decoded_qual.push(b.saturating_sub(33));
        }
        *record_buf.quality_scores_mut() = decoded_qual.into();
    }
}

/// Create an unmapped record to store FASTQ info.
fn create_unmapped_record(record: &fastq::Record, is_paired: bool, is_first: bool) -> RecordBuf {
    let mut unmapped = RecordBuf::default();

    *unmapped.name_mut() = Some(record.name().to_vec().into());
    *unmapped.sequence_mut() = record.sequence().to_vec().into();

    // Set MAPQ=0 for unmapped reads (SAM convention)
    if let Ok(mq) = noodles::sam::alignment::record::MappingQuality::try_from(0u8) {
        *unmapped.mapping_quality_mut() = Some(mq);
    }

    // Decode FASTQ ASCII (Phred+33) to raw Phred score
    let qual = record.quality_scores();
    let mut decoded_qual = Vec::with_capacity(qual.len());
    for &b in qual {
        decoded_qual.push(b.saturating_sub(33));
    }
    *unmapped.quality_scores_mut() = decoded_qual.into();

    let mut flags = Flags::UNMAPPED;
    
    if is_paired {
        flags |= Flags::SEGMENTED;
        if is_first {
            flags |= Flags::FIRST_SEGMENT;
        } else {
            flags |= Flags::LAST_SEGMENT;
        }
    }
    *unmapped.flags_mut() = flags;

    unmapped
}

/// Set CIGAR string from minimap2 alignment.
/// Note: minimap2's CIGAR doesn't include soft clips for unaligned query regions.
/// We add them based on query_start and query_end (both 0-based).
fn set_cigar(
    record_buf: &mut RecordBuf,
    mapping: &minimap2::Mapping,
    alignment: &minimap2::Alignment,
    seq_len: usize,
    is_reverse: bool,
) {
    if let Some(cigar_ops) = &alignment.cigar {
        // Pre-allocate capacity: soft clips (0-2) + cigar operations
        let estimated_capacity = cigar_ops.len() + 2;
        let mut ops: Vec<Op> = Vec::with_capacity(estimated_capacity);

        // Calculate leading and trailing soft clip sizes
        let query_len = mapping.query_len.map(|l| l.get() as i32).unwrap_or(seq_len as i32);

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
        // Test with long-read scATAC-seq data
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

    #[test]
    #[ignore]
    fn test_alignment_pair() -> Result<()> {
        use std::io::BufRead;

        // Test configuration
        let ref_path = "/data/Public/genome/GRCh38/GRCh38.primary_assembly.genome.fa.gz";
        let read1_path = "/data2/xurui/projects/archive/minimap2-rs/data/sr_pair/SRR891272_1.subset20.fastq";
        let read2_path = "/data2/xurui/projects/archive/minimap2-rs/data/sr_pair/SRR891272_2.subset20.fastq";
        let cli_sam_path = "/data2/xurui/projects/archive/minimap2-rs/data/sr_pair/alignment/SRR891272_cli.subset20.sam"; // generated by minimap2 -ax sr ref.fa read1.fastq read2.fastq
        let output_sam_path = "/data2/xurui/projects/archive/minimap2-rs/data/sr_pair/alignment/SRR891272_rs.subset20.sam";

        // Initialize aligner with "sr" preset
        let opts = Minimap2Opts::new(PathBuf::from(ref_path))
            .with_preset(minimap2::Preset::Sr);
        let mut aligner = Minimap2Aligner::new(opts)?;
        let header = aligner.get_header().clone();

        // Read paired FASTQ files
        let mut reader1 = fastq::io::Reader::new(BufReader::new(File::open(read1_path)?));
        let mut reader2 = fastq::io::Reader::new(BufReader::new(File::open(read2_path)?));

        let records1: Vec<_> = reader1.records().collect::<std::result::Result<_, _>>()?;
        let records2: Vec<_> = reader2.records().collect::<std::result::Result<_, _>>()?;
        assert_eq!(records1.len(), records2.len(), "R1 and R2 must have same length");

        // Align and collect only primary alignments
        let mut primary_alignments = Vec::new();
        for (r1, r2) in records1.iter().zip(records2.iter()) {
            let (ali1, ali2) = aligner.align_read_pair(r1, r2)?;

            for aln in ali1.into_iter().chain(ali2.into_iter()) {
                if !aln.flags().is_secondary() && !aln.flags().is_supplementary() {
                    primary_alignments.push(aln);
                }
            }
        }

        // Write our output SAM
        let output_file = File::create(output_sam_path)?;
        let mut writer = sam::io::Writer::new(output_file);
        writer.write_header(&header)?;
        for aln in &primary_alignments {
            writer.write_alignment_record(&header, aln)?;
        }
        println!("Wrote {} primary alignments to {}", primary_alignments.len(), output_sam_path);

        // Helper: extract columns 1-9 from SAM file (primary alignments only)
        let extract_cols = |path: &str| -> Result<Vec<String>> {
            let file = File::open(path)?;
            let mut lines = Vec::new();
            for line in BufReader::new(file).lines() {
                let line = line?;
                if line.starts_with('@') { continue; }
                let fields: Vec<&str> = line.split('\t').collect();
                if fields.len() < 9 { continue; }
                let flag: u16 = fields[1].parse().unwrap_or(0);
                if (flag & 0x100) != 0 || (flag & 0x800) != 0 { continue; }
                lines.push(fields[0..9].join("\t"));
            }
            Ok(lines)
        };

        // Compare last 40 lines of columns 1-9
        let our_lines = extract_cols(output_sam_path)?;
        let cli_lines = extract_cols(cli_sam_path)?;

        let our_tail: Vec<_> = our_lines.iter().rev().take(40).rev().cloned().collect();
        let cli_tail: Vec<_> = cli_lines.iter().rev().take(40).rev().cloned().collect();

        println!("\n========== COMPARISON (last 40 lines, columns 1-9) ==========");
        let mut mismatches = 0;
        let compare_count = our_tail.len().min(cli_tail.len());

        for i in 0..compare_count {
            if our_tail[i] != cli_tail[i] {
                mismatches += 1;
                println!("Line {} mismatch:", i + 1);
                println!("  Ours: {}", our_tail[i]);
                println!("  CLI:  {}", cli_tail[i]);
            }
        }

        println!("\nCompared {} lines, {} mismatches", compare_count, mismatches);
        assert_eq!(mismatches, 0, "Found {} mismatches in SAM output", mismatches);

        Ok(())
    }
}
