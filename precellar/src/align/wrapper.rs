//! Wrapper for community minimap2 crate to provide Minimap2Aligner interface

use anyhow::Result;
use noodles::sam::{self, alignment::record_buf::RecordBuf};
use noodles::fastq;
use std::sync::Arc;
use std::path::PathBuf;

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
        // Build the aligner using community minimap2 API
        let mut builder = minimap2::Aligner::builder();
        
        // Set preset if provided, otherwise use map-ont for long reads
        builder = match opts.preset.as_ref() {
            Some(preset) => builder.preset(preset.clone()),
            None => builder.map_ont(),
        };

        // Load the index
        let aligner = builder
            .with_cigar()
            .with_index(&opts.index_path, None)?;

        // Build SAM header from minimap2 index
        let header = Self::build_header(&aligner);

        Ok(Self {
            aligner,
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
        
        // Map the sequence
        let mappings = self.aligner.map(
            seq,
            false,  // is_rev
            false,  // print_cigar_in_seq
            None,   // max_frag_len
            None,   // extra_flags
            Some(name),
        )?;

        // Convert minimap2 mappings to RecordBuf
        let mut results = Vec::new();
        
        for mapping in mappings {
            let mut record_buf = RecordBuf::default();
            
            // Set name
            record_buf.name_mut().clear();
            record_buf.name_mut().extend_from_slice(name);
            
            // Set sequence and quality
            record_buf.sequence_mut().clear();
            record_buf.sequence_mut().extend(
                seq.iter().map(|&b| b.try_into().unwrap_or(noodles::core::Position::MIN))
            );
            
            if let Some(qual) = record.quality_scores() {
                record_buf.quality_scores_mut().clear();
                record_buf.quality_scores_mut().extend(
                    qual.iter().map(|&q| q.try_into().unwrap_or(0))
                );
            }
            
            // Set mapping information
            if let Some(target_name) = mapping.target_name {
                let ref_id = self.header
                    .reference_sequences()
                    .get_index_of(target_name.as_bytes())
                    .map(|i| i as i32)
                    .unwrap_or(-1);
                    
                if ref_id >= 0 {
                    record_buf.set_reference_sequence_id(Some(ref_id));
                }
            }
            
            record_buf.set_alignment_start(
                noodles::core::Position::try_from(mapping.target_start as usize + 1).ok()
            );
            
            record_buf.set_mapping_quality(
                noodles::sam::alignment::record::MappingQuality::try_from(mapping.mapq).ok()
            );
            
            // Set flags
            let mut flags = noodles::sam::alignment::record::Flags::empty();
            if mapping.strand == minimap2::Strand::Reverse {
                flags |= noodles::sam::alignment::record::Flags::REVERSE_COMPLEMENTED;
            }
            if mapping.is_primary {
                // Primary alignment (no flag set)
            } else {
                flags |= noodles::sam::alignment::record::Flags::SECONDARY;
            }
            record_buf.set_flags(flags);
            
            // Parse and set CIGAR
            if let Some(cigar_str) = mapping.cigar {
                if let Ok(cigar) = cigar_str.parse() {
                    *record_buf.cigar_mut() = cigar;
                }
            }
            
            // Add alignment score and other tags
            let data = record_buf.data_mut();
            data.insert(
                noodles::sam::alignment::record::data::field::tag::Tag::ALIGNMENT_SCORE,
                noodles::sam::alignment::record_buf::data::field::value::Value::Int32(mapping.alignment_score as i32),
            );
            
            results.push(record_buf);
        }
        
        Ok(results)
    }

    fn build_header(aligner: &minimap2::Aligner<minimap2::Built>) -> sam::Header {
        let mut header = sam::Header::default();
        
        // Add reference sequences from index
        for (name, length) in aligner.get_index_names().iter().zip(aligner.get_index_lengths()) {
            use noodles::sam::header::record::value::map::ReferenceSequence;
            use noodles::core::Position;
            
            let ref_seq = ReferenceSequence::new(
                name.clone(),
                Position::try_from(*length as usize).unwrap()
            );
            
            header.reference_sequences_mut()
                .insert(name.as_bytes().to_vec(), ref_seq);
        }
        
        header
    }
}