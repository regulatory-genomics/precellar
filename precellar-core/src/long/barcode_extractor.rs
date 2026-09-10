use anyhow::Result;
use indexmap::{IndexMap, IndexSet};
use std::collections::HashSet;
use std::ops::Range;
use std::sync::{Arc, RwLock};

use seqspec::region::{LibSpec, Region};
use seqspec::Modality;

use super::{
    barcode_index::BarcodeIndex,
    collect_end_regions, collect_target_flanks,
    sequence_aligner::{
        fitting_alignment_with_span, CompositeAlignmentResult, CompositePattern, FittingAligner,
    },
    EndRegions,
};

/// Final long-read barcode extraction result
#[derive(Debug, Clone)]
pub struct LongReadBarcodeResult {
    /// Extracted barcode sequence (standardized from whitelist)
    pub barcode: Option<Vec<u8>>,
    /// UMI sequence and quality in designed-library orientation
    pub umi: Option<noodles_fastq::Record>,
    /// Extraction confidence (0.0 if no barcode extracted)
    pub confidence: f64,
    /// Whether the end segments were reverse-complemented during extraction
    pub is_reverse_complemented: bool,
    /// Trim length for the 5' side in designed-library orientation
    pub five_prime_trim: usize,
    /// Trim length for the 3' side in designed-library orientation
    pub three_prime_trim: usize,
}

impl LongReadBarcodeResult {
    /// Whether extraction was successful
    pub fn is_success(&self) -> bool {
        self.barcode.is_some()
    }
}

/// Per-read statistics from the long-read barcode extraction process, used for QC accumulation.
#[derive(Debug, Clone, Default)]
pub struct LrBarcodeExtractionStats {
    /// Whether the chosen orientation was forward (false) or reverse (true)
    pub is_reverse: bool,
    /// Whether composite alignment passed the quality gate
    pub composite_pass: bool,
    /// Per barcode-group consensus results: (group_name, intersection_hit)
    pub consensus_results: Vec<(String, bool)>,
}

/// Successfully extracted barcode with metadata (may contain multiple tied-best candidates)
#[derive(Debug, Clone)]
pub struct ExtractedBarcode {
    pub region_id: String,
    pub region_name: String, // For grouping same-barcode regions across ends
    pub barcodes: Vec<Vec<u8>>, // All tied-best candidates from whitelist
    pub confidence: f64,     // confidence = 1.0 - (min_edit_distance / barcode_length)
    /// Observed span in normalized end-segment coordinates.
    ///
    /// When `barcodes` contains tied-best candidates, this span is traced only
    /// from `barcodes[0]`. Multi-end consensus may later choose another barcode,
    /// so the span is not guaranteed to belong to the final combined barcode.
    pub traced_span: Range<usize>,
}

/// Evidence for read orientation determination using composite alignment
#[derive(Debug, Clone)]
pub struct OrientationEvidence {
    /// Total composite alignment score (sum of 5' and 3' scores)
    pub total_score: i32,
    /// Average match rate of fixed regions
    pub avg_fixed_match_rate: f64,
    /// Number of fixed regions with match rate >= threshold
    pub num_good_fixed_regions: usize,
    /// Total number of fixed regions expected
    pub total_fixed_regions: usize,
}

impl OrientationEvidence {
    /// Check if evidence meets conservative thresholds for accepting orientation
    /// Conservative thresholds: ≥2 anchors, ≥0.8 confidence
    pub fn meets_threshold(&self) -> bool {
        self.num_good_fixed_regions >= 2 && self.avg_fixed_match_rate >= 0.8
    }

    /// Build from composite alignment results of both ends
    pub fn from_composite(
        five_prime: &Option<CompositeAlignmentResult>,
        three_prime: &Option<CompositeAlignmentResult>,
    ) -> Self {
        let mut total_score = 0i32;
        let mut fixed_rates: Vec<(f64, usize)> = Vec::new(); // (match_rate, region_length)

        for result in [five_prime, three_prime].into_iter().flatten() {
            total_score += result.score;
            for mapping in &result.region_mappings {
                if !mapping.is_spacer {
                    let len = mapping.read_end.saturating_sub(mapping.read_start);
                    fixed_rates.push((mapping.match_rate, len));
                }
            }
        }

        let num_good = fixed_rates.iter().filter(|&&(r, _)| r >= 0.8).count();
        let total_len: usize = fixed_rates.iter().map(|(_, len)| len).sum();
        let avg_rate = if total_len == 0 {
            0.0
        } else {
            fixed_rates
                .iter()
                .map(|(r, len)| r * *len as f64)
                .sum::<f64>()
                / total_len as f64
        };

        OrientationEvidence {
            total_score,
            avg_fixed_match_rate: avg_rate,
            num_good_fixed_regions: num_good,
            total_fixed_regions: fixed_rates.len(),
        }
    }
}

/// End segment with precomputed composite alignment for reuse during orientation detection.
#[derive(Debug)]
struct EndSegmentWithAlignment {
    /// The cut segment sequence (empty if sequence was too short)
    sequence: Vec<u8>,
    /// Precomputed composite alignment result
    composite_result: Option<CompositeAlignmentResult>,
}

impl EndSegmentWithAlignment {
    /// Create an empty segment (used when source sequence is too short)
    fn empty() -> Self {
        Self {
            sequence: Vec::new(),
            composite_result: None,
        }
    }

    /// Check if this segment is empty (no sequence was extracted)
    fn is_empty(&self) -> bool {
        self.sequence.is_empty()
    }
}

/// Find the best barcode match using k-mer indexed search
fn find_best_barcode_match(
    candidate_seq: &[u8],
    normalized_window: &[u8],
    barcode_region: &Arc<RwLock<Region>>,
    whitelist_indices: &IndexMap<String, BarcodeIndex>,
    reverse_barcode_for_tracing: bool,
) -> Option<ExtractedBarcode> {
    let region_guard = barcode_region.read().unwrap();
    let region_id = region_guard.region_id.clone();
    let region_name = region_guard.name.clone();
    drop(region_guard);

    // Get whitelist index for this region
    let index = whitelist_indices.get(&region_id)?;

    // Find all tied-best matches using k-mer voting + fitting alignment
    let (matched_barcodes, confidence) = index.find_best_match(candidate_seq)?;
    let reverse_barcode;
    let barcode_for_tracing = if reverse_barcode_for_tracing {
        reverse_barcode = seqspec::utils::rev_compl(&matched_barcodes[0]);
        reverse_barcode.as_slice()
    } else {
        matched_barcodes[0].as_slice()
    };
    let (traced_distance, traced_span) =
        fitting_alignment_with_span(barcode_for_tracing, normalized_window);
    debug_assert_eq!(
        confidence,
        1.0 - traced_distance as f64 / matched_barcodes[0].len().max(1) as f64
    );

    Some(ExtractedBarcode {
        region_id,
        region_name,
        barcodes: matched_barcodes,
        confidence,
        traced_span,
    })
}

/// Locate the barcode extraction window using a strict topological hierarchy
/// based on aligned fixed region boundaries.
///
/// Hierarchy:
/// 1. Left-Adjacent Anchoring: the immediate left neighbor is a fixed region →
///    start from its read_end, extend right by up to 1.2× barcode_len, capped by
///    the next right fixed region's read_start.
/// 2. Right-Adjacent Anchoring: the immediate right neighbor is a fixed region →
///    anchor to its read_start, extract left by up to 1.2× barcode_len.
/// 3. Sandwiched (Non-Adjacent): barcode between two non-adjacent fixed regions →
///    extract the entire gap between them.
/// 4. Terminal / Edge: barcode at the EndSegment edge → extract from segment boundary
///    to nearest anchor.
fn locate_barcode_extraction_window(
    barcode_region: &Arc<RwLock<Region>>,
    end_regions: &EndRegions,
    composite_result: &CompositeAlignmentResult,
    sequence_length: usize,
) -> Option<(usize, usize)> {
    let barcode_length = {
        let r = barcode_region.read().unwrap();
        r.max_len as usize
    };
    let target_length = ((barcode_length as f64) * 1.2).ceil() as usize;

    // Search window must be above a minimum length threshold to avoid spurious extractions
    let min_length = ((barcode_length as f64) * 0.8).floor() as usize;

    // Find barcode position in forward-order end regions.
    let barcode_pos = end_regions
        .regions
        .iter()
        .position(|r| Arc::ptr_eq(r, barcode_region))?;

    // Helper: look up a region's aligned coordinates from composite result
    let get_fixed_coords = |region: &Arc<RwLock<Region>>| -> Option<(usize, usize)> {
        composite_result
            .region_mappings
            .iter()
            .find(|m| Arc::ptr_eq(&m.region, region) && !m.is_spacer)
            .map(|m| (m.read_start, m.read_end))
    };

    // Check immediate left neighbor for structural adjacency.
    let left_adjacent = if barcode_pos > 0 {
        let left_region = &end_regions.regions[barcode_pos - 1];
        if left_region.read().unwrap().sequence_type.is_fixed() {
            get_fixed_coords(left_region)
        } else {
            None
        }
    } else {
        None
    };

    // Check immediate right neighbor for structural adjacency.
    let right_adjacent = if barcode_pos + 1 < end_regions.regions.len() {
        let right_region = &end_regions.regions[barcode_pos + 1];
        if right_region.read().unwrap().sequence_type.is_fixed() {
            get_fixed_coords(right_region)
        } else {
            None
        }
    } else {
        None
    };

    // Find nearest fixed region on each side (not necessarily adjacent).
    let nearest_left_fixed = (0..barcode_pos).rev().find_map(|i| {
        let region = &end_regions.regions[i];
        if region.read().unwrap().sequence_type.is_fixed() {
            get_fixed_coords(region)
        } else {
            None
        }
    });

    let nearest_right_fixed = ((barcode_pos + 1)..end_regions.regions.len()).find_map(|i| {
        let region = &end_regions.regions[i];
        if region.read().unwrap().sequence_type.is_fixed() {
            get_fixed_coords(region)
        } else {
            None
        }
    });

    let (start, end) = if let Some((_, left_end)) = left_adjacent {
        // Case 1: Left-Adjacent Anchoring
        let max_end = left_end + target_length;
        let capped_end = if let Some((right_start, _)) = nearest_right_fixed {
            max_end.min(right_start)
        } else {
            max_end.min(sequence_length)
        };
        (left_end, capped_end)
    } else if let Some((right_start, _)) = right_adjacent {
        // Case 2: Right-Adjacent Anchoring
        let min_start = right_start.saturating_sub(target_length);
        let capped_start = if let Some((_, left_end)) = nearest_left_fixed {
            min_start.max(left_end)
        } else {
            min_start
        };
        (capped_start, right_start)
    } else {
        match (nearest_left_fixed, nearest_right_fixed) {
            // Case 3: Sandwiched (Non-Adjacent) Anchoring
            (Some((_, left_end)), Some((right_start, _))) => (left_end, right_start),
            // Case 4: Terminal / Edge Extraction
            (None, Some((right_start, _))) => (0, right_start),
            (Some((_, left_end)), None) => (left_end, sequence_length),
            (None, None) => {
                return None;
            }
        }
    };

    let actual_length = end.saturating_sub(start);
    if actual_length >= min_length && end <= sequence_length {
        Some((start, end))
    } else {
        None
    }
}

/// Main barcode extractor for long reads.
/// Acts as an orchestrator, delegating tasks to specialized components.
#[derive(Debug)]
pub struct BarcodeExtractor {
    five_prime_regions: EndRegions,
    three_prime_regions: EndRegions,
    five_prime_trim_regions: Vec<Arc<RwLock<Region>>>, // Flanking non-target regions for 5' end
    three_prime_trim_regions: Vec<Arc<RwLock<Region>>>, // Flanking non-target regions for 3' end
    /// Prebuilt composite patterns (built once, reused for all reads)
    five_prime_composite: CompositePattern,
    three_prime_composite: CompositePattern,
    /// Barcode whitelist indices: region_id -> BarcodeIndex
    whitelist_indices: IndexMap<String, BarcodeIndex>,
    /// Region IDs whose extracted candidate should be reverse-complemented before whitelist matching
    rc_regions: HashSet<String>,
    /// All barcode group names that must be resolved for a valid barcode (deduplicated, first-seen order)
    required_barcode_groups: IndexSet<String>,
    /// UMI region declared by the selected modality, if any.
    umi_region: Option<Arc<RwLock<Region>>>,
    /// Direct region neighbor on the UMI's 5' side, when usable as an anchor.
    umi_five_prime_anchor: Option<Arc<RwLock<Region>>>,
    /// Direct region neighbor on the UMI's 3' side, when usable as an anchor.
    umi_three_prime_anchor: Option<Arc<RwLock<Region>>>,
}

impl BarcodeExtractor {
    /// Minimum average fixed-region match rate for barcode extraction
    const MIN_MATCH_RATE: f64 = 0.7;

    /// Create a new barcode extractor for the given library_spec, modality, and whitelists.
    pub fn new(
        lib_spec: &LibSpec,
        modality: &Modality,
        whitelists: IndexMap<String, IndexSet<Vec<u8>>>,
    ) -> Result<Self> {
        // Check that all barcode regions have non-empty whitelists and build indices
        let mut whitelist_indices = IndexMap::new();
        for (region_id, whitelist) in &whitelists {
            if whitelist.is_empty() {
                anyhow::bail!(
                    "Barcode region '{}' does not have a whitelist. Long-read processing requires all barcode regions to have whitelists.",
                    region_id
                );
            }
            // Build k-mer index for fast matching
            whitelist_indices.insert(region_id.clone(), BarcodeIndex::new(whitelist));
        }

        let (five_prime_regions, three_prime_regions) = collect_end_regions(lib_spec, modality)?;
        let (five_prime_trim_regions, three_prime_trim_regions) =
            collect_target_flanks(lib_spec, modality)?;
        let five_prime_composite = five_prime_regions.build_composite_pattern();
        let three_prime_composite = three_prime_regions.build_composite_pattern();

        let modality_region = lib_spec
            .get_modality(modality)
            .ok_or_else(|| anyhow::anyhow!("Cannot find specified modality: {:?}", modality))?;
        let modality_guard = modality_region.read().unwrap();
        let umi_positions: Vec<_> = modality_guard
            .subregions
            .iter()
            .enumerate()
            .filter_map(|(idx, region)| region.read().unwrap().region_type.is_umi().then_some(idx))
            .collect();
        anyhow::ensure!(
            umi_positions.len() <= 1,
            "long-read modality {:?} contains multiple UMI regions; only one is supported",
            modality
        );

        let umi_region = umi_positions
            .first()
            .map(|&idx| modality_guard.subregions[idx].clone());
        let is_umi_anchor = |region: &&Arc<RwLock<Region>>| {
            let region = region.read().unwrap();
            region.sequence_type.is_fixed() || region.region_type.is_barcode()
        };
        let umi_five_prime_anchor = umi_positions.first().and_then(|&idx| {
            idx.checked_sub(1)
                .and_then(|idx| modality_guard.subregions.get(idx))
                .filter(is_umi_anchor)
                .cloned()
        });
        let umi_three_prime_anchor = umi_positions.first().and_then(|&idx| {
            modality_guard
                .subregions
                .get(idx + 1)
                .filter(is_umi_anchor)
                .cloned()
        });
        drop(modality_guard);

        if let Some(umi) = &umi_region {
            if umi_five_prime_anchor.is_none() && umi_three_prime_anchor.is_none() {
                log::warn!(
                    "UMI region '{}' has no directly adjacent fixed/barcode anchor; long-read UMI extraction is disabled for modality {:?}",
                    umi.read().unwrap().region_id,
                    modality,
                );
            }
        }

        // Validate minimum total fixed sequence length (12bp) only for ends with barcodes
        const MIN_FIXED_LEN: usize = 12;
        if five_prime_regions.has_barcode {
            let five_prime_fixed_len = five_prime_composite.total_fixed_len();
            if five_prime_fixed_len < MIN_FIXED_LEN {
                anyhow::bail!(
                    "5' end has insufficient total fixed sequence length ({} bp < {} bp minimum) for composite alignment",
                    five_prime_fixed_len, MIN_FIXED_LEN
                );
            }
        }
        if three_prime_regions.has_barcode {
            let three_prime_fixed_len = three_prime_composite.total_fixed_len();
            if three_prime_fixed_len < MIN_FIXED_LEN {
                anyhow::bail!(
                    "3' end has insufficient total fixed sequence length ({} bp < {} bp minimum) for composite alignment",
                    three_prime_fixed_len, MIN_FIXED_LEN
                );
            }
        }

        // Collect region IDs that require reverse-complementing before whitelist matching
        // and collect all barcode group names (deduplicated)
        let mut rc_regions = HashSet::new();
        let mut required_barcode_groups = IndexSet::new();
        for region in five_prime_regions
            .regions
            .iter()
            .chain(three_prime_regions.regions.iter())
        {
            let guard = region.read().unwrap();
            if guard.region_type.is_barcode() {
                required_barcode_groups.insert(guard.name.clone());
                if let Some(onlist) = &guard.onlist {
                    if onlist.rc {
                        rc_regions.insert(guard.region_id.clone());
                    }
                }
            }
        }

        Ok(Self {
            five_prime_regions,
            three_prime_regions,
            five_prime_trim_regions,
            three_prime_trim_regions,
            five_prime_composite,
            three_prime_composite,
            whitelist_indices,
            rc_regions,
            required_barcode_groups,
            umi_region,
            umi_five_prime_anchor,
            umi_three_prime_anchor,
        })
    }

    /// Whether this layout has a UMI with at least one supported direct anchor.
    pub fn expects_umi(&self) -> bool {
        self.umi_region.is_some()
            && (self.umi_five_prime_anchor.is_some() || self.umi_three_prime_anchor.is_some())
    }

    fn should_analyze_end(&self, end_regions: &EndRegions) -> bool {
        end_regions.has_barcode
            || [
                self.umi_five_prime_anchor.as_ref(),
                self.umi_three_prime_anchor.as_ref(),
            ]
            .into_iter()
            .flatten()
            .any(|anchor| {
                end_regions
                    .regions
                    .iter()
                    .any(|region| Arc::ptr_eq(region, anchor))
            })
    }

    /// Extract barcode from a FASTQ record with automatic orientation detection.
    ///
    /// This method uses an adaptive forward-first strategy with anchor reuse:
    /// 1. Sample end windows and find anchors for forward orientation
    /// 2. If forward meets threshold, extract barcodes using cached anchors (fast path)
    /// 3. Otherwise, sample and RC end windows, find anchors for reverse orientation
    /// 4. Compare evidence, choose better orientation, extract barcodes using cached anchors
    ///
    /// IMPORTANT: Only end segments are processed and RC'd, not the full read.
    /// Anchors found during orientation detection are reused for barcode extraction.
    pub fn extract_barcode(
        &self,
        record: &noodles_fastq::Record,
    ) -> Result<(LongReadBarcodeResult, LrBarcodeExtractionStats)> {
        let sequence = record.sequence();
        let quality = record.quality_scores();
        anyhow::ensure!(
            sequence.len() == quality.len(),
            "long-read FASTQ sequence and quality lengths differ"
        );

        // Step 1: Sample end windows and composite-align for forward orientation
        let forward_5p = if self.should_analyze_end(&self.five_prime_regions) {
            self.sample_and_analyze_end_window(
                sequence,
                quality,
                &self.five_prime_regions,
                &self.five_prime_composite,
                false,
            )?
        } else {
            EndSegmentWithAlignment::empty()
        };
        let forward_3p = if self.should_analyze_end(&self.three_prime_regions) {
            self.sample_and_analyze_end_window(
                sequence,
                quality,
                &self.three_prime_regions,
                &self.three_prime_composite,
                false,
            )?
        } else {
            EndSegmentWithAlignment::empty()
        };

        let forward_evidence = OrientationEvidence::from_composite(
            &forward_5p.composite_result,
            &forward_3p.composite_result,
        );

        // Step 2: If forward meets threshold, use it directly (fast path)
        if forward_evidence.meets_threshold() {
            let (result, consensus_results) =
                self.extract_barcodes_from_segments(record, forward_5p, forward_3p, false)?;
            let stats = LrBarcodeExtractionStats {
                is_reverse: false,
                composite_pass: true,
                consensus_results,
            };
            return Ok((result, stats));
        }

        // Step 3: Forward didn't meet threshold, try reverse orientation
        // Only Rc end segments, not full read.
        let reverse_5p = if self.should_analyze_end(&self.five_prime_regions) {
            self.sample_and_analyze_end_window(
                sequence,
                quality,
                &self.five_prime_regions,
                &self.five_prime_composite,
                true,
            )?
        } else {
            EndSegmentWithAlignment::empty()
        };
        let reverse_3p = if self.should_analyze_end(&self.three_prime_regions) {
            self.sample_and_analyze_end_window(
                sequence,
                quality,
                &self.three_prime_regions,
                &self.three_prime_composite,
                true,
            )?
        } else {
            EndSegmentWithAlignment::empty()
        };

        let reverse_evidence = OrientationEvidence::from_composite(
            &reverse_5p.composite_result,
            &reverse_3p.composite_result,
        );

        // Step 4: Compare and choose better orientation
        let should_rc = self.compare_evidence(&forward_evidence, &reverse_evidence);
        let chosen_evidence = if should_rc {
            &reverse_evidence
        } else {
            &forward_evidence
        };

        // Step 5: Gate on composite alignment quality before barcode extraction
        if chosen_evidence.avg_fixed_match_rate < Self::MIN_MATCH_RATE {
            #[cfg(debug_assertions)]
            {
                let read_name = std::str::from_utf8(record.name()).unwrap_or("<invalid>");
                eprintln!(
                    "[DEBUG] Read '{}': low composite alignment quality \
                    (chosen avg_fixed_match_rate: {:.3}, threshold: {:.1}; \
                    forward: {:.3}, reverse: {:.3})",
                    read_name,
                    chosen_evidence.avg_fixed_match_rate,
                    Self::MIN_MATCH_RATE,
                    forward_evidence.avg_fixed_match_rate,
                    reverse_evidence.avg_fixed_match_rate,
                );
            }
            return Ok((
                LongReadBarcodeResult {
                    barcode: None,
                    umi: None,
                    confidence: 0.0,
                    is_reverse_complemented: should_rc,
                    five_prime_trim: if should_rc {
                        self.calculate_trim_length(
                            &reverse_5p,
                            &self.five_prime_regions,
                            &self.five_prime_trim_regions,
                        )
                    } else {
                        self.calculate_trim_length(
                            &forward_5p,
                            &self.five_prime_regions,
                            &self.five_prime_trim_regions,
                        )
                    },
                    three_prime_trim: if should_rc {
                        self.calculate_trim_length(
                            &reverse_3p,
                            &self.three_prime_regions,
                            &self.three_prime_trim_regions,
                        )
                    } else {
                        self.calculate_trim_length(
                            &forward_3p,
                            &self.three_prime_regions,
                            &self.three_prime_trim_regions,
                        )
                    },
                },
                LrBarcodeExtractionStats {
                    is_reverse: should_rc,
                    composite_pass: false,
                    consensus_results: Vec::new(),
                },
            ));
        }

        // Step 6: Extract barcodes using cached anchors from the chosen orientation
        let (result, consensus_results) = if should_rc {
            self.extract_barcodes_from_segments(record, reverse_5p, reverse_3p, true)?
        } else {
            self.extract_barcodes_from_segments(record, forward_5p, forward_3p, false)?
        };
        let stats = LrBarcodeExtractionStats {
            is_reverse: should_rc,
            composite_pass: true,
            consensus_results,
        };
        Ok((result, stats))
    }

    /// Sample an end window and find fixed-region alignment in one pass.
    /// Returns the sampled window with precomputed anchors for reuse.
    /// Returns an empty window if the FASTQ sequence is too short.
    fn sample_and_analyze_end_window(
        &self,
        sequence: &[u8],
        _quality: &[u8],
        end_regions: &EndRegions,
        composite: &CompositePattern,
        should_rc: bool,
    ) -> Result<EndSegmentWithAlignment> {
        let window_length = end_regions.calculate_cut_length();

        if sequence.len() < window_length {
            return Ok(EndSegmentWithAlignment::empty());
        }

        // Keep end segments in forward order of the designed library structure.
        // In this coordinate system, left = 5' side and right = 3' side.
        let seq_segment = match (end_regions.end_type, should_rc) {
            (super::EndType::FivePrime, false) | (super::EndType::ThreePrime, true) => {
                &sequence[..window_length]
            }
            (super::EndType::ThreePrime, false) | (super::EndType::FivePrime, true) => {
                let start = sequence.len() - window_length;
                &sequence[start..]
            }
        };

        let final_seq = if should_rc {
            seqspec::utils::rev_compl(seq_segment)
        } else {
            seq_segment.to_vec()
        };

        // Perform composite alignment
        let composite_result = FittingAligner::align_composite(&final_seq, composite)?;

        Ok(EndSegmentWithAlignment {
            sequence: final_seq,
            composite_result: Some(composite_result),
        })
    }

    /// Extract barcodes from precomputed segments with composite alignment.
    /// Returns (LongReadBarcodeResult, consensus_results) where consensus_results
    /// contains (group_name, intersection_hit) for each multi-end barcode group.
    fn extract_barcodes_from_segments(
        &self,
        record: &noodles_fastq::Record,
        five_prime: EndSegmentWithAlignment,
        three_prime: EndSegmentWithAlignment,
        is_reverse_complemented: bool,
    ) -> Result<(LongReadBarcodeResult, Vec<(String, bool)>)> {
        let mut extracted_barcodes = Vec::new();
        let five_prime_trim = self.calculate_trim_length(
            &five_prime,
            &self.five_prime_regions,
            &self.five_prime_trim_regions,
        );
        let three_prime_trim = self.calculate_trim_length(
            &three_prime,
            &self.three_prime_regions,
            &self.three_prime_trim_regions,
        );

        if !five_prime.is_empty() {
            let barcode_results =
                self.extract_barcodes_from_composite(&five_prime, &self.five_prime_regions)?;
            extracted_barcodes.extend(barcode_results);
        }

        if !three_prime.is_empty() {
            let barcode_results =
                self.extract_barcodes_from_composite(&three_prime, &self.three_prime_regions)?;
            extracted_barcodes.extend(barcode_results);
        }

        let umi = self.extract_umi(
            record,
            &five_prime,
            &three_prime,
            &extracted_barcodes,
            is_reverse_complemented,
        );
        let (mut result, consensus_results) = self.combine_barcodes(extracted_barcodes)?;
        result.umi = umi;
        result.is_reverse_complemented = is_reverse_complemented;
        result.five_prime_trim = five_prime_trim;
        result.three_prime_trim = three_prime_trim;
        Ok((result, consensus_results))
    }

    fn calculate_trim_length(
        &self,
        segment: &EndSegmentWithAlignment,
        end_regions: &EndRegions,
        trim_regions: &[Arc<RwLock<Region>>],
    ) -> usize {
        if trim_regions.is_empty() {
            return 0;
        }

        let fallback = Self::sum_average_lengths(trim_regions);
        let anchor = match Self::find_trim_anchor(end_regions) {
            Some(anchor) => anchor,
            None => return fallback,
        };
        let composite_result = match &segment.composite_result {
            Some(result) => result,
            None => return fallback,
        };
        let mapping = match composite_result
            .region_mappings
            .iter()
            .find(|m| Arc::ptr_eq(&m.region, &anchor) && !m.is_spacer)
        {
            Some(mapping) => mapping,
            None => return fallback,
        };
        let base = match end_regions.end_type {
            super::EndType::FivePrime => mapping.read_end,
            super::EndType::ThreePrime => segment.sequence.len().saturating_sub(mapping.read_start),
        };
        let offset = Self::sum_average_lengths(Self::trim_offset_regions(
            trim_regions,
            &anchor,
            end_regions.end_type,
        ));
        base + offset
    }

    fn find_trim_anchor(end_regions: &EndRegions) -> Option<Arc<RwLock<Region>>> {
        match end_regions.end_type {
            super::EndType::FivePrime => end_regions
                .regions
                .iter()
                .rfind(|region| region.read().unwrap().sequence_type.is_fixed())
                .cloned(),
            super::EndType::ThreePrime => end_regions
                .regions
                .iter()
                .find(|region| region.read().unwrap().sequence_type.is_fixed())
                .cloned(),
        }
    }

    fn trim_offset_regions<'a>(
        trim_regions: &'a [Arc<RwLock<Region>>],
        anchor: &Arc<RwLock<Region>>,
        end_type: super::EndType,
    ) -> &'a [Arc<RwLock<Region>>] {
        match trim_regions
            .iter()
            .position(|region| Arc::ptr_eq(region, anchor))
        {
            Some(anchor_idx) => match end_type {
                super::EndType::FivePrime => &trim_regions[anchor_idx + 1..],
                super::EndType::ThreePrime => &trim_regions[..anchor_idx],
            },
            None => &[],
        }
    }

    fn sum_average_lengths(regions: &[Arc<RwLock<Region>>]) -> usize {
        regions
            .iter()
            .map(|region| {
                let region = region.read().unwrap();
                ((region.min_len + region.max_len) / 2) as usize
            })
            .sum()
    }

    /// Extract barcodes using composite alignment result.
    /// Uses topological hierarchy to locate barcode extraction windows
    /// based on aligned fixed region boundaries.
    fn extract_barcodes_from_composite(
        &self,
        segment: &EndSegmentWithAlignment,
        end_regions: &EndRegions,
    ) -> Result<Vec<ExtractedBarcode>> {
        let composite_result = match &segment.composite_result {
            Some(r) => r,
            None => return Ok(Vec::new()),
        };

        let barcode_regions = end_regions.get_barcode_regions();
        let mut extracted_barcodes = Vec::new();

        for barcode_region in &barcode_regions {
            let window = locate_barcode_extraction_window(
                barcode_region,
                end_regions,
                composite_result,
                segment.sequence.len(),
            );

            if let Some((start, end)) = window {
                let raw_seq = &segment.sequence[start..end];
                let rc_seq;
                let should_rc = self
                    .rc_regions
                    .contains(&barcode_region.read().unwrap().region_id);
                let candidate_seq = if should_rc {
                    rc_seq = seqspec::utils::rev_compl(raw_seq);
                    rc_seq.as_slice()
                } else {
                    raw_seq
                };

                if let Some(mut matched) = find_best_barcode_match(
                    candidate_seq,
                    raw_seq,
                    barcode_region,
                    &self.whitelist_indices,
                    should_rc,
                ) {
                    matched.traced_span =
                        start + matched.traced_span.start..start + matched.traced_span.end;
                    extracted_barcodes.push(matched);
                }
            }
        }

        Ok(extracted_barcodes)
    }

    fn extract_umi(
        &self,
        record: &noodles_fastq::Record,
        five_prime: &EndSegmentWithAlignment,
        three_prime: &EndSegmentWithAlignment,
        barcodes: &[ExtractedBarcode],
        is_reverse_complemented: bool,
    ) -> Option<noodles_fastq::Record> {
        let umi_region = self.umi_region.as_ref()?;
        if !self.expects_umi() {
            return None;
        }
        let umi_length = {
            let umi = umi_region.read().unwrap();
            (umi.min_len as usize + umi.max_len as usize) / 2
        };

        // Fixed anchors precede barcode anchors; within each type the 5' side
        // precedes the 3' side. Once an observed anchor is found, extraction
        // failure for that anchor does not change the topology preference.
        for fixed_anchor in [true, false] {
            for (anchor, is_five_prime_neighbor) in [
                (self.umi_five_prime_anchor.as_ref(), true),
                (self.umi_three_prime_anchor.as_ref(), false),
            ] {
                let Some(anchor) = anchor else {
                    continue;
                };
                let anchor_matches_type = {
                    let anchor = anchor.read().unwrap();
                    if fixed_anchor {
                        anchor.sequence_type.is_fixed()
                    } else {
                        anchor.region_type.is_barcode()
                    }
                };
                if !anchor_matches_type {
                    continue;
                }

                let Some((anchor_span, end_type)) =
                    self.observed_anchor_span(anchor, five_prime, three_prime, barcodes)
                else {
                    continue;
                };
                let normalized_umi_span = if is_five_prime_neighbor {
                    anchor_span.end..anchor_span.end.checked_add(umi_length)?
                } else {
                    anchor_span.start.checked_sub(umi_length)?..anchor_span.start
                };
                let segment = match end_type {
                    super::EndType::FivePrime => five_prime,
                    super::EndType::ThreePrime => three_prime,
                };
                let raw_span = Self::normalized_to_raw_span(
                    normalized_umi_span,
                    end_type,
                    is_reverse_complemented,
                    record.sequence().len(),
                    segment.sequence.len(),
                )?;

                let mut sequence = record.sequence().get(raw_span.clone())?.to_vec();
                let mut quality = record.quality_scores().get(raw_span)?.to_vec();
                sequence.make_ascii_uppercase();
                if !sequence
                    .iter()
                    .all(|base| matches!(base, b'A' | b'C' | b'G' | b'T'))
                {
                    return None;
                }
                if is_reverse_complemented {
                    sequence = seqspec::utils::rev_compl(&sequence);
                    quality.reverse();
                }
                return Some(noodles_fastq::Record::new(
                    record.definition().clone(),
                    sequence,
                    quality,
                ));
            }
        }

        None
    }

    fn observed_anchor_span(
        &self,
        anchor: &Arc<RwLock<Region>>,
        five_prime: &EndSegmentWithAlignment,
        three_prime: &EndSegmentWithAlignment,
        barcodes: &[ExtractedBarcode],
    ) -> Option<(Range<usize>, super::EndType)> {
        for (end_regions, segment) in [
            (&self.five_prime_regions, five_prime),
            (&self.three_prime_regions, three_prime),
        ] {
            if !end_regions
                .regions
                .iter()
                .any(|region| Arc::ptr_eq(region, anchor))
            {
                continue;
            }

            let anchor_guard = anchor.read().unwrap();
            let span = if anchor_guard.region_type.is_barcode() {
                barcodes
                    .iter()
                    .find(|barcode| barcode.region_id == anchor_guard.region_id)
                    .map(|barcode| barcode.traced_span.clone())
            } else {
                segment.composite_result.as_ref().and_then(|result| {
                    result
                        .region_mappings
                        .iter()
                        .find(|mapping| Arc::ptr_eq(&mapping.region, anchor) && !mapping.is_spacer)
                        .map(|mapping| mapping.read_start..mapping.read_end)
                })
            };
            if let Some(span) = span {
                return Some((span, end_regions.end_type));
            }
        }

        None
    }

    fn normalized_to_raw_span(
        normalized_span: Range<usize>,
        end_type: super::EndType,
        is_reverse_complemented: bool,
        read_length: usize,
        segment_length: usize,
    ) -> Option<Range<usize>> {
        if segment_length > read_length {
            return None;
        }
        let uses_prefix = matches!(
            (end_type, is_reverse_complemented),
            (super::EndType::FivePrime, false) | (super::EndType::ThreePrime, true)
        );
        let raw_segment = if uses_prefix {
            0..segment_length
        } else {
            read_length - segment_length..read_length
        };
        let raw_span = if is_reverse_complemented {
            raw_segment.end.checked_sub(normalized_span.end)?
                ..raw_segment.end.checked_sub(normalized_span.start)?
        } else {
            raw_segment.start.checked_add(normalized_span.start)?
                ..raw_segment.start.checked_add(normalized_span.end)?
        };

        (raw_span.start <= raw_span.end && raw_span.end <= read_length).then_some(raw_span)
    }

    /// Get whitelist indices reference for external access
    pub fn whitelist_indices(&self) -> &IndexMap<String, BarcodeIndex> {
        &self.whitelist_indices
    }

    /// Get 5' end regions for external access
    pub fn five_prime_regions(&self) -> &EndRegions {
        &self.five_prime_regions
    }

    /// Get 3' end regions for external access
    pub fn three_prime_regions(&self) -> &EndRegions {
        &self.three_prime_regions
    }

    /// Compare evidence and choose which orientation is better
    /// Returns true if reverse-complement should be used
    fn compare_evidence(
        &self,
        forward: &OrientationEvidence,
        reverse: &OrientationEvidence,
    ) -> bool {
        // Primary: more good fixed regions found
        if reverse.num_good_fixed_regions > forward.num_good_fixed_regions {
            return true;
        }
        if forward.num_good_fixed_regions > reverse.num_good_fixed_regions {
            return false;
        }
        // Secondary: higher total alignment score
        reverse.total_score > forward.total_score
    }

    /// Combine multiple extracted barcodes into final result.
    ///
    /// Groups by `region_name` to handle same-barcode regions appearing at both ends.
    /// For each group:
    /// - Single entry: pick first candidate from its tied-best list
    /// - Multiple entries (same barcode at both ends):
    ///   1. Intersect candidate sets
    ///   2. Non-empty intersection: pick first from intersection
    ///   3. Empty intersection: pick from the entry with highest confidence
    fn combine_barcodes(
        &self,
        barcodes: Vec<ExtractedBarcode>,
    ) -> Result<(LongReadBarcodeResult, Vec<(String, bool)>)> {
        if barcodes.is_empty() {
            return Ok((
                LongReadBarcodeResult {
                    barcode: None,
                    umi: None,
                    confidence: 0.0,
                    is_reverse_complemented: false,
                    five_prime_trim: 0,
                    three_prime_trim: 0,
                },
                Vec::new(),
            ));
        }

        // Group by region_name, preserving first-seen order
        let mut groups: IndexMap<String, Vec<&ExtractedBarcode>> = IndexMap::new();
        for bc in &barcodes {
            groups.entry(bc.region_name.clone()).or_default().push(bc);
        }

        // All required barcode groups must be resolved
        if !self
            .required_barcode_groups
            .iter()
            .all(|name| groups.contains_key(name))
        {
            return Ok((
                LongReadBarcodeResult {
                    barcode: None,
                    umi: None,
                    confidence: 0.0,
                    is_reverse_complemented: false,
                    five_prime_trim: 0,
                    three_prime_trim: 0,
                },
                Vec::new(),
            ));
        }

        // Resolve each group to a single barcode
        let mut combined_barcode = Vec::new();
        let mut total_confidence = 0.0;
        let mut num_resolved = 0usize;
        let mut consensus_results = Vec::new();

        for (region_name, entries) in &groups {
            let (chosen_barcode, chosen_confidence) = if entries.len() == 1 {
                // Single entry: pick first candidate
                (entries[0].barcodes[0].clone(), entries[0].confidence)
            } else {
                // Multiple entries: try intersection
                let (bc, conf, hit) = Self::resolve_multi_end_barcode(entries);
                consensus_results.push((region_name.clone(), hit));
                (bc, conf)
            };

            combined_barcode.extend(&chosen_barcode);
            total_confidence += chosen_confidence;
            num_resolved += 1;
        }

        let average_confidence = if num_resolved > 0 {
            total_confidence / num_resolved as f64
        } else {
            0.0
        };

        Ok((
            LongReadBarcodeResult {
                barcode: Some(combined_barcode),
                umi: None,
                confidence: average_confidence,
                is_reverse_complemented: false,
                five_prime_trim: 0,
                three_prime_trim: 0,
            },
            consensus_results,
        ))
    }

    /// Resolve a barcode region that appears at multiple ends.
    ///
    /// Strategy:
    /// 1. Intersect candidate sets from all entries
    /// 2. Non-empty intersection → pick first from intersection
    /// 3. Empty intersection → pick from entry with highest confidence;
    ///    if tied, pool all candidates and pick first
    fn resolve_multi_end_barcode(entries: &[&ExtractedBarcode]) -> (Vec<u8>, f64, bool) {
        use std::collections::HashSet;

        // Build intersection of candidate sets
        let mut intersection: HashSet<Vec<u8>> = entries[0].barcodes.iter().cloned().collect();
        for entry in &entries[1..] {
            let entry_set: HashSet<Vec<u8>> = entry.barcodes.iter().cloned().collect();
            intersection = intersection.intersection(&entry_set).cloned().collect();
        }

        if !intersection.is_empty() {
            // Pick the first one that appears in the highest-confidence entry's candidate list
            let best_entry = entries
                .iter()
                .max_by(|a, b| a.confidence.partial_cmp(&b.confidence).unwrap())
                .unwrap();
            let chosen = best_entry
                .barcodes
                .iter()
                .find(|bc| intersection.contains(*bc))
                .unwrap_or_else(|| intersection.iter().next().unwrap());
            return (chosen.clone(), best_entry.confidence, true);
        }

        // Empty intersection: pick from the highest-confidence entry
        let best_entry = entries
            .iter()
            .max_by(|a, b| a.confidence.partial_cmp(&b.confidence).unwrap())
            .unwrap();
        (best_entry.barcodes[0].clone(), best_entry.confidence, false)
    }
}

#[cfg(test)]
mod tests {
    use super::super::sequence_aligner::{fitting_alignment_distance, RegionMapping};
    use super::*;
    use indexmap::IndexMap;
    use seqspec::{RegionType, SequenceType};

    fn create_test_region(
        id: &str,
        region_type: RegionType,
        sequence_type: SequenceType,
        min_len: u32,
        max_len: u32,
    ) -> Arc<RwLock<Region>> {
        Arc::new(RwLock::new(Region {
            region_id: id.to_string(),
            region_type,
            name: id.to_string(),
            sequence_type,
            sequence: String::new(),
            min_len,
            max_len,
            onlist: None,
            subregions: vec![],
        }))
    }

    fn make_fixed_mapping(region: &Arc<RwLock<Region>>, start: usize, end: usize) -> RegionMapping {
        RegionMapping {
            region: region.clone(),
            read_start: start,
            read_end: end,
            is_spacer: false,
            match_rate: 0.95,
        }
    }

    fn make_spacer_mapping(
        region: &Arc<RwLock<Region>>,
        start: usize,
        end: usize,
    ) -> RegionMapping {
        RegionMapping {
            region: region.clone(),
            read_start: start,
            read_end: end,
            is_spacer: true,
            match_rate: 0.0,
        }
    }

    #[test]
    fn test_sample_and_analyze_end_window_keeps_forward_order() {
        let five_prime_barcode =
            create_test_region("bc5", RegionType::Barcode, SequenceType::Onlist, 4, 4);
        let three_prime_barcode =
            create_test_region("bc3", RegionType::Barcode, SequenceType::Onlist, 5, 5);

        let mut five_prime_regions = EndRegions::new(super::super::EndType::FivePrime);
        five_prime_regions.add_region(five_prime_barcode);

        let mut three_prime_regions = EndRegions::new(super::super::EndType::ThreePrime);
        three_prime_regions.add_region(three_prime_barcode);

        let extractor = BarcodeExtractor {
            five_prime_regions: five_prime_regions.clone(),
            three_prime_regions: three_prime_regions.clone(),
            five_prime_trim_regions: Vec::new(),
            three_prime_trim_regions: Vec::new(),
            five_prime_composite: CompositePattern {
                pattern: Vec::new(),
                spans: Vec::new(),
                pos_to_span: Vec::new(),
            },
            three_prime_composite: CompositePattern {
                pattern: Vec::new(),
                spans: Vec::new(),
                pos_to_span: Vec::new(),
            },
            whitelist_indices: IndexMap::new(),
            rc_regions: HashSet::new(),
            required_barcode_groups: IndexSet::new(),
            umi_region: None,
            umi_five_prime_anchor: None,
            umi_three_prime_anchor: None,
        };

        let sequence = b"AGTCAAAACCCCGTTA";

        let forward_5p = extractor
            .sample_and_analyze_end_window(
                sequence,
                &[],
                &five_prime_regions,
                &extractor.five_prime_composite,
                false,
            )
            .unwrap();
        let forward_3p = extractor
            .sample_and_analyze_end_window(
                sequence,
                &[],
                &three_prime_regions,
                &extractor.three_prime_composite,
                false,
            )
            .unwrap();
        let reverse_5p = extractor
            .sample_and_analyze_end_window(
                sequence,
                &[],
                &five_prime_regions,
                &extractor.five_prime_composite,
                true,
            )
            .unwrap();
        let reverse_3p = extractor
            .sample_and_analyze_end_window(
                sequence,
                &[],
                &three_prime_regions,
                &extractor.three_prime_composite,
                true,
            )
            .unwrap();

        // 5‘ Barcode length is 4bp, so the sampled window is ceil(4 * 1.15) = 5bp.
        // 3’ Barcode length is 5bp, so the sampled window is ceil(5 * 1.15) = 6bp.
        assert_eq!(forward_5p.sequence, b"AGTCA");
        assert_eq!(forward_3p.sequence, b"CCGTTA");
        assert_eq!(reverse_5p.sequence, b"TAACG");
        assert_eq!(reverse_3p.sequence, b"TTGACT");
    }

    #[test]
    fn test_calculate_five_prime_trim_from_anchor() {
        let primer = create_test_region(
            "primer",
            RegionType::IlluminaP5,
            SequenceType::Random,
            10,
            10,
        );
        let fixed = create_test_region("fixed5", RegionType::Linker, SequenceType::Fixed, 6, 6);
        let barcode = create_test_region("bc5", RegionType::Barcode, SequenceType::Onlist, 8, 8);
        let umi = create_test_region("umi5", RegionType::Umi, SequenceType::Random, 4, 4);

        let mut end_regions = EndRegions::new(super::super::EndType::FivePrime);
        end_regions.add_region(primer.clone());
        end_regions.add_region(fixed.clone());
        end_regions.add_region(barcode.clone());
        end_regions.add_region(umi.clone());

        let extractor = BarcodeExtractor {
            five_prime_regions: end_regions.clone(),
            three_prime_regions: EndRegions::new(super::super::EndType::ThreePrime),
            five_prime_trim_regions: vec![
                primer.clone(),
                fixed.clone(),
                barcode.clone(),
                umi.clone(),
            ],
            three_prime_trim_regions: Vec::new(),
            five_prime_composite: CompositePattern {
                pattern: Vec::new(),
                spans: Vec::new(),
                pos_to_span: Vec::new(),
            },
            three_prime_composite: CompositePattern {
                pattern: Vec::new(),
                spans: Vec::new(),
                pos_to_span: Vec::new(),
            },
            whitelist_indices: IndexMap::new(),
            rc_regions: HashSet::new(),
            required_barcode_groups: IndexSet::new(),
            umi_region: None,
            umi_five_prime_anchor: None,
            umi_three_prime_anchor: None,
        };
        let segment = EndSegmentWithAlignment {
            sequence: vec![b'A'; 40],
            composite_result: Some(CompositeAlignmentResult {
                score: 40,
                region_mappings: vec![
                    make_fixed_mapping(&fixed, 12, 18),
                    make_spacer_mapping(&barcode, 18, 25),
                    make_spacer_mapping(&umi, 25, 28),
                ],
            }),
        };

        assert_eq!(
            extractor.calculate_trim_length(
                &segment,
                &end_regions,
                &extractor.five_prime_trim_regions
            ),
            30
        );
    }

    #[test]
    fn test_calculate_three_prime_trim_from_anchor() {
        let umi = create_test_region("umi3", RegionType::Umi, SequenceType::Random, 4, 4);
        let barcode = create_test_region("bc3", RegionType::Barcode, SequenceType::Onlist, 8, 8);
        let fixed = create_test_region("fixed3", RegionType::Linker, SequenceType::Fixed, 6, 6);
        let adapter = create_test_region("adapter3", RegionType::Umi, SequenceType::Random, 10, 10);

        // The length of sampled end segment is (8 + 6 + 10) * 1.15 = 27.6 -> 28bp.
        let mut end_regions = EndRegions::new(super::super::EndType::ThreePrime);
        end_regions.add_region(barcode.clone());
        end_regions.add_region(fixed.clone());
        end_regions.add_region(adapter.clone());

        let extractor = BarcodeExtractor {
            five_prime_regions: EndRegions::new(super::super::EndType::FivePrime),
            three_prime_regions: end_regions.clone(),
            five_prime_trim_regions: Vec::new(),
            three_prime_trim_regions: vec![
                umi.clone(),
                barcode.clone(),
                fixed.clone(),
                adapter.clone(),
            ],
            five_prime_composite: CompositePattern {
                pattern: Vec::new(),
                spans: Vec::new(),
                pos_to_span: Vec::new(),
            },
            three_prime_composite: CompositePattern {
                pattern: Vec::new(),
                spans: Vec::new(),
                pos_to_span: Vec::new(),
            },
            whitelist_indices: IndexMap::new(),
            rc_regions: HashSet::new(),
            required_barcode_groups: IndexSet::new(),
            umi_region: None,
            umi_five_prime_anchor: None,
            umi_three_prime_anchor: None,
        };
        let segment = EndSegmentWithAlignment {
            sequence: vec![b'A'; 28],
            composite_result: Some(CompositeAlignmentResult {
                score: 40,
                region_mappings: vec![
                    make_spacer_mapping(&barcode, 3, 12),
                    make_fixed_mapping(&fixed, 12, 18),
                    make_spacer_mapping(&adapter, 18, 28),
                ],
            }),
        };

        assert_eq!(
            extractor.calculate_trim_length(
                &segment,
                &end_regions,
                &extractor.three_prime_trim_regions
            ),
            28
        );
    }

    #[test]
    fn test_calculate_trim_falls_back_to_average_lengths() {
        let primer = create_test_region(
            "primer",
            RegionType::IlluminaP5,
            SequenceType::Random,
            10,
            10,
        );
        let fixed = create_test_region("fixed5", RegionType::Linker, SequenceType::Fixed, 6, 6);
        let barcode = create_test_region("bc5", RegionType::Barcode, SequenceType::Onlist, 8, 8);

        let mut end_regions = EndRegions::new(super::super::EndType::FivePrime);
        end_regions.add_region(primer.clone());
        end_regions.add_region(fixed.clone());
        end_regions.add_region(barcode.clone());

        let extractor = BarcodeExtractor {
            five_prime_regions: end_regions.clone(),
            three_prime_regions: EndRegions::new(super::super::EndType::ThreePrime),
            five_prime_trim_regions: vec![primer.clone(), fixed.clone(), barcode.clone()],
            three_prime_trim_regions: Vec::new(),
            five_prime_composite: CompositePattern {
                pattern: Vec::new(),
                spans: Vec::new(),
                pos_to_span: Vec::new(),
            },
            three_prime_composite: CompositePattern {
                pattern: Vec::new(),
                spans: Vec::new(),
                pos_to_span: Vec::new(),
            },
            whitelist_indices: IndexMap::new(),
            rc_regions: HashSet::new(),
            required_barcode_groups: IndexSet::new(),
            umi_region: None,
            umi_five_prime_anchor: None,
            umi_three_prime_anchor: None,
        };
        let segment = EndSegmentWithAlignment {
            sequence: vec![b'A'; 20],
            composite_result: Some(CompositeAlignmentResult {
                score: 10,
                region_mappings: vec![make_spacer_mapping(&barcode, 8, 16)],
            }),
        };

        assert_eq!(
            extractor.calculate_trim_length(
                &segment,
                &end_regions,
                &extractor.five_prime_trim_regions
            ),
            24
        );
    }

    #[test]
    fn test_locate_barcode_case1_left_adjacent() {
        // Layout: fixed1(16bp) -> barcode(10bp) -> fixed2(16bp)
        // Case 1: left neighbor (fixed1) is structurally adjacent
        // Window: [fixed1.read_end, fixed1.read_end + 1.2*10], capped by fixed2.read_start
        let fixed1 = create_test_region("fixed1", RegionType::Linker, SequenceType::Fixed, 16, 16);
        let barcode = create_test_region("bc1", RegionType::Barcode, SequenceType::Onlist, 10, 10);
        let fixed2 = create_test_region("fixed2", RegionType::Linker, SequenceType::Fixed, 16, 16);

        let mut end_regions = EndRegions::new(super::super::EndType::FivePrime);
        end_regions.add_region(fixed1.clone());
        end_regions.add_region(barcode.clone());
        end_regions.add_region(fixed2.clone());

        let composite_result = CompositeAlignmentResult {
            score: 50,
            region_mappings: vec![
                make_fixed_mapping(&fixed1, 0, 17),
                make_spacer_mapping(&barcode, 17, 26),
                make_fixed_mapping(&fixed2, 26, 42),
            ],
        };

        let window =
            locate_barcode_extraction_window(&barcode, &end_regions, &composite_result, 100);
        let (start, end) = window.unwrap();

        // start = fixed1.read_end = 17
        // max_end = 17 + ceil(10*1.2) = 17 + 12 = 29, capped by fixed2.read_start = 26
        assert_eq!(start, 17);
        assert_eq!(end, 26);
    }

    #[test]
    fn test_locate_barcode_case2_right_adjacent() {
        // Layout: random -> barcode(10bp) -> fixed1(16bp)
        // Case 2: right neighbor (fixed1) is structurally adjacent
        // Window: [fixed1.read_start - 1.2*10, fixed1.read_start]
        let random = create_test_region(
            "random",
            RegionType::IlluminaP5,
            SequenceType::Random,
            20,
            20,
        );
        let barcode = create_test_region("bc1", RegionType::Barcode, SequenceType::Onlist, 10, 10);
        let fixed1 = create_test_region("fixed1", RegionType::Linker, SequenceType::Fixed, 16, 16);

        let mut end_regions = EndRegions::new(super::super::EndType::FivePrime);
        end_regions.add_region(random.clone());
        end_regions.add_region(barcode.clone());
        end_regions.add_region(fixed1.clone());

        // Composite only contains fixed1 (single fixed region)
        let composite_result = CompositeAlignmentResult {
            score: 30,
            region_mappings: vec![make_fixed_mapping(&fixed1, 30, 46)],
        };

        let window =
            locate_barcode_extraction_window(&barcode, &end_regions, &composite_result, 100);
        let (start, end) = window.unwrap();

        // end = fixed1.read_start = 30
        // min_start = 30 - 12 = 18, no left fixed to cap
        assert_eq!(start, 18);
        assert_eq!(end, 30);
    }

    #[test]
    fn test_locate_barcode_case3_sandwiched_non_adjacent() {
        // Layout: fixed1(16bp) -> umi(4bp) -> barcode(10bp) -> umi2(4bp) -> fixed2(16bp)
        // Case 3: neither immediate neighbor is fixed, but both sides have fixed regions
        // Window: entire gap [fixed1.read_end, fixed2.read_start]
        let fixed1 = create_test_region("fixed1", RegionType::Linker, SequenceType::Fixed, 16, 16);
        let umi1 = create_test_region("umi1", RegionType::Umi, SequenceType::Random, 4, 4);
        let barcode = create_test_region("bc1", RegionType::Barcode, SequenceType::Onlist, 10, 10);
        let umi2 = create_test_region("umi2", RegionType::Umi, SequenceType::Random, 4, 4);
        let fixed2 = create_test_region("fixed2", RegionType::Linker, SequenceType::Fixed, 16, 16);

        let mut end_regions = EndRegions::new(super::super::EndType::FivePrime);
        end_regions.add_region(fixed1.clone());
        end_regions.add_region(umi1.clone());
        end_regions.add_region(barcode.clone());
        end_regions.add_region(umi2.clone());
        end_regions.add_region(fixed2.clone());

        let composite_result = CompositeAlignmentResult {
            score: 50,
            region_mappings: vec![
                make_fixed_mapping(&fixed1, 0, 16),
                make_spacer_mapping(&umi1, 16, 20),
                make_spacer_mapping(&barcode, 20, 30),
                make_spacer_mapping(&umi2, 30, 34),
                make_fixed_mapping(&fixed2, 34, 50),
            ],
        };

        let window =
            locate_barcode_extraction_window(&barcode, &end_regions, &composite_result, 100);
        let (start, end) = window.unwrap();

        // Entire gap between fixed regions
        assert_eq!(start, 16); // fixed1.read_end
        assert_eq!(end, 34); // fixed2.read_start (note: NOT 50, which is read_end)
    }

    #[test]
    fn test_locate_barcode_case4_terminal_edge() {
        // Layout: barcode(10bp) -> umi(4bp) -> fixed1(16bp)
        // Case 4: barcode at left edge, no fixed region on the left side
        // Window: [0, fixed1.read_start] (but fixed1 is not adjacent, umi in between)
        let barcode = create_test_region("bc1", RegionType::Barcode, SequenceType::Onlist, 10, 10);
        let umi = create_test_region("umi1", RegionType::Umi, SequenceType::Random, 4, 4);
        let fixed1 = create_test_region("fixed1", RegionType::Linker, SequenceType::Fixed, 16, 16);

        let mut end_regions = EndRegions::new(super::super::EndType::FivePrime);
        end_regions.add_region(barcode.clone());
        end_regions.add_region(umi.clone());
        end_regions.add_region(fixed1.clone());

        let composite_result = CompositeAlignmentResult {
            score: 30,
            region_mappings: vec![make_fixed_mapping(&fixed1, 14, 30)],
        };

        let window =
            locate_barcode_extraction_window(&barcode, &end_regions, &composite_result, 50);
        let (start, end) = window.unwrap();

        // Terminal: from segment start (0) to nearest right fixed region's read_start
        assert_eq!(start, 0);
        assert_eq!(end, 14); // fixed1.read_start
    }

    #[test]
    fn test_locate_barcode_window_too_small() {
        // Layout: fixed1(16bp) -> barcode(10bp) -> fixed2(16bp)
        // But fixed regions are so close that window < 0.8 * barcode_len
        let fixed1 = create_test_region("fixed1", RegionType::Linker, SequenceType::Fixed, 16, 16);
        let barcode = create_test_region("bc1", RegionType::Barcode, SequenceType::Onlist, 10, 10);
        let fixed2 = create_test_region("fixed2", RegionType::Linker, SequenceType::Fixed, 16, 16);

        let mut end_regions = EndRegions::new(super::super::EndType::FivePrime);
        end_regions.add_region(fixed1.clone());
        end_regions.add_region(barcode.clone());
        end_regions.add_region(fixed2.clone());

        // Gap between fixed regions is only 5bp < floor(10*0.8) = 8bp
        let composite_result = CompositeAlignmentResult {
            score: 50,
            region_mappings: vec![
                make_fixed_mapping(&fixed1, 0, 16),
                make_spacer_mapping(&barcode, 16, 21),
                make_fixed_mapping(&fixed2, 21, 37),
            ],
        };

        let window =
            locate_barcode_extraction_window(&barcode, &end_regions, &composite_result, 100);
        assert!(window.is_none(), "Window should be rejected when too small");
    }

    #[test]
    fn test_orientation_evidence_thresholds() {
        // Case 1: Meets all thresholds
        let evidence = OrientationEvidence {
            total_score: 100,
            avg_fixed_match_rate: 0.90,
            num_good_fixed_regions: 2,
            total_fixed_regions: 2,
        };
        assert!(evidence.meets_threshold());

        // Case 2: Below match rate threshold
        let evidence = OrientationEvidence {
            total_score: 50,
            avg_fixed_match_rate: 0.60,
            num_good_fixed_regions: 0,
            total_fixed_regions: 2,
        };
        assert!(!evidence.meets_threshold());

        // Case 3: No good fixed regions
        let evidence = OrientationEvidence {
            total_score: 30,
            avg_fixed_match_rate: 0.50,
            num_good_fixed_regions: 0,
            total_fixed_regions: 1,
        };
        assert!(!evidence.meets_threshold());

        // Case 4: Edge case - exactly at threshold
        let evidence = OrientationEvidence {
            total_score: 80,
            avg_fixed_match_rate: 0.80,
            num_good_fixed_regions: 2,
            total_fixed_regions: 2,
        };
        assert!(evidence.meets_threshold());
    }

    #[test]
    fn test_fitting_alignment_distance() {
        // Test with exact matches
        assert_eq!(fitting_alignment_distance(b"ATCG", b"ATCG"), 0);

        // Test with partial matches (short sequence found in long sequence)
        assert_eq!(fitting_alignment_distance(b"ATC", b"ATCA"), 0);
        assert_eq!(fitting_alignment_distance(b"ATCG", b"GGATCGCC"), 0);

        // Test with empty sequences
        assert_eq!(fitting_alignment_distance(b"", b"ATCG"), 0);
        assert_eq!(fitting_alignment_distance(b"ATCG", b""), 4);
    }

    #[test]
    fn test_reverse_complement_barcode_span_uses_normalized_window_coordinates() {
        let barcode =
            create_test_region("barcode", RegionType::Barcode, SequenceType::Onlist, 4, 4);
        let whitelist = IndexSet::from_iter([b"AAGC".to_vec()]);
        let indices = IndexMap::from_iter([("barcode".to_string(), BarcodeIndex::new(&whitelist))]);
        let normalized_window = b"GCTTGCTT";
        let matching_window = seqspec::utils::rev_compl(normalized_window);

        let matched = find_best_barcode_match(
            &matching_window,
            normalized_window,
            &barcode,
            &indices,
            true,
        )
        .unwrap();

        assert_eq!(matched.barcodes[0], b"AAGC");
        assert_eq!(matched.traced_span, 4..8);
    }

    #[test]
    fn test_extract_umi_from_barcode_anchor_in_both_orientations() {
        use noodles_fastq::record::Definition;

        let fixed = create_test_region("fixed", RegionType::Linker, SequenceType::Fixed, 16, 16);
        fixed.write().unwrap().sequence = "ACGTACGTACGTACGT".to_string();
        let barcode =
            create_test_region("barcode", RegionType::Barcode, SequenceType::Onlist, 10, 10);
        let umi = create_test_region("umi", RegionType::Umi, SequenceType::Random, 6, 6);
        let target = create_test_region("target", RegionType::Cdna, SequenceType::Random, 50, 500);
        let modality_region = Region {
            region_id: "rna".to_string(),
            region_type: RegionType::Modality(Modality::RNA),
            name: "RNA".to_string(),
            sequence_type: SequenceType::Joined,
            sequence: String::new(),
            min_len: 0,
            max_len: 0,
            onlist: None,
            subregions: vec![fixed, barcode, umi, target],
        };
        let lib_spec = LibSpec::new(vec![modality_region]).unwrap();
        let whitelist = IndexSet::from_iter([b"AACCGGTTAA".to_vec()]);
        let extractor = BarcodeExtractor::new(
            &lib_spec,
            &Modality::RNA,
            IndexMap::from_iter([("barcode".to_string(), whitelist)]),
        )
        .unwrap();
        assert!(extractor.expects_umi());

        let sequence = [
            b"ACGTACGTACGTACGT".as_slice(),
            b"AACCGGTTAA".as_slice(),
            b"aagtcc".as_slice(),
            vec![b'G'; 80].as_slice(),
        ]
        .concat();
        let mut quality = vec![b'I'; sequence.len()];
        quality[26..32].copy_from_slice(b"#$%&'(");
        let forward = noodles_fastq::Record::new(
            Definition::new("forward", ""),
            sequence.clone(),
            quality.clone(),
        );
        let (result, _) = extractor.extract_barcode(&forward).unwrap();
        assert!(!result.is_reverse_complemented);
        let extracted_umi = result.umi.unwrap();
        assert_eq!(extracted_umi.sequence(), b"AAGTCC");
        assert_eq!(extracted_umi.quality_scores(), b"#$%&'(");

        let uppercase_sequence = sequence.to_ascii_uppercase();
        let reverse = noodles_fastq::Record::new(
            Definition::new("reverse", ""),
            seqspec::utils::rev_compl(&uppercase_sequence),
            quality.iter().rev().copied().collect::<Vec<_>>(),
        );
        let (result, _) = extractor.extract_barcode(&reverse).unwrap();
        assert!(result.is_reverse_complemented);
        let extracted_umi = result.umi.unwrap();
        assert_eq!(extracted_umi.sequence(), b"AAGTCC");
        assert_eq!(extracted_umi.quality_scores(), b"#$%&'(");
    }

    #[test]
    fn test_umi_anchor_prefers_fixed_and_falls_back_to_barcode() {
        use noodles_fastq::record::Definition;

        let outer_fixed =
            create_test_region("outer", RegionType::Linker, SequenceType::Fixed, 16, 16);
        outer_fixed.write().unwrap().sequence = "ACGTACGTACGTACGT".to_string();
        let barcode =
            create_test_region("barcode", RegionType::Barcode, SequenceType::Onlist, 10, 10);
        let umi = create_test_region("umi", RegionType::Umi, SequenceType::Random, 4, 4);
        let inner_fixed =
            create_test_region("inner", RegionType::Linker, SequenceType::Fixed, 16, 16);
        inner_fixed.write().unwrap().sequence = "TGCATGCATGCATGCA".to_string();
        let target = create_test_region("target", RegionType::Cdna, SequenceType::Random, 50, 500);
        let modality_region = Region {
            region_id: "rna".to_string(),
            region_type: RegionType::Modality(Modality::RNA),
            name: "RNA".to_string(),
            sequence_type: SequenceType::Joined,
            sequence: String::new(),
            min_len: 0,
            max_len: 0,
            onlist: None,
            subregions: vec![
                outer_fixed.clone(),
                barcode,
                umi,
                inner_fixed.clone(),
                target,
            ],
        };
        let lib_spec = LibSpec::new(vec![modality_region]).unwrap();
        let whitelist = IndexSet::from_iter([b"AACCGGTTAA".to_vec()]);
        let extractor = BarcodeExtractor::new(
            &lib_spec,
            &Modality::RNA,
            IndexMap::from_iter([("barcode".to_string(), whitelist)]),
        )
        .unwrap();
        let sequence = [
            b"ACGTACGTACGTACGT".as_slice(),
            b"AACCGGTTAA".as_slice(),
            b"TGCA".as_slice(),
            b"TGCATGCATGCATGCA".as_slice(),
            vec![b'G'; 50].as_slice(),
        ]
        .concat();
        let record = noodles_fastq::Record::new(
            Definition::new("read", ""),
            sequence.clone(),
            vec![b'I'; sequence.len()],
        );
        let mut barcode = ExtractedBarcode {
            region_id: "barcode".to_string(),
            region_name: "barcode".to_string(),
            barcodes: vec![b"AACCGGTTAA".to_vec()],
            confidence: 1.0,
            // Deliberately offset: the fixed anchor must win while observed.
            traced_span: 15..25,
        };
        let mut five_prime = EndSegmentWithAlignment {
            sequence: sequence[..50].to_vec(),
            composite_result: Some(CompositeAlignmentResult {
                score: 50,
                region_mappings: vec![
                    make_fixed_mapping(&outer_fixed, 0, 16),
                    make_fixed_mapping(&inner_fixed, 30, 46),
                ],
            }),
        };
        let three_prime = EndSegmentWithAlignment::empty();

        let extracted = extractor
            .extract_umi(
                &record,
                &five_prime,
                &three_prime,
                &[barcode.clone()],
                false,
            )
            .unwrap();
        assert_eq!(extracted.sequence(), b"TGCA");

        // With the preferred 3' fixed anchor unavailable on this read, extraction
        // falls back to the observed 5' barcode span.
        five_prime
            .composite_result
            .as_mut()
            .unwrap()
            .region_mappings
            .retain(|mapping| !Arc::ptr_eq(&mapping.region, &inner_fixed));
        barcode.traced_span = 16..26;
        let extracted = extractor
            .extract_umi(&record, &five_prime, &three_prime, &[barcode], false)
            .unwrap();
        assert_eq!(extracted.sequence(), b"TGCA");
    }

    #[test]
    fn test_unsupported_umi_topology_disables_extraction() {
        let fixed = create_test_region("fixed", RegionType::Linker, SequenceType::Fixed, 16, 16);
        fixed.write().unwrap().sequence = "ACGTACGTACGTACGT".to_string();
        let barcode =
            create_test_region("barcode", RegionType::Barcode, SequenceType::Onlist, 10, 10);
        let spacer = create_test_region("spacer", RegionType::Named, SequenceType::Random, 4, 4);
        let umi = create_test_region("umi", RegionType::Umi, SequenceType::Random, 4, 4);
        let target = create_test_region("target", RegionType::Cdna, SequenceType::Random, 50, 500);
        let modality_region = Region {
            region_id: "rna".to_string(),
            region_type: RegionType::Modality(Modality::RNA),
            name: "RNA".to_string(),
            sequence_type: SequenceType::Joined,
            sequence: String::new(),
            min_len: 0,
            max_len: 0,
            onlist: None,
            subregions: vec![fixed, barcode, spacer, umi, target],
        };
        let lib_spec = LibSpec::new(vec![modality_region]).unwrap();
        let whitelist = IndexSet::from_iter([b"AACCGGTTAA".to_vec()]);
        let extractor = BarcodeExtractor::new(
            &lib_spec,
            &Modality::RNA,
            IndexMap::from_iter([("barcode".to_string(), whitelist)]),
        )
        .unwrap();

        assert!(!extractor.expects_umi());
        assert!(extractor.umi_region.is_some());
    }

    #[test]
    fn test_multiple_umi_regions_are_rejected() {
        let fixed = create_test_region("fixed", RegionType::Linker, SequenceType::Fixed, 16, 16);
        fixed.write().unwrap().sequence = "ACGTACGTACGTACGT".to_string();
        let barcode =
            create_test_region("barcode", RegionType::Barcode, SequenceType::Onlist, 10, 10);
        let umi1 = create_test_region("umi1", RegionType::Umi, SequenceType::Random, 4, 4);
        let umi2 = create_test_region("umi2", RegionType::Umi, SequenceType::Random, 4, 4);
        let target = create_test_region("target", RegionType::Cdna, SequenceType::Random, 50, 500);
        let modality_region = Region {
            region_id: "rna".to_string(),
            region_type: RegionType::Modality(Modality::RNA),
            name: "RNA".to_string(),
            sequence_type: SequenceType::Joined,
            sequence: String::new(),
            min_len: 0,
            max_len: 0,
            onlist: None,
            subregions: vec![fixed, barcode, umi1, umi2, target],
        };
        let lib_spec = LibSpec::new(vec![modality_region]).unwrap();
        let whitelist = IndexSet::from_iter([b"AACCGGTTAA".to_vec()]);
        let error = BarcodeExtractor::new(
            &lib_spec,
            &Modality::RNA,
            IndexMap::from_iter([("barcode".to_string(), whitelist)]),
        )
        .unwrap_err();

        assert!(error.to_string().contains("multiple UMI regions"));
    }
}
