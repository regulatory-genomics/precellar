# Long-Read Strand Orientation

This page describes the current record-wise strand orientation logic used during long-read barcode extraction.

## Why Orientation Is Record-Wise

Long-read assays cannot assume that every FASTQ record arrives in the same orientation as the designed library structure. In `seqspec`, `Read.strand` describes whether a FASTQ read file is expected to have the same orientation as the library structure or the opposite orientation. For long-read preprocessing, `Read.strand` is typically set to `unstranded`, so the YAML-defined library structure is treated as the reference layout, while the implementation determines the orientation of each FASTQ record independently rather than assigning one global orientation to the whole read set.

## Core Principle

Orientation is inferred from fixed-region evidence near the read ends, not from transcript alignment.

This record-wise orientation step is independent of `ChemistryStrandedness`. `ChemistryStrandedness` indicates whether library preparation preserves the strand specificity of the original RNA molecule and, for paired-end data, defines the mapping geometry of the read pair. By contrast, the long-read logic described here only decides whether an individual FASTQ record should be interpreted in forward or reverse-complement orientation relative to the YAML-defined end-region layout for barcode extraction.

The algorithm compares two hypotheses for the same record:

1. forward orientation relative to the YAML-defined library layout
2. reverse-complement orientation relative to that same layout

Only end segments are used for this comparison. The full read is not reverse-complemented just to decide orientation.

## Core Algorithm

The current implementation uses an adaptive forward-first workflow.

1. Cut the 5' and 3' end windows in forward orientation.
2. Run composite alignment on those windows.
3. Summarize the resulting fixed-region evidence.
4. If the forward evidence already meets the early-accept threshold, accept forward immediately.
5. Otherwise, repeat the same process on reverse-complemented end segments.
6. Compare forward and reverse evidence and choose the better orientation.
7. Gate barcode extraction on the chosen orientation's alignment quality.

This avoids the cost of evaluating both orientations on every read when forward evidence is already unambiguous.

The comparison between forward and reverse is based on `OrientationEvidence` derived from composite-alignment results.

For both ends together, the evidence tracks:

- total composite alignment score
- average fixed-region match rate
- number of fixed regions with match rate at least `0.8`
- total number of fixed regions observed

Only non-spacer regions contribute to fixed-region evidence.

Forward orientation is accepted immediately when the evidence meets a conservative threshold:

- at least 2 good fixed anchors
- average fixed-region match rate at least `0.8`

If forward orientation does not meet the early threshold, the extractor computes reverse evidence and then compares the two candidates.

Selection order is:

1. Prefer the orientation with more good fixed regions.
2. If tied, prefer the orientation with the higher total composite score.

The chosen orientation is recorded as `is_reverse_complemented` in the final long-read barcode result.

Choosing an orientation is not enough by itself. The chosen orientation must still satisfy a minimum composite quality requirement before barcode extraction is allowed to continue.

- minimum average fixed-region match rate: `0.7`

If the chosen orientation falls below this threshold:

- no barcode is returned
- orientation metadata is still preserved
- composite-alignment results are still returned when available

This prevents weak anchor evidence from driving spurious barcode extraction.

## End-Segment Handling

The orientation logic operates on the same bounded end segments used by barcode extraction.

- The segment length for each end is `alignment_window_length()` from the corresponding `EndRegions` collection.
- For reverse orientation, only the selected end segment is reverse-complemented.
- For the 3' end in non-reverse-complement mode, the segment is reversed before alignment so the segment ordering matches the composite pattern's outer-to-inner layout.
- After alignment, region coordinates are converted back into the original segment coordinate system.

This design keeps orientation detection localized and allows the same cached anchors to be reused immediately by barcode extraction.

## Relationship To Trimming

The same composite-alignment results used for orientation are also reused to estimate trimming boundaries.

- For the 5' end, trimming anchors to the innermost fixed region's `read_end`.
- For the 3' end, trimming anchors to the innermost fixed region's `read_start` after coordinate normalization.
- Heuristic trailing-region length and intermediate-gap estimates are added on top of the aligned anchor position.

This reuse keeps orientation, barcode extraction, and target trimming consistent with the same end-region interpretation.
