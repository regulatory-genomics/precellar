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

## End-Window Sampling And Coordinate Handling

The orientation logic uses the same bounded windows that are later reused for barcode extraction.

- Window length for each end is `ceil(sum(max_len of end regions) * 1.15)`.
- End windows are normalized into forward order of the designed library structure before alignment.
- In that coordinate system, `left = 5' side` and `right = 3' side`.
- If `should_rc` is true, only the bounded end window is reverse-complemented.

This keeps orientation detection localized and allows the same cached anchors to be reused immediately by barcode extraction.

## Reuse For Trimming

The same composite-alignment results are reused to compute precise trimming boundaries, so that orientation detection, barcode extraction, and trimming all share a single end-region interpretation.

### Trim Regions vs End Regions

Trimming uses a broader set of flanking regions than barcode extraction. The extractor collects two region sets from the library specification:

- **End regions** contain only the regions from the first (or last) fixed/barcode anchor to the target boundary, used for composite alignment and barcode extraction.
- **Trim regions** contain all non-target regions on each side of the target, including any outer regions (e.g., primers, adapters) that fall outside the end-region window.

This distinction matters because trimming must remove the entire non-target portion of the read, not just the barcode-adjacent portion.

### Anchor-Based Precise Trimming

For each end, the extractor computes the trim length as follows:

1. Find the **trim anchor** from the end regions:
   - For the 5' end: the last (rightmost) fixed region in the end-region collection.
   - For the 3' end: the first (leftmost) fixed region in the end-region collection.
2. Look up the anchor's aligned position in the composite alignment result.
3. Compute a **base length** from the anchor position:
   - For the 5' end: the anchor's aligned `read_end`.
   - For the 3' end: `segment_length - anchor.read_start`.
4. Compute an **offset**: sum of average lengths (`(min_len + max_len) / 2`) of trim regions that lie between the anchor and the target.
5. Total trim length = base + offset.

### Fallback

If no trim anchor is found in the composite alignment (e.g., because the anchor region was not matched, or no composite alignment was performed for that end), the extractor falls back to the sum of average lengths of all trim regions for that end.

### Orientation-Aware Application

Trim lengths are reported as `five_prime_trim` and `three_prime_trim` in designed-library orientation. When applying trims to the raw read:

- If the read is in forward orientation: `head_trim = five_prime_trim`, `tail_trim = three_prime_trim`.
- If the read is reverse-complemented: `head_trim = three_prime_trim`, `tail_trim = five_prime_trim`.

The target sequence is then `original_seq[head_trim .. len - tail_trim]`, without reverse-complementing the raw read itself.
