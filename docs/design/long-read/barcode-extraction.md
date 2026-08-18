# Long-Read Barcode Extraction

This page describes the current long-read barcode extraction design implemented in `precellar-core/src/long/`.

## Scope

The long-read barcode workflow is designed for assays where barcode regions are embedded near one or both read ends and must be recovered from noisy long-read sequence. The algorithm uses library structure from `seqspec`, composite alignment against fixed regions, and whitelist-guided barcode resolution.

## Inputs

- A `seqspec` library specification for the selected modality
- Per-barcode-region whitelists
- A long-read FASTQ record

## End Region Discovery

The extractor first derives a 5' and 3' `EndRegions` collection from the library specification.

- Region order always follows the forward order of the designed library structure (`5' -> 3'`).
- In end-segment coordinates, `left = 5' side` and `right = 3' side`.
- The 5' collection keeps all regions from the library 5' side through the last fixed/barcode before the target.
- The 3' collection keeps all regions from the first fixed/barcode after the target through the library 3' side.
- Ends without barcode regions are skipped by the extractor.

This means the barcode workflow only reasons over the structured end regions that can anchor extraction, not the full library layout.

## Trim Region Discovery

In addition to end regions, the extractor collects **trim regions** for each side of the target. Trim regions include all non-target regions on each side, not just the anchor-adjacent subset used for barcode extraction. For example, given a layout:

```
primer -> barcode -> linker -> cDNA(target) -> UMI -> adapter
```

- 5' end regions: `[barcode, linker]` (from first anchor to target)
- 5' trim regions: `[primer, barcode, linker]` (all non-target regions on the 5' side)
- 3' end regions: `[adapter]` (from first anchor after target to end, if adapter is fixed)
- 3' trim regions: `[UMI, adapter]` (all non-target regions on the 3' side)

The trim regions are used to compute how much sequence to remove from each end of the raw read before alignment. The trimming algorithm and its fallback behavior are described in [Strand Orientation — Reuse For Trimming](strand-orientation.md#reuse-for-trimming).

## Composite Alignment Pattern

Each end builds one reusable composite pattern.

- Fixed regions contribute their real sequence.
- Non-fixed regions between fixed anchors contribute `N` spacers.
- Spacer length is the average of `min_len` and `max_len` for that region.
- The total fixed sequence length must be at least 12 bp for any end that contains a barcode.

Composite alignment uses one fitting alignment instead of aligning fixed regions independently. This reduces the need to reconcile multiple local alignments and preserves the expected topological order of anchor regions.

The composite scoring model is:

- fixed-base match: `+2`
- fixed-base mismatch: `-1`
- `N` spacer match to any base: `+1`
- gap open: `-1`
- gap extend: `-1`

This scoring intentionally tolerates indels and variable-length regions that are common in long-read data.

## End-Window Sampling

Barcode extraction never aligns the full read for anchor detection.

- Each end cuts only a bounded segment.
- Segment length is `sum(max_len of end regions) * 1.15`, rounded up.
- The 15% buffer allows moderate indel drift during long-read sequencing.

Each end segment is normalized into the same forward order as its `EndRegions` collection before composite alignment. The extractor only reverse-complements bounded end windows when evaluating the reverse-orientation hypothesis.

## Barcode Extraction Workflow

Once the chosen orientation is known, barcode extraction proceeds per end.

1. Reuse the cached composite alignment for that end.
2. Enumerate barcode regions in the `EndRegions` collection.
3. Locate an extraction window for each barcode region using aligned fixed-region boundaries.
4. Match the extracted candidate sequence against the whitelist index.
5. Combine resolved barcode pieces across ends into one final barcode.

The barcode extractor only proceeds when the chosen orientation has average fixed-region match rate at least `0.7`.

## Extraction Window Hierarchy

Barcode windows are not taken from fixed offsets. They are derived from alignment anchors using a strict topological hierarchy.

### 1. Left-Adjacent Anchoring

If the barcode's immediate left neighbor is a fixed region:

- start at that fixed region's aligned `read_end`
- extend right by up to `1.2 * barcode_length`
- cap the end at the next right fixed region if present

### 2. Right-Adjacent Anchoring

If the barcode's immediate right neighbor is a fixed region:

- end at that fixed region's aligned `read_start`
- extend left by up to `1.2 * barcode_length`
- cap the start at the nearest left fixed region if present

### 3. Sandwiched Non-Adjacent Anchoring

If neither immediate neighbor is fixed, but fixed regions exist on both sides:

- use the full gap between the nearest left and right fixed anchors

### 4. Terminal Edge Extraction

If the barcode sits at the segment boundary:

- extract between the segment edge and the nearest fixed anchor

### Length Gate

The candidate window must still be long enough to be plausible.

- target window length is `ceil(1.2 * barcode_length)`
- minimum accepted window length is `floor(0.8 * barcode_length)`

Windows shorter than the minimum are rejected.

## Whitelist Matching

Barcode identification uses a two-stage search.

### Stage 1: K-mer Voting

Each whitelist is indexed once using overlapping 6-mers.

- key: encoded 6-mer
- value: barcode IDs containing that 6-mer

For each extracted candidate sequence:

- generate overlapping 6-mers
- count votes for whitelist barcodes sharing those 6-mers
- keep the top 50 candidates, including ties at the cutoff

If the whitelist has 50 or fewer barcodes, or the candidate is shorter than 6 bp, the algorithm skips voting and falls back to direct comparison.

### Stage 2: Fitting-Alignment Distance

The extractor then computes fitting-alignment distance only on the retained candidates.

- keep all tied-best barcodes
- compute confidence as `1 - min_edit_distance / barcode_length`
- require confidence at least `0.7`

This returns an `ExtractedBarcode` containing the region ID, region name, tied-best barcode candidates, and confidence.

## Multi-End Barcode Resolution

Barcode extraction does not immediately return a single resolved barcode for each physical barcode region. Instead, each physical barcode region produces an entry containing:

- its unique `region_id`
- its logical barcode label `region_name`
- a tied-best candidate barcode list
- a confidence score for that entry

The final barcode may contain contributions from multiple barcode regions, including the same logical barcode observed at both ends. To resolve these into one output barcode, the extractor performs a consensus-style resolution across entries that share the same logical barcode name.

Resolved barcode entries are grouped by `region_name`.

- Single entry: choose the first tied-best candidate.
- Multiple entries for the same logical barcode:
  - intersect candidate sets
  - if the intersection is non-empty, choose a consensus barcode from the shared candidates
  - otherwise choose the first candidate from the highest-confidence entry

The final output barcode is the concatenation of the resolved groups in first-seen order. Final confidence is the average confidence across resolved groups.

## Failure Modes

Long-read barcode extraction returns no barcode when any of the following blocks the workflow:

- the read is shorter than the required end-segment length, computed as `sum(max_len of end regions) * 1.15`; the extra 15% buffer allows moderate indel drift in long-read sequencing
- the chosen orientation has insufficient composite match quality, meaning the selected orientation's average fixed-region match rate is below `0.7`
- no valid extraction window can be located, for example because no usable anchor topology is found or the candidate window is shorter than `floor(0.8 * barcode_length)`
- whitelist matching fails, meaning no tied-best barcode candidate reaches confidence `0.7`, where confidence is `1 - min_edit_distance / barcode_length`

In those cases, the extractor still preserves orientation metadata and computes trim lengths when composite-alignment results are available. Even when no barcode is returned, the trim values allow the pipeline to remove non-target regions from the read. For read-length failures, an end may remain empty because no composite alignment is attempted for that end; the trim length then falls back to the sum of average region lengths.
