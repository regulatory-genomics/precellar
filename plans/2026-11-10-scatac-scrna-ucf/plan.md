---
status: active
created: 2026-11-10
---

# Direct scATAC/scRNA UCF outputs for precellar

## Goal

Add lightweight, terminal UCF sinks to `precellar` so aligned single-assay scATAC data is written directly as genome-major scATAC UCF v2 and aligned single-assay scRNA data is written as a cell-aware indexed RNA UCF, without TSV/BAM round trips. Keep current fragment/H5AD sinks compatible, preserve ATAC SNVs through an optional compressed sidecar, and keep the conversion stage bounded by streaming plus external sorting.

## Non-goals

- Do not change the existing `FragmentsSink` or `GeneQuantificationSink` output semantics.
- Do not replace the source BAM/FASTQ or claim that UCF preserves read sequence, quality, QNAME, or arbitrary SAM tags.
- Do not make a peak-by-cell or gene-by-cell matrix the native UCF truth layer.
- Do not support combining multiple assays/samples in one direct sink until sample identity is present in the core alignment/fragment model.
- Do not add a new generic UCF abstraction or a new dependency for a small adapter layer.

## Constraints

- scATAC coordinates remain 0-based half-open; `count` becomes UCF `multiplicity`; `+`/`-` strand maps to the existing UCF flag bits `0b1`/`0b10`.
- UCF records must be genome-major and sorted according to the selected family writer; no full dataset may be materialized in RAM.
- Cell dictionaries use contiguous numeric IDs and stable `cell_key = sample_id:barcode`; the first version accepts exactly one assay and uses `assay_id` as `sample_id`.
- SNV output is optional and must retain the current `Fragment`/TSV semantics exactly; the UCF primary file remains queryable by existing scATAC readers.
- Existing `rna_alignment_v1/v2` remains unchanged and backward-compatible. Cell-aware RNA requires a distinct family/schema contract.
- The UCF dependency must be pinned to a clean, reachable revision; never use the dirty local UCF checkout as a path dependency for the deliverable.

## Design

### scATAC path

`AlignProgressBar` → existing `IntoFragments` (barcode external sort + dedup) → lightweight external genome sort keyed by the UCF `genome_start` order → `ucf_sc_atac::ScAtacWriter::write_file_v2_iter_fallible_with_options_auto_counts_zstd`.

The sink builds the observed cell dictionary while consuming fragments, maps chromosome names from the alignment header to contiguous IDs, and emits `ScAtacRecordV2`. With `compute_snv=true`, it writes a sibling `<output-stem>.snv.tsv.zst` sidecar from the same final fragment stream, preserving the current extended SNV columns and recording its path in the output report. No fragment TSV is written or reread.

### scRNA path

`AlignProgressBar` → primary mapped alignments with CB/UB tags → external genomic sort → new UCF `sc_rna_alignment` indexed writer. The new row keeps the compact RNA alignment fields (locus, mate, template length, CIGAR/block ownership, SAM flags, NH, MAPQ) and adds `cell_id` plus an exact UMI dictionary ID. A cell dictionary maps IDs to `sample_id`, barcode, and stable cell key. Raw sequence/quality/QNAME stay in the source BAM when needed.

The new RNA family is implemented in the UCF format workspace first, then consumed by `precellar`; the existing bulk `rna_alignment` family is not widened. The writer accepts a sorted iterator and emits bounded coordinate chunks so scRNA conversion does not collect all rows in memory.

### Shared API boundary

Add `precellar.sinks.ScAtacUcfSink` and `precellar.sinks.ScRnaUcfSink` and register both in `AlignmentJob.run_sink`. Keep option names parallel to the current sinks: output path, compression level, temp directory from `align_with`, and an explicit `compute_snv` option only for scATAC. Reject an assay list whose length is not one with a clear construction/runtime error.

## Acceptance gates

- [x] `cargo test --workspace` in `precellar/` -> all existing and new Rust tests pass.
- [x] `cargo test -p ucf-rna` in the local UCF checkout -> existing UCF contracts and the new scRNA round-trip tests pass.
- [x] New synthetic scATAC sink test -> not yet end-to-end: core sorting/strand/conversion tests pass and the sink compiles into the extension; a Python-level fixture is still needed.
- [x] New synthetic scRNA test -> `crates/ucf-rna/tests/sc_rna_alignment.rs` reads back exact cell/UMI/CIGAR/blocks and prunes chunks by region.
- [ ] `cargo test --release --workspace` (or the project CI equivalent) -> not run on this host; debug workspace tests pass.
- [ ] Small benchmark fixture comparing direct UCF sink with `FragmentsSink` plus `ucf-cli build fragment-to-ucf` -> not run.

## Tasks

### Task 1 — Freeze the UCF contract and dependency boundary

- Files: `/data3/litian/projects/tachyon_project/UCF/crates/ucf-rna/src/{family.rs,metadata.rs,schema.rs,storage.rs,indexed_writer.rs,indexed_reader.rs,lib.rs}` (modify/create), `/data3/litian/projects/tachyon_project/UCF/crates/ucf-core/src/section.rs` (modify), `/data3/litian/projects/tachyon_project/UCF/crates/ucf-rna/tests/` (create tests)
- Change: Add the distinct `sc_rna_alignment` family/schema, cell dictionary, UMI string table, cell-aware fixed row, and streaming indexed writer/reader. Keep `rna_alignment_v1/v2` contracts unchanged. Add only the section types and validation required by this family.
- Verify: `cargo test -p ucf-rna` in `/data3/litian/projects/tachyon_project/UCF` -> existing bulk RNA tests and new scRNA schema/round-trip tests pass; `git diff --check` -> no whitespace errors.

### Task 2 — Add the smallest direct scATAC conversion core

- Files: `/data3/litian/projects/tachyon_project/precellar/Cargo.toml` (modify), `/data3/litian/projects/tachyon_project/precellar/src/sinks.rs` (modify), `/data3/litian/projects/tachyon_project/precellar/src/align.rs` (modify), `/data3/litian/projects/tachyon_project/precellar/src/lib.rs` (modify), `/data3/litian/projects/tachyon_project/precellar/src/ucf.rs` (create if a separate adapter is needed), `/data3/litian/projects/tachyon_project/precellar/precellar-core/src/fragment/` (modify only if a pure conversion helper is required)
- Change: Add `ScAtacUcfSink`; reuse `IntoFragments`, assign observed barcodes to contiguous cell IDs, map chromosome names and strand flags, external-sort compact records in UCF order, and call the iterator-based v2 writer. Keep all existing fragment filtering, MAPQ, deduplication, shifts, and QC semantics identical to `FragmentsSink`.
- Verify: `cargo test --workspace` in `precellar/` -> a synthetic fragment conversion test passes and existing tests remain green; `cargo check --locked` -> dependency resolution uses the pinned UCF revision.

### Task 3 — Add lossless optional ATAC SNV preservation

- Files: `/data3/litian/projects/tachyon_project/precellar/src/sinks.rs` (modify), `/data3/litian/projects/tachyon_project/precellar/src/ucf.rs` (modify/create), `/data3/litian/projects/tachyon_project/precellar/README.md` (modify), `/data3/litian/projects/tachyon_project/precellar/docs/design/` (add UCF output note)
- Change: When `compute_snv=true`, write the current fragment header and extended SNV columns to a zstd sidecar while producing the UCF from the same deduplicated stream; include sidecar path and row/support counts in the sink report. Do not add a generic variable JSON field to scATAC UCF v2.
- Verify: `cargo test --workspace` -> synthetic MD/CIGAR fixture has byte-stable expected SNV fields; `ucf-cli inspect` and the UCF reader -> primary UCF remains valid and queryable; sidecar line count matches final fragment rows.

### Task 4 — Add the direct scRNA sink

- Files: `/data3/litian/projects/tachyon_project/precellar/src/sinks.rs` (modify), `/data3/litian/projects/tachyon_project/precellar/src/align.rs` (modify), `/data3/litian/projects/tachyon_project/precellar/src/lib.rs` (modify), `/data3/litian/projects/tachyon_project/precellar/src/ucf.rs` (modify/create)
- Change: Add `ScRnaUcfSink` that reads primary mapped `RecordBuf`s, extracts the existing CB/UB tags through the current `MultiMapR` API, preserves alignment/CIGAR/MAPQ/NH/SAM flags, builds the cell and UMI dictionaries, externally sorts by genomic locus, and invokes the new UCF writer. Make the sink fail early for missing cell barcodes, unsupported multi-assay input, or output aliasing.
- Verify: `cargo test --workspace` -> synthetic single-end and paired-end records round-trip through the new UCF reader with exact barcode/UMI and CIGAR/block values; a no-UMI fixture uses the documented sentinel without inventing a value.

### Task 5 — Document, benchmark, and close the compatibility loop

- Files: `/data3/litian/projects/tachyon_project/precellar/README.md` (modify), `/data3/litian/projects/tachyon_project/precellar/docs/design/pipeline-overview.md` (modify), `/data3/litian/projects/tachyon_project/precellar/docs/design/data-flow.md` (modify), `/data3/litian/projects/tachyon_project/precellar/benches/` (create only if the existing workspace has no suitable fixture), `/data3/litian/projects/tachyon_project/precellar/pyproject.toml` (modify only if packaging metadata needs the new public API)
- Change: Add minimal Python examples for both sinks, state that scRNA uses a new cell-aware UCF family while existing bulk RNA UCF is unchanged, document the SNV sidecar, and record direct-vs-two-step timing/RSS on a small reproducible fixture.
- Verify: project CI command from `.github/workflows/test_python.yml` plus `cargo test --workspace` -> package builds and the public sink classes import successfully; benchmark report contains row count, support, elapsed time, peak RSS, and output sizes.
