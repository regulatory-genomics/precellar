# Findings

## Verified facts

- `precellar` is a Rust/PyO3 workspace. The terminal output extension point is `src/align.rs::AlignmentSink`; a job is consumed exactly once and the sink receives a streaming `AlignProgressBar`.
- Current outputs are split by modality: `FragmentsSink` converts aligned ATAC records through `precellar-core::fragment::IntoFragments`, while `GeneQuantificationSink` creates an AnnData H5AD through transcriptome annotation and UMI counting.
- `IntoFragments` already performs an external sort by cell barcode, deduplicates per barcode/UMI/alignment fingerprint, then returns genome-sorted fragments within each barcode. `Fragment` contains coordinates, barcode, multiplicity, strand, and optional `ExtendedFields` with SNV summaries.
- `FragmentsSink(compute_snv=true)` currently preserves SNVs only through `Fragment`'s TSV `Display` implementation. The current scATAC UCF record has coordinates, `cell_id`, `multiplicity`, and `flags`; it has no variable-length SNV field.
- UCF `ucf-sc-atac` exposes an iterator-based v2 writer, requires contiguous cell/chromosome dictionaries, and validates `genome_start` order as `(chromosome, genomic block, start, cell, end, flags)`. The CLI already maps `+` to flag bit `0b1` and `-` to `0b10`.
- UCF `ucf-rna` currently stores **bulk** BAM alignments. Its 40-byte row has no cell barcode, UMI, or SNV, and the current writer materializes a `RnaDatasetV1` before writing. Its documented unsupported list explicitly includes single-cell barcode/UMI counting.
- The UCF workspace is a separate repository. The local checkout has uncommitted CLI/browser changes, so `precellar` must not depend on an uncommitted path checkout. The clean UCF baseline currently inspected is commit `3632149c34e27b6b0fe0a0ea90bd977627a42cef`; the eventual dependency must be pinned to a clean, reachable commit.
- The current `AlignmentInfo`/`Fragment` types do not carry assay/sample identity. A first direct sink therefore has to be single-assay-per-output, or the core fragment model must be widened. The lightweight choice is to reject multi-assay jobs and use `assay_id` as the sample ID.

## Decisions

| Decision | Rationale | Date |
| --- | --- | --- |
| Add type-specific `ScAtacUcfSink` and `ScRnaUcfSink`; do not create a generic modality-switching `UcfSink`. | The two UCF families have different contracts and data models; separate sinks keep the Python API and Rust branches small. | 2026-11-10 |
| Write scATAC UCF directly from the post-dedup fragment stream, with one external genome sort and no TSV round trip. | Reuses the existing deduplication semantics and removes an avoidable write/parse/compress cycle. | 2026-11-10 |
| Preserve ATAC SNVs in an optional compressed sidecar with explicit provenance, rather than changing the stable scATAC v2 row contract in the first refactor. | SNVs are variable-length and currently exist only in the fragment extension; a new UCF payload section would require coordinated writer, reader, validator, and browser changes. | 2026-11-10 |
| Define a separate `sc_rna_alignment` UCF contract instead of reusing `rna_alignment_v1/v2`. | Reusing bulk RNA UCF would silently discard cell barcode and UMI semantics, which is unacceptable for scRNA. | 2026-11-10 |
| Keep raw sequence/quality/QNAME/arbitrary SAM tags outside the native UCF truth layer. | Existing RNA UCF is alignment-oriented and intentionally compact; full read fidelity remains the BAM/source artifact. | 2026-11-10 |
| Start with one assay per UCF output; fail early for multi-assay jobs. | The current fragment/alignment model has no sample field, and adding one would be a larger cross-cutting refactor than requested. | 2026-11-10 |
| `precellar` consumes UCF crates by relative path (`../UCF/crates/...`). | This workspace already links `ucf-sc-atac`/`ucf-rna` this way from `tachyon-browser` and `tachyon`; the UCF repo's clean commits are on a local Gitea branch that GitHub does not have. | 2026-11-10 |
| scRNA UMI comes from the `RX` tag (`Tag::UMI_SEQUENCE`), not `UB`. | `precellar-core::align::aligners` writes the UMI into `Tag::UMI_SEQUENCE`, and `MultiMap::umi()` reads the same tag. | 2026-11-10 |
| The scRNA row carries `block_count`, computed from CIGAR at conversion time. | The UCF writer validates block/CIGAR totals against metadata, and the sink must know them before streaming. | 2026-11-10 |
| Real ATAC smoke test passes with BWA-MEM2. | 2,500 real read pairs from `/data3/litian/macrophage/data/ATAC/ATAC_test/` yielded 1,136 unique fragments; direct scATAC UCF validation passed and the SNV sidecar matched the fragment TSV multiset. | 2026-11-10 |
| The local ATAC FASTQ is not a representative barcode-valid scATAC library under the DNBelab template. | The template reported 0.16%/0.12% exact barcode matches; the successful sink run used a generic no-whitelist assay and one derived 8/22-base index-read file. | 2026-11-10 |
| MiniBWA is unsafe for this real FASTQ/assembly combination. | The run terminated in `ext/minibwa/format.c:187` with `write_sam_cigar` assertion failure; BWA-MEM2 completed the same reference/data run. | 2026-11-10 |
| New UCF sinks inherit a temp-directory panic from `IntoFragments`. | Passing a nonexistent `temp_dir` reached `precellar-core/src/fragment/mod.rs:457`, where `ExternalSorterBuilder::build().unwrap()` panicked with `NotFound`; the directory must be created or the sorter must return an error. | 2026-11-10 |
| Standalone precellar packaging is not yet possible from the pushed branch. | `Cargo.toml` uses sibling path dependencies for all three UCF crates; a clone of `precellar` without `../UCF` cannot resolve them. | 2026-11-10 |
| scRNA conversion currently drops mapped records without CB silently. | `row_from_record` returns `Ok(None)` for missing barcode and `ScRnaUcfSink` ignores that case; no drop counter or failure policy is exposed. | 2026-11-10 |
| The scRNA sink iterates all `MultiMapR` alignments while the converter comment and plan promise primary mapped alignments. | `src/ucf.rs:468-483` calls `record.iter()`, while `row_from_record` only filters `is_unmapped`; secondary/supplementary rows can enter the UCF. | 2026-11-10 |
| Direct output paths are not atomic. | Both sinks pass the final output path directly to UCF writers; iterator errors or sidecar write errors can leave a plausible partial file at the requested path. | 2026-11-10 |
| Sort memory is not controlled by the public chunk option. | Both new sinks hard-code the external-sort input chunk size to `10_000_000`; `chunk_rows` only controls UCF output chunks. | 2026-11-10 |
| UCF row IDs are not range-validated at the writer boundary. | `ScRnaWriter` validates offsets/counts and metadata totals but does not check `cell_id < cells.len()` or a non-sentinel `umi_id < umis.len()`. | 2026-11-10 |

## Rejected options

- Convert `fragments.tsv.zst` or `BAM` after the sink finishes. This adds a second full scan and loses the direct-stream performance advantage.
- Reuse `rna_alignment_v1/v2` for scRNA. It produces a valid bulk RNA UCF but loses cell/UMI identity and therefore does not satisfy the requested scRNA conversion.
- Store a peak-by-cell or gene-by-cell matrix as the UCF truth. It would discard fragment/alignment-level information and contradict the existing fragment-first UCF direction.
- Add a generic free-form JSON/SNV blob to every UCF row. It would defeat the compact indexed layout and make readers pay for optional information on every query.

## Open questions

- Confirm that “scRNA → UCF” means raw genome alignments with cell/UMI identity, not only the existing gene-count H5AD. The plan assumes raw alignment UCF because it preserves the most information and matches UCF's RNA alignment model.
- Confirm whether a native `sc_rna_alignment` family is acceptable as a coordinated change in the separate UCF repository. It is required for cell-aware, queryable scRNA UCF; a precellar-only change cannot provide that contract without data loss.
