---
status: in_progress
current_task: 5
updated: 2026-11-10
---

# Execution State

## Now

Implementation complete on branches `UCF/feat/sc-rna-alignment` and
`precellar/feat/scatac-scrna-ucf`. Remaining from `plan.md`: release build and
the direct-vs-two-step benchmark fixture.

## Blockers

- None for implementation. Both repos build and test on this host.

## Log

### 2026-11-10

- Did: created the plan directory.
- Did: inspected precellar sink lifecycle, ATAC fragment/SNV path, transcriptome alignment data, UCF scATAC writer, and UCF bulk RNA contract.
- Did: recorded the minimal design: direct scATAC sink plus optional SNV sidecar; distinct cell-aware scRNA UCF family; no changes to existing sinks.
- Did (UCF): added `sc_rna_alignment` family, 48-byte cell-aware row, cell dictionary, UMI table, streaming chunked writer, chunk-pruned reader, section types, and `crates/ucf-rna/tests/sc_rna_alignment.rs`. Bulk `rna_alignment_v1/v2` untouched.
- Did (precellar-core): added `ucf` module with a shared external sorter, scATAC row/dictionary helpers, and scRNA conversion from `RecordBuf` (CB/RX tags, CIGAR codes, block count, NH).
- Did (precellar): added `ScAtacUcfSink` (direct scATAC UCF v2 + optional zstd SNV sidecar) and `ScRnaUcfSink` (cell-aware scRNA UCF), single-assay guard, and sink registration.
- Verified: `cargo test -p ucf-rna` -> 18 tests pass; `cargo test --workspace` in precellar -> 12 + 111 + 5 + 19 tests pass; `cargo check --workspace` has no new warnings.
- Verified: built the debug cdylib and imported it with Python 3.12; both sink classes construct and reject invalid options.
- Next: run the release build and the direct-vs-two-step benchmark, then commit.
