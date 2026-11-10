---
status: in_progress
current_task: 5
updated: 2026-11-10
---

# Execution State

## Now

Implementation is pushed on branches `UCF/feat/sc-rna-alignment` and
`precellar/feat/scatac-scrna-ucf`. Real scATAC FASTQ validation passed through
BWA-MEM2 and the direct UCF sink; the real scRNA smoke test is running in the
background.

## Blockers

- The UCF branch keeps unrelated pre-existing CLI/table-query edits uncommitted.
- MiniBWA aborts on the real ATAC FASTQ with a C assertion; BWA-MEM2 succeeds.
- A nonexistent user-supplied `temp_dir` causes an existing fragment sorter
  `unwrap()` panic instead of a returned error.

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
- Verified: real ATAC FASTQ (`2,500` read pairs) with hg38 and BWA-MEM2 produced
  `1,136` unique fragments through both `FragmentsSink` and `ScAtacUcfSink`.
  `ucf-cli inspect` and `ucf-cli validate` passed for the direct UCF; the SNV
  sidecar contained `1,136` data rows and its row multiset matched the fragment
  TSV exactly.
- Observed: the DNBelab barcode template gave only `0.08%` valid barcodes on
  the local ATAC FASTQ, so that run is a real alignment/sink smoke test rather
  than a representative single-cell barcode benchmark. A generic no-whitelist
  assay with a derived index-read file was used for the successful UCF run.
- Next: finish the real scRNA smoke test, record the review findings, and use
  `plans/2026-11-10-scatac-scrna-ucf-refactor/` for the next implementation
  phase.
