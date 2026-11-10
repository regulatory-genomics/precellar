---
status: done
current_task: 6
updated: 2026-11-10
---

# Execution State

## Now

All six refactor tasks are implemented, verified on real FASTQ data, and pushed.

## Blockers

- CI needs an SSH deploy key (or an `insteadOf` rewrite) to fetch the private
  `regulatory-genomics/UCF` git dependency.
- MiniBWA still aborts on the real ATAC FASTQ; that is an upstream aligner
  issue and is recorded rather than worked around.

## Log

### 2026-11-10

- Recorded review findings from the pushed branches.
- Confirmed real ATAC direct-UCF success with BWA-MEM2: 2,500 read pairs,
  1,136 unique fragments, valid UCF, and sidecar multiset equality.
- Recorded MiniBWA assertion failure and nonexistent-temp-directory panic as
  operational risks.
- Created this next-stage hardening plan.
- Did (Task 1): replaced sibling path dependencies with a pinned git revision
  (`e333d0bd37d316d2d8b89dabe5b993fa35030078`) on
  `ssh://git@github.com/regulatory-genomics/UCF.git`; pushed the UCF branch to
  GitHub; committed `Cargo.lock`; a copy of the tree without `../UCF` resolves
  with `cargo metadata --locked` and type-checks with `cargo check --locked`.
- Did (Task 2): `classify_record` returns a disposition; the sink applies
  `alignments=primary|all` and `on_missing_barcode=drop|error` and reports every
  exclusion class.
- Did (Task 3): `prepare_temp_dir` creates/validates the sorter directory at the
  job boundary; sorter errors propagate as `Result` items instead of panics;
  scRNA UCF builds stage to a temp file and rename; the SNV sidecar uses a
  `NamedTempFile` and persists only after a successful write.
- Did (Task 4): UCF writer and reader reject out-of-range chromosome/cell/UMI
  references; `BuildInfo` now carries epoch seconds and the input FASTQ paths.
- Did (Task 5): both sinks expose `sort_chunk_size` separately from
  `chunk_rows` and report it.
- Did (Task 6): added `scripts/verify_ucf_sinks.py` and a real-fixture reader
  test.
- Verified: real ATAC (2,500 pairs) -> 1,136 rows, `ucf-cli validate` passed,
  sidecar multiset equal; real 10x RNA v3 (2,500 reads) -> 364 stored rows, 297
  cells, 330 UMIs, 2,131 unmapped, 5 missing-CB counted, `ucf-cli validate`
  passed, `ScRnaReader` round-trip passed.
- Verified: `cargo test --workspace` (precellar) and `cargo test -p ucf-rna`
  (UCF) pass; deleting `temp_dir` no longer panics.
