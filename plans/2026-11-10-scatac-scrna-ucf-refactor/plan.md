---
status: active
created: 2026-11-10
---

# Harden direct scATAC/scRNA UCF conversion

## Goal

Turn the pushed direct UCF sinks into a release-safe path: standalone-buildable,
loss-aware, atomic on failure, explicit about primary/secondary alignment policy,
and measurable on real FASTQ data without hiding dropped records or relying on
unbounded sorter defaults.

## Non-goals

- Do not widen the stable scATAC v2 row with arbitrary SNV payloads in this phase.
- Do not add multi-assay merging until sample identity is represented in the core
  alignment/fragment model.
- Do not make MiniBWA silently recover from a native assertion; either fix the
  upstream aligner or document and reject the unsafe input/aligner combination.
- Do not replace the source BAM/FASTQ as the lossless evidence artifact.

## Constraints

- Preserve the successful real-data semantics: scATAC direct UCF rows and the
  optional SNV sidecar must remain multiset-equivalent to `FragmentsSink`.
- Keep UCF rows genome-major and externally sorted without materializing the
  complete dataset in RAM.
- A failed conversion must not leave a file at the requested final path that
  looks valid.
- Every intentional drop must be represented in the report, and missing CB
  behavior must be an explicit option/policy rather than an implicit discard.
- The precellar package must build from its own clone using a clean, reachable
  UCF revision and a committed lockfile.

## Tasks

### Task 1 — Make the dependency boundary release-safe

- Files: `Cargo.toml`, `Cargo.lock`, packaging/CI metadata if required, and the
  UCF repository release/remote metadata.
- Change: replace sibling path-only UCF dependencies with a clean reachable
  revision strategy that works in a standalone precellar checkout; regenerate
  and commit the lockfile; verify that a checkout without `../UCF` can build.
- Verify: clean-clone or copied-tree `cargo check --locked`; Python extension
  import using the project CI Python version.

### Task 2 — Define explicit scRNA record-retention policy

- Files: `precellar-core/src/ucf/sc_rna.rs`, `src/ucf.rs`, tests under
  `precellar-core/tests/`.
- Change: filter primary mapped alignments according to the selected policy;
  count secondary/supplementary/unmapped/missing-CB records separately; make
  missing-CB behavior fail-fast or report-and-drop through a documented sink
  option. Keep the no-UMI sentinel explicit.
- Verify: synthetic records with each SAM flag class and missing CB produce the
  expected report counts and exact UCF row count.

### Task 3 — Add atomic output and temp-directory handling

- Files: `src/ucf.rs`, `precellar-core/src/fragment/mod.rs` or its sorter helper,
  and integration tests.
- Change: create missing temporary directories or return a normal error;
  write UCF and SNV sidecar to unique temporary paths in the destination
  directory; rename only after all sections, sidecars, and validation succeed;
  clean temporary files on failure.
- Verify: nonexistent `temp_dir`, injected writer failure, and sidecar write
  failure leave no final artifact and return a Python exception rather than a
  panic.

### Task 4 — Validate the UCF boundary and preserve provenance

- Files: `UCF/crates/ucf-rna/src/{sc_schema.rs,sc_writer.rs,sc_reader.rs}`,
  `precellar/src/ucf.rs`, and UCF tests.
- Change: validate cell/UMI IDs, chromosome IDs, row counts, schema/payload
  contract, and chunk references at write/read boundaries; populate
  `BuildInfo` with the actual source paths, timestamp, and relevant sink/options
  instead of placeholder values.
- Verify: malformed IDs/contracts are rejected; a real output round-trips with
  non-empty provenance fields and exact dictionary references.

### Task 5 — Make sorter resource policy explicit

- Files: `src/ucf.rs`, `precellar-core/src/ucf/mod.rs`, sink constructors/docs,
  and benchmark scripts.
- Change: separate input-sort chunk size from output UCF chunk size, expose or
  derive it from the existing job chunk/temp configuration, and record it in
  the output report. Avoid hard-coded `10_000_000` input records.
- Verify: the same real fixture runs at two sorter sizes; row counts and
  sidecar multisets match while peak RSS and elapsed time are recorded.

### Task 6 — Close real-data integration coverage

- Files: `tests/` or a repository-approved integration-test harness, plus
  `plans/2026-11-10-scatac-scrna-ucf-refactor/` run records.
- Change: retain a small real FASTQ smoke fixture manifest or an ignored local
  fixture recipe for the 10x scRNA and ATAC examples; test direct UCF against
  the existing fragment/BAM path; include UCF CLI validation and row/dictionary
  equality checks. Record the MiniBWA abort separately from BWA-MEM2 success.
- Verify: real ATAC and real scRNA runs complete with nonzero mapped rows,
  exact counts, valid UCF files, explicit drop metrics, timing, peak RSS, and
  output sizes.

## Acceptance gates

- [x] Standalone precellar clone builds with `cargo check --locked` (copied
  tree with no sibling `UCF` checkout).
- [x] Existing UCF and precellar tests remain green.
- [x] Missing-CB and non-primary policies are tested and reported.
- [x] Failed writes are atomic and panic-free.
- [x] UCF writer/reader reject invalid dictionary references and preserve
  provenance.
- [x] Sort resource policy is configurable and reported.
- [x] Real ATAC and scRNA FASTQ smoke runs pass with nonzero output rows.

## Outcome

| Task | Result |
| --- | --- |
| 1 | Git dependency pinned to `e333d0b`; `Cargo.lock` committed; standalone `cargo check --locked` passes. |
| 2 | `alignments` and `on_missing_barcode` options plus per-class record counts. |
| 3 | Temp-dir creation, `Result`-based sorter, staged scRNA UCF, atomic SNV sidecar. |
| 4 | Range validation on write/read, real `BuildInfo`. |
| 5 | `sort_chunk_size` on both sinks, reported in QC. |
| 6 | Reproducible verification script and real-fixture reader test. |
