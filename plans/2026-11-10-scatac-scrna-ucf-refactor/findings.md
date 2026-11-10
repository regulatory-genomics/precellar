# Findings

## Review findings

1. **Release blocker — sibling path dependencies.** `Cargo.toml` points at
   `../UCF/crates/...`; the pushed precellar branch is not standalone-buildable.
2. **Data-loss policy is implicit.** Missing CB returns `None` and is silently
   ignored by the scRNA sink.
3. **Alignment multiplicity is ambiguous.** The sink iterates all `MultiMapR`
   alignments although the converter documents primary-only behavior.
4. **Failure is not atomic.** Final UCF and sidecar paths are opened directly,
   so a late error can leave a plausible partial artifact.
5. **Resource control is incomplete.** Input external-sort chunks are fixed at
   `10_000_000` and are unrelated to the public `chunk_rows` option.
6. **Boundary validation is incomplete.** ScRNA row dictionary IDs and some
   schema/provenance fields are not validated or populated strictly enough.

## Real-data evidence

- ATAC input: `/data3/litian/macrophage/data/ATAC/ATAC_test/`, 2,500 paired
  reads, 150 bp.
- Reference: `/data3/litian/database/genome/hg38.fa` with a local BWA-MEM2
  index.
- BWA-MEM2 direct sink result: 1,136 unique fragments; UCF CLI inspect and
  validate passed; SNV sidecar had 1,136 data rows and matched the fragment TSV
  data-row multiset exactly.
- The DNBelab barcode template was not representative for this input. A
  generic no-whitelist assay plus a derived index-read file was used for the
  successful sink smoke test, so this does not validate realistic barcode
  diversity.
- MiniBWA on the same real input terminated with a native C assertion in
  `ext/minibwa/format.c:187`; BWA-MEM2 completed successfully.

## Refactor results

- The upstream UCF repository is private, so the dependency is a pinned SSH git
  revision rather than HTTPS. Standalone resolution needs an SSH key in CI.
- A copied tree with no sibling `UCF` checkout resolves with
  `cargo metadata --locked` and type-checks with `cargo check --locked`.
- Deleting `temp_dir` before a run now succeeds: the directory is created at the
  job boundary and the sorter returns errors instead of panicking.
- Real ATAC after hardening: 1,136 rows, `ucf-cli validate` passed, SNV sidecar
  multiset equal to `FragmentsSink`.
- Real 10x RNA v3 after hardening: 364 stored rows, 297 cells, 330 UMIs,
  2,131 unmapped, 5 missing-CB, 0 secondary/supplementary, `ucf-cli validate`
  passed, and `ScRnaReader` decoded exactly 364 rows with valid dictionary
  references.
- Default `sort_chunk_size` is 1,000,000 records; measured peak RSS on these
  smoke fixtures is dominated by the index/aligner, not the sorter, so a larger
  fixture is still needed to calibrate a production default.
