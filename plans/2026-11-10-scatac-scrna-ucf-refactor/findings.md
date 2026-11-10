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
