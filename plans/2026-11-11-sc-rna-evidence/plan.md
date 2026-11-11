# sc_rna_evidence_v1 sidecar

## Goal

Add an optional per-row mismatch/indel sidecar to `ScRnaUcfSink` so a single
scRNA alignment pass yields both the splicing-capable UCF and the base-level
evidence needed for SNV / ASE / RNA-editing work.

## Contract

- Sidecar name: `<output>.evidence.tsv.zst` (overridable).
- One file per UCF, produced in the same pass and published atomically with it.
- Row alignment: `row_index` is the 0-based ordinal of the stored alignment in
  UCF order, i.e. the order the writer streams rows after external sorting.
  Rows with no events are omitted, so the file is sparse but strictly ordered.
- Columns: `row_index, read_pos, ref_offset, op, ref_base, read_base, base_qual`.
  - `read_pos`: 0-based position in the original read sequence (soft clips included).
  - `ref_offset`: 0-based offset from the alignment start.
  - `op`: `X` substitution, `I` insertion, `D` deletion.
  - `ref_base`/`read_base`: uppercase base, `.` when not applicable.
  - `base_qual`: Phred score of the read base, `0` for deletions.
- Source of truth is the SAM `MD:Z` tag plus CIGAR and read sequence/quality.
  Rows without `MD` still emit insertions/deletions and are counted in metrics.

## Non-goals

- No genotype calls or allele thresholds; the sidecar stores raw evidence.
- No change to `sc_rna_alignment_v1`; the UCF row stays 48 bytes.

## Acceptance gates

- [ ] Unit tests cover substitution, insertion, deletion, soft-clip offset,
      explicit `X` CIGAR, and rows without `MD`.
- [ ] `ScRnaUcfSink(compute_snv=True)` writes both files and reports event counts.
- [ ] Failure of either write rolls back the pair.
- [ ] Real sample: sidecar row indices are strictly increasing and < row_count,
      and events match the SAM `MD` of the same records.
- [ ] `cargo test --workspace` passes.

## Tasks

1. Expose the existing MD parser for reuse.
2. Add `RnaEvent` + `rna_events()` to `precellar-core`, gated by a policy flag.
3. Thread `events` through `ScRnaRow` and the external sort.
4. Add sidecar writing, staged publication, and metrics to `ScRnaUcfSink`.
5. Verify on one Mannens sample, then run the 23-sample array.
