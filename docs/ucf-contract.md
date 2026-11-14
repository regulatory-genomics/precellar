# Producer contract for UCF consumers

This checkout writes the UCF payloads consumed by Tachyon and Tachyon Browser.
The UCF dependency is pinned in `Cargo.toml` to the repository revision used by
the current serving stack; a producer update must update that revision and rerun
the downstream UCF fixtures.

## scATAC

- family: `sc_atac`
- schema: `sc_atac_fragment_v2`
- coordinates: 0-based half-open `[start, end)`
- `cell_id` is the physical payload row id; it is not a barcode or project row
- `cell_key` is the stable identity `sample_id:barcode`
- `multiplicity` is retained as weighted fragment support and must not be
  expanded by serving code unless an explicit compatibility API asks for it
- cell dictionary order is the payload row order
- the writer rejects invalid spans and preserves nonzero multiplicity

## scRNA

The scRNA UCF sink preserves the declared cell key, UMI, alignment blocks,
CIGAR operations, primary/secondary flags, and explicit retention counts. The
v2 storage layout changes physical ordering only; it does not change these
logical fields.

## Verification

Run the producer's focused tests and the UCF reader compatibility tests before
publishing a new producer revision. The test in `src/ucf.rs` pins the schema,
coordinate, identity, and multiplicity fields that downstream code relies on.
