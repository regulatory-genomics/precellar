# Core Algorithms

This document describes the core algorithmic responsibilities in `precellar` at a design level.

## Barcode And UMI Extraction

The pipeline derives cell barcodes, sample barcodes, and UMIs from read segments defined by the assay description.

Key expectations:

- The assay definition identifies where each tag is located in the sequencing reads.
- Tag extraction is deterministic for a given assay configuration.
- Extracted tags are carried forward as structured metadata for later filtering, alignment, and quantification.

## Modality-Aware Alignment

The alignment step depends on the target modality.

- RNA workflows produce gene-oriented quantification outputs.
- ATAC workflows produce genomic fragment-oriented outputs.
- The same assay abstraction should support different aligners without changing the user-facing setup pattern more than necessary.

## Output Materialization

After alignment, the pipeline converts intermediate records into the requested output type.

- Gene expression workflows emit count-matrix style outputs.
- Chromatin workflows emit fragment-style outputs.
- QC summaries are generated alongside the main artifact.

## Invariants

- Read-role interpretation must stay consistent with the assay specification.
- Output type must match the selected modality and processing path.
- The preprocessing pipeline should expose a predictable interface even when assay layouts differ.
