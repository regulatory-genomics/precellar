# Pipeline Overview

`precellar` preprocesses single-cell genomics data from raw sequencing reads and assay metadata into downstream analysis outputs such as gene count matrices or fragment files.

## Goals

- Provide a unified preprocessing interface across multiple assay designs.
- Separate assay description from execution logic through `seqspec`-driven configuration.
- Produce outputs that are consistent with the target modality and aligner.

## High-Level Stages

1. Load assay metadata and read structure from a `seqspec` description.
2. Associate FASTQ inputs with the expected read roles for each modality.
3. Parse biological tags such as barcodes and UMIs from configured read segments.
4. Route reads through the appropriate alignment workflow for the selected modality.
5. Aggregate aligned records into modality-specific output artifacts.
6. Report quality-control metrics alongside the main output.

## Design Principles

- Assay-specific structure should live in data, not hard-coded branches.
- Modality-specific workflows should share a common entrypoint where possible.
- Outputs should preserve enough metadata for downstream reproducibility.

## Scope

This design documentation focuses on stable intended behavior. Experimental ideas and temporary notes should stay outside the public design docs until they are mature enough to describe as supported behavior.
