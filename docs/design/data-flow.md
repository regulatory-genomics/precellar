# Data Flow

This page outlines the major data transitions in the preprocessing workflow.

## Inputs

- Assay specification
- Sequencing reads
- Reference index for the selected aligner
- Runtime configuration such as thread count and output type

## Processing Flow

1. The assay specification defines the expected read layout and logical read names.
2. Input FASTQ files are attached to those logical reads.
3. The pipeline extracts barcodes, UMIs, and modality-specific sequence content.
4. Reads are aligned against the configured reference.
5. Aligned records are transformed into quantification or fragment outputs.
6. QC metrics are summarized for the final report.

## UCF Sequence

```text
aligned records
  -> fragment dedup (ATAC) or primary alignments (RNA)
  -> external genome sort
  -> UCF writer (bounded chunks)
```

No intermediate fragment TSV is produced. The SNV sidecar, when requested, is
emitted from the same deduplicated fragment stream in UCF order.

## Outputs

- Primary modality-specific output file
- Quality-control summary object or report
- Optional UCF companion files (ATAC SNV sidecar)

## Notes For Contributors

When algorithm behavior changes, update these design docs to capture the intended stable behavior rather than transient implementation details.
