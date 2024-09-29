# Single-cell genomics preprocessing package

This tool is an automated pipeline for preprocessing single-cell genomics data.
It is designed to take raw data (fastq files) from a variety of single-cell genomics
platforms and a seqspec file as input, and output a count matrix (RNA) or a fragment file (ATAC)
for downstream analysis.

The seqspec files for common platforms can be found here: https://github.com/IGVF/seqspec.