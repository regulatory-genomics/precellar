# Universal Single-Cell Genomics Preprocessing package

![PyPI](https://img.shields.io/pypi/v/precellar)
![PyPI - Downloads](https://img.shields.io/pypi/dm/precellar)
![Continuous integration](https://github.com/regulatory-genomics/precellar/workflows/test-python-package/badge.svg)
![GitHub Repo stars](https://img.shields.io/github/stars/regulatory-genomics/precellar?style=social)

This tool is an automated pipeline for preprocessing single-cell genomics data.
It is designed to take raw data (fastq files) from a variety of single-cell genomics
platforms and a seqspec file as input, and output a count matrix (RNA) or a fragment file (ATAC)
for downstream analysis. The seqspec files for common platforms can be found here: https://github.com/IGVF/seqspec.

In addition to FASTQ, precellar can consume lossless SeqArc archives. SeqArc
archives containing R1/R2 and index reads are declared as grouped seqspec
inputs; see [SeqArc input](seqarc.md) for the supported YAML format and archive
requirements.

## Installation

### Stable version

```
pip install precellar
```

### Development version

```
pip install 'git+https://github.com/regulatory-genomics/precellar.git'
```

CI builds also need a repository secret named `UCF_REPO_TOKEN` with read access
to the private UCF repository. The Rust setup action rewrites the pinned SSH
dependency to authenticated HTTPS before Cargo resolves it.

## Examples

Alignment uses the chainable `FastqPipeline -> AlignmentJob -> terminal output`
API. See [Alignment Pipeline API](docs/alignment_pipeline.md) for complete
configuration, middleware, output, QC, and migration documentation.

Each example dataset below contains approximately 2.5 million fastq records.

> [!NOTE]
> You need to **change the paths to the reference genome** in the examples below.
> The reference genome should be downloaded and indexed before running the examples.
> STAR genome index for human and mouse can be downloaded from [here](https://www.10xgenomics.com/support/software/cell-ranger/downloads#reference-downloads).
> BWA-MEM2, MiniBWA, and Minimap2 indexes can be built by their aligner constructors when `fasta` is provided.

### UCF output

`ScAtacUcfSink` and `ScRnaUcfSink` write UCF directly from the alignment stream.
There is no intermediate fragment TSV or BAM, so the conversion costs one extra
external sort rather than a full write/read/compress cycle.

```python
# scATAC: genome-major fragments in scATAC UCF v2
atac_qc = (
    precellar.FastqPipeline(atac_assay)
    .align_with(precellar.aligners.BwaMem2("/data/Public/BWA_MEM2_index/GRCh38"))
    .run_sink(
        precellar.sinks.ScAtacUcfSink(
            "sample.sc-atac.v2.ucf",
            compute_snv=True,   # optional: writes sample.sc-atac.v2.ucf.snv.tsv.zst
        )
    )
)

# scRNA: cell-aware alignments in scRNA UCF (family sc_rna_alignment)
rna_qc = (
    precellar.FastqPipeline(rna_assay, modality="rna")
    .align_with(precellar.aligners.Star("STAR_reference/refdata-gex-GRCh38-2024-A"))
    .run_sink(precellar.sinks.ScRnaUcfSink("sample.sc-rna.ucf"))
)
```

Both sinks take **exactly one assay per output**; the assay id becomes the UCF
sample id. Cell dictionaries contain observed barcodes only. Smart-seq and
other UMI-free libraries keep an explicit "no UMI" sentinel in the scRNA payload.

Outputs are staged in the destination directory and moved into place only after
every section and sidecar is written, so a failed run leaves no file at the
requested path. A missing `temp_dir` is created instead of failing.

`ScRnaUcfSink` makes record retention explicit:

| Option | Default | Meaning |
|---|---|---|
| `alignments` | `"primary"` | Store primary alignments only, or `"all"` mapped alignments. |
| `on_missing_barcode` | `"drop"` | Skip and count mapped records without `CB`, or `"error"`. |
| `chunk_rows` | `65536` | Rows per independently decodable UCF chunk. |
| `sort_chunk_size` | `1000000` | Records buffered in memory before an external-sort run spills to disk. |
| `storage` | `"v1"` | `"v1"` writes `sc_rna_alignment_v1` (fixed 48-byte rows). `"v2"` writes chunked columns with packed UMIs. |
| `block_size_bp` | `8388608` | v2 genomic block used by the layout and the start-offset column. |
| `cell_tile_size` | `4096` | v2 cell tile; rows store a `u16` rank inside the tile. Must not exceed 65536. |
| `has_mate_columns` | `True` | v2 emits mate columns. Set `False` for single-end data; a paired row then fails instead of dropping its mate. |
| `has_template_length_column` | `True` | v2 emits a template-length column; same rule as `has_mate_columns`. |

`storage="v2"` writes the same logical information (exact CIGAR, aligned blocks,
flags, NH, MAPQ, cell barcode, exact UMI, optional mate/template). It changes only
the physical layout:

- rows are sorted by `(chromosome, genomic block, cell tile, start, cell rank, UMI)`;
- block and CIGAR offsets are derived from per-row counts instead of stored;
- canonical 12 bp A/C/G/T UMIs are stored as 24-bit codes, with a small exception
  table for other sequences; the reader returns the exact original string;
- aligned blocks keep a one-bit gap type (`D` versus `N`) in the block encoding.

The SNV evidence sidecar is written in the same order as the alignment file for
both storage versions, so `row_index` stays a valid link.

The QC report counts `records_stored`, `records_unmapped`,
`records_secondary`, `records_supplementary`, `records_missing_barcode`,
`records_missing_reference`, and `records_empty_cigar`, so no record disappears
without a number.

`ScAtacUcfSink` takes the same `chunk_rows`/`sort_chunk_size` pair.

`ScAtacUcfSink` keeps the standard scATAC UCF contract. Unless
`compute_snv=True`, no SNV work is done. When enabled, SNVs are written to a
zstd TSV sidecar with the same columns as `FragmentsSink(compute_snv=True)`, so
existing UCF readers do not pay for variable-length SNV data.

`ScRnaUcfSink(compute_snv=True)` writes a row-aligned `sc_rna_evidence_v1`
sidecar (`<output>.evidence.tsv.zst`) in the same pass. Each line is one
substitution, insertion, or deletion:

```
row_index  read_pos  ref_offset  op  ref_base  read_base  base_qual
```

`row_index` is the 0-based ordinal of the stored alignment in UCF order, so the
sidecar links to the UCF without widening its 48-byte row. `read_pos` is the
position in the original read, `ref_offset` is relative to the alignment start,
and `op` is `X`, `I`, or `D`. Rows with no events are omitted. Evidence comes
from `MD:Z` plus CIGAR, sequence, and base qualities; `min_base_qual` can drop
low-quality events at write time. This is raw evidence, not a genotype call.

### Gene Expression

`GeneQuantificationSink` normally reads transcript annotation from the
configured STAR aligner. To quantify RNA alignments produced by another
aligner, such as long-read Minimap2, provide a STAR index containing
`geneInfo.tab`, `transcriptInfo.tab`, `exonInfo.tab`, `chrNameLength.txt`
and `chrStart.txt` explicitly:

Build the aligner and STAR indexes from the same reference FASTA, and ensure the
GTF uses matching chromosome names. Reference names and lengths are validated
before reads are processed, so mismatches such as `chr1` versus `1` fail
immediately instead of inflating intergenic reads.

```python
rna_qc = (
    precellar.FastqPipeline(assay, modality="rna")
    .align_with(
        precellar.aligners.Minimap2(
            "GRCh38.splice_hq.mmi",
            preset="splice:hq",
            build_if_missing=False,
        )
    )
    .run_sink(
        precellar.sinks.GeneQuantificationSink(
            "gene_matrix.h5ad",
            strandedness="reverse",
            transcriptome_index="GRCh38_STAR_index",
        )
    )
)
```

<details>
<summary>10x scRNA-seq v3</summary>

```python
import precellar

assay = precellar.Assay('https://raw.githubusercontent.com/regulatory-genomics/precellar/refs/heads/main/seqspec_templates/10x_rna_v3.yaml')

data = precellar.examples.txg_rna_v3()
assay.add_illumina_reads('rna')
assay.update_read('rna-R1', fastq=data['R1'])
assay.update_read('rna-R2', fastq=data['R2'])

rna_qc = (
    precellar.FastqPipeline(assay)
    .align_with(
        precellar.aligners.Star("STAR_reference/refdata-gex-GRCm39-2024-A"),
        num_threads=8,
    )
    .run_sink(precellar.sinks.GeneQuantificationSink("gene_matrix.h5ad"))
)
print(rna_qc)
```

</details>


<details>
<summary>sci-RNA-seq3</summary>

```python
import precellar

assay = precellar.Assay('https://raw.githubusercontent.com/regulatory-genomics/precellar/refs/heads/main/seqspec_templates/sci_rna_seq3.yaml')

data = precellar.examples.sci_rna_seq3()
assay.update_read('R1', fastq=data['R1'])
assay.update_read('R2', fastq=data['R2'])

rna_qc = (
    precellar.FastqPipeline(assay, modality="rna")
    .align_with(
        precellar.aligners.Star("STAR_reference/refdata-gex-GRCm39-2024-A"),
        num_threads=8,
    )
    .run_sink(precellar.sinks.GeneQuantificationSink("gene_matrix.h5ad"))
)
print(rna_qc)
```

</details>


<details>
<summary>MARS-seq</summary>

```python
import precellar

assay = precellar.Assay('https://raw.githubusercontent.com/regulatory-genomics/precellar/refs/heads/main/seqspec_templates/mars_seq.yaml')

data = precellar.examples.mars_seq()
assay.add_illumina_reads(modality='rna')
assay.update_read('rna-R1', fastq=data['R1'])
assay.update_read('rna-R2', fastq=data['R2'])

rna_qc = (
    precellar.FastqPipeline(assay, modality="rna")
    .align_with(
        precellar.aligners.Star("STAR_reference/refdata-gex-GRCm39-2024-A"),
        num_threads=8,
    )
    .run_sink(precellar.sinks.GeneQuantificationSink("gene_matrix.h5ad"))
)
print(rna_qc)
```

</details>

<details>
<summary>MGI DNBelab C4 scRNA-seq v1</summary>

```python
import precellar

assay = precellar.Assay('https://raw.githubusercontent.com/regulatory-genomics/precellar/refs/heads/main/seqspec_templates/dnbelabc4_rna_v1_dark.yaml')

data = precellar.examples.dnbelabc4_rna_v1()
assay.update_read('R1', fastq=data['R1'])
assay.update_read('R2', fastq=data['R2'])

rna_qc = (
    precellar.FastqPipeline(assay, modality="rna")
    .align_with(
        precellar.aligners.Star("STAR_reference/refdata-gex-GRCh38-2024-A"),
        num_threads=8,
    )
    .run_sink(precellar.sinks.GeneQuantificationSink("gene_matrix.h5ad"))
)
print(rna_qc)
```

</details>

### Chromatin accessibility and protein-DNA interactions

<details>
<summary>10x scATAC-seq</summary>

```python
import precellar

assay = precellar.Assay('https://raw.githubusercontent.com/regulatory-genomics/precellar/refs/heads/main/seqspec_templates/10x_atac.yaml')

data = precellar.examples.txg_atac()
assay.add_illumina_reads('atac', forward_strand_workflow=True)
assay.update_read('atac-I2', fastq=data['I2'])
assay.update_read('atac-R1', fastq=data['R1'])
assay.update_read('atac-R2', fastq=data['R2'])
atac_qc = (
    precellar.FastqPipeline(assay)
    .align_with(
        precellar.aligners.BwaMem2("/data/Public/BWA_MEM2_index/GRCh38"),
        num_threads=8,
    )
    .run_sink(precellar.sinks.FragmentsSink("fragments.tsv.zst"))
)
print(atac_qc)
```

</details>

<details>
<summary>MGI DNBelab C4 scATAC-seq v1</summary>

```python
import precellar

assay = precellar.Assay('https://raw.githubusercontent.com/regulatory-genomics/precellar/refs/heads/main/seqspec_templates/dnbelabc4_atac_v1.yaml')

data = precellar.examples.dnbelabc4_atac_v1()
assay.update_read('R1', fastq=data['R1'])
assay.update_read('R2', fastq=data['R2'])

qc = (
    precellar.FastqPipeline(assay)
    .align_with(
        precellar.aligners.BwaMem2("/data/Public/BWA_MEM2_index/GRCh38"),
        num_threads=8,
    )
    .run_sink(precellar.sinks.FragmentsSink("fragments.tsv.zst"))
)
print(qc)
```

</details>


<details>
<summary>dscATAC-seq</summary>

```python
import precellar

assay = precellar.Assay('https://raw.githubusercontent.com/regulatory-genomics/precellar/refs/heads/main/seqspec_templates/dscATAC.yaml')

data = precellar.examples.dsc_atac()
assay.update_read('R1', fastq=data['R1'])
assay.update_read('R2', fastq=data['R2'])

atac_qc = (
    precellar.FastqPipeline(assay, modality="atac")
    .align_with(
        precellar.aligners.BwaMem2("/data/Public/BWA_MEM2_index/GRCm39"),
        num_threads=8,
    )
    .run_sink(precellar.sinks.FragmentsSink("fragments.tsv.zst"))
)
print(atac_qc)
```

</details>


<details>
<summary>scifi-ATAC-seq</summary>

```python
import precellar

assay = precellar.Assay('https://raw.githubusercontent.com/regulatory-genomics/precellar/refs/heads/main/seqspec_templates/scifi_atac.yaml')

data = precellar.examples.scifi_atac()
assay.update_read('I2', fastq=data['I2'])
assay.update_read('R1', fastq=data['R1'])
assay.update_read('R2', fastq=data['R2'])

atac_qc = (
    precellar.FastqPipeline(assay, modality="atac")
    .align_with(
        precellar.aligners.BwaMem2("BWA_MEM2_index/Zea_mays"),
        num_threads=8,
    )
    .run_sink(precellar.sinks.FragmentsSink("fragments.tsv.zst"))
)
print(atac_qc)
```

</details>


### Multi-Omics

<details>
<summary>10x single-cell multiome (Gene expression + ATAC)</summary>

```python
import precellar

assay = precellar.Assay('https://raw.githubusercontent.com/regulatory-genomics/precellar/refs/heads/main/seqspec_templates/10x_rna_atac.yaml')

data = precellar.examples.txg_multiome()
assay.add_illumina_reads('rna')
assay.update_read('rna-R1', fastq=data['rna-R1'])
assay.update_read('rna-R2', fastq=data['rna-R2'])

assay.add_illumina_reads('atac', forward_strand_workflow=True)
assay.update_read('atac-I2', fastq=data['atac-I2'])
assay.update_read('atac-R1', fastq=data['atac-R1'])
assay.update_read('atac-R2', fastq=data['atac-R2'])

rna_qc = (
    precellar.FastqPipeline(assay, modality="rna")
    .align_with(
        precellar.aligners.Star("STAR_reference/refdata-gex-GRCm39-2024-A"),
        num_threads=8,
    )
    .run_sink(precellar.sinks.GeneQuantificationSink("gene_matrix.h5ad"))
)
print(rna_qc)

atac_qc = (
    precellar.FastqPipeline(assay, modality="atac")
    .align_with(
        precellar.aligners.BwaMem2("/data/Public/BWA_MEM2_index/GRCm39"),
        num_threads=8,
    )
    .run_sink(precellar.sinks.FragmentsSink("fragments.tsv.zst"))
)
print(atac_qc)
```

</details>


<details>
<summary>SHARE-seq</summary>

```python
import precellar

assay = precellar.Assay('https://raw.githubusercontent.com/regulatory-genomics/precellar/refs/heads/main/seqspec_templates/share_seq.yaml')

data = precellar.examples.share_seq()
assay.update_read('rna-I1', fastq=data['rna-I1'])
assay.update_read('rna-R1', fastq=data['rna-R1'])
assay.update_read('rna-R2', fastq=data['rna-R2'], min_len=10, max_len=10)

assay.update_read('atac-I1', fastq=data['atac-I1'])
assay.update_read('atac-R1', fastq=data['atac-R1'])
assay.update_read('atac-R2', fastq=data['atac-R2'])

rna_qc = (
    precellar.FastqPipeline(assay, modality="rna")
    .align_with(
        precellar.aligners.Star("/data/Public/STAR_reference/refdata-gex-GRCh38-2024-A/star/"),
        num_threads=8,
    )
    .run_sink(precellar.sinks.GeneQuantificationSink("gene_matrix.h5ad"))
)
print(rna_qc)

atac_qc = (
    precellar.FastqPipeline(assay, modality="atac")
    .align_with(
        precellar.aligners.BwaMem2("/data/Public/BWA_MEM2_index/GRCh38"),
        num_threads=8,
    )
    .run_sink(precellar.sinks.FragmentsSink("fragments.tsv.zst"))
)
print(atac_qc)
```

</details>


<details>
<summary>SNARE-seq</summary>

```python
import precellar

assay = precellar.Assay('https://raw.githubusercontent.com/regulatory-genomics/precellar/refs/heads/main/seqspec_templates/snare_seq.yaml')

data = precellar.examples.snare_seq()
assay.update_read('rna-R1', fastq=data['rna-R1'])
assay.update_read('rna-R2', fastq=data['rna-R2'])

assay.update_read('atac-I1', fastq=data['atac-I1'])
assay.update_read('atac-R1', fastq=data['atac-R1'])
assay.update_read('atac-R2', fastq=data['atac-R2'])

rna_qc = (
    precellar.FastqPipeline(assay, modality="rna")
    .align_with(
        precellar.aligners.Star("/data/Public/STAR_reference/GRCm39/"),
        num_threads=8,
    )
    .run_sink(precellar.sinks.GeneQuantificationSink("gene_matrix.h5ad"))
)
print(rna_qc)

atac_qc = (
    precellar.FastqPipeline(assay, modality="atac")
    .align_with(
        precellar.aligners.BwaMem2("/data/Public/BWA_MEM2_index/GRCm39"),
        num_threads=8,
    )
    .run_sink(precellar.sinks.FragmentsSink("fragments.tsv.zst"))
)
print(atac_qc)
```

</details>


<details>
<summary>Droplet Paired-Tag</summary>

```python
import precellar

assay = precellar.Assay('https://raw.githubusercontent.com/regulatory-genomics/precellar/refs/heads/main/seqspec_templates/droplet_paired_tag.yaml')

data = precellar.examples.droplet_paired_tag()
assay.add_illumina_reads('rna')
assay.update_read('rna-R1', fastq=data['rna-R1'])
assay.update_read('rna-R2', fastq=data['rna-R2'])

assay.add_illumina_reads('atac')
assay.update_read('atac-I2', fastq=data['atac-I2'])
assay.update_read('atac-R1', fastq=data['atac-R1'])
assay.update_read('atac-R2', fastq=data['atac-R2'])

rna_qc = (
    precellar.FastqPipeline(assay, modality="rna")
    .align_with(
        precellar.aligners.Star("STAR_reference/refdata-gex-GRCm39-2024-A"),
        num_threads=8,
    )
    .run_sink(precellar.sinks.GeneQuantificationSink("gene_matrix.h5ad"))
)
print(rna_qc)

atac_qc = (
    precellar.FastqPipeline(assay, modality="atac")
    .align_with(
        precellar.aligners.BwaMem2("/data/Public/BWA_MEM2_index/GRCm39"),
        num_threads=8,
    )
    .run_sink(precellar.sinks.FragmentsSink("fragments.tsv.zst"))
)
print(atac_qc)
```

</details>
