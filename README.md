# Single-cell genomics preprocessing package

![PyPI](https://img.shields.io/pypi/v/precellar)
![PyPI - Downloads](https://img.shields.io/pypi/dm/precellar)
![Continuous integration](https://github.com/regulatory-genomics/precellar/workflows/test-python-package/badge.svg)
![GitHub Repo stars](https://img.shields.io/github/stars/regulatory-genomics/precellar?style=social)

This tool is an automated pipeline for preprocessing single-cell genomics data.
It is designed to take raw data (fastq files) from a variety of single-cell genomics
platforms and a seqspec file as input, and output a count matrix (RNA) or a fragment file (ATAC)
for downstream analysis. The seqspec files for common platforms can be found here: https://github.com/IGVF/seqspec.

## Installation

### Stable version

```
pip install precellar
```

### Development version

```
pip install 'git+https://github.com/regulatory-genomics/precellar.git#egg=precellar&subdirectory=python'
```

## Examples

### 10X scATAC-seq

```python
import precellar

assay = precellar.Assay('https://raw.githubusercontent.com/regulatory-genomics/precellar/refs/heads/main/seqspec_templates/10x_atac.yaml')
assay.add_illumina_reads('atac')
assay.update_read('atac-R1', fastq='R1.fastq.gz')
assay.update_read('atac-I2', fastq='R2.fastq.gz')
assay.update_read('atac-R2', fastq='R3.fastq.gz')
qc = precellar.align(
    assay,
    "/data/Public/BWA_MEM2_index/GRCh38",
    modality="atac",
    output='fragments.tsv.zst',
    output_type='fragment',
    num_threads=32,
)
print(qc)
```

### 10X scRNA-seq

```python
import precellar

assay = precellar.Assay('https://raw.githubusercontent.com/regulatory-genomics/precellar/refs/heads/main/seqspec_templates/10x_rna_v3.yaml')
assay.add_illumina_reads('rna')
assay.update_read('rna-R1', fastq='R1.fastq.gz')
assay.update_read('rna-R2', fastq='R2.fastq.gz')
qc = precellar.align(
    assay,
    "/data/STAR_reference/star_2.7.1",
    modality="rna",
    output="gene_matrix.h5ad",
    output_type="gene_quantification",
    num_threads=32,
)
print(qc)
```

### 10X single-cell multiome (Gene expression + ATAC)

```python
import precellar

assay = precellar.Assay('https://raw.githubusercontent.com/regulatory-genomics/precellar/refs/heads/main/seqspec_templates/10x_rna_atac.yaml')

assay.add_illumina_reads('rna')
assay.update_read('rna-R1', fastq='gex_R1.fastq.gz')
assay.update_read('rna-R2', fastq='gex_R2.fastq.gz')

assay.add_illumina_reads('atac')
assay.update_read('atac-R1', fastq='atac_R1.fastq.gz')
assay.update_read('atac-I2', fastq='atac_R2.fastq.gz')
assay.update_read('atac-R2', fastq='atac_R3.fastq.gz')

rna_qc = precellar.align(
    assay,
    "/data/STAR_reference/star_2.7.1",
    modality="rna",
    output="gene_matrix.h5ad",
    output_type="gene_quantification",
    num_threads=32,
)
atac_qc = precellar.align(
    assay,
    "/data/Public/BWA_MEM2_index/GRCh38",
    modality="atac",
    output='fragments.tsv.zst',
    output_type='fragment',
    num_threads=32,
)
```

For more information, please refer to the documentation: https://lab.kaizhang.org/precellar/.