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

For more information, please refer to the documentation: https://lab.kaizhang.org/precellar/.