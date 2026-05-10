# SeqSpec Authoring Guide

This document describes how to write a seqspec YAML file for `precellar`. It covers fields that the pipeline reads, lists all enum values with their semantics, and highlights constraints that will cause errors at runtime.

## 1. File Overview

A seqspec YAML describes a sequencing assay. It has two main sections:

- **`sequence_spec`** — defines the sequencing reads (FASTQ files), their strand orientation, and which part of the library each read sequences.
- **`library_spec`** — defines the designed library structure: the physical arrangement of all regions (adapters, barcodes, UMIs, target inserts, etc.) in 5'→3' forward order.

The link between the two is the `primer_id` field on each Read, which points to a region in `library_spec` where sequencing begins.

Minimal skeleton:

```yaml
!Assay
seqspec_version: 0.3.0
assay_id: my_assay
name: My Assay
doi: ""
date: 2025-01-01
description: Example assay
modalities:
- rna
lib_struct: ""
library_protocol: Custom
library_kit: Custom
sequence_protocol: Illumina
sequence_kit: Custom
chemistry_strandedness: reverse
sequence_spec:
- !Read
  read_id: R1
  name: Read 1
  modality: rna
  primer_id: truseq_read1      # Must match a region_id in library_spec
  min_len: 28
  max_len: 28
  strand: pos
  files:
  - file_id: R1.fastq.gz
    filename: R1.fastq.gz
    filetype: fastq
    filesize: 0
    url: /path/to/R1.fastq.gz
    urltype: local
    md5: ""
library_spec:
- !Region
  region_id: rna
  region_type: rna
  name: RNA
  sequence_type: joined
  sequence: ""
  min_len: 0
  max_len: 0
  onlist: null
  regions:
  - !Region
    region_id: truseq_read1
    region_type: truseq_read1
    ...
  - !Region
    region_id: barcode
    region_type: barcode
    ...
  - !Region
    region_id: cdna
    region_type: cdna
    ...
```

## 2. Top-Level Fields

All fields below are required by the YAML format (serde deserialization will fail if they are missing). However, only the fields marked with **used** are actually read by `precellar`'s processing logic. The rest are metadata — any non-empty string value will satisfy the parser.

| Field | Used by precellar | Description |
|-------|:-----------------:|-------------|
| `seqspec_version` | no | Format version. Use `0.3.0`. |
| `assay_id` | no | Unique identifier for this assay configuration. |
| `name` | no | Human-readable assay name. |
| `doi` | no | Publication DOI or empty string. |
| `date` | no | Date string (e.g., `2025-01-01`). |
| `description` | no | Free-text description. |
| `modalities` | **yes** | List of modality names present in this assay. |
| `lib_struct` | no | URL to library structure diagram or empty string. |
| `library_protocol` | no | Protocol name string. |
| `library_kit` | no | Kit name string. |
| `sequence_protocol` | no | Sequencing platform (e.g., `Illumina`, `Oxford Nanopore`). |
| `sequence_kit` | no | Sequencing kit name string. |
| `chemistry_strandedness` | **yes** | See below. Optional in YAML, but **required at runtime for RNA modality**. |
| `sequence_spec` | **yes** | List of Read entries. |
| `library_spec` | **yes** | List of top-level modality Region entries. |

### `chemistry_strandedness`

Determines the strandedness relationship between R1 and the original RNA molecule. Used during alignment to set transcript strand orientation.

| Value | Meaning |
|-------|---------|
| `forward` | R1 matches the strand of the original RNA molecule. |
| `reverse` | R1 is the reverse complement of the original RNA molecule. Most 10x assays use this. |
| `unstranded` | The assay captures both strands equally, or strandedness is not applicable. |

### `modalities`

Each entry must be one of:

| Value | Description |
|-------|-------------|
| `rna` | RNA / gene expression |
| `atac` | Chromatin accessibility |
| `dna` | Genomic DNA |
| `protein` | Protein (e.g., CITE-seq) |
| `tag` | Sample tags / hashtags |
| `crispr` | CRISPR guide RNA |

Each modality listed here must have a corresponding top-level Region in `library_spec` with a matching `region_type`.

## 3. `library_spec` — Library Structure

### Structure Rules

1. **Flat child regions (no deep nesting).** The top-level Region is the modality wrapper. Its `regions` list contains the direct children — the actual structural regions. Children must not have their own `regions` (i.e., no further nesting). Violating this causes a runtime error.

2. **Forward order.** Regions are listed in 5'→3' forward order of the designed library structure, regardless of how reads are sequenced.

3. **One modality per top-level Region.** For multi-modal assays (e.g., RNA + ATAC), list multiple top-level Regions.

4. **One target region recommended.** Each modality should have exactly one target region (`cdna` or `gdna`). This is not validated at load time, but the long-read processing logic uses the target to separate 5' and 3' end regions, so multiple or missing targets will produce incorrect results.

### Region Fields

```yaml
- !Region
  region_id: barcode1          # Unique identifier (referenced by primer_id, onlist, etc.)
  region_type: barcode          # See §4 for all valid types
  name: Cell Barcode            # Human-readable name (used for consensus grouping in long reads)
  sequence_type: onlist         # See §5 for all valid types
  sequence: NNNNNNNNNNNNNNNN    # Known sequence for fixed types; N's or X's for others
  min_len: 16                   # Minimum length in bp
  max_len: 16                   # Maximum length in bp
  onlist: !Onlist               # Whitelist configuration (see §6); null if not applicable
    file_id: whitelist.txt
    ...
  regions: null                 # Must be null for leaf regions
```

**`region_id` must be unique** within the entire `library_spec` (across all modalities). It is used as the key for whitelist index lookup and `primer_id` resolution. Duplicate `region_id` values cause a hard error at load time (`"Duplicate region id: ..."`).

**`sequence` must be written in 5'→3' forward orientation**, matching the region order in `library_spec`.

**`name` is used for barcode consensus grouping** in long-read processing. When the same barcode appears at both ends of the library, give both copies the same `name` but different `region_id` values (e.g., `inner_barcode_5p` / `inner_barcode_3p` with `name: "inner barcode"` on both). The extractor groups extracted barcodes by `name` and runs multi-end consensus resolution on groups with more than one entry.

**`min_len` and `max_len` impact on processing:**

- **Short reads:** Used by `SegmentInfo.split()` to determine how many bases to extract for each segment. Fixed-length regions (`min_len == max_len`) are extracted at exact positions; variable-length regions use pattern matching against neighboring fixed anchors.
- **Long reads:** The sum of `max_len` across all end regions (× 1.15 buffer) determines the end-window size for barcode extraction. For non-fixed regions in composite alignment, spacer length is `(min_len + max_len) / 2`.

### Top-Level Modality Region

The top-level Region must use:

```yaml
- !Region
  region_id: rna               # Any unique id
  region_type: rna              # Must match a modality value
  name: RNA                     # Human-readable
  sequence_type: joined         # Must be "joined" for top-level regions
  sequence: ""                  # Can be empty
  min_len: 0
  max_len: 0
  onlist: null
  regions:                      # List of child regions
  - ...
```

## 4. `region_type` Reference

All valid values and their YAML serialization names:

### Target Regions

These define the region of interest (insert). `precellar` uses them to determine the boundary between barcode/structural regions and the sequenced content.

| YAML value | Description |
|------------|-------------|
| `cdna` | Complementary DNA — used in RNA assays. |
| `gdna` | Genomic DNA — used in ATAC / DNA assays. |

Each modality should have exactly one target region. For long reads, the target separates the 5' end regions from the 3' end regions.

### Barcode Regions

| YAML value | Description |
|------------|-------------|
| `barcode` | Cell barcode, sample barcode, or combinatorial index. The only type recognized by `is_barcode()`. Must have `sequence_type: onlist` and an `onlist` configuration with a whitelist to enable barcode correction (short reads) or whitelist matching (long reads). |

### UMI Region

| YAML value | Description |
|------------|-------------|
| `umi` | Unique Molecular Identifier. Typically `sequence_type: random`. |

### Sequencing Primer Regions

These are valid targets for `primer_id` in `sequence_spec`. They mark where a read begins sequencing.

| YAML value | Description |
|------------|-------------|
| `custom_primer` | Generic sequencing primer. |
| `truseq_read1` | Illumina TruSeq Read 1 primer. |
| `truseq_read2` | Illumina TruSeq Read 2 primer. |
| `nextera_read1` | Nextera Read 1 primer. |
| `nextera_read2` | Nextera Read 2 primer. |
| `illumina_p5` | Illumina P5 adapter (also valid as primer). |
| `illumina_p7` | Illumina P7 adapter (also valid as primer). |

**Note:** `primer_id` in a Read should reference a region whose `region_type` is one of these. A warning is logged if it does not, but processing will still proceed.

### Structural / Adapter Regions

These regions provide structural context. For long reads, `fixed` sequence types among these serve as alignment anchors.

| YAML value | Description |
|------------|-------------|
| `linker` | Linker sequence between functional regions. |
| `named` | Generic named region (e.g., adapter with a custom name). |
| `s5` | Nextera S5 adapter. |
| `s7` | Nextera S7 adapter. |
| `me1` | Mosaic End 1. |
| `me2` | Mosaic End 2. |
| `index5` | i5 index. |
| `index7` | i7 index. |
| `poly_a` | Poly-A tail. |
| `poly_t` | Poly-T region. |
| `poly_g` | Poly-G region. |
| `poly_c` | Poly-C region. |
| `methyl` | Methylation region. |
| `hic` | Hi-C specific region. |
| `fastq` | FASTQ format region. |
| `fastq_link` | Link to FASTQ data. |

### Modality Types (top-level only)

Used only for the top-level modality wrapper region:

`rna`, `atac`, `dna`, `protein`, `tag`, `crispr`

## 5. `sequence_type` Reference

| Value | Meaning | `sequence` field | Use case |
|-------|---------|-----------------|----------|
| `fixed` | Known, invariant sequence. | Must contain the actual nucleotide sequence. | Adapters, linkers, primers. Used as alignment anchors in long-read composite alignment. |
| `random` | Unknown sequence. | Typically `N`'s or `X`'s (not used by pipeline). | UMIs, primers with unknown sequence, poly-T. |
| `onlist` | Sequence is one of a known set (whitelist). | Typically `N`'s. | Barcodes. Requires `onlist` field with whitelist file. |
| `joined` | Composite region made of child regions. | Can be empty. | **Only for the top-level modality region.** Must have `regions` list. |

**Long-read constraint:** In each end-region collection that contains a barcode, the total length of all `fixed` regions must be at least **12 bp** for composite alignment to work. This is enforced at runtime.

## 6. `onlist` — Whitelist Configuration

Required for `barcode` regions with `sequence_type: onlist`. Defines which whitelist file to use for barcode matching.

```yaml
onlist: !Onlist
  file_id: 3M-february-2018.txt.gz
  filename: 3M-february-2018.txt.gz
  filetype: txt
  filesize: 0
  url: /path/to/3M-february-2018.txt.gz
  urltype: local
  md5: ""
  location: local
  reverse_complement: false     # Optional, defaults to false
```

| Field | Required | Description |
|-------|----------|-------------|
| `file_id` | yes | Identifier for this file. |
| `filename` | yes | File name. |
| `filetype` | yes | File type (typically `txt`). |
| `filesize` | yes | File size in bytes (can be `0`). |
| `url` | yes | Path or URL to the whitelist file. |
| `urltype` | yes | One of: `local`, `ftp`, `http`, `https`. |
| `location` | optional | One of: `local`, `remote`. |
| `md5` | yes | MD5 checksum or empty string. |
| `reverse_complement` | optional | Default `false`. When `true`, the extracted barcode candidate is reverse-complemented before matching against this whitelist. |

### Whitelist File Format

The whitelist file must contain **one barcode sequence per line**, with no headers, indices, or extra columns. Each line is read as-is and used as a whitelist entry.

```
AAGAAAGTTGTCGGTGTCTTTGTG
TCGATTCCGTTTGTAGTCGTCTGT
GAGTCTTGTGTCCCAGTTACCAGG
```

Files with additional columns (e.g., `1\tAAGAAAGTT...`) will cause errors because the entire line — including the index and tab character — is treated as the barcode sequence. If your whitelist has extra columns, preprocess it with `cut -f2 input.txt > whitelist.txt` or equivalent.

### When to Use `reverse_complement: true`

In some library designs (e.g., plate-based scATAC-seq with long-read sequencing), the same cell barcode appears at both ends of the library, but the 3' copy is the reverse complement of the 5' copy in the forward library orientation. Setting `reverse_complement: true` on the 3' barcode's onlist ensures that:

1. Both ends match against the same whitelist entries.
2. Multi-end barcode consensus resolution correctly merges the two observations.

**Short-read behavior:** The `reverse_complement` flag is XOR'd with the Read's `strand` direction to determine whether to RC the extracted barcode before whitelist comparison.

**Long-read behavior:** The extracted candidate sequence is reverse-complemented before whitelist matching, independent of per-read orientation detection, as it was normalized to forward strand before.

## 7. `sequence_spec` — Read Definitions

Each entry in `sequence_spec` describes one sequencing read (one FASTQ file).

```yaml
- !Read
  read_id: R1
  name: Read 1
  modality: rna
  primer_id: truseq_read1
  min_len: 28
  max_len: 28
  strand: pos
  files:
  - file_id: sample_R1.fastq.gz
    filename: sample_R1.fastq.gz
    filetype: fastq
    filesize: 0
    url: /path/to/sample_R1.fastq.gz
    urltype: local
    md5: ""
```

| Field | Required | Description |
|-------|----------|-------------|
| `read_id` | yes | Unique identifier for this read. |
| `name` | optional | Human-readable name. |
| `modality` | yes | Which modality this read belongs to (`rna`, `atac`, etc.). |
| `primer_id` | yes | `region_id` of the sequencing primer region in `library_spec` where this read starts. Must be a direct child of the modality region and must have a sequencing primer `region_type`. |
| `min_len` | yes | Minimum read length in bp. Records shorter than this are discarded. |
| `max_len` | yes | Maximum read length in bp. For long reads, use `2147483647` (`i32::MAX`, meaning no upper limit). |
| `strand` | yes | Orientation of the read relative to the library structure. See below. |
| `files` | optional | List of FASTQ file entries. Can be left empty and set later via the Python API (`assay.update_read(read_id, fastq=...)`). |

### `strand`

| Value | Meaning | When to use |
|-------|---------|-------------|
| `pos` | Read sequence matches the forward (5'→3') library orientation. | Most R1 reads in short-read Illumina assays. |
| `neg` | Read sequence is the reverse complement of the forward library orientation. | R2 reads in paired-end Illumina assays. |
| `unstranded` | Read orientation varies per record. | **Long-read assays** (e.g., Oxford Nanopore). Each record's orientation is detected dynamically. |

**Constraint:** Short-read assays must not use `unstranded`. This is validated at runtime.

### How `primer_id` Links Reads to the Library

The `primer_id` tells `precellar` where in the library structure this read begins sequencing:

1. Find the region with matching `region_id` in the modality's child regions.
2. Skip that primer region.
3. Extract all subsequent regions as segments of the read.
4. If `strand: neg`, iterate in reverse (from the 3' end backward).

**Example:** Given this library structure:

```
illumina_p5 → truseq_read1 → barcode → umi → cdna → truseq_read2 → illumina_p7
```

- A Read with `primer_id: truseq_read1` and `strand: pos` will produce segments: `[barcode, umi, cdna, truseq_read2, illumina_p7]`
- A Read with `primer_id: truseq_read2` and `strand: neg` will iterate in reverse from `truseq_read2`, producing segments: `[cdna, umi, barcode, ...]` (reversed)

The `max_len` of the Read truncates the segment list — only segments that fit within the read length are kept.

## 8. Short-Read vs Long-Read Differences

The processing path is determined by the `strand` field and read lengths.

### Short-Read Processing

- **`strand: pos` or `strand: neg`** — short-read path.
- Segments are split from the read using fixed-length extraction and pattern matching (KMP algorithm for fixed-sequence anchors).
- Barcode correction uses quality scores and mismatch-based probabilistic correction.
- Each barcode region must have fixed, known positions relative to the primer.

### Long-Read Processing

- **`strand: unstranded`** — long-read path.
- Typically only one Read entry in `sequence_spec` (the full-length nanopore read).
- `primer_id` points to the outermost region on one end.
- `precellar` derives **end regions** from `library_spec`:
  - **5' end regions:** All regions from the library start through the last fixed/barcode region before the target.
  - **3' end regions:** All regions from the first fixed/barcode region after the target through the library end.
- End regions are used for composite alignment-based barcode extraction.

**Key requirements for long-read library_spec:**

1. **Fixed anchors are critical.** Regions with `sequence_type: fixed` serve as alignment anchors. Their `sequence` field must contain the actual nucleotide sequence. The total fixed sequence per end must be ≥ 12 bp.
2. **All barcode regions must have whitelists.** Unlike short reads, long-read barcode extraction requires every barcode region to have a non-empty `onlist`.
3. **`min_len` / `max_len` accuracy matters.** They determine the end-window sampling size and spacer lengths in composite alignment.
4. **Target region placement.** The target (`cdna` or `gdna`) separates 5' and 3' processing. Ensure it is correctly positioned in the region order.

## 9. Multi-Modal Assays

For assays that produce multiple modalities (e.g., RNA + ATAC), list all modalities at the top level and provide separate library and read entries for each.

```yaml
modalities:
- rna
- atac

sequence_spec:
- !Read
  read_id: rna-R1
  modality: rna
  primer_id: rna-truseq_read1
  ...
- !Read
  read_id: atac-R1
  modality: atac
  primer_id: atac-nextera_read1
  ...

library_spec:
- !Region
  region_id: rna
  region_type: rna
  sequence_type: joined
  regions:
  - region_id: rna-truseq_read1
    ...
  - region_id: rna-barcode
    ...
  - region_id: rna-cdna
    ...

- !Region
  region_id: atac
  region_type: atac
  sequence_type: joined
  regions:
  - region_id: atac-nextera_read1
    ...
  - region_id: atac-barcode
    ...
  - region_id: atac-gdna
    ...
```

**Tips:**
- Use prefixed `region_id` values (e.g., `rna-barcode`, `atac-barcode`) to avoid collisions.
- Each Read's `modality` field determines which library branch it is associated with.
- `chemistry_strandedness` is only used by the RNA alignment path. It has no effect on non-RNA modalities (ATAC, DNA, etc.).

## 10. Common Patterns

### 10x Chromium scRNA-seq (Short-Read)

```yaml
chemistry_strandedness: reverse
library_spec:
- !Region
  region_id: rna
  region_type: rna
  sequence_type: joined
  regions:
  - region_id: illumina_p5
    region_type: illumina_p5
    sequence_type: random
    min_len: 29
    max_len: 29
  - region_id: truseq_read1
    region_type: truseq_read1
    sequence_type: fixed
    sequence: TCTTTCCCTACACGACGCTCTTCCGATCT   # Actual primer sequence
    min_len: 29
    max_len: 29
  - region_id: barcode
    region_type: barcode
    sequence_type: onlist
    sequence: NNNNNNNNNNNNNNNN
    min_len: 16
    max_len: 16
    onlist: !Onlist
      file_id: 3M-february-2018.txt.gz
      filename: 3M-february-2018.txt.gz
      filetype: txt
      filesize: 0
      url: /path/to/whitelist.txt.gz
      urltype: local
      md5: ""
      location: local
  - region_id: umi
    region_type: umi
    sequence_type: random
    sequence: NNNNNNNNNNNN
    min_len: 12
    max_len: 12
  - region_id: cdna
    region_type: cdna
    sequence_type: random
    min_len: 1
    max_len: 1000
```

R1 reads barcode + UMI (`primer_id: truseq_read1`, `strand: pos`).
R2 reads cDNA in reverse (`primer_id: truseq_read2`, `strand: neg`).

### scNanoATAC-seq (Long-Read, Dual Barcode at Both Ends)

This example demonstrates multi-barcode extraction with consensus resolution and `reverse_complement` onlist handling.

Library structure (5'→3'):
```
start_linker → outer_barcode → linker → inner_barcode → adapter → gDNA → adapter → inner_barcode(RC) → linker → outer_barcode(RC) → end_linker
```

```yaml
sequence_spec:
- !Read
  read_id: R1
  modality: atac
  primer_id: start_linker
  min_len: 100                    # Low threshold
  max_len: 2147483647
  strand: unstranded              # Per-record orientation detection

library_spec:
- !Region
  region_id: atac
  region_type: atac
  sequence_type: joined
  regions:
  - region_id: start_linker
    region_type: custom_primer
    sequence_type: fixed
    sequence: ATCT
    min_len: 4
    max_len: 4
  - region_id: outer_barcode_5p
    region_type: barcode
    name: outer barcode           # Same name on both ends → triggers consensus
    sequence_type: onlist
    sequence: NNNNNNNNNNNNNNNNNNNNNNNN
    min_len: 24
    max_len: 24
    onlist: !Onlist
      file_id: 96_barcode.txt
      url: /path/to/96_barcode.txt
      urltype: local
      ...
  - region_id: linker_5p
    region_type: linker
    sequence_type: fixed
    sequence: CTACACGACGCTCTTCCGATCT
    min_len: 22
    max_len: 22
  - region_id: inner_barcode_5p
    region_type: barcode
    name: inner barcode           # Same name on both ends → triggers consensus
    sequence_type: onlist
    min_len: 24
    max_len: 24
    onlist: !Onlist
      file_id: 96_barcode.txt
      url: /path/to/96_barcode.txt
      urltype: local
      ...
  - region_id: adapter_5p
    region_type: named
    sequence_type: fixed
    sequence: TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG
    min_len: 33
    max_len: 33
  - region_id: gdna
    region_type: gdna
    sequence_type: random
    min_len: 1000
    max_len: 2147483647
  - region_id: adapter_3p
    region_type: named
    sequence_type: fixed
    sequence: CTGTCTCTTATACACATCTGACGCTGCCGACGA
    min_len: 33
    max_len: 33
  - region_id: inner_barcode_3p
    region_type: barcode
    name: inner barcode           # Same name as 5' copy
    sequence_type: onlist
    min_len: 24
    max_len: 24
    onlist: !Onlist
      file_id: 96_barcode.txt
      url: /path/to/96_barcode.txt
      urltype: local
      reverse_complement: true    # 3' copy is RC of whitelist
      ...
  - region_id: linker_3p
    region_type: linker
    sequence_type: fixed
    sequence: AGATCGGAAGAGCGTCGTGTAG
    min_len: 22
    max_len: 22
  - region_id: outer_barcode_3p
    region_type: barcode
    name: outer barcode           # Same name as 5' copy
    sequence_type: onlist
    min_len: 24
    max_len: 24
    onlist: !Onlist
      file_id: 96_barcode.txt
      url: /path/to/96_barcode.txt
      urltype: local
      reverse_complement: true    # 3' copy is RC of whitelist
      ...
  - region_id: end_linker
    region_type: linker
    sequence_type: fixed
    sequence: AGAT
    min_len: 4
    max_len: 4
```

Key points:
- `outer_barcode_5p` and `outer_barcode_3p` share `name: "outer barcode"` → multi-end consensus resolution.
- `inner_barcode_5p` and `inner_barcode_3p` share `name: "inner barcode"` → multi-end consensus resolution.
- 3' barcodes set `reverse_complement: true` because they appear as RC of the whitelist in forward library orientation.
- All barcode groups must be resolved for a valid barcode output (CB tag = outer_barcode + inner_barcode = 48bp).

## 11. Validation Checklist

Use this checklist to verify a seqspec YAML before running `precellar`. Each item is annotated with its enforcement level:

- **error** — causes a runtime error if violated.
- **warn** — logs a warning but processing continues.
- **not checked** — precellar does not validate this; incorrect values will silently produce wrong results.

### Structure

- [ ] Child regions of the modality region must not have their own `regions` sub-list. (**error** — `validate_structure()`)
- [ ] Top-level modality regions use `sequence_type: joined`. (**not checked** — but required for segment extraction to work)
- [ ] Each modality listed in `modalities` has a corresponding top-level Region in `library_spec`. (**not checked**)
- [ ] Regions within each modality are listed in 5'→3' forward order. (**not checked** — incorrect order silently breaks extraction)
- [ ] One target region (`cdna` or `gdna`) per modality. (**not checked** — but long-read end-region derivation depends on it)

### Read–Library Linking

- [ ] Every Read's `primer_id` matches a `region_id` that exists as a direct child of the modality region. (**error** — `verify()`)
- [ ] The referenced primer region has a sequencing primer `region_type`. (**warn** — `update_read()`)
- [ ] Every Read's `modality` matches one of the entries in the top-level `modalities` list. (**not checked**)

### Strand

- [ ] Short-read assays use `strand: pos` or `strand: neg` — never `unstranded`. (**error** — `validate_strands()`)
- [ ] Long-read assays use `strand: unstranded`. (**not checked** — but using `pos`/`neg` will skip per-read orientation detection)

### Barcodes

- [ ] Every `barcode` region has `sequence_type: onlist`. (**not checked** — but barcode correction will not work without it)
- [ ] Every `barcode` region has an `onlist` entry with a valid whitelist file path/URL. (**not checked** at load time)
- [ ] **Long-read only:** All barcode regions must have non-empty whitelists. (**error** — `BarcodeExtractor::new()`)
- [ ] `reverse_complement: true` is set only on barcode regions whose extracted sequence (in forward library orientation) is the RC of the whitelist entries. (**not checked** — incorrect setting silently breaks matching)

### Long-Read Specific

- [ ] Total fixed-sequence length ≥ 12 bp in each end-region collection that contains a barcode. (**error** — `BarcodeExtractor::new()`)
- [ ] Fixed regions have their actual nucleotide sequence in the `sequence` field, written in 5'→3' forward orientation. (**not checked** — wrong sequences or orientation silently break composite alignment)
- [ ] `min_len` and `max_len` are accurate — they determine end-window size and composite alignment spacer lengths. (**not checked**)

### RNA Modality

- [ ] `chemistry_strandedness` is set (`forward`, `reverse`, or `unstranded`). (**error** at alignment time for RNA modality)

### Files

- [ ] All `url` paths point to existing files (FASTQ, whitelist). (**error** — at file open time)
- [ ] `urltype` matches the URL scheme (`local` for filesystem paths, `https` for remote URLs). (**not checked** — wrong type causes file-not-found errors)
