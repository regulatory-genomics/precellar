#!/usr/bin/env python3
"""End-to-end check of the direct UCF sinks against the fragment sink.

Runs one real scATAC and one real scRNA FASTQ set through precellar, then
compares the direct UCF output with the legacy path and validates the result
with ``ucf-cli``. Inputs are supplied on the command line because real FASTQ
data must not live in the repository.

Example::

    python scripts/verify_ucf_sinks.py \
        --ucf-cli /path/to/ucf-cli \
        --index /data/index/GRCh38 --aligner bwa-mem2 \
        --atac-assay seqspec_templates/generic_atac.yaml \
        --atac-r1 ATAC_R1.fastq.gz --atac-r2 ATAC_R2.fastq.gz \
        --atac-i1 ATAC_I1.fastq.gz \
        --scrna-assay seqspec_templates/10x_rna_v3.yaml \
        --scrna-r1 RNA_R1.fastq.gz --scrna-r2 RNA_R2.fastq.gz \
        --scrna-i1 RNA_I1.fastq.gz \
        --output-dir /tmp/ucf-sink-smoke
"""

from __future__ import annotations

import argparse
import gzip
import json
import re
import subprocess
import sys
from pathlib import Path

import precellar


def fastq_records(path: Path) -> int:
    with gzip.open(path, "rt") as handle:
        return sum(1 for _ in handle) // 4


def fragment_rows(path: Path) -> list[str]:
    """Read a zstd-compressed TSV and return its data rows."""
    import zstandard

    with path.open("rb") as handle, zstandard.ZstdDecompressor().stream_reader(handle) as reader:
        text = reader.read().decode()
    return [line for line in text.splitlines() if line and not line.startswith("#")]


def atac_assay(options: argparse.Namespace):
    assay = precellar.Assay(str(options.atac_assay))
    assay.update_read(
        "R1", modality="atac", primer_id="atac-read1", is_reverse=False,
        fastq=str(options.atac_r1),
    )
    assay.update_read(
        "R2", modality="atac", primer_id="atac-read2", is_reverse=True,
        fastq=str(options.atac_r2),
    )
    assay.update_read(
        "I1", modality="atac", primer_id="atac-read2", is_reverse=False,
        fastq=str(options.atac_i1),
    )
    return assay


def scrna_assay(options: argparse.Namespace):
    assay = precellar.Assay(str(options.scrna_assay))
    assay.update_read(
        "rna-R1", modality="rna", primer_id="rna-truseq_read1", is_reverse=False,
        fastq=str(options.scrna_r1),
    )
    assay.update_read(
        "rna-R2", modality="rna", primer_id="rna-truseq_read2", is_reverse=True,
        fastq=str(options.scrna_r2),
    )
    assay.update_read(
        "rna-I1", modality="rna", primer_id="rna-truseq_read2", is_reverse=False,
        fastq=str(options.scrna_i1),
    )
    return assay


def aligner(options: argparse.Namespace):
    if options.aligner == "bwa-mem2":
        return precellar.aligners.BwaMem2(str(options.index), build_if_missing=False)
    if options.aligner == "minibwa":
        return precellar.aligners.MiniBwa(
            str(options.index), build_if_missing=False, preset="sr"
        )
    raise SystemExit(f"unsupported aligner: {options.aligner}")


def run(assay, modality: str, sink, options: argparse.Namespace):
    return (
        precellar.FastqPipeline(assay, modality=modality)
        .align_with(
            aligner(options),
            num_threads=options.threads,
            chunk_size=options.chunk_size,
            temp_dir=str(options.output_dir / "tmp"),
        )
        .run_sink(sink)
    )


def ucf_cli(options: argparse.Namespace, command: str, path: Path) -> str:
    return subprocess.run(
        [str(options.ucf_cli), command, str(path)],
        check=True,
        capture_output=True,
        text=True,
    ).stdout


def inspected_row_count(text: str, path: Path) -> int:
    match = re.search(r"row_count=(\d+)", text)
    if match is None:
        raise RuntimeError(f"ucf-cli inspect did not report row_count for {path}")
    return int(match.group(1))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ucf-cli", type=Path, required=True)
    parser.add_argument("--index", type=Path, required=True)
    parser.add_argument("--aligner", choices=["bwa-mem2", "minibwa"], default="bwa-mem2")
    parser.add_argument("--atac-assay", type=Path, required=True)
    parser.add_argument("--atac-r1", type=Path, required=True)
    parser.add_argument("--atac-r2", type=Path, required=True)
    parser.add_argument("--atac-i1", type=Path, required=True)
    parser.add_argument("--scrna-assay", type=Path)
    parser.add_argument("--scrna-r1", type=Path)
    parser.add_argument("--scrna-r2", type=Path)
    parser.add_argument("--scrna-i1", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--chunk-size", type=int, default=1000)
    parser.add_argument("--sort-chunk-size", type=int, default=1_000_000)
    options = parser.parse_args()
    try:
        import zstandard  # noqa: F401
    except ImportError:
        parser.error("zstandard is required for the sidecar multiset check")
    options.output_dir.mkdir(parents=True, exist_ok=True)

    report: dict[str, object] = {"command": sys.argv}

    atac_dir = options.output_dir / "atac"
    atac_dir.mkdir(exist_ok=True)
    report["atac"] = run(
        atac_assay(options),
        "atac",
        precellar.sinks.FragmentsSink(str(atac_dir / "fragments.tsv.zst"), compute_snv=True),
        options,
    )
    report["atac_ucf"] = run(
        atac_assay(options),
        "atac",
        precellar.sinks.ScAtacUcfSink(
            str(atac_dir / "direct.ucf"),
            compute_snv=True,
            sort_chunk_size=options.sort_chunk_size,
        ),
        options,
    )

    inspect = ucf_cli(options, "inspect", atac_dir / "direct.ucf")
    report["atac_rows"] = inspected_row_count(inspect, atac_dir / "direct.ucf")
    fragments = fragment_rows(atac_dir / "fragments.tsv.zst")
    sidecar = fragment_rows(atac_dir / "direct.ucf.snv.tsv.zst")
    report["atac_fragment_rows"] = len(fragments)
    report["atac_sidecar_rows"] = len(sidecar)
    report["atac_sidecar_matches_fragments"] = sorted(fragments) == sorted(sidecar)
    if not report["atac_sidecar_matches_fragments"]:
        raise RuntimeError("SNV sidecar rows do not match FragmentsSink rows")
    ucf_cli(options, "validate", atac_dir / "direct.ucf")

    if options.scrna_assay is not None:
        if not options.scrna_r1 or not options.scrna_r2 or not options.scrna_i1:
            parser.error("--scrna-assay requires --scrna-r1, --scrna-r2, and --scrna-i1")
        rna_dir = options.output_dir / "scrna"
        rna_dir.mkdir(exist_ok=True)
        report["scrna_ucf"] = run(
            scrna_assay(options),
            "rna",
            precellar.sinks.ScRnaUcfSink(
                str(rna_dir / "direct.ucf"),
                sort_chunk_size=options.sort_chunk_size,
            ),
            options,
        )
        inspect = ucf_cli(options, "inspect", rna_dir / "direct.ucf")
        report["scrna_rows"] = inspected_row_count(inspect, rna_dir / "direct.ucf")
        ucf_cli(options, "validate", rna_dir / "direct.ucf")
        report["scrna_fastq_reads"] = fastq_records(options.scrna_r2)

    print(json.dumps(report, indent=2, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
