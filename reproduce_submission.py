"""Replay a reviewed two-SNV hypothesis against local inputs; no discovery claim."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import platform
import re
import sys
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    os.umask(0o077)
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("hypothesis", "vcf", "fasta", "gtf", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error("Output must not already exist")
    hypothesis = json.loads(args.hypothesis.read_text())
    variants = hypothesis["primary_pair"]
    if len(variants) != 2 or any(len(v[2]) != 1 or len(v[3]) != 1 for v in variants):
        raise ValueError("This replay supports one reviewed pair of SNVs only")
    if variants[0][0] != variants[1][0]:
        raise ValueError("This replay requires a single chromosome")
    matched = {}
    sample = None
    with gzip.open(args.vcf, "rt") as stream:
        for line in stream:
            if line.startswith("#CHROM"):
                samples = line.rstrip().split("\t")[9:]
                if samples != [hypothesis["proband_id"]]:
                    raise ValueError("VCF sample identity mismatch or multiple samples")
                sample = samples[0]
            if line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            for chrom, pos, ref, alt in variants:
                if (fields[0].removeprefix("chr"), int(fields[1]), fields[3]) != (
                    chrom.removeprefix("chr"),
                    pos,
                    ref,
                ) or alt not in fields[4].split(","):
                    continue
                key = f"{chrom}:{pos}:{ref}>{alt}"
                if key in matched:
                    raise ValueError("Duplicate candidate key in VCF")
                call = dict(zip(fields[8].split(":"), fields[9].split(":"), strict=True))
                alt_index = fields[4].split(",").index(alt) + 1
                if str(alt_index) not in re.split(r"[/|]", call["GT"]):
                    raise ValueError("Selected ALT is not carried")
                matched[key] = {"alt_index": alt_index, "FILTER": fields[6], **call}
    if len(matched) != 2 or sample is None:
        raise ValueError("Expected two unique carried alleles and one sample")
    with gzip.open(args.fasta, "rt") as stream:
        header = next(stream).strip().removeprefix(">").split()[0]
        if header.removeprefix("chr") != variants[0][0].removeprefix("chr"):
            raise ValueError("Reference chromosome mismatch")
        sequence = "".join(line.strip() for line in stream).upper()
    intervals = []
    with gzip.open(args.gtf, "rt") as stream:
        for line in stream:
            if f'transcript_id "{hypothesis["transcript"]}"' not in line:
                continue
            fields = line.split("\t")
            if fields[2] == "CDS":
                if fields[6] != "+" or fields[0] != variants[0][0]:
                    raise ValueError("This reviewed replay expects a plus-strand transcript")
                intervals.append((int(fields[3]), int(fields[4])))
    intervals.sort()
    coding = "".join(sequence[start - 1 : end] for start, end in intervals)
    consequences = []
    for index, (chrom, pos, ref, alt) in enumerate(variants):
        if sequence[pos - 1] != ref:
            raise ValueError("Reference allele mismatch")
        containing = [(start, end) for start, end in intervals if start <= pos <= end]
        if len(containing) != 1:
            raise ValueError("Expected unique coding overlap")
        offset = sum(end - start + 1 for start, end in intervals if end < pos)
        offset += pos - containing[0][0]
        codon = coding[offset // 3 * 3 : offset // 3 * 3 + 3]
        mutant = codon[: offset % 3] + alt + codon[offset % 3 + 1 :]
        if (offset + 1, codon, mutant) != (
            hypothesis["expected_cds_positions"][index],
            hypothesis["expected_reference_codons"][index],
            hypothesis["expected_alternate_codons"][index],
        ):
            raise ValueError("Frozen coding consequence mismatch")
        consequences.append(
            {
                "chrom": chrom,
                "pos": pos,
                "cds_pos": offset + 1,
                "residue": offset // 3 + 1,
                "ref_codon": codon,
                "alt_codon": mutant,
            }
        )
    args.output.mkdir(mode=0o700)
    fields = [
        "proband_id",
        "chrom_1",
        "pos_1",
        "ref_1",
        "alt_1",
        "chrom_2",
        "pos_2",
        "ref_2",
        "alt_2",
        "epcr",
        "finding_type",
        "notes",
    ]
    result_csv = args.output / "reproduced.csv"
    with result_csv.open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(fields)
        writer.writerow(
            [
                sample,
                *variants[0],
                *variants[1],
                f"{hypothesis['epcr']:.2f}",
                "primary",
                hypothesis["notes"],
            ]
        )
    receipt = {
        "status": "PASS",
        "scope": "frozen-hypothesis replay, not causal validation",
        "command": sys.argv,
        "python": platform.python_version(),
        "platform": platform.platform(),
        "selected_calls": matched,
        "consequences": consequences,
        "inputs": {
            name: {"path": str(getattr(args, name)), "sha256": sha256(getattr(args, name))}
            for name in ("hypothesis", "vcf", "fasta", "gtf")
        },
        "csv_sha256": sha256(result_csv),
    }
    (args.output / "receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(
        json.dumps(
            {
                "status": "PASS",
                "selected_carried_alleles": len(matched),
                "reference_and_codon_checks": 2,
                "csv_sha256": receipt["csv_sha256"],
            }
        )
    )


if __name__ == "__main__":
    main()
