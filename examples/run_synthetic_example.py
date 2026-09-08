"""Exercise source/submission identity mapping using entirely invented data."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import os
import subprocess
import sys
from pathlib import Path


def main() -> None:
    os.umask(0o077)
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(mode=0o700, parents=True, exist_ok=False)
    source = Path(__file__).resolve().parent
    hypothesis = json.loads((source / "hypothesis.json").read_text())
    inputs = args.output / "inputs"
    inputs.mkdir()
    vcf = inputs / "synthetic.vcf.gz"
    fasta = inputs / "synthetic.fa.gz"
    gtf = inputs / "synthetic.gtf.gz"
    with gzip.open(vcf, "wt") as stream:
        stream.write(
            "##fileformat=VCFv4.2\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSOURCE_SAMPLE\n"
            "chr1\t4\t.\tA\tG\t99\tPASS\t.\tGT:AD:DP\t0/1:10,10:20\n"
            "chr1\t8\t.\tC\tT\t99\tPASS\t.\tGT:AD:DP\t0/1:12,12:24\n"
        )
    with gzip.open(fasta, "wt") as stream:
        stream.write(">chr1\nATGAAAGCCTAA\n")
    with gzip.open(gtf, "wt") as stream:
        stream.write(
            "chr1\tsynthetic\tCDS\t1\t9\t.\t+\t0\t"
            'gene_id "SYNTHETIC_GENE"; transcript_id "SYNTHETIC_TRANSCRIPT.1";\n'
        )
    checks = []
    cases = [
        ("mapped", hypothesis, "EXAMPLE_CASE", 0),
        (
            "same_identifier",
            {k: v for k, v in hypothesis.items() if k != "vcf_sample_id"}
            | {"proband_id": "SOURCE_SAMPLE"},
            "SOURCE_SAMPLE",
            0,
        ),
        (
            "wrong_source_rejected",
            hypothesis | {"vcf_sample_id": "WRONG_SOURCE"},
            None,
            1,
        ),
    ]
    for name, case, expected_identifier, expected_exit in cases:
        case_file = inputs / f"{name}.json"
        case_file.write_text(json.dumps(case, indent=2) + "\n")
        output = args.output / name
        command = [
            sys.executable,
            str(source.parent / "reproduce_submission.py"),
            "--hypothesis",
            str(case_file),
            "--vcf",
            str(vcf),
            "--fasta",
            str(fasta),
            "--gtf",
            str(gtf),
            "--output",
            str(output),
        ]
        run = subprocess.run(command, capture_output=True, text=True, check=False)
        if run.returncode != expected_exit:
            raise RuntimeError(f"Unexpected exit for {name}: {run.stderr}")
        if expected_exit == 0:
            with (output / "reproduced.csv").open() as stream:
                rows = list(csv.DictReader(stream))
            if len(rows) != 1 or rows[0]["proband_id"] != expected_identifier:
                raise RuntimeError("Source-to-submission identity mapping failed")
        elif output.exists() or "VCF sample identity mismatch" not in run.stderr:
            raise RuntimeError("Wrong sample was not rejected before output")
        checks.append(
            {
                "name": name,
                "command": command,
                "exit_code": run.returncode,
                "expected_exit_code": expected_exit,
                "status": "PASS",
                "stdout": run.stdout,
                "stderr": run.stderr,
            }
        )
    receipt = {
        "scope": "Synthetic software checks only; no real genome or clinical claim",
        "status": "PASS",
        "checks": checks,
    }
    (args.output / "receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print("PASS: mapped identifier, same-identifier compatibility, wrong-source rejection")


if __name__ == "__main__":
    main()
