"""Validate the public MVA Track 1 CSV contract without echoing case content.

Schema validation is not biological review or authorization to submit.
"""

from __future__ import annotations

import argparse
import csv
import io
import json
from decimal import Decimal, InvalidOperation
from pathlib import Path
from typing import Any, TextIO

REQUIRED = (
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
)


def validate(stream: TextIO) -> dict[str, Any]:
    reader = csv.DictReader(stream)
    fields = reader.fieldnames or []
    errors: list[str] = []
    if len(fields) != len(set(fields)):
        errors.append("duplicate column names")
    if set(fields) not in (set(REQUIRED), set(REQUIRED) | {"notes"}):
        errors.append("expected required Track 1 columns and optional notes only")
    count = 0
    seen: set[tuple[str, ...]] = set()
    for count, row in enumerate(reader, 1):
        if count > 10:
            errors.append("more than ten rows")
            break
        prefix = f"row {count}: "
        if None in row or any(value is None for value in row.values()):
            errors.append(prefix + "malformed CSV row")
            continue
        row = {key: value.strip() for key, value in row.items()}
        for key in ("proband_id", "chrom_1", "pos_1", "ref_1", "alt_1"):
            if not row.get(key):
                errors.append(prefix + f"missing {key}")
        for suffix in ("1", "2"):
            allele = [row.get(f"{key}_{suffix}", "") for key in ("chrom", "pos", "ref", "alt")]
            if suffix == "2" and not any(allele):
                continue
            if not all(allele):
                errors.append(prefix + f"incomplete allele {suffix}")
                continue
            if not allele[1].isdigit() or int(allele[1]) < 1:
                errors.append(prefix + f"position {suffix} must be a positive integer")
            if allele[2] == allele[3]:
                errors.append(prefix + f"identical reference and alternate allele {suffix}")
        try:
            score = Decimal(row.get("epcr", ""))
            valid = score.is_finite() and Decimal(0) < score <= Decimal(1)
        except InvalidOperation:
            valid = False
        if not valid:
            errors.append(prefix + "epcr must be finite and in (0, 1]")
        if row.get("finding_type") not in {"primary", "secondary"}:
            errors.append(prefix + "finding_type must be primary or secondary")
        first = tuple(row.get(f"{key}_1", "") for key in ("chrom", "pos", "ref", "alt"))
        second = tuple(row.get(f"{key}_2", "") for key in ("chrom", "pos", "ref", "alt"))
        if all(second) and first == second:
            errors.append(prefix + "paired alleles must differ")
        normalized_pair = sorted([first, second]) if any(second) else [first]
        identity = (
            row.get("proband_id", ""),
            *(part for allele in normalized_pair for part in allele),
        )
        if identity in seen:
            errors.append(prefix + "duplicate finding")
        seen.add(identity)
    if count == 0:
        errors.append("at least one finding is required")
    return {
        "status": "FAIL" if errors else "PASS",
        "row_count": count,
        "errors": errors,
        "limitations": [
            "Does not verify GRCh38 reference alleles, genotype, phase or causal evidence",
            "Does not calibrate EPCR or check report and repository requirements",
            "Does not submit or authorize submission",
        ],
    }


def smoke_test() -> None:
    header = ",".join(REQUIRED) + "\n"
    good = "synthetic,1,10,A,G,,,,,0.5,primary\n"
    assert validate(io.StringIO(header + good))["status"] == "PASS"
    for bad in (
        good.replace("0.5", "NaN"),
        good.replace("0.5", "0"),
        good.replace(",10,", ",0,"),
        good.replace(",,,,,", ",1,,,,"),
    ):
        assert validate(io.StringIO(header + bad))["status"] == "FAIL"
    assert validate(io.StringIO(header + good * 11))["status"] == "FAIL"
    assert validate(io.StringIO(header))["status"] == "FAIL"
    same_pair = "synthetic,1,10,A,G,1,10,A,G,0.5,primary\n"
    pair = "synthetic,1,10,A,G,1,20,C,T,0.5,primary\n"
    reverse = "synthetic,1,20,C,T,1,10,A,G,0.5,primary\n"
    assert validate(io.StringIO(header + same_pair))["status"] == "FAIL"
    assert validate(io.StringIO(header + pair + reverse))["status"] == "FAIL"
    print("PASS: synthetic CSV schema, finite EPCR, positions, pair completeness and row limit")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csv_path", nargs="?", type=Path)
    parser.add_argument("--smoke-test", action="store_true")
    args = parser.parse_args()
    if args.smoke_test:
        smoke_test()
        return
    if args.csv_path is None:
        parser.error("csv_path is required")
    with args.csv_path.open(newline="", encoding="utf-8-sig") as stream:
        result = validate(stream)
    print(json.dumps(result, indent=2))
    raise SystemExit(0 if result["status"] == "PASS" else 1)


if __name__ == "__main__":
    main()
