"""Stream a whole public Ensembl 1000G VCF against a private coding/splice subset.

Exact allele matches only. Missingness never implies rarity. No network requests.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path

POPULATIONS = ("AFR", "AMR", "EAS", "EUR", "SAS")


def key(chrom: str, pos: str, ref: str, alt: str) -> tuple[str, int, str, str]:
    chrom = chrom.removeprefix("chr")
    return ("MT" if chrom == "M" else chrom, int(pos), ref, alt)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1048576), b""):
            h.update(block)
    return h.hexdigest()


def run(target: Path, reference: Path, output: Path) -> dict[str, int]:
    output.mkdir(parents=True, exist_ok=False, mode=0o700)
    targets = set()
    with target.open() as stream:
        for row in csv.DictReader(stream, delimiter="\t"):
            flags = row.get("location_flags", "")
            if row.get("allele_carriage") != "carried":
                continue
            if not any(
                flag in flags for flag in ("CDS", "splice_acceptor_2bp", "splice_donor_2bp")
            ):
                continue
            targets.add(key(row["chrom"], row["pos"], row["ref"], row["alt"]))
            if len(targets) > 1000000:
                raise ValueError("target cap exceeded; partition target file locally")
    matches: dict[tuple[str, int, str, str], list[str]] = {}
    conflicts = set()
    reference_ids: dict[tuple[str, int, str, str], set[str]] = {}
    counts = {
        "target_alleles": len(targets),
        "reference_records": 0,
        "matched_alleles": 0,
        "malformed_frequency_rows": 0,
        "duplicate_reference_allele_rows": 0,
    }
    outfile = output / "coding_splice_population.csv"
    with gzip.open(reference, "rt") as source, outfile.open("w", newline="") as dest:
        writer = csv.writer(dest)
        writer.writerow(
            [
                "chrom",
                "pos",
                "ref",
                "alt",
                "population_match",
                *[f"AF_{p}" for p in POPULATIONS],
                "reference_ids",
            ]
        )
        for line in source:
            if line.startswith("#"):
                continue
            cells = line.rstrip("\n").split("\t")
            counts["reference_records"] += 1
            alts = cells[4].split(",")
            candidates = [(i, key(cells[0], cells[1], cells[3], alt)) for i, alt in enumerate(alts)]
            selected = [(i, allele) for i, allele in candidates if allele in targets]
            if not selected:
                continue
            info = dict(item.split("=", 1) for item in cells[7].split(";") if "=" in item)
            for i, allele in selected:
                values = []
                for pop in POPULATIONS:
                    raw = info.get(pop, ".").split(",")
                    value = raw[i] if len(raw) == len(alts) else "."
                    if value != ".":
                        try:
                            if not 0 <= float(value) <= 1:
                                raise ValueError
                        except ValueError:
                            counts["malformed_frequency_rows"] += 1
                            value = "."
                    values.append(value)
                if allele in matches:
                    counts["duplicate_reference_allele_rows"] += 1
                    if matches[allele] != values:
                        conflicts.add(allele)
                else:
                    matches[allele] = values
                reference_ids.setdefault(allele, set()).add(cells[2])
        for allele in sorted(targets):
            if allele in conflicts:
                status, values = "conflicting_reference_records", ["."] * len(POPULATIONS)
            elif allele in matches:
                status, values = "exact_allele_match", matches[allele]
            else:
                status, values = "unknown_not_found", ["."] * len(POPULATIONS)
            writer.writerow(
                [*allele, status, *values, ";".join(sorted(reference_ids.get(allele, set())))]
            )
    counts["matched_alleles"] = len(matches)
    counts["conflicting_reference_alleles"] = len(conflicts)
    counts["unknown_alleles"] = len(targets - matches.keys())
    receipt = {
        **counts,
        "inputs": [{"path": str(p), "sha256": sha256(p)} for p in (target, reference)],
        "method_sha256": sha256(Path(__file__)),
        "output_sha256": sha256(outfile),
        "limitations": [
            "Exact representation only; no left normalization",
            "Only fully called carried CDS/splice phenotype-gene subset",
            "1000G ascertainment/sample-size/ancestry limits; absent is unknown",
            "No caller AF or global minor-allele MAF substituted for ALT frequency",
            "HPO presence/negation and causal interpretation unreviewed",
        ],
    }
    (output / "receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
    for path in output.iterdir():
        path.chmod(0o600)
    return counts


def smoke() -> None:
    import tempfile

    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        target = root / "target.tsv"
        target.write_text(
            "chrom\tpos\tref\talt\tallele_carriage\tlocation_flags\nchr1\t10\tA\tT\tcarried\tCDS\n1\t20\tC\tG\tcarried\tsplice_donor_2bp\n1\t30\tC\tG\tuncarried\tCDS\n"
        )
        reference = root / "reference.vcf.gz"
        with gzip.open(reference, "wt") as stream:
            stream.write(
                "##fileformat=VCFv4.2\n1\t10\t.\tA\tG,T\t.\t.\tAFR=0.2,0.7;EUR=0.1,0.8;MAF=0.01\n"
            )
        result = run(target, reference, root / "out")
        assert result["target_alleles"] == 2 and result["matched_alleles"] == 1
        with (root / "out" / "coding_splice_population.csv").open() as stream:
            rows = list(csv.DictReader(stream))
        assert rows[0]["AF_AFR"] == "0.7" and rows[0]["AF_EUR"] == "0.8"
        assert all(None not in row for row in rows) and rows[0]["reference_ids"] == "."
        assert rows[1]["population_match"] == "unknown_not_found"
    print("PASS: multiallelic Number=A mapping, unknown missingness, carriage and location subset")


def main() -> None:
    os.umask(0o077)
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target", type=Path)
    parser.add_argument("--reference", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--smoke-test", action="store_true")
    args = parser.parse_args()
    if args.smoke_test:
        smoke()
    elif args.target and args.reference and args.output:
        print(json.dumps(run(args.target, args.reference, args.output)))
    else:
        parser.error("target, reference and output required")


if __name__ == "__main__":
    main()
