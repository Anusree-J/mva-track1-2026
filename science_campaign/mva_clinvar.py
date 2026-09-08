"""Private, bounded-memory exact-allele ClinVar join; stdout is aggregate only.

No patient queries are sent to any API. Whole public ClinVar VCF is downloaded
separately. Exact representation matches are not normalized-equivalence matches.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import platform
import re
import sqlite3
import sys
from collections import Counter
from datetime import UTC, datetime
from pathlib import Path
from typing import TextIO

FIELDS = [
    "chrom",
    "pos",
    "ref",
    "alt",
    "allele_index",
    "genotype",
    "alt_dosage",
    "phase_status",
    "phase_set",
    "phase_id",
    "physical_phase_genotype",
    "genotype_completeness",
    "filter",
    "qual",
    "depth",
    "allele_depths",
    "genotype_quality",
    "clinvar_id",
    "geneinfo",
    "clnsig",
    "clnrevstat",
    "clndn",
    "clndisdb",
    "clnsigconf",
    "clnalleleid",
    "origin",
    "mc",
    "af_exac",
    "af_esp",
    "af_tgp",
    "classification_bucket",
    "case_multiallelic",
    "clinvar_multiallelic",
]
INFO_FIELDS = [
    "GENEINFO",
    "CLNSIG",
    "CLNREVSTAT",
    "CLNDN",
    "CLNDISDB",
    "CLNSIGCONF",
    "ALLELEID",
    "ORIGIN",
    "MC",
    "AF_EXAC",
    "AF_ESP",
    "AF_TGP",
]


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for data in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(data)
    return h.hexdigest()


def vcf_open(path: Path) -> TextIO:
    return gzip.open(path, "rt") if path.suffix == ".gz" else path.open()


def chrom_key(value: str) -> str:
    value = value.removeprefix("chr")
    return "MT" if value == "M" else value


def info_map(value: str) -> dict[str, str]:
    return dict(item.split("=", 1) for item in value.split(";") if "=" in item)


def bucket(value: str) -> str:
    value = value.lower()
    if "conflict" in value:
        return "conflicting"
    if "pathogenic" in value:
        return "pathogenic_or_likely_pathogenic"
    if "uncertain" in value:
        return "uncertain"
    if "benign" in value:
        return "benign_or_likely_benign"
    return "other"


def build_reference(reference: Path, db: sqlite3.Connection) -> dict[str, object]:
    db.execute("DROP TABLE IF EXISTS annotation")
    db.execute(
        "CREATE TABLE annotation (chrom TEXT, pos INTEGER, ref TEXT, alt TEXT, payload TEXT)"
    )
    count = 0
    multi = 0
    metadata: dict[str, object] = {}
    batch = []
    with vcf_open(reference) as fh:
        for line in fh:
            if line.startswith("##fileDate=") or line.startswith("##reference="):
                key, value = line[2:].rstrip().split("=", 1)
                metadata[key] = value
            if line.startswith("#"):
                continue
            cells = line.rstrip("\n").split("\t")
            annotation = info_map(cells[7])
            alts = cells[4].split(",")
            multi += len(alts) > 1
            # Standard ClinVar is biallelic. If that changes, preserve raw values
            # and flag rather than pretending allele-indexed INFO was resolved.
            payload = json.dumps(
                [cells[2], *[annotation.get(k, "") for k in INFO_FIELDS], len(alts) > 1]
            )
            for alt in alts:
                batch.append((chrom_key(cells[0]), int(cells[1]), cells[3], alt, payload))
                count += 1
            if len(batch) >= 10000:
                db.executemany("INSERT INTO annotation VALUES (?,?,?,?,?)", batch)
                batch.clear()
    if batch:
        db.executemany("INSERT INTO annotation VALUES (?,?,?,?,?)", batch)
    db.execute("CREATE INDEX allele_key ON annotation (chrom,pos,ref,alt)")
    db.commit()
    metadata.update(allele_records=count, multiallelic_records=multi)
    return metadata


def join_case(case: Path, db: sqlite3.Connection, output: Path) -> dict[str, object]:
    counts: Counter[str] = Counter()
    categories: Counter[str] = Counter()
    with vcf_open(case) as source, output.open("w", newline="") as dest:
        writer = csv.writer(dest)
        writer.writerow(FIELDS)
        for line in source:
            if line.startswith("#CHROM") and len(line.rstrip().split("\t")) != 10:
                raise ValueError("Expected exactly one case sample")
            if line.startswith("#"):
                continue
            cells = line.rstrip("\n").split("\t")
            counts["case_records"] += 1
            fmt = dict(zip(cells[8].split(":"), cells[9].split(":"), strict=False))
            gt = fmt.get("GT", ".")
            alleles = re.split(r"[/|]", gt)
            if "." in alleles:
                counts["partial_or_missing_gt_records"] += 1
            alts = cells[4].split(",")
            for allele_index, alt in enumerate(alts, 1):
                counts["case_alt_alleles"] += 1
                dosage = alleles.count(str(allele_index))
                if not dosage:
                    counts["alt_not_in_called_genotype"] += 1
                    continue
                counts["called_alt_alleles"] += 1
                matches = db.execute(
                    "SELECT payload FROM annotation WHERE chrom=? AND pos=? AND ref=? AND alt=?",
                    (chrom_key(cells[0]), int(cells[1]), cells[3], alt),
                ).fetchall()
                if matches:
                    counts["matched_called_alt_alleles"] += 1
                for (payload,) in matches:
                    cv = json.loads(payload)
                    classification = bucket(cv[2])
                    categories[classification] += 1
                    if any(v not in ("", ".") for v in cv[10:13]):
                        counts["annotation_rows_with_any_legacy_population_af"] += 1
                    writer.writerow(
                        [
                            chrom_key(cells[0]),
                            cells[1],
                            cells[3],
                            alt,
                            allele_index,
                            gt,
                            dosage,
                            "phased" if "|" in gt else "unphased",
                            fmt.get("PS", ""),
                            fmt.get("PID", ""),
                            fmt.get("PGT", ""),
                            "partial_or_missing" if "." in alleles else "complete",
                            cells[6],
                            cells[5],
                            fmt.get("DP", ""),
                            fmt.get("AD", ""),
                            fmt.get("GQ", ""),
                            *cv[:-1],
                            classification,
                            len(alts) > 1,
                            cv[-1],
                        ]
                    )
                    counts["annotation_rows"] += 1
    return {**counts, "classification_counts": dict(categories)}


def run(case: Path, reference: Path, output_dir: Path, source_url: str) -> dict[str, object]:
    os.umask(0o077)
    output_dir.mkdir(parents=True, exist_ok=False, mode=0o700)
    output_dir.chmod(0o700)
    start = datetime.now(UTC).isoformat()
    db = sqlite3.connect(output_dir / "clinvar-index.sqlite")
    db.execute("PRAGMA cache_size=-32768")
    db.execute("PRAGMA temp_store=FILE")
    reference_metadata = build_reference(reference, db)
    aggregates = join_case(case, db, output_dir / "clinvar-exact-matches.csv")
    db.close()
    receipt: dict[str, object] = {
        "started_utc": start,
        "finished_utc": datetime.now(UTC).isoformat(),
        "python": sys.version,
        "platform": platform.platform(),
        "command": sys.argv,
        "script_sha256": digest(Path(__file__)),
        "case_input": {"path": str(case), "sha256": digest(case), "bytes": case.stat().st_size},
        "public_reference": {
            "path": str(reference),
            "sha256": digest(reference),
            "bytes": reference.stat().st_size,
            "url": source_url,
            **reference_metadata,
        },
        "output": {
            "path": str(output_dir / "clinvar-exact-matches.csv"),
            "sha256": digest(output_dir / "clinvar-exact-matches.csv"),
            "bytes": (output_dir / "clinvar-exact-matches.csv").stat().st_size,
        },
        "aggregates": aggregates,
        "limitations": [
            "Exact GRCh38 CHROM/POS/REF/ALT representation only; no left "
            "normalization, liftover or haplotype matching.",
            "Multiallelic case alleles included only when present in GT; "
            "original genotype and allele index retained.",
            "Multiallelic reference INFO remains raw and flagged; "
            "allele-specific assignment requires separate review.",
            "Classifications, conflicts and review status are distinct; "
            "pathogenic classification is not patient causality.",
            "Case AF is not population frequency and is unused. Legacy ClinVar "
            "population AF fields are optional and incomplete.",
            "No parent genotypes, segregation, clinical interpretation or "
            "phenotype prioritization performed.",
            "Unmatched variants, structural variants and reference-normalization "
            "differences remain unresolved.",
        ],
    }
    (output_dir / "receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(
        json.dumps(
            {
                "status": "complete",
                "aggregates": aggregates,
                "reference_metadata": reference_metadata,
            }
        )
    )
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", type=Path, required=True)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--source-url", default="https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz"
    )
    args = parser.parse_args()
    run(args.case, args.reference, args.output_dir, args.source_url)


if __name__ == "__main__":
    main()
