"""Local condition-level HPO literal overlap; does not assert case phenotypes."""

from __future__ import annotations

import argparse
import csv
import json
import os
import platform
import re
import sys
from collections import Counter, defaultdict
from datetime import UTC, datetime
from pathlib import Path
from urllib.parse import unquote

from science_campaign.mva_clinvar import digest

TARGET_BUCKETS = {"pathogenic_or_likely_pathogenic", "conflicting", "uncertain"}
FIELDS = [
    "chrom",
    "pos",
    "ref",
    "alt",
    "clinvar_id",
    "classification_bucket",
    "clnsig",
    "clnrevstat",
    "clnsigconf",
    "condition_group_index",
    "condition_id",
    "condition_name_raw",
    "condition_group_raw",
    "all_clndn_raw",
    "all_clndisdb_raw",
    "condition_group_count",
    "condition_identifier_count",
    "assertion_condition_ambiguity",
    "condition_lookup_status",
    "hpoa_disease_name",
    "exact_literal_overlap_count",
    "exact_literal_overlap_ids",
    "negated_disease_literal_overlap_ids",
    "explicit_inheritance_hpo_ids",
    "negated_inheritance_hpo_ids",
    "inheritance_status",
    "case_phenotype_qualifier_status",
    "support_interpretation",
    "hpoa_annotation_references",
]


def read_hpoa(path: Path) -> tuple[dict[str, list[dict[str, str]]], dict[str, object]]:
    records: dict[str, list[dict[str, str]]] = defaultdict(list)
    aspects: Counter[str] = Counter()
    with path.open() as source:
        for row in csv.DictReader(
            (line for line in source if not line.startswith("#")), delimiter="\t"
        ):
            records[row["database_id"]].append(row)
            aspects[row["aspect"]] += 1
    return dict(records), {"aspect_counts": dict(aspects), "disease_count": len(records)}


def identifiers(group: str) -> list[str]:
    found = re.findall(r"(?:^|,)(OMIM|Orphanet|ORPHA):(?:ORPHA:)?(\d+)(?=,|$)", unquote(group))
    return sorted({f"{'OMIM' if db == 'OMIM' else 'ORPHA'}:{number}" for db, number in found})


def annotate(
    candidate: dict[str, str], literals: set[str], database: dict[str, list[dict[str, str]]]
) -> list[dict[str, str | int]]:
    groups = candidate.get("clndisdb", "").split("|")
    names = candidate.get("clndn", "").split("|")
    results = []
    for index, group in enumerate(groups):
        ids = identifiers(group)
        ambiguous = len(groups) > 1 or len(ids) > 1 or len(names) != len(groups)
        for condition_id in ids or [""]:
            annotation = database.get(condition_id, [])
            positive = {
                r["hpo_id"] for r in annotation if r["aspect"] == "P" and not r["qualifier"]
            }
            negative = {
                r["hpo_id"] for r in annotation if r["aspect"] == "P" and r["qualifier"] == "NOT"
            }
            inherited = {
                r["hpo_id"] for r in annotation if r["aspect"] == "I" and not r["qualifier"]
            }
            neg_inherited = {
                r["hpo_id"] for r in annotation if r["aspect"] == "I" and r["qualifier"] == "NOT"
            }
            overlap = positive & literals
            matched = [r for r in annotation if r["hpo_id"] in overlap or r["aspect"] == "I"]
            result: dict[str, str | int] = {k: candidate.get(k, "") for k in FIELDS[:9]}
            result.update(
                {
                    "condition_group_index": index + 1,
                    "condition_id": condition_id,
                    "condition_name_raw": names[index] if len(names) == len(groups) else "",
                    "condition_group_raw": group,
                    "all_clndn_raw": candidate.get("clndn", ""),
                    "all_clndisdb_raw": candidate.get("clndisdb", ""),
                    "condition_group_count": len(groups),
                    "condition_identifier_count": len(ids),
                    "assertion_condition_ambiguity": str(ambiguous).lower(),
                    "condition_lookup_status": "matched"
                    if annotation
                    else "unmatched_identifier"
                    if ids
                    else "no_supported_identifier",
                    "hpoa_disease_name": "|".join(sorted({r["disease_name"] for r in annotation})),
                    "exact_literal_overlap_count": len(overlap),
                    "exact_literal_overlap_ids": "|".join(sorted(overlap)),
                    "negated_disease_literal_overlap_ids": "|".join(sorted(negative & literals)),
                    "explicit_inheritance_hpo_ids": "|".join(sorted(inherited)),
                    "negated_inheritance_hpo_ids": "|".join(sorted(neg_inherited)),
                    "inheritance_status": "explicit_public_annotation"
                    if inherited
                    else "unavailable",
                    "case_phenotype_qualifier_status": "unreviewed_literal_only",
                    "support_interpretation": "literal_overlap_only_not_clinical_support",
                    "hpoa_annotation_references": "|".join(
                        sorted({r["reference"] for r in matched})
                    ),
                }
            )
            results.append(result)
    return results


def run(clinvar: Path, literals_path: Path, hpoa: Path, output_dir: Path) -> dict[str, object]:
    os.umask(0o077)
    output_dir.mkdir(parents=True, exist_ok=False, mode=0o700)
    started = datetime.now(UTC).isoformat()
    literals_raw = json.loads(literals_path.read_text())
    if not isinstance(literals_raw, list) or not all(
        re.fullmatch(r"HP:\d{7}", x) for x in literals_raw
    ):
        raise ValueError("Expected a JSON array of literal HPO IDs")
    literals = set(literals_raw)
    database, metadata = read_hpoa(hpoa)
    counts: Counter[str] = Counter()
    output = output_dir / "condition-evidence.csv"
    with clinvar.open() as source, output.open("w", newline="") as dest:
        writer = csv.DictWriter(dest, fieldnames=FIELDS)
        writer.writeheader()
        for candidate in csv.DictReader(source):
            if candidate["classification_bucket"] not in TARGET_BUCKETS:
                continue
            counts["input_target_alleles"] += 1
            rows = annotate(candidate, literals, database)
            counts["alleles_with_any_matched_condition"] += any(
                r["condition_lookup_status"] == "matched" for r in rows
            )
            counts["alleles_with_any_exact_literal_overlap"] += any(
                r["exact_literal_overlap_count"] for r in rows
            )
            for row in rows:
                writer.writerow(row)
                counts["condition_rows"] += 1
                counts[str(row["condition_lookup_status"])] += 1
                counts["rows_with_explicit_public_inheritance"] += (
                    row["inheritance_status"] == "explicit_public_annotation"
                )
    receipt: dict[str, object] = {
        "started_utc": started,
        "finished_utc": datetime.now(UTC).isoformat(),
        "command": sys.argv,
        "python": sys.version,
        "platform": platform.platform(),
        "script_sha256": digest(Path(__file__)),
        "digest_dependency_sha256": digest(Path(__file__).with_name("mva_clinvar.py")),
        "inputs": [
            {"path": str(p), "sha256": digest(p), "bytes": p.stat().st_size}
            for p in (clinvar, literals_path, hpoa)
        ],
        "output": {"path": str(output), "sha256": digest(output), "bytes": output.stat().st_size},
        "aggregates": dict(counts),
        "public_hpoa_metadata": metadata,
        "limitations": [
            "Case HPO literals have unreviewed presence/absence qualifiers; "
            "overlap is not positive clinical evidence.",
            "Exact condition and term matching only; no ontology ancestor "
            "expansion or disease-name fuzzy match.",
            "OMIM/Orphanet identifiers only. Unmatched and unsupported IDs "
            "remain in raw condition fields.",
            "Multiple assertion conditions or identifier aliases are marked "
            "ambiguous and never merged for overlap.",
            "Inheritance terms are explicit HPOA aspect I annotations, not "
            "inferred case inheritance fit.",
            "No inference of segregation, phase, de novo status, penetrance, or causality.",
        ],
    }
    (output_dir / "receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(
        json.dumps(
            {"status": "complete", "aggregates": dict(counts), "public_hpoa_metadata": metadata}
        )
    )
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--clinvar", type=Path, required=True)
    parser.add_argument("--hpo-literals", type=Path, required=True)
    parser.add_argument("--hpoa", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    run(args.clinvar, args.hpo_literals, args.hpoa, args.output_dir)


if __name__ == "__main__":
    main()
