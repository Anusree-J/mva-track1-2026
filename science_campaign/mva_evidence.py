"""Local MVA evidence consolidation; review priorities are not causal scores."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import sys
from collections import Counter, defaultdict
from datetime import UTC, datetime
from pathlib import Path
from typing import Any


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1048576), b""):
            h.update(chunk)
    return h.hexdigest()


def key(row: dict[str, Any]) -> tuple[str, int, str, str]:
    chrom = row["chrom"].removeprefix("chr")
    return ("MT" if chrom == "M" else chrom, int(row["pos"]), row["ref"], row["alt"])


def review_level(value: str) -> int:
    value = value.lower()
    if "practice_guideline" in value:
        return 4
    if "reviewed_by_expert_panel" in value:
        return 3
    if "multiple_submitters" in value and "no_conflicts" in value:
        return 2
    if "criteria_provided" in value:
        return 1
    return 0


def numeric_af(value: str) -> float | None:
    try:
        af = float(value)
    except (ValueError, TypeError):
        return None
    return af if 0 <= af <= 1 else None


def consolidate(
    targets: list[dict[str, Any]], clinical: list[dict[str, Any]], population: list[dict[str, Any]]
) -> dict[str, Any]:
    """Keep every source row and expose conservative, explicit review ordering."""
    categories = {r.get("classification_bucket", "") for r in clinical}
    pathogenic = "pathogenic_or_likely_pathogenic" in categories
    benign = "benign_or_likely_benign" in categories
    conflict = (
        "conflicting" in categories
        or (pathogenic and benign)
        or any(r.get("clnsigconf", "") not in {"", "."} for r in clinical)
    )
    supported = max(
        (
            review_level(r.get("clnrevstat", ""))
            for r in clinical
            if r.get("classification_bucket") == "pathogenic_or_likely_pathogenic"
        ),
        default=0,
    )
    af_evidence = []
    excluded_legacy = 0
    for row in clinical:
        for source in ["af_exac", "af_esp", "af_tgp"]:
            value = numeric_af(row.get(source, ""))
            if value is None:
                continue
            if row.get("clinvar_multiallelic", "").lower() == "true":
                excluded_legacy += 1
                continue
            af_evidence.append(
                {
                    "source": "ClinVar_legacy_" + source,
                    "value": value,
                    "match": "exact_allele",
                    "clinvar_id": row.get("clinvar_id", ""),
                }
            )
    for row in population:
        if row.get("population_match") != "exact_allele_match":
            continue
        for column, text in row.items():
            if column.startswith("AF_") and (value := numeric_af(text)) is not None:
                af_evidence.append(
                    {
                        "source": "population_resource_" + column,
                        "value": value,
                        "match": "exact_allele",
                    }
                )
    maximum_af = max((a["value"] for a in af_evidence), default=None)
    frequency = (
        "unknown"
        if maximum_af is None
        else "common_gt_1pct"
        if maximum_af > 0.01
        else "observed_AF_le_1pct_not_proof_of_rarity"
    )
    common = frequency == "common_gt_1pct"
    hpo_count = max((int(t.get("exact_hpo_count", 0)) for t in targets), default=0)
    if (pathogenic or conflict) and common:
        group = "05_assertion_frequency_discordance_review"
    elif conflict:
        group = "40_conflicting_assertions_review"
    elif pathogenic and targets and supported >= 2:
        group = "10_supported_pathogenic_assertion_with_unreviewed_HPO_overlap"
    elif pathogenic and targets:
        group = "20_pathogenic_assertion_with_unreviewed_HPO_overlap"
    elif pathogenic:
        group = "30_pathogenic_assertion_outside_HPO_coding_splice_panel"
    elif common:
        group = "90_common_allele_low_priority"
    elif benign:
        group = "80_benign_assertion_low_priority"
    elif "uncertain" in categories:
        group = "60_uncertain_assertion_coding_splice_review"
    else:
        group = "50_no_pathogenic_assertion_coding_splice_review"
    quality_flags = set()
    for row in targets:
        if row.get("filter") != "PASS":
            quality_flags.add("non_PASS_retained")
        if row.get("allele_carriage") != "carried":
            quality_flags.add("incomplete_or_unconfirmed_carriage")
    for row in clinical:
        if row.get("filter") != "PASS":
            quality_flags.add("non_PASS_retained")
        if row.get("genotype_completeness") != "complete":
            quality_flags.add("incomplete_GT_retained")
    genotype_rows = [
        {
            "GT": r.get("genotypes", ""),
            "AD": r.get("AD", ""),
            "DP": r.get("DP", ""),
            "GQ": r.get("GQ", ""),
            "FILTER": r.get("filter", ""),
            "QUAL": r.get("qual", ""),
            "PS": r.get("PS", ""),
            "PID": r.get("PID", ""),
            "PGT": r.get("PGT", ""),
            "allele_carriage": r.get("allele_carriage", ""),
        }
        for r in targets
    ]
    genotype_rows += [
        {
            "GT": r.get("genotype", ""),
            "AD": r.get("allele_depths", ""),
            "DP": r.get("depth", ""),
            "GQ": r.get("genotype_quality", ""),
            "FILTER": r.get("filter", ""),
            "QUAL": r.get("qual", ""),
            "PS": r.get("phase_set", ""),
            "PID": r.get("phase_id", ""),
            "PGT": r.get("physical_phase_genotype", ""),
            "genotype_completeness": r.get("genotype_completeness", ""),
        }
        for r in clinical
    ]
    for row in genotype_rows:
        for field, floor in [("DP", 10), ("GQ", 20)]:
            try:
                value = float(row.get(field, ""))
            except ValueError:
                quality_flags.add(field + "_missing")
            else:
                if value < floor:
                    quality_flags.add(field + "_below_" + str(floor))
    missing = [
        "normalized_reference_validation",
        "HPO_presence_negation_review",
        "transcript_consequence_prediction",
        "inheritance_segregation",
        "phase_trans_configuration",
        "orthogonal_functional_evidence",
    ]
    if not clinical:
        missing.append("exact_ClinVar_match_not_found_not_benign")
    if maximum_af is None:
        missing.append("population_frequency_unknown_not_rare")
    return {
        "review_priority_group": group,
        "maximum_pathogenic_review_level": supported,
        "exact_hpo_count_unreviewed": hpo_count,
        "HPO_qualifier": "unknown",
        "frequency_status": frequency,
        "maximum_observed_matching_AF": maximum_af,
        "frequency_evidence": af_evidence,
        "excluded_ambiguous_legacy_AF_count": excluded_legacy,
        "clinical_categories": sorted(categories),
        "assertion_conflict": conflict,
        "gene_symbols": sorted({r["gene_symbol"] for r in targets}),
        "location_flags": sorted(
            {flag for r in targets for flag in r.get("location_flags", "").split(";") if flag}
        ),
        "variant_classes": sorted({r.get("class", "") for r in targets}),
        "quality_flags": sorted(quality_flags),
        "genotype_evidence": genotype_rows,
        "phenotype_panel_rows": targets,
        "clinvar_assertion_rows": clinical,
        "population_rows": population,
        "missing_evidence": missing,
        "outside_panel_review": not targets,
    }


def run(
    target: Path, clinvar: Path, output: Path, population: Path | None = None
) -> dict[str, Any]:
    output.mkdir(parents=True, exist_ok=False, mode=0o700)
    panel: dict[Any, Any] = defaultdict(list)
    with target.open() as stream:
        for row in csv.DictReader(stream, delimiter="\t"):
            if row.get("allele_carriage") != "carried":
                raise ValueError("Target input must contain carried alleles only")
            panel[key(row)].append(row)
    assertions: dict[Any, Any] = defaultdict(list)
    outside = set()
    clinvar_total = 0
    with clinvar.open() as stream:
        for row in csv.DictReader(stream):
            clinvar_total += 1
            allele = key(row)
            carried = row.get("alt_dosage", "")
            if carried and int(carried) <= 0:
                raise ValueError("ClinVar row is not a carried ALT")
            if allele in panel or row.get("classification_bucket") in {
                "pathogenic_or_likely_pathogenic",
                "conflicting",
            }:
                assertions[allele].append(row)
                if allele not in panel:
                    outside.add(allele)
    # Keep other assertions for every outside allele selected on any P/LP/conflict row.
    if outside:
        assertions = defaultdict(list)
        with clinvar.open() as stream:
            for row in csv.DictReader(stream):
                if (allele := key(row)) in panel or allele in outside:
                    assertions[allele].append(row)
    populations: dict[Any, Any] = defaultdict(list)
    if population:
        with population.open() as stream:
            for row in csv.DictReader(stream):
                populations[key(row)].append(row)
    records = []
    for allele in panel.keys() | outside:
        record = consolidate(panel[allele], assertions[allele], populations[allele])
        record.update(dict(zip(["chrom", "pos", "ref", "alt"], allele, strict=True)))
        records.append(record)
    records.sort(
        key=lambda r: (
            r["review_priority_group"],
            -r["maximum_pathogenic_review_level"],
            -r["exact_hpo_count_unreviewed"],
            r["chrom"],
            r["pos"],
            r["ref"],
            r["alt"],
        )
    )
    columns = [
        "review_order",
        "review_priority_group",
        "chrom",
        "pos",
        "ref",
        "alt",
        "gene_symbols",
        "variant_classes",
        "location_flags",
        "exact_hpo_count_unreviewed",
        "HPO_qualifier",
        "clinical_categories",
        "maximum_pathogenic_review_level",
        "assertion_conflict",
        "frequency_status",
        "maximum_observed_matching_AF",
        "quality_flags",
        "outside_panel_review",
        "missing_evidence",
        "genotype_evidence",
        "clinvar_assertion_rows",
        "frequency_evidence",
        "phenotype_panel_rows",
    ]
    with (
        (output / "evidence-ledger.jsonl").open("w") as ledger,
        (output / "review-priority.csv").open("w") as csv_file,
        (output / "outside-panel-review.csv").open("w") as outside_file,
    ):
        writer = csv.DictWriter(csv_file, fieldnames=columns)
        outside_writer = csv.DictWriter(outside_file, fieldnames=columns)
        writer.writeheader()
        outside_writer.writeheader()
        for rank, record in enumerate(records, 1):
            record["review_order"] = rank
            ledger.write(json.dumps(record, sort_keys=True) + "\n")
            row = {
                k: json.dumps(record[k], sort_keys=True)
                if isinstance(record[k], dict | list)
                else record[k]
                for k in columns
            }
            writer.writerow(row)
            if record["outside_panel_review"]:
                outside_writer.writerow(row)
    summary = {
        "status": "local_evidence_consolidated_not_submission_ready",
        "unique_panel_alleles": len(panel) - len(outside),
        "outside_panel_alleles": len(outside),
        "total_review_alleles": len(records),
        "clinvar_input_rows": clinvar_total,
        "review_group_counts": dict(Counter(r["review_priority_group"] for r in records)),
        "frequency_status_counts": dict(Counter(r["frequency_status"] for r in records)),
        "panel_alleles_with_exact_ClinVar_match": sum(
            bool(r["clinvar_assertion_rows"]) for r in records if not r["outside_panel_review"]
        ),
        "quality_flag_counts": dict(Counter(f for r in records for f in r["quality_flags"])),
        "population_join_supplied": population is not None,
        "interpretation": "Heuristic review order; not a disease or causal ranking",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "command": sys.argv,
        "input_receipts": [
            {"path": str(p), "bytes": p.stat().st_size, "sha256": digest(p)}
            for p in [target, clinvar, *([population] if population else [])]
        ],
        "source_sha256": digest(Path(__file__)),
    }
    (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    (output / "method_snapshot.py").write_bytes(Path(__file__).read_bytes())
    receipt = {
        "status": "PASS",
        "files": [
            {"name": p.name, "bytes": p.stat().st_size, "sha256": digest(p)}
            for p in sorted(output.iterdir())
        ],
    }
    (output / "receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
    for p in output.iterdir():
        p.chmod(0o600)
    return {k: v for k, v in summary.items() if k not in {"command", "input_receipts"}}


def smoke() -> None:
    target = {
        "gene_symbol": "SYNTH",
        "exact_hpo_count": "2",
        "class": "indel",
        "location_flags": "CDS",
        "allele_carriage": "carried",
        "filter": "PASS",
        "DP": "30",
        "GQ": "50",
    }
    benign = {
        "classification_bucket": "benign_or_likely_benign",
        "af_tgp": "0.02",
        "clinvar_multiallelic": "False",
        "genotype_completeness": "complete",
    }
    pathogenic = {
        "classification_bucket": "pathogenic_or_likely_pathogenic",
        "clnrevstat": "reviewed_by_expert_panel",
        "af_tgp": ".",
        "genotype_completeness": "complete",
    }
    assert consolidate([target], [benign], [])["review_priority_group"].startswith("90_")
    assert consolidate([target], [pathogenic], [])["review_priority_group"].startswith("10_")
    assert consolidate([], [pathogenic], [])["review_priority_group"].startswith("30_")
    assert consolidate([target], [pathogenic, benign], [])["assertion_conflict"]
    assert consolidate([target], [], [])["frequency_status"] == "unknown"
    ambiguous = {**benign, "clinvar_multiallelic": "True"}
    assert consolidate([target], [ambiguous], [])["frequency_status"] == "unknown"
    assert (
        consolidate(
            [target], [pathogenic], [{"population_match": "unknown_not_found", "AF_EUR": "0.5"}]
        )["frequency_status"]
        == "unknown"
    )
    assert consolidate(
        [target], [pathogenic], [{"population_match": "exact_allele_match", "AF_EUR": "0.5"}]
    )["review_priority_group"].startswith("05_")
    print("PASS: synthetic priorities, conflicts, unknown AF and exact-allele AF gating")


def main() -> None:
    os.umask(0o077)
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target", type=Path)
    parser.add_argument("--clinvar", type=Path)
    parser.add_argument("--population", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--smoke-test", action="store_true")
    args = parser.parse_args()
    if args.smoke_test:
        smoke()
    elif args.target and args.clinvar and args.output:
        print(json.dumps(run(args.target, args.clinvar, args.output, args.population)))
    else:
        parser.error("target, clinvar and output required")


if __name__ == "__main__":
    main()
