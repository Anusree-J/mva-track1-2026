"""Local conservative genotype/condition worksheet; no causal inference."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1048576), b""):
            h.update(chunk)
    return h.hexdigest()


def key(row: dict[str, Any]) -> tuple[str, int, str, str]:
    chrom = row["chrom"].removeprefix("chr")
    return ("MT" if chrom == "M" else chrom, int(row["pos"]), row["ref"], row["alt"])


def genes(value: str) -> list[str]:
    return sorted({g.split(":")[0] for g in value.split("|") if g and g != "."})


def genotype_state(gt: str, alt_index: int) -> str:
    alleles = re.split(r"[/|]", gt)
    if "." in alleles or not all(a.isdigit() for a in alleles):
        return "partial_or_missing"
    dosage = alleles.count(str(alt_index))
    if not dosage:
        return "uncarried"
    if len(alleles) == 1:
        return "haploid_alt_ploidy_context_unresolved"
    if len(alleles) != 2:
        return "non_diploid_context_unresolved"
    return "homozygous_alt" if dosage == 2 else "heterozygous_alt"


def phase_pair(first: dict[str, Any], second: dict[str, Any]) -> str:
    if key(first)[0] != key(second)[0]:
        return "different_chromosomes_not_phase_comparable"
    gt1, gt2 = first["genotype"], second["genotype"]
    idx1, idx2 = int(first["allele_index"]), int(second["allele_index"])
    if (
        genotype_state(gt1, idx1) != "heterozygous_alt"
        or genotype_state(gt2, idx2) != "heterozygous_alt"
    ):
        return "not_two_complete_diploid_heterozygotes"
    if "|" not in gt1 or "|" not in gt2:
        return "unphased_not_trans_evidence"
    ps1, ps2 = first.get("phase_set", ""), second.get("phase_set", "")
    if not ps1 or ps1 == "." or ps1 != ps2:
        return "phase_set_missing_or_different_unresolved"
    hap1, hap2 = gt1.split("|").index(str(idx1)), gt2.split("|").index(str(idx2))
    if hap1 == hap2:
        return "same_PS_same_haplotype_cis_supported_not_independent_biallelic_hit"
    return "same_PS_opposite_haplotypes_supported_not_parental_origin_or_causality"


def inheritance_notes(state: str, terms: set[str]) -> list[str]:
    notes = []
    if "HP:0000007" in terms:
        if state == "homozygous_alt":
            notes.append("AR_homozygous_zygosity_compatible_not_causal_confirmation")
        elif state == "heterozygous_alt":
            notes.append("AR_second_functional_allele_and_phase_unresolved")
        else:
            notes.append("AR_genotype_or_ploidy_context_unresolved")
    if "HP:0000006" in terms:
        if state == "heterozygous_alt":
            notes.append("AD_heterozygous_zygosity_compatible_not_de_novo_proof")
        else:
            notes.append("AD_allele_specific_zygosity_and_ploidy_review_required")
    if terms - {"HP:0000006", "HP:0000007"}:
        notes.append("other_inheritance_annotations_preserved_not_interpreted")
    if not terms:
        notes.append("condition_inheritance_not_resolved")
    return notes


def run(
    clinvar: Path, conditions: Path, targets: Path, gene_panel: Path, output: Path
) -> dict[str, Any]:
    output.mkdir(parents=True, exist_ok=False, mode=0o700)
    assertions = []
    with clinvar.open() as f:
        for row in csv.DictReader(f):
            if row["classification_bucket"] == "pathogenic_or_likely_pathogenic":
                assertions.append(row)
    condition_rows = defaultdict(list)
    with conditions.open() as f:
        for row in csv.DictReader(f):
            condition_rows[(key(row), row["clinvar_id"])].append(row)
    by_gene = defaultdict(list)
    with targets.open() as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if row.get("allele_carriage") != "carried":
                raise ValueError("Counterpart target set must contain carried alleles only")
            by_gene[row["gene_symbol"]].append(row)
    with gene_panel.open() as f:
        panel_genes = {r["gene_symbol"] for r in csv.DictReader(f, delimiter="\t")}
    worksheets: list[dict[str, Any]] = []
    pairs: list[dict[str, Any]] = []
    for assertion in assertions:
        allele = key(assertion)
        state = genotype_state(assertion["genotype"], int(assertion["allele_index"]))
        gene_names = genes(assertion.get("geneinfo", ""))
        matched_conditions = condition_rows[(allele, assertion["clinvar_id"])]
        condition_notes = []
        for row in matched_conditions:
            terms = set(re.findall(r"HP:\d{7}", row.get("explicit_inheritance_hpo_ids", "")))
            negated = set(re.findall(r"HP:\d{7}", row.get("negated_inheritance_hpo_ids", "")))
            contradictory = terms & negated
            note = inheritance_notes(state, terms - negated)
            if contradictory:
                note.append("positive_negative_condition_inheritance_conflict")
            condition_notes.append(
                {
                    "condition_source_row": row,
                    "genotype_screen": note,
                    "HPO_phenotype_qualifier": "unknown",
                    "not_a_variant_condition_assignment": True,
                }
            )
        coverage = {}
        pair_count = 0
        seen = set()
        for gene in gene_names:
            if gene not in panel_genes:
                coverage[gene] = "outside_HPO_gene_panel_counterpart_search_not_covered"
                continue
            coverage[gene] = "targeted_coding_splice_only_no_exhaustive_second_allele_search"
            for counterpart in by_gene[gene]:
                counterpart_key = key(counterpart)
                if counterpart_key == allele or (gene, counterpart_key) in seen:
                    continue
                seen.add((gene, counterpart_key))
                second = {
                    "chrom": counterpart["chrom"],
                    "pos": counterpart["pos"],
                    "ref": counterpart["ref"],
                    "alt": counterpart["alt"],
                    "genotype": counterpart["genotypes"],
                    "allele_index": counterpart["alt_index"],
                    "phase_set": counterpart.get("PS", ""),
                }
                pair = {
                    "index_allele": dict(zip(["chrom", "pos", "ref", "alt"], allele, strict=True)),
                    "clinvar_id": assertion["clinvar_id"],
                    "gene_symbol": gene,
                    "index_assertion_row": assertion,
                    "counterpart_target_row": counterpart,
                    "counterpart_zygosity": genotype_state(
                        second["genotype"], int(second["allele_index"])
                    ),
                    "phase_screen": phase_pair(assertion, second),
                    "interpretation": "counterpart_for_review_not_compound_heterozygous_causality",
                    "counterpart_functional_impact": "unresolved",
                    "counterpart_population_frequency": "not_joined_in_this_worksheet",
                }
                pairs.append(pair)
                pair_count += 1
        worksheets.append(
            {
                "allele": dict(zip(["chrom", "pos", "ref", "alt"], allele, strict=True)),
                "clinvar_id": assertion["clinvar_id"],
                "gene_symbols": gene_names,
                "observed_genotype_state": state,
                "assertion_source_row": assertion,
                "condition_rows_with_screens": condition_notes,
                "counterpart_search_coverage": coverage,
                "counterpart_pair_rows": pair_count,
                "absence_interpretation": "no_observed_counterpart_is_not_absence_of_second_allele",
                "limits": [
                    "No de novo proof, segregation or parent-of-origin inference",
                    "No compound heterozygous diagnosis",
                    "Gene membership does not establish shared condition mechanism",
                    "Reference normalization and allele consequence remain unresolved",
                    "HPO literal presence/negation and condition ambiguity unreviewed",
                    "CNV, deep intronic and outside-panel variants not searched here",
                ],
            }
        )
    with (output / "inheritance-worksheet.jsonl").open("w") as f:
        for row in worksheets:
            f.write(json.dumps(row, sort_keys=True) + "\n")
    with (output / "same-gene-counterpart-review.jsonl").open("w") as f:
        for row in pairs:
            f.write(json.dumps(row, sort_keys=True) + "\n")
    summary = {
        "status": "conservative_local_screen_not_inheritance_diagnosis",
        "P_LP_assertion_rows": len(worksheets),
        "unique_P_LP_alleles": len({key(r) for r in assertions}),
        "genotype_state_counts": dict(Counter(r["observed_genotype_state"] for r in worksheets)),
        "assertion_rows_with_condition_evidence": sum(
            bool(r["condition_rows_with_screens"]) for r in worksheets
        ),
        "condition_screen_counts": dict(
            Counter(
                n
                for r in worksheets
                for c in r["condition_rows_with_screens"]
                for n in c["genotype_screen"]
            )
        ),
        "counterpart_pair_rows": len(pairs),
        "phase_screen_counts": dict(Counter(r["phase_screen"] for r in pairs)),
        "gene_search_coverage_counts": dict(
            Counter(v for r in worksheets for v in r["counterpart_search_coverage"].values())
        ),
        "no_parental_or_trans_claim_from_unphased_GT": True,
        "command": sys.argv,
        "source_sha256": digest(Path(__file__)),
        "inputs": [
            {"path": str(p), "bytes": p.stat().st_size, "sha256": digest(p)}
            for p in [clinvar, conditions, targets, gene_panel]
        ],
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
    return {k: v for k, v in summary.items() if k not in {"command", "inputs"}}


def smoke() -> None:
    assert genotype_state("1/1", 1) == "homozygous_alt"
    assert genotype_state("0/1", 1) == "heterozygous_alt"
    assert genotype_state("1/.", 1) == "partial_or_missing"
    assert genotype_state("1", 1) == "haploid_alt_ploidy_context_unresolved"
    assert genotype_state("0/1", 2) == "uncarried"
    one = {
        "chrom": "1",
        "pos": "10",
        "ref": "A",
        "alt": "G",
        "genotype": "1|0",
        "allele_index": "1",
        "phase_set": "123",
    }
    two = {**one, "pos": "20"}
    assert phase_pair(one, two).startswith("same_PS_same_haplotype")
    assert phase_pair(one, {**two, "genotype": "0|1"}).startswith("same_PS_opposite")
    assert phase_pair(one, {**two, "genotype": "0/1"}) == "unphased_not_trans_evidence"
    assert phase_pair(one, {**two, "phase_set": "456"}).startswith("phase_set_missing")
    assert phase_pair(one, {**two, "genotype": "1|."}).startswith("not_two_complete")
    assert inheritance_notes("heterozygous_alt", {"HP:0000007"}) == [
        "AR_second_functional_allele_and_phase_unresolved"
    ]
    assert "not_de_novo" in inheritance_notes("heterozygous_alt", {"HP:0000006"})[0]
    print("PASS: synthetic homo/het/partial/haploid and same-PS cis/opposite/unphased controls")


def main() -> None:
    os.umask(0o077)
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ["clinvar", "conditions", "targets", "gene-panel", "output"]:
        parser.add_argument("--" + name, type=Path)
    parser.add_argument("--smoke-test", action="store_true")
    args = parser.parse_args()
    if args.smoke_test:
        smoke()
    elif all([args.clinvar, args.conditions, args.targets, args.gene_panel, args.output]):
        print(
            json.dumps(
                run(args.clinvar, args.conditions, args.targets, args.gene_panel, args.output)
            )
        )
    else:
        parser.error("clinvar, conditions, targets, gene-panel and output required")


if __name__ == "__main__":
    main()
