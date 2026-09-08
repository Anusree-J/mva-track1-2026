"""Private representation-conditional CDS length-change screen, not pathogenicity."""

import argparse
import csv
import gzip
import hashlib
import json
import os
import re
from collections import Counter, defaultdict
from pathlib import Path

os.umask(0o077)
MARGIN = 3


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for c in iter(lambda: f.read(1048576), b""):
            h.update(c)
    return h.hexdigest()


def key(row):
    c = row["chrom"].removeprefix("chr")
    return ("MT" if c == "M" else c, int(row["pos"]), row["ref"], row["alt"])


def changed_span(pos, ref, alt):
    if not re.fullmatch("[ACGTacgt]+", ref) or not re.fullmatch("[ACGTacgt]+", alt):
        return None
    ref, alt = ref.upper(), alt.upper()
    prefix = 0
    while prefix < min(len(ref), len(alt)) and ref[prefix] == alt[prefix]:
        prefix += 1
    changed_ref, changed_alt = ref[prefix:], alt[prefix:]
    while changed_ref and changed_alt and changed_ref[-1] == changed_alt[-1]:
        changed_ref, changed_alt = changed_ref[:-1], changed_alt[:-1]
    start = pos + prefix
    if not changed_ref and not changed_alt:
        return None
    return {
        "changed_start": start,
        "changed_reference_length": len(changed_ref),
        "changed_alternate_length": len(changed_alt),
        "net_length_change": len(changed_alt) - len(changed_ref),
        "common_prefix_trimmed": prefix,
        "checked_left_flank": start - 1,
        "checked_right_flank": start + len(changed_ref) if changed_ref else start,
    }


def segment_effect(span, left, right):
    if span is None:
        return "unknown_non_sequence_or_no_change"
    # Both unchanged flanks must also be strictly internal with a 3-base buffer.
    if span["checked_left_flank"] < left + MARGIN or span["checked_right_flank"] > right - MARGIN:
        return "unknown_boundary_or_outside_same_CDS_segment"
    if span["net_length_change"] % 3:
        return "conditional_predicted_frameshift"
    return "conditional_in_frame_length_change"


def smoke():
    insertion = changed_span(20, "A", "AT")
    deletion = changed_span(20, "AT", "A")
    assert insertion["checked_left_flank"] == 20 and insertion["checked_right_flank"] == 21
    assert deletion["changed_start"] == 21 and deletion["checked_right_flank"] == 22
    assert segment_effect(insertion, 10, 40) == "conditional_predicted_frameshift"
    assert segment_effect(deletion, 10, 40) == "conditional_predicted_frameshift"
    assert (
        segment_effect(changed_span(20, "A", "ATTT"), 10, 40)
        == "conditional_in_frame_length_change"
    )
    assert segment_effect(changed_span(10, "AT", "A"), 10, 40).startswith("unknown_boundary")
    assert segment_effect(changed_span(38, "ATTTT", "A"), 10, 40).startswith("unknown_boundary")
    assert changed_span(20, "AAC", "ATC")["changed_reference_length"] == 1
    print("PASS: synthetic insertion/deletion, modulo3, boundaries and allele trimming")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("targets", "gtf", "evidence", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    ROOT, TARGET, GTF, EVIDENCE = args.output, args.targets, args.gtf, args.evidence
    ROOT.mkdir(mode=0o700, parents=True, exist_ok=True)
    for name in ["conditional-coding-indel-worksheet.jsonl", "summary.json", "receipt.json"]:
        if (ROOT / name).exists():
            raise FileExistsError("Frozen output exists; choose a fresh script/output directory")
    smoke()
    targets = defaultdict(list)
    with TARGET.open() as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if row["class"] == "indel":
                assert row["allele_carriage"] == "carried"
                targets[key(row)].append(row)
    genes = {r["gene_symbol"] for rows in targets.values() for r in rows}
    cds = defaultdict(list)
    with gzip.open(GTF, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if cols[2] != "CDS":
                continue
            attrs = dict(re.findall(r'(\w+) "([^"]+)"', cols[8]))
            gene = attrs.get("gene_name", "")
            if gene not in genes:
                continue
            chrom = cols[0].removeprefix("chr")
            chrom = "MT" if chrom == "M" else chrom
            cds[(chrom, gene)].append(
                {
                    "start": int(cols[3]),
                    "end": int(cols[4]),
                    "transcript_id": attrs["transcript_id"],
                    "strand": cols[6],
                    "CDS_phase": cols[7],
                }
            )
    annotations = {}
    with EVIDENCE.open() as f:
        for line in f:
            row = json.loads(line)
            if key(row) in targets:
                annotations[key(row)] = row
    records = []
    for allele, rows in sorted(targets.items()):
        span = changed_span(allele[1], allele[2], allele[3])
        effects = []
        overlapping_segments = 0
        for gene in sorted({r["gene_symbol"] for r in rows}):
            for segment in cds[(allele[0], gene)]:
                if span is None:
                    continue
                if (
                    segment["end"] < span["checked_left_flank"]
                    or segment["start"] > span["checked_right_flank"]
                ):
                    continue
                overlapping_segments += 1
                effect = segment_effect(span, segment["start"], segment["end"])
                effects.append({"gene_symbol": gene, **segment, "effect": effect})
        labels = {r["effect"] for r in effects}
        status = (
            "conditional_predicted_frameshift"
            if "conditional_predicted_frameshift" in labels
            else "conditional_in_frame_length_change"
            if "conditional_in_frame_length_change" in labels
            else "unknown_boundary_or_no_clear_same_CDS_segment"
        )
        a = annotations.get(allele, {})
        records.append(
            {
                "allele": dict(zip(["chrom", "pos", "ref", "alt"], allele, strict=True)),
                "status": status,
                "status_scope": "any_overlapping_transcript_segment",
                "changed_span": span,
                "transcript_segment_effects": effects,
                "overlapping_CDS_segment_count": overlapping_segments,
                "original_target_rows": rows,
                "frequency_status": a.get("frequency_status", "unknown"),
                "frequency_evidence": a.get("frequency_evidence", []),
                "clinical_categories": a.get("clinical_categories", []),
                "clinvar_assertion_rows": a.get("clinvar_assertion_rows", []),
                "quality_flags": a.get("quality_flags", []),
                "reference_base_validation": "unknown_no_FASTA",
                "left_normalization": "not_performed",
                "NMD_or_protein_effect": "unknown",
                "clinical_pathogenicity_from_this_screen": False,
            }
        )
    output = ROOT / "conditional-coding-indel-worksheet.jsonl"
    output.write_text("".join(json.dumps(r, sort_keys=True) + "\n" for r in records))
    summary = {
        "target_rows": sum(map(len, targets.values())),
        "unique_carried_indels": len(records),
        "status_counts": dict(Counter(r["status"] for r in records)),
        "frequency_status_counts": dict(Counter(r["frequency_status"] for r in records)),
        "conditional_frameshift_frequency_counts": dict(
            Counter(
                r["frequency_status"]
                for r in records
                if r["status"] == "conditional_predicted_frameshift"
            )
        ),
        "conditional_frameshift_without_exact_ClinVar_count": sum(
            r["status"] == "conditional_predicted_frameshift" and not r["clinvar_assertion_rows"]
            for r in records
        ),
        "CDS_margin_bases": MARGIN,
        "scope": "HPO-targeted carried coding/splice indels only",
        "causal_or_clinical_interpretation": "none",
    }
    (ROOT / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    receipt = {
        "status": "PASS",
        "sources": [{"path": str(p), "sha256": digest(p)} for p in [TARGET, GTF, EVIDENCE]],
        "method_sha256": digest(Path(__file__)),
        "outputs": [
            {"name": p.name, "bytes": p.stat().st_size, "sha256": digest(p)}
            for p in [output, ROOT / "summary.json"]
        ],
    }
    (ROOT / "receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
    for p in ROOT.iterdir():
        p.chmod(0o600)
    print(json.dumps(summary))


if __name__ == "__main__":
    main()
