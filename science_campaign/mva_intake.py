"""Private, local-only streaming MVA structural intake (no causal ranking).

Callset AF is deliberately never interpreted as population frequency. Outputs
contain sensitive case data and must remain outside Git and hosted-model context.
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
import struct
import sys
import time
import xml.etree.ElementTree as ET
import zipfile
from collections import Counter
from pathlib import Path
from typing import Any


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def carriage(gt: str, alt_index: int) -> str:
    alleles = re.split(r"[/|]", gt)
    if str(alt_index) in alleles:
        return "carried_partial" if "." in alleles else "carried"
    return "missing_or_partial" if "." in alleles else "uncarried"


def variant_class(ref: str, alt: str) -> str:
    if alt == "*":
        return "spanning_deletion"
    if alt.startswith("<") or "[" in alt or "]" in alt:
        return "symbolic_or_breakend"
    if len(ref) == len(alt) == 1:
        return "SNV"
    if len(ref) == len(alt):
        return "MNV"
    return "indel"


def intake(source: Path, output: Path) -> dict[str, Any]:
    output.mkdir(mode=0o700, parents=True, exist_ok=False)
    os.chmod(output, 0o700)
    vcfs = list(source.glob("*.vcf.gz"))
    docs = list(source.glob("*.docx"))
    if len(vcfs) != 1 or len(docs) != 1:
        raise ValueError("Expected exactly one VCF and one phenotype DOCX")
    vcf, doc = vcfs[0], docs[0]
    index = Path(str(vcf) + ".tbi")
    start = time.monotonic()
    summary: dict[str, Any] = {
        "schema_version": 1,
        "record_count": 0,
        "alternate_allele_count": 0,
        "sample_count": None,
        "info_fields": [],
        "format_fields": [],
        "annotation_available": False,
        "candidate_ranking_status": "blocked_missing_functional_population_annotations",
        "reference_build": "unresolved",
        "reference_build_validation": "header_only",
        "population_frequency_available": False,
        "callset_AF_is_not_population_AF": True,
    }
    counts: dict[str, Counter[str]] = {
        key: Counter()
        for key in [
            "classes",
            "filters",
            "contigs",
            "genotype_states",
            "missing_format_fields",
            "quality_bins",
        ]
    }
    info_fields: set[str] = set()
    formats: set[str] = set()
    header_build = False
    with gzip.open(index, "rb") as stream:
        magic = stream.read(4)
        if magic != b"TBI\x01":
            raise ValueError("Invalid TBI magic")
        header = struct.unpack("<8i", stream.read(32))
        names = stream.read(header[7]).split(b"\x00")
        summary["index_reference_count"] = header[0]
        index_names = {n.decode() for n in names if n}
    with (output / "variant_evidence_groundwork.tsv").open("w") as out:
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "record_number",
                "chrom",
                "pos",
                "ref",
                "alt",
                "alt_index",
                "class",
                "filter",
                "qual",
                "genotypes",
                "DP",
                "GQ",
                "AD",
                "PS",
                "PID",
                "PGT",
                "allele_carriage",
                "tier",
                "reason",
            ]
        )
        with gzip.open(vcf, "rt") as stream:
            for line in stream:
                if line.startswith("##"):
                    if line.startswith("##reference=") and "GRCh38" in line:
                        header_build = True
                    match = re.match(r"##(INFO|FORMAT)=<ID=([^,>]+)", line)
                    if match:
                        (info_fields if match[1] == "INFO" else formats).add(match[2])
                    continue
                if line.startswith("#CHROM"):
                    summary["sample_count"] = max(0, len(line.rstrip("\n").split("\t")) - 9)
                    continue
                if line.startswith("#"):
                    continue
                columns = line.rstrip("\n").split("\t")
                if len(columns) < 8:
                    raise ValueError("Invalid VCF record structure")
                chrom, pos, _, ref, alts, qual, filt = columns[:7]
                int(pos)
                summary["record_count"] += 1
                counts["contigs"][chrom] += 1
                counts["filters"][filt] += 1
                fmt = columns[8].split(":") if len(columns) > 8 else []
                samples = [dict(zip(fmt, s.split(":"), strict=False)) for s in columns[9:]]
                genotypes = [s.get("GT", ".") for s in samples]
                for sample, gt in zip(samples, genotypes, strict=True):
                    alleles = re.split(r"[/|]", gt)
                    state = (
                        "missing_or_partial"
                        if "." in alleles
                        else "hom_ref"
                        if set(alleles) == {"0"}
                        else "hom_alt_or_haploid_alt"
                        if len(set(alleles)) == 1
                        else "heterozygous"
                    )
                    counts["genotype_states"][state] += 1
                    for field in ["GT", "DP", "GQ", "AD"]:
                        if sample.get(field, ".") in ["", "."]:
                            counts["missing_format_fields"][field] += 1
                if qual == ".":
                    counts["quality_bins"]["missing"] += 1
                else:
                    counts["quality_bins"]["below30" if float(qual) < 30 else "atleast30"] += 1
                for alt_index, alt in enumerate(alts.split(","), 1):
                    kind = variant_class(ref, alt)
                    counts["classes"][kind] += 1
                    summary["alternate_allele_count"] += 1
                    # All classes retained. This is an evidence queue, not a causal score.
                    writer.writerow(
                        [
                            summary["record_count"],
                            chrom,
                            pos,
                            ref,
                            alt,
                            alt_index,
                            kind,
                            filt,
                            qual,
                            ";".join(genotypes),
                            ";".join(s.get("DP", ".") for s in samples),
                            ";".join(s.get("GQ", ".") for s in samples),
                            ";".join(s.get("AD", ".") for s in samples),
                            ";".join(s.get("PS", ".") for s in samples),
                            ";".join(s.get("PID", ".") for s in samples),
                            ";".join(s.get("PGT", ".") for s in samples),
                            ";".join(carriage(gt, alt_index) for gt in genotypes),
                            "UNRANKED",
                            "requires_allele_normalization_and_annotation",
                        ]
                    )
    with zipfile.ZipFile(doc) as archive:
        root = ET.fromstring(archive.read("word/document.xml"))
    ns = {"w": "http://schemas.openxmlformats.org/wordprocessingml/2006/main"}
    paragraphs = [
        "".join(t.text or "" for t in p.findall(".//w:t", ns)) for p in root.findall(".//w:p", ns)
    ]
    text = "\n".join(paragraphs)
    hpo = sorted(set(re.findall(r"HP:\d{7}", text)))
    (output / "phenotype-private.txt").write_text(text)
    (output / "hpo-explicit-private.json").write_text(json.dumps(hpo, indent=2) + "\n")
    summary.update({key: dict(value) for key, value in counts.items()})
    summary.update(
        {
            "info_fields": sorted(info_fields),
            "format_fields": sorted(formats),
            "reference_build": "GRCh38" if header_build else "unresolved",
            "index_covers_observed_contig_names": set(counts["contigs"]) <= index_names,
            "index_validation": "header_names_only_not_random_access",
            "phenotype_paragraph_count": len(paragraphs),
            "phenotype_character_count": len(text),
            "explicit_hpo_count": len(hpo),
            "runtime_seconds": round(time.monotonic() - start, 3),
            "environment": {"python": sys.version, "platform": platform.platform()},
            "command": sys.argv,
            "inputs": [
                {"path": str(p), "bytes": p.stat().st_size, "sha256": sha256(p)}
                for p in [vcf, index, doc]
            ],
            "method_sha256": sha256(Path(__file__)),
            "limitations": [
                "No reference-sequence validation or left normalization",
                "No population-frequency/functional annotation",
                "No inheritance inference from a single sample",
                "Explicit HPO extraction only, no phenotype interpretation",
                "No copy-number inference from absence of symbolic alleles",
            ],
        }
    )
    (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    receipt = {
        "status": "PASS",
        "files": [
            {"name": p.name, "bytes": p.stat().st_size, "sha256": sha256(p)}
            for p in sorted(output.iterdir())
        ],
    }
    (output / "receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
    for path in output.iterdir():
        os.chmod(path, 0o600)
    return summary


def smoke_test() -> None:
    import tempfile

    with tempfile.TemporaryDirectory() as tmp:
        base = Path(tmp)
        source = base / "input"
        source.mkdir()
        with gzip.open(source / "synthetic.vcf.gz", "wt") as out:
            out.write("##fileformat=VCFv4.2\n##reference=synthetic_GRCh38\n")
            out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsynthetic\n")
            out.write("1\t10\t.\tA\tG,AT\t40\tPASS\t.\tGT:DP:GQ\t1/2:20:50\n")
            out.write("1\t20\t.\tC\t<DEL>\t.\tLowQual\t.\tGT\t./.\n")
        with gzip.open(source / "synthetic.vcf.gz.tbi", "wb") as out:
            out.write(b"TBI\x01" + struct.pack("<8i", 1, 2, 1, 2, 0, 35, 0, 2) + b"1\x00")
        with zipfile.ZipFile(source / "synthetic.docx", "w") as out:
            out.writestr(
                "word/document.xml",
                '<w:document xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main">'
                "<w:body><w:p><w:r><w:t>HP:0000001 synthetic</w:t>"
                "</w:r></w:p></w:body></w:document>",
            )
        result = intake(source, base / "output")
        assert result["record_count"] == 2 and result["alternate_allele_count"] == 3
        assert result["classes"] == {"SNV": 1, "indel": 1, "symbolic_or_breakend": 1}
        assert result["genotype_states"] == {"heterozygous": 1, "missing_or_partial": 1}
        assert result["explicit_hpo_count"] == 1
        assert all((p.stat().st_mode & 0o777) == 0o600 for p in (base / "output").iterdir())
    assert carriage("0/1", 2) == "uncarried"
    assert carriage("1/.", 1) == "carried_partial"
    assert carriage("./.", 1) == "missing_or_partial"
    print("PASS: synthetic multiallelic, indel, symbolic, missingness, HPO and permissions")


def annotate(intake_dir: Path, resources: Path, output: Path) -> dict[str, Any]:
    """Exact HPO matches and gene-span overlaps; neither implies causality."""
    from bisect import bisect_right
    from collections import defaultdict

    output.mkdir(mode=0o700, parents=True, exist_ok=False)
    terms = set(json.loads((intake_dir / "hpo-explicit-private.json").read_text()))
    matches: dict[str, set[str]] = defaultdict(set)
    diseases: dict[str, set[str]] = defaultdict(set)
    ids: dict[str, set[str]] = defaultdict(set)
    resource = resources / "genes_to_phenotype.txt"
    with resource.open() as stream:
        for row in csv.DictReader(stream, delimiter="\t"):
            if row["hpo_id"] in terms:
                symbol = row["gene_symbol"]
                matches[symbol].add(row["hpo_id"])
                diseases[symbol].add(row["disease_id"])
                ids[symbol].add(row["ncbi_gene_id"])
    ordered = sorted(matches, key=lambda gene: (-len(matches[gene]), gene))
    with (output / "hpo_exact_gene_priority.tsv").open("w") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "rank",
                "gene_symbol",
                "ncbi_gene_ids",
                "exact_hpo_count",
                "hpo_ids",
                "disease_ids",
                "interpretation",
            ]
        )
        for rank, gene in enumerate(ordered, 1):
            writer.writerow(
                [
                    rank,
                    gene,
                    ";".join(sorted(ids[gene])),
                    len(matches[gene]),
                    ";".join(sorted(matches[gene])),
                    ";".join(sorted(diseases[gene])),
                    "exact_match_only_presence_and_negation_not_reviewed",
                ]
            )
    intervals: dict[str, list[tuple[int, int, str, str]]] = defaultdict(list)
    features = defaultdict(list)
    transcript_exons = defaultdict(list)
    gtf = resources / "gencode.v50.primary_assembly.annotation.gtf.gz"
    with gzip.open(gtf, "rt") as stream:
        for line in stream:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if cols[2] not in {"gene", "exon", "CDS"}:
                continue
            attrs = dict(re.findall(r'(\w+) "([^"]+)"', cols[8]))
            gene = attrs.get("gene_name", "")
            if gene in matches:
                chrom = cols[0].removeprefix("chr")
                chrom = "MT" if chrom == "M" else chrom
                left, right = int(cols[3]), int(cols[4])
                if cols[2] == "gene":
                    intervals[chrom].append((left, right, gene, attrs["gene_id"]))
                else:
                    transcript = attrs["transcript_id"]
                    features[(chrom, gene)].append((left, right, cols[2], transcript, cols[6]))
                    if cols[2] == "exon":
                        transcript_exons[(chrom, gene, transcript, cols[6])].append((left, right))
    for (chrom, gene, transcript, strand), exons in transcript_exons.items():
        ordered_exons = sorted(set(exons))
        for index, (left, right) in enumerate(ordered_exons):
            if index > 0:
                label = "splice_acceptor_2bp" if strand == "+" else "splice_donor_2bp"
                features[(chrom, gene)].append((left - 2, left - 1, label, transcript, strand))
            if index < len(ordered_exons) - 1:
                label = "splice_donor_2bp" if strand == "+" else "splice_acceptor_2bp"
                features[(chrom, gene)].append((right + 1, right + 2, label, transcript, strand))
    feature_starts, feature_ends = {}, {}
    for key, feature_spans in features.items():
        feature_spans.sort()
        feature_starts[key] = [feature_span[0] for feature_span in feature_spans]
        peak = 0
        ends = []
        for feature_span in feature_spans:
            peak = max(peak, feature_span[1])
            ends.append(peak)
        feature_ends[key] = ends
    starts = {}
    prefix_ends = {}
    for chrom, spans in intervals.items():
        spans.sort()
        starts[chrom] = [span[0] for span in spans]
        peak = 0
        ends = []
        for span in spans:
            peak = max(peak, span[1])
            ends.append(peak)
        prefix_ends[chrom] = ends
    location_flags: Counter[str] = Counter()
    coding_classes: Counter[str] = Counter()
    overlap_count = 0
    overlap_classes: Counter[str] = Counter()
    genes_with_variants = set()
    with (
        (intake_dir / "variant_evidence_groundwork.tsv").open() as stream,
        (output / "phenotype_gene_span_variant_overlaps.tsv").open("w") as out,
    ):
        reader = csv.DictReader(stream, delimiter="\t")
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        columns = list(reader.fieldnames or [])
        writer.writerow(
            columns
            + [
                "gene_symbol",
                "ensembl_gene_id",
                "exact_hpo_count",
                "overlap_evidence",
                "location_flags",
                "feature_transcript_strand",
            ]
        )
        for row in reader:
            chrom = row["chrom"].removeprefix("chr")
            chrom = "MT" if chrom == "M" else chrom
            if chrom not in intervals:
                continue
            begin = int(row["pos"])
            end = begin + len(row["ref"]) - 1
            spans = intervals[chrom]
            idx = bisect_right(starts[chrom], end) - 1
            while idx >= 0 and prefix_ends[chrom][idx] >= begin:
                left, right, gene, gene_id = spans[idx]
                if right >= begin and left <= end:
                    flags, details = set(), set()
                    key = (chrom, gene)
                    if key in features:
                        feature_idx = bisect_right(feature_starts[key], end) - 1
                        while feature_idx >= 0 and feature_ends[key][feature_idx] >= begin:
                            fleft, fright, flag, transcript, strand = features[key][feature_idx]
                            if fright >= begin and fleft <= end:
                                flags.add(flag)
                                details.add(f"{flag}|{transcript}|{strand}")
                            feature_idx -= 1
                    for flag in flags:
                        location_flags[flag] += 1
                    if "CDS" in flags:
                        coding_classes[row["class"]] += 1
                    writer.writerow(
                        [row[col] for col in columns]
                        + [
                            gene,
                            gene_id,
                            len(matches[gene]),
                            "reference_span_overlap_not_functional_consequence",
                            ";".join(sorted(flags)),
                            ";".join(sorted(details)),
                        ]
                    )
                    overlap_count += 1
                    overlap_classes[row["class"]] += 1
                    genes_with_variants.add(gene)
                idx -= 1
    matched_terms = set().union(*matches.values()) if matches else set()
    genes_in_gtf = {s[2] for spans in intervals.values() for s in spans}
    result = {
        "method": "exact_HPO_then_gene_reference_span_overlap_v1",
        "hpo_term_count": len(terms),
        "matched_hpo_term_count": len(matched_terms),
        "matched_gene_count": len(matches),
        "matched_genes_in_gtf": len(genes_in_gtf),
        "matched_genes_missing_gtf": len(set(matches) - genes_in_gtf),
        "genes_with_overlaps": len(genes_with_variants),
        "allele_gene_overlap_rows": overlap_count,
        "overlap_classes": dict(overlap_classes),
        "location_flag_rows": dict(location_flags),
        "CDS_overlap_classes": dict(coding_classes),
        "gene_exact_match_count_histogram": dict(Counter(map(len, matches.values()))),
        "source_sha256": sha256(Path(__file__)),
        "command": sys.argv,
        "inputs": [
            {"path": str(p), "sha256": sha256(p)}
            for p in [
                resource,
                gtf,
                intake_dir / "hpo-explicit-private.json",
                intake_dir / "variant_evidence_groundwork.tsv",
            ]
        ],
        "limitations": [
            "No HPO semantic similarity or ontology propagation",
            "Explicit HPO terms not reviewed for presence/negation",
            "Gene-name exact join; unmatched aliases retained as missing",
            "Reference-span overlap is not a molecular consequence",
            "All quality filters retained; no population rarity assessment",
            "SV overlap uses REF span only, not END/breakend resolution",
            "No causal ranking or clinical conclusion",
        ],
    }
    (output / "summary.json").write_text(json.dumps(result, indent=2) + "\n")
    receipt = {
        "status": "PASS",
        "files": [
            {"name": p.name, "bytes": p.stat().st_size, "sha256": sha256(p)}
            for p in sorted(output.iterdir())
        ],
    }
    (output / "receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
    for p in output.iterdir():
        p.chmod(0o600)
    return result


def main() -> None:
    os.umask(0o077)
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--smoke-test", action="store_true")
    parser.add_argument("--intake-dir", type=Path)
    parser.add_argument("--resources", type=Path)
    args = parser.parse_args()
    if args.smoke_test:
        smoke_test()
        return
    if args.intake_dir and args.resources and args.output:
        result = annotate(args.intake_dir, args.resources, args.output)
        print(
            json.dumps(
                {k: v for k, v in result.items() if k not in {"inputs", "command", "limitations"}}
            )
        )
        return
    if not args.source or not args.output:
        parser.error("--source and --output are required")
    result = intake(args.source, args.output)
    print(
        json.dumps(
            {
                k: result[k]
                for k in [
                    "record_count",
                    "alternate_allele_count",
                    "classes",
                    "sample_count",
                    "reference_build",
                    "explicit_hpo_count",
                    "candidate_ranking_status",
                    "runtime_seconds",
                ]
            }
        )
    )


if __name__ == "__main__":
    main()
