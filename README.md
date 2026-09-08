# Local verification-first variant review tools

Prototype research code. Python 3.12.11 standard library; no runtime third-party dependency, network request, credential, patient genome or final patient prediction is bundled here. Keep inputs and generated evidence in a private external directory. No trained model or random seed is involved. Interpretation and EPCR are manual research judgments; code reproduces annotation and a frozen reviewed hypothesis, not biological causality.

`science_campaign/` contains the seven analysis/validation modules. `select_targets.py` makes the carried coding/two-base-splice selection explicit. `conditional_indels.py` is the reviewed, overwrite-guarded indel screen adapted only to accept input/output arguments instead of personal absolute paths. Historical evidence remains attributed to its original frozen source. `reproduce_submission.py` accepts an externally supplied frozen two-SNV hypothesis, verifies its VCF sample/carriage and plus-strand transcript codons, and writes a CSV and hashed receipt. Its intentionally narrow scope is stated in its CLI and errors.

From this code directory, a local pipeline can be run with separately authorized organizer files in SOURCE and whole public resources in RESOURCES. Use fresh output paths. Public reference URLs, filenames, releases and SHA256 values are in `public-references.json`. RESOURCES must contain `genes_to_phenotype.txt`, `phenotype.hpoa` and `gencode.v50.primary_assembly.annotation.gtf.gz`. HPO URLs are rolling resources: a future download may differ, so record and verify its hash rather than silently treating it as the frozen snapshot. Organizer-controlled input identifiers and hashes remain in private receipts. No patient-specific query or download is performed by this code.

```bash
python -m science_campaign.mva_intake --source "$SOURCE" --output "$OUT/intake"
python -m science_campaign.mva_intake --intake-dir "$OUT/intake" --resources "$RESOURCES" --output "$OUT/annotation"
python select_targets.py --input "$OUT/annotation/phenotype_gene_span_variant_overlaps.tsv" --output "$OUT/annotation/coding_splice_carried_targets.tsv"
python -m science_campaign.mva_clinvar --case "$VCF" --reference "$CLINVAR" --output-dir "$OUT/clinvar"
python -m science_campaign.mva_population --target "$OUT/annotation/coding_splice_carried_targets.tsv" --reference "$POPULATION" --output "$OUT/population"
python -m science_campaign.mva_evidence --target "$OUT/annotation/coding_splice_carried_targets.tsv" --clinvar "$OUT/clinvar/clinvar-exact-matches.csv" --population "$OUT/population/coding_splice_population.csv" --output "$OUT/evidence"
python -m science_campaign.mva_condition --clinvar "$OUT/clinvar/clinvar-exact-matches.csv" --hpo-literals "$OUT/intake/hpo-explicit-private.json" --hpoa "$RESOURCES/phenotype.hpoa" --output-dir "$OUT/conditions"
python -m science_campaign.mva_inheritance --clinvar "$OUT/clinvar/clinvar-exact-matches.csv" --conditions "$OUT/conditions/condition-evidence.csv" --targets "$OUT/annotation/coding_splice_carried_targets.tsv" --gene-panel "$OUT/annotation/hpo_exact_gene_priority.tsv" --output "$OUT/inheritance"
python conditional_indels.py --targets "$OUT/annotation/coding_splice_carried_targets.tsv" --gtf "$RESOURCES/gencode.v50.primary_assembly.annotation.gtf.gz" --evidence "$OUT/evidence/evidence-ledger.jsonl" --output "$OUT/indels"
```

Environment variables represent local paths: SOURCE is the three-file dataset directory, RESOURCES contains the named HPO/GTF files, OUT is a fresh private working directory, VCF is the organizer VCF, CLINVAR is the complete GRCh38 release VCF and POPULATION the complete Ensembl115 population VCF. No downloading is performed by this sequence. Post-feedback verification on 8 September 2026 executed all nine commands above in order on the original published code and pinned local inputs under Python 3.12.11/macOS ARM64. Every command exited 0; combined wall time was 187.350 seconds. Selected downstream table row multisets matched the historical outputs, preserving duplicate counts; intake and annotation scientific summaries also matched. Runtime, command/path metadata and the previously documented indel scope label are separately accounted for. This demonstrates reproducibility on those inputs, not independent rediscovery or causal validation. See `METHODS-POST-FEEDBACK.md`; controlled inputs, full outputs and private receipts are not in this repository.

Final replay accepts `--hypothesis`, `--vcf`, `--fasta`, `--gtf`, `--output`. The hypothesis JSON specifies proband_id, primary_pair, transcript, expected_cds_positions, expected_reference_codons, expected_alternate_codons, epcr and notes. Optional vcf_sample_id records a source sample name when the submission identifier differs; absent it, proband_id is also the expected VCF sample. Supply only a locally reviewed hypothesis. The replay outputs do not establish segregation, variant pathogenicity or diagnosis. The CSV validator checks format and values, not causal truth.

Limits: exact allele matching without comprehensive left normalization, location/length consequence screens rather than full clinical annotation, literal HPO target selection before manual contextual review, incomplete frequency and outside-panel coverage, no BAM/trio/CNV/SV/noncoding functional analysis, and no clinical probability calibration. Tests used targeted synthetic and actual replay checks rather than a full production regression suite. Analysis files and outputs must not be interpreted as permission to redistribute controlled data.

## Synthetic replay example

This entirely invented 12-base sequence is a software demonstration. It is not a real GRCh38 region or a patient prediction. Run from this repository using a fresh local output directory:

```bash
python examples/run_synthetic_example.py --output /tmp/mva-synthetic-example
```

Observed output on Python 3.12.11:

```text
PASS: mapped identifier, same-identifier compatibility, wrong-source rejection
```

The runner creates tiny local compressed FASTA/VCF/GTF inputs, uses `examples/hypothesis.json`, executes the actual replay CLI, verifies different source/submission identifiers, confirms backward-compatible equal identifiers, and checks that a wrong source sample is rejected before output creation. Expected rejection exits 1 and is recorded as a passing negative control. It requires no network or patient data. The expected codons and metadata are declared rather than inferred by an LLM.

For a separately reviewed local two-SNV hypothesis, the complete replay form is:

```bash
python reproduce_submission.py --hypothesis "$HYPOTHESIS" --vcf "$VCF" --fasta "$CHROMOSOME_FASTA" --gtf "$GTF" --output "$REPLAY_OUT"
python -m science_campaign.mva_submission "$REPLAY_OUT/reproduced.csv"
```

All variables are local paths. The first command verifies sample carriage/reference/codons and reproduces the frozen row; the second checks CSV syntax. Neither command discovers a cause or confirms trans phase.
