# Local verification-first variant review tools

Prototype research code. Python 3.12.11 standard library; no runtime third-party dependency, network request, credential, patient genome or final patient prediction is bundled here. Keep inputs and generated evidence in a private external directory. No trained model or random seed is involved. Interpretation and EPCR are manual research judgments; code reproduces annotation and a frozen reviewed hypothesis, not biological causality.

`science_campaign/` contains the seven analysis/validation modules. `select_targets.py` makes the carried coding/two-base-splice selection explicit. `conditional_indels.py` is the reviewed, overwrite-guarded indel screen adapted only to accept input/output arguments instead of personal absolute paths. Historical evidence remains attributed to its original frozen source. `reproduce_submission.py` accepts an externally supplied frozen two-SNV hypothesis, verifies its VCF sample/carriage and plus-strand transcript codons, and writes a CSV and hashed receipt. Its intentionally narrow scope is stated in its CLI and errors.

From this code directory, a local pipeline can be run with separately authorized organizer files in SOURCE and whole public resources in RESOURCES. Use fresh output paths. The same reference snapshots used for the submitted analysis are identified by SHA256 in the privately supplied report packet; avoid silently replacing them with current releases.

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

Environment variables represent local paths: SOURCE is the three-file dataset directory, RESOURCES contains the named HPO/GTF files, OUT is a fresh private working directory, VCF is the organizer VCF, CLINVAR is the complete GRCh38 release VCF and POPULATION the complete Ensembl115 population VCF. No downloading is performed by this sequence. The earlier component executions are documented in original receipts; this assembled shell sequence has not been claimed as a newly executed whole-pipeline run.

Final replay accepts `--hypothesis`, `--vcf`, `--fasta`, `--gtf`, `--output`. The hypothesis JSON specifies proband_id, primary_pair, transcript, expected_cds_positions, expected_reference_codons, expected_alternate_codons, epcr and notes. Supply only a locally reviewed hypothesis. The replay outputs do not establish segregation, variant pathogenicity or diagnosis. The CSV validator checks format and values, not causal truth.

Limits: exact allele matching without comprehensive left normalization, location/length consequence screens rather than full clinical annotation, literal HPO target selection before manual contextual review, incomplete frequency and outside-panel coverage, no BAM/trio/CNV/SV/noncoding functional analysis, and no clinical probability calibration. Tests used targeted synthetic and actual replay checks rather than a full production regression suite. Analysis files and outputs must not be interpreted as permission to redistribute controlled data.
