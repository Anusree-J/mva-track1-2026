# Track 1 methods supplement — 8 September 2026

Team: **anusreejj**. Approach: **Model 1, local annotation followed by explicit manual curation**. This supplement was prepared after the first scored entry. It adds verification and documentation; it does not change the submitted prediction, ordering, EPCR or biological interpretation. The pre-feedback report and subsequent identifier-correction report remain preserved privately. No additional scored submission was used for this methods work.

## Method abstract

We use a local, standard-library Python workflow to turn a single-sample variant callset and an organizer-supplied phenotype document into an auditable review queue. Intake retains ALT-specific carriage, genotype quality, variant class and available phase fields. Whole public HPO resources define a literal-phenotype gene set, and GENCODE supplies strand-aware coding/exon/splice-location annotations. Exact-allele ClinVar joins retain conflicting assertions and condition labels. Population annotations use the matching alternate allele, preserving unknown frequency rather than substituting callset AF or treating missingness as rarity.

Condition-specific review, inheritance screening and conservative indel geometry expose cases where generic gene overlap or an apparent frameshift is insufficient. Downstream manual review interprets phenotype context, disease mechanism, transcript relevance, phase and alternative explanations. The final prediction and subjective EPCR are curated judgments; deterministic replay regenerates the frozen row only after validating its source sample, carried alleles, reference bases and transcript codons. Source sample identity is explicitly separated from the organizer's submission identifier.

Strengths are traceability, preservation of contradictory/missing evidence, and targeted checks that prevent allele, transcript and identifier conflation. A post-feedback execution of the exact published nine-command sequence completed in 187.350 seconds with matching scientific outputs on the pinned local inputs. The synthetic example independently exercises identifier mapping and wrong-source rejection without patient data.

Limitations include exact matching without comprehensive normalization, phenotype-targeted interpretation, geometric rather than comprehensive consequence annotation, incomplete population coverage and manual causal ranking. Segregation, functional damage, structural/noncoding mechanisms and clinical utility were not comprehensively established. A single-case competition match does not measure cross-case sensitivity, calibrated probability or clinical validity. This is an auditable research workflow, not a newly trained predictive model or an asserted novel diagnostic algorithm.

## Automated output and manual review

The annotation tables and replayed CSV are generated automatically. **Selection of the final hypothesis and EPCR underwent downstream manual review and curation**, aided by Codex and a fresh independent reasoning review. Manual work checked narrative positivity/chronology, correlated features, condition-specific inheritance, exact allele assertions, transcript relevance and alternative explanations. A reviewer challenged the initial confidence estimate, and the frozen prediction used the lower agreed estimate. Scoring feedback did not drive a new biological selection.

The workflow can represent proposed pairs, but unphased heterozygous calls do not establish compound heterozygosity. Secondary/incidental findings were considered during review and omitted from the compact primary-focused prediction. Full details remain in the report submitted to the organizer. All 211 targeted indels underwent computational screening and class-level inspection; selected plausible alternatives received deeper transcript/mechanism review. Equal clinical depth for every indel and exhaustive genome interpretation are not claimed.

## Public and proprietary data

The competition supplied controlled-access genomic and phenotype data. Auxiliary biological evidence came from public HPO/GENCODE, ClinVar, population references, GRCh38 sequence and primary literature. **No proprietary auxiliary biological dataset was used.** Access controls on organizer data are retained; this statement does not relabel those data as unrestricted public material. `public-references.json` identifies the actual public files with URLs, versions and hashes. The official source is the [SageBio competition dataset](https://huggingface.co/datasets/SageBio/mva-hackathon-2026-data), introduced on the [Synapse project syn76251147](https://www.synapse.org/Synapse:syn76251147/wiki/642892). No DOI or prescribed formal citation was present in the public Synapse landing-page text inspected on 8 September; those links identify the source without inventing a citation.

## Reproducibility, runtime and cost

The exact published nine-command sequence was executed after scoring using unchanged source under **Python 3.12.11, macOS 26.6.2 ARM64** and eight existing input snapshots. Combined measured command wall time: **187.350 seconds**. All commands exited 0. No new reference transfer, paid compute job or model training was part of this replay. Downloads, earlier analysis and manual reasoning are excluded from this timing. Peak memory for this assembled run was not measured, so no new memory benchmark is claimed.

| Stage | Seconds |
|---|---:|
| Intake | 22.241 |
| HPO/GENCODE annotation | 33.924 |
| Carried target selection | 0.561 |
| ClinVar exact joins | 59.505 |
| Population overlay | 56.053 |
| Evidence consolidation | 0.600 |
| Condition comparison | 0.562 |
| Inheritance screening | 0.163 |
| Conditional indel screen | 13.741 |

Scientific row multisets matched earlier outputs for 3,979 target rows, 47,987 ClinVar rows, 3,978 population rows, 4,107 evidence rows, 1,404 condition rows, seven inheritance records, six counterpart pairs and 211 indel records. Duplicate counts were preserved. The previously documented indel scope-label addition was checked separately; every original indel evidence field matched. This is local reproducibility evidence, not a second independent case or blinded rediscovery.

The replay incurred **no separately billed cloud/API execution charge**; local hardware, electricity, prior public-resource transfers and the existing Pro subscription are not treated as free. Total project monetary cost and total manual/AI effort were not reliably metered, so a complete dollar or labor estimate is unavailable. No new GPU or raw-read processing is required for this demonstrated workflow. Broader cohort scalability remains untested.

## AI disclosure and post-feedback boundary

Provider: OpenAI, Codex desktop under ChatGPT **Pro**. The user explicitly attested that active Codex training sharing, including full-environment sharing if shown, was OFF before minimal case-content reasoning. This is user attestation, not an independently inspected Codex settings panel. The fresh scientific reviewer used GPT-6-Astra with high reasoning effort. No patient-specific AlphaGenome query or unrelated-provider case dispatch contributed to the prediction.

The original source-sample/submission-identifier error was rejected before scoring, then corrected transparently after explicit approval. The accepted entry received the official automated match result; subsequent methods verification did not change its biological fields or EPCR. Public entrant filenames were incidentally visible during a later request to locate our leaderboard row; no entrant answer was opened or used. Earlier prior-exposure and non-blinded-discovery disclosures remain in the frozen organizer report. The five remaining attempts have not been used for this supplement.

## Acknowledgement

This work was made possible through the Hackathon, organized by Sage Bionetworks in partnership with the MVA Society, Hugging Face, and BEACON (The Benchmarking, Evaluation, and Assessment Consortium for Science), with prize sponsorship from AWS and Anthropic. We are deeply grateful to the child and their family who generously contributed their data and their story to advance research into this rare disease. We acknowledge their trust in making this Hackathon possible.
