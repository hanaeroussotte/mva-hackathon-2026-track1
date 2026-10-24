# REPO_SCOPE.md — what goes in the GitHub repository, and what is held back

**Status: the repository is NOT public-ready, and this document does not claim
it is.** One question must be answered by the organisers first (§3).

## 1. What the two sets of rules each require

| source | requirement |
|---|---|
| submission tab (`tabs/submit_track1.py`, retrieved 2026-10-24) | repository may stay private during the hackathon, **must be public once it ends** so judges can review; must contain documented reproducible code and configuration; submission files and methods optional but recommended |
| gated data-use agreement (`SageBio/mva-hackathon-2026-data`) | "I will not release or otherwise grant data access to anyone"; delete all data after the hackathon; notify `cnb@sagebase.org`; submissions are published CC-BY 4.0 |

These are compatible for code and incompatible for nothing — **provided the
repository contains code rather than the child's data.** The rules ask for
reproducible code, not for the genome.

## 2. The scoped package: `mva_track1_repo_scoped.zip`

Contains **only** material whose shareability follows directly from the rules:

| included | why it is permissible |
|---|---|
| `src/*.py` (all 11 scripts) | explicitly required: "all scripts and configuration files needed to reproduce your results" |
| `environment.yml` | configuration, explicitly required |
| `work/evaluation.py` | the organisers' own evaluator, already public on the challenge Space |
| `work/bench_cases.json` | benchmark case definitions — planted ClinVar variants, no proband genotypes |
| `hanaeroussotte_track1_report.md` | the submission report; methods are optional-but-recommended, and the submitted variants it names are themselves CC-BY 4.0 per the agreement |
| `submissions/hanaeroussotte_mechanism-gated_FINAL.csv` | the submission, explicitly CC-BY 4.0 |
| `docs/SCORING_SPECIFICATION.md`, `docs/BENCHMARK_EXPLAINED.md`, `docs/SUBMISSION_CHECKLIST.md`, `docs/CONTRIBUTIONS.md`, `docs/CHANGELOG.md`, this file | methods and process documentation containing no individual-level genotype data beyond the submitted variants |
| benchmark result tables: `variant_recovery_summary.tsv`, `variant_recovery_cases.tsv`, `failed_case_trace.tsv`, `baseline_comparison.tsv`, `competing_allele_availability.tsv`, `missense_policy_on_new_benchmark.tsv`, `pair_aggregation_sensitivity.tsv` | synthetic cases; no proband genotypes |
| `README.md`, `docs/CHECKSUMS_REPO.txt` | navigation and integrity |

## 3. Held back pending clarification — the exact question to ask

**Question for the organisers:**
*"The submission rules require a public repository with reproducible code, while
the data-use agreement prohibits releasing the Data. Does 'the Data' extend to
genotype-scale derived outputs — for example a full ranked candidate table
listing hundreds of the proband's variant positions with allelic depths — or is
the prohibition limited to the supplied sequence files and their direct copies?
If derived tables must be excluded, may the public repository contain code and
methods only?"*

Until that is answered, these are **excluded** from the scoped package:

| file | what it contains |
|---|---|
| `work/cds_vep.tsv`, `work/cds_genotypes.tsv`, `work/cds_phase.tsv` | 29,907 proband coding calls with genotypes — closest to raw data; excluded under any reading |
| `work/candidates_pre_af.tsv` | 861 rows of proband variants with AD/DP/GQ |
| `docs/final_ranked_hypotheses.tsv` | 448 hypotheses naming proband variant positions |
| `docs/af_lookup_status.tsv` | 841 proband variants |
| `docs/transcript_selection_suspects.tsv` | 275 proband variant positions |
| `docs/candidate_count_reconciliation.tsv` | 21 proband variants with genotypes |
| `docs/bub1b_allele_evidence.tsv`, `docs/allele_evidence_table.tsv` | per-allele detail for the submitted variants with read depths |
| `docs/kmer_counts.tsv`, `docs/phase_states.tsv`, `docs/upstream_scan.tsv`, `docs/locus_spliceai.tsv` | proband read-level and locus data |
| `docs/phenotype_match_details.tsv`, `docs/phenotype_feature_map.tsv`, `config/proband_hpo.txt` | the eight HPO terms, term by term |
| `docs/METHODS_FINAL.md`, `docs/CASE_EXPLAINED.md`, `docs/BUB1B_EVIDENCE.md`, `docs/WRITING_HANDOVER.md`, `docs/QUESTIONS_I_SHOULD_BE_ABLE_TO_ANSWER.md`, `docs/historical/*`, `STATUS.md` | prose carrying individual-level genotype and phenotype detail |
| `bench_inputs.tar.gz`, `mva_track1_package.zip` | the complete private handover archives |

**Note the consequence honestly:** with these excluded, the repository supports
the **ranking rerun** only if `work/candidates_pre_af.tsv` is restored. As
shipped, the scoped package lets a reviewer read and audit every line of code and
every benchmark result, but not re-execute the proband ranking. That trade is
deliberate and reversible the moment the question in §3 is answered.

## 4. Repository link

`https://github.com/hanaeroussotte/mva-hackathon-2026-track1`

**STATUS: UNVERIFIED — not checked in this pass.** `github.com` and `api.github.com` were
unreachable from this environment: the local proxy accepted the connection and
then timed out after 25–28 s on three consecutive attempts, while other
allowlisted domains returned HTTP 200. The repository was created and pushed
earlier in this project and was recorded then as **private**; that is a
previously recorded state, not a check executed now. **Confirm the link resolves
and decide its visibility yourself before entering it on the submission form.**

No repository was published, made public, or modified in this pass.
