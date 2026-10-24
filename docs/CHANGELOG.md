# CHANGELOG — final correction pass, 2026-10-23

Two passes are recorded. Nothing in the prediction changed in either.
**No prediction row was altered: variants, ordering, EPCR and finding types are
identical to the earlier submission, and both files score identically on all
five evaluator scenarios.**

## Pass 2 (this pass) — corrections

| # | change | status before | status now |
|---|---|---|---|
| 1 | F5/COL17A1/RNF31 annotation mismatch | attributed to the **benchmark's** truth selection | **retracted.** Demonstrated to be a defect in `pick_tc()` (`02_vep_annotate.py`): it takes the first MANE-tagged transcript, so a neighbouring gene's MANE transcript displaces a coding consequence in the intended gene. VEP's own `most_severe_consequence` is the coding term in all four. Not fixed. |
| 2 | affected-case count | "four alleles, three absent genes" | **four alleles across three cases**: F5 = one incomplete pair; COL17A1 + RNF31 = **two** absent genes. All three kept in the denominators; no sensitivity analysis offered. |
| 3 | scope of the annotation defect | not measured | **275 of 29,907 proband coding calls (0.92%)** are candidate instances; the three submitted alleles are unaffected. |
| 4 | 583 vs 562 candidates | attributed to the missing `lowconf_regions.bed` | **retracted** — a missing −1.5 penalty cannot reduce a written count, and both figures come from the same run. Demonstrated: **15** hets beyond the top-4 `comphet` cap + **6** unhandled `GT 1/2` genotypes. |
| 5 | "variant evidence earns its place in the shortlist" | asserted from the 0/45 evidence-only arm | **removed.** Every arm shares that shortlist and none varies it. |
| 6 | phase | "unresolvable from short reads at any depth" | "**unresolved in the analysed reads and marker set**"; the limit is fragment length against gap size on this library. |
| 7 | *cis* | "the diagnosis is wrong entirely" | "would **undermine this proposed two-allele explanation**" — not exclusion of the gene. |
| 8 | p.Asn1002Lys transcript fate | "49 nt downstream → escapes NMD" | **no premature stop is introduced**; it is a single amino-acid substitution, so no decay rule applies. Its last-exon position is simply where it is. |
| 9 | Suijkerbuijk 2010 substitutions | five described collectively as destabilising | **corrected from the primary text.** Reduced abundance: I909T, R727C (~2× turnover), L1012P. **Y155C and R550Q did *not* reduce protein levels**; LAP-R550Q rescued both functions. R814H appears only in the patient-line panel. Residue 1002 untested. |
| 10 | regulatory variation | "not assessable" alongside a reported 60 kb scan | **reconciled by class.** Searched: small variants in the 60 kb upstream window (85 PASS, 52 het, no CADD filter) and locus SpliceAI. Inadequately assessed: structural variants, mobile-element insertions, single-exon CNV. |
| 11 | allele 2 frequency | "`absent` is materially correct" | "the measured frequency **receives the same rarity bonus** as a verified absence" — unknown, absent and observed ultra-rare are different evidence states. |
| 12 | neutral-frequency result | 36/45 and 37/45 both appeared | **36/45** everywhere, with the one changed case named (**COQ7**, rank 1 → 3). Primary and sensitivity arms separated. |
| 13 | frequency conflation | implied repaired | **disclosed as unrepaired**: the primary configuration still treats a failed lookup as an absence. |
| 14 | "most prioritisers rank genes" | asserted | **removed** — no comparison against other tools was run. |
| 15 | submission identity | STATUS carried the old checksum against the final filename; PROVENANCE said the original stands | both corrected; FINAL and HISTORICAL designated explicitly, checksums regenerated from the files. |
| 16 | reproduction commands | `04_shortlist.py` shown with no arguments | corrected to the 5 and 6 required positional arguments, **executed in a clean directory**, log at `docs/REPRO_LOG.txt`; `environment.yml` added. |
| 17 | data-sharing statement | "patient data is not in this repository" | replaced: raw sequence excluded, **but several reports and derived files contain individual-level genotype or phenotype**; three-way categorisation in `PROVENANCE.json`; official terms must be confirmed before publishing. |
| 18 | references | none | 20-entry reference list; those not read in full are marked. |

## Pass 1 — corrections

| # | change |
|---|---|
| 19 | "No *cis* pair ever proposed / 13 absent" — **false**. 11 of 15 proposed at ranks 16–24; 2 absent. Only "none in the top ten" survives. |
| 20 | "Zero wrong partners" as evidence about partner selection — **withdrawn**; 80 of 80 injected genes had zero competing alleles. |
| 21 | The current 80 benchmark cases are a **new benchmark, not a reproduction**; historical 68/80, 70/80 and 47/80 superseded for variant-level claims. |
| 22 | Gene recovery (39/45) and exact-configuration recovery (37/45) reported separately throughout. |
| 23 | Exact recovery **does** improve with depth in `single_allele` (17 → 19 at rank 3), not in complete-truth strata. |
| 24 | `mean` vs `min`: same endpoints, **not identical rankings** — 25 of 80 scores differ, 2 ranks change. |
| 25 | Ranking granularity is **one hypothesis per gene**, not per gene/model. |
| 26 | *LZTR1*: ranking hypothesis (`incomplete_recessive` → Noonan) separated from reporting rationale (dominant schwannomatosis). The old note had these backwards. |
| 27 | "7 of 8 features explained" removed; EPCR labelled an ordering value. |

## Defects found and deliberately left unfixed

Re-annotating or re-shortlisting would change the candidate set after the
results were seen. Both are documented instead:

1. `pick_tc()` first-MANE transcript selection (items 1–3).
2. Unhandled `GT 1/2` in `04_shortlist.py` (item 4).
