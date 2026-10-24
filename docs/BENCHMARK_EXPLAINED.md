# BENCHMARK_EXPLAINED.md

What the benchmark is, exactly how it was built, what each number means, and
where it is weak. Read this before quoting any benchmark figure.

**There are two benchmarks in this project's history. They are not comparable.**

| | historical | **current** |
|---|---|---|
| when | earlier in the project | 2026-10-23 |
| ClinVar release | 2026-10-04 | **2026-10-19** |
| injected alleles recorded? | **no** — gene ranks only | **yes** |
| gene overlap with the other | — | **4 of 80** |
| headline quoted | 68/80, later 70/80 vs 47/80 | see below |

Because the old runner never saved the injected alleles, case identity cannot
be verified and **the current set is a new benchmark, not a reproduction.** No
paired comparison with the historical figures is valid, and the historical
figures are superseded for all variant-level claims. The one advantage: these
cases were never used to tune anything.

---

## 1. How the cases were built

Generator: `src/09_spikein_benchmark.py` plus the case builder recorded in
`handoff/bench_cases.json`. Deterministic, **seed 20261007**.

**Source of truth labels.** Each case plants one or two variants of known
identity into this proband's own callset; the truth label is exactly those
planted variants. Nothing is inferred from an external answer key.

**Eligibility for a truth variant** — all four conditions:

1. ClinVar **Pathogenic** or **Likely pathogenic**;
2. a coding molecular consequence (missense, nonsense, frameshift, start/stop
   lost, splice donor/acceptor);
3. population frequency **< 1e-3** in every one of AF_ESP, AF_EXAC, AF_TGP;
4. single-nucleotide, with a gene assignment.

192,241 such variants across 5,193 genes
(`work/clinvar_cds_plp_af.tsv`).

**Eligibility for a gene** — it must have a recessive disease (HP:0000007)
with **≥ 8 annotated phenotype terms**, and must *not* be (a) an MVA gene,
(b) a gene appearing in the real submission, or (c) **a gene already in the
proband's own candidate set**. Pools: 2,171 genes with ≥2 eligible variants,
2,327 with ≥1.

**Exclusion (c) is the benchmark's most consequential design decision** and is
revisited in §5.

**Strata**, 80 cases total:

| stratum | n | construction | genotype |
|---|---:|---|---|
| `comphet_trans` | 30 | two eligible variants in the gene, **≥1,000 bp apart**, no phase tags | both 0/1 |
| `comphet_cis` | 15 | two eligible variants **≤60 bp apart**, given a shared `PID` and `0\|1` `PGT` so GATK reads them as *cis* | both 0/1 |
| `single_allele` | 20 | one eligible variant | 0/1 |
| `hom` | 15 | one eligible variant | 1/1 |

**Phenotypes are annotation-derived.** Each case's eight HPO terms are
**sampled from the phenotype annotation of the disease being planted**. This
is the benchmark's second major limitation (§5).

**Genotype quality is resampled from the real data.** AD/DP/GQ for each planted
variant is drawn from the proband's own 18,741 heterozygous or 10,960
homozygous coding calls, so quality metrics are realistic rather than ideal.

**What was masked.** All arms reported here run `ARMS=masked`: the planted
variants are **removed from the ClinVar evidence file**, so the pipeline cannot
see the Pathogenic assertion that made them eligible. Without masking the task
would be trivial.

**What is shared across cases.** The entire background: all 80 cases use this
same proband's 29,907 coding calls. Cases are independent at gene and disease
level (80 distinct genes, 80 distinct diseases) but **not** at background
level.

---

## 2. Endpoint definitions and denominators

These were conflated in earlier reporting. They are distinct:

| endpoint | definition | denominator |
|---|---|---|
| **correct gene @ k** | the gene appears in the top *k* ranked hypotheses | cases in the stratum |
| **exact configuration @ k** | a top-*k* hypothesis names *exactly* the truth allele set | cases in the stratum |
| **one correct allele** | the gene's own hypothesis is a strict subset of the truth | cases in the stratum |
| **correct allele + wrong partner** | a hypothesis contains ≥1 truth allele **and** ≥1 non-truth allele | cases in the stratum |
| **absent from the full ranking** | the gene appears nowhere in the ranked output | cases in the stratum |
| **outside the top ten** | the gene is ranked, but at rank > 10 | **not the same as absent** |

The last two were previously conflated — an error that produced a false claim
about the cis controls (§4).

**These four strata measure four different things and must not be summed:**

* `comphet_trans` and `hom` have a **complete causal truth** → exact-recovery
  is meaningful, and the official evaluator can score them.
* `single_allele` measures retrieval of an **incomplete hypothesis** — a lead
  for follow-up, not a solved diagnosis.
* `comphet_cis` has **no valid causal answer**; its endpoint is *not* promoting
  a non-causal pair. Feeding it to an evaluator that requires a true answer
  would be meaningless.

---

## 3. Current results

Source: `docs/variant_recovery_summary.tsv`, `docs/variant_recovery_cases.tsv`,
`docs/variant_recovery_evaluator.tsv`.

| stratum | n | correct gene @1 | **exact config @1** | @10 | one allele only | absent from full ranking | mean rank points | mean F-max |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| comphet_trans | 30 | 26 | **24 (80%)** | 24 | 5 | 1 | 85.17 | 0.883 |
| hom | 15 | 13 | **13 (87%)** | 13 | 0 | 2 | 86.67 | 0.867 |
| single_allele | 20 | 17 | 17 | 19 | 0 | 1 | — | — |
| comphet_cis | 15 | 2 | **0** | 0 | 2 | 2 | — | — |

**Complete-truth total: 37 of 45.** Wilson 95% CI 69–91%. Rank points and
F-max come from the challenge's own `evaluation.py`, applied to the pipeline's
hypothesis ranking (**not** the submission row shape, which adds backstop
single-allele rows — see §6).

**Exact recovery improves with rank depth in the `single_allele` stratum
(17 → 19 at rank 3) and not in the complete-truth strata.** The ranking retains
**one hypothesis per gene** — the best-scoring configuration for that gene
(`src/05_rank_and_submit.py`, the `max(options, …)` selection) — so deeper ranks
add *other genes*, never a better configuration for the same gene. Depth helps
only where the correct configuration was already formed but ranked below first,
which is what happens in the two `single_allele` cases.

### Baselines

All three use the **same cases, candidate pool, annotations, inheritance gate
and truth definitions**. They are **internal ablations, not independent
external tools** — the full method's shortlist is upstream of all of them.

| arm | rule | exact @1 (45) | gene @1 (45) | challenge *BUB1B* rank |
|---|---|---:|---:|---:|
| B1 evidence+rarity only | `total = e_var` — phenotype ignored | **0/45** | 0/45 | 3 |
| B2 phenotype only | `total = e_pheno` — variant evidence ignored, inheritance gate retained | **37/45** | **42/45** | 1 |
| FULL | `total = e_var × e_pheno` | **37/45** | 39/45 | 1 |

Paired exact McNemar, FULL vs B2: exact config @1 **37 vs 37, zero discordant
pairs, p = 1.000**; correct gene @1 **39 vs 42, p = 0.250** with all three
discordant cases favouring B2.

**On this benchmark the variant-evidence multiplier adds nothing to exact
recovery.** B1's 0/45 shows only that variant evidence and rarity alone cannot
order these candidates.

Two things this result does **not** establish: it says nothing about the
shortlist's value, because every arm — including both baselines — is built on
the same `04_shortlist.py` output and none of them varies it; and
**non-significance does not prove the methods are equivalent** — 45 cases cannot
exclude a modest difference.

### Missense-evidence policy

| policy | exact @1 (45) | gene @1 (45) |
|---|---:|---:|
| `once` (current default) | 37/45 | 39/45 |
| `both` (pre-correction) | **39/45** | **41/45** |
| `none` | 26/45 | 33/45 |

`both` recovers two cases `once` misses (**HMGCR, PTRHD1**), exactly those the
failure trace predicts (§4). p = 0.50 — not significant, so `both` is **not**
claimed to be better. But the double-counting correction has a **measured cost
with an understood mechanism**, and that is worth stating. `none` being far
worse confirms missense evidence matters somewhere in the pipeline.

### Frequency-handling sensitivity

`rarity_points()` awards `absent_gnomAD` (+1.5) to a verified absence **and** to
a failed lookup. In this run **287 of 841 (34%)** frequency lookups were
unresolved (`docs/af_lookup_status.tsv`).

| unresolved treated as | exact @1 (45) | challenge *BUB1B* rank |
|---|---:|---:|
| absent, +1.5 (current) | 37/45 | 1 |
| unknown, 0.0 | 36/45 | 1 |
| common, −6.0 (worst case) | 17/45 | 2 |

`absent` is the primary (submitted) configuration; the other two are sensitivity
arms. **The neutral result is 36/45, not 37/45.** The single case that changes
is **COQ7** (`hom` stratum): its truth allele `16:19075775 T>A` has an
unresolved lookup, so `var_score` falls 3.0 → 1.5 and the gene moves from rank 1
(0.6988) to rank 3 (0.4512) — still recovered, no longer at rank 1. The `common`
arm is an implausible bound rather than a scenario, since the planted variants
are rare by construction. **The primary configuration still conflates unknown
with absent**; that is disclosed here, not repaired. The baseline comparison is
unchanged under `AF_UNKNOWN=unknown` (B2 37/45, FULL 37/45).

---

## 4. Failure analysis

Source: `docs/failed_case_trace.tsv`. Eight complete-truth cases did not return
the exact configuration at rank 1 — **six** in `comphet_trans`, **two** in `hom`.

| failure mode | n | cases |
|---|---:|---|
| **an allele was lost during filtering** → only one allele reached the candidate set, so `comphet` (needing two rows) was unavailable and `incomplete_recessive` was selected | **5** | HMGCR, POPDC3, PTRHD1, SGSH (ranks 1, 1, 2, 2), F5 (rank 7) |
| **the gene never entered the candidate set** → both alleles filtered out | **3** | SOD1, TPRKB, NOS3 |
| a complete pair was selected but ranked below first | 0 | — |
| a wrong partner was chosen | **0** | — (but see §5 — no wrong partner was ever available) |

**One mechanism explains eight of the nine dropped alleles.** They are missense
variants whose AlphaMissense call is *ambiguous* (+0.4) or *likely_benign*
(−0.5). Under `MISSENSE_EVIDENCE=once` CADD is **discarded for missense**, so a
CADD of 21–28 contributes nothing and the allele falls below the 1.0 evidence
floor. The ninth, F5 `1:169586232:G:C`, was stored as a
`downstream_gene_variant` (MODIFIER) **because of a defect in this pipeline's
transcript selection**, not by design: VEP's own MANE transcript for *F5*
(ENST00000367797, `NM_000130.5`) calls it `stop_gained`, but `pick_tc()` in
`02_vep_annotate.py` takes the first MANE-tagged transcript in the list, which
at this overlapping locus belongs to the neighbouring gene *SELP*. The same
defect affects both COL17A1 alleles (assigned to *SLK*) and the RNF31 allele
(assigned to *PSME2*) — **four alleles across three cases**: one incomplete
pair (F5) and two absent genes (COL17A1, RNF31). All three stay in the reported
denominators. In the proband's own data, 275 of 29,907 annotated coding calls
(0.92%) carry a selected upstream/downstream consequence and are candidate
instances of the same defect; the three submitted alleles are not among them.
See METHODS §5.6.

These are ClinVar Pathogenic/Likely-pathogenic variants by construction, yet
AlphaMissense disagrees with that assertion. The policy therefore discards the
predictor that would have retained them in favour of one that contradicts
ClinVar. This is a real cost of a principled correction, and it is also the
reason `both` recovers HMGCR and PTRHD1.

**It bears directly on the challenge case.** *BUB1B* p.Asn1002Lys survives the
evidence floor *only* because AlphaMissense calls it `likely_pathogenic`
(+1.5). Had the call been `ambiguous`, this allele would have been dropped
exactly as those nine were — and the *BUB1B* hypothesis would have been
`incomplete_recessive` on the nonsense allele alone.

---

## 5. Limitations

1. **No competing alleles existed.** In **80 of 80** cases the injected gene
   contained **zero** other candidate alleles
   (`docs/competing_allele_availability.tsv`) — a direct consequence of
   excluding genes already in the proband's candidate set. So
   "**zero correct-allele-with-wrong-partner**" carries **no information about
   partner selection**: the method was never offered a wrong partner. An
   earlier version of this report presented that zero as evidence the pair
   score does not pad a correct allele. **That claim is withdrawn.**
2. **Phenotypes are derived from the answer.** Sampling the eight terms from
   the planted disease's own annotation makes the phenotype match far cleaner
   than reality. The real case has three specific matches of eight and one
   feature (nephrocalcinosis) probably unrelated to the diagnosis.
3. **One shared genomic background** across all 80 cases.
4. **Selected variant classes only.** Coding single-nucleotide variants with an
   existing ClinVar pathogenic assertion. Excluded: indels, structural
   variants, deep-intronic and regulatory variants, mobile-element insertions,
   and — importantly — **any variant absent from ClinVar**, the class to which
   the real p.Asn1002Lys belongs.
5. **Phase is never tested.** Trans cases are *constructed* trans and assumed
   so; there is no read-based phasing step anywhere in the benchmark.
6. **Cis controls are demoted, not suppressed** (§6).
7. **Test-set reuse.** These 80 cases are untouched by any tuning. The
   *historical* 80 were used to choose the `MISSENSE_EVIDENCE` default and are
   therefore development data; their figures must not be read as held-out.

---

## 6. Two reporting corrections

**Cis controls.** An earlier report claimed "no cis pair was ever proposed" and
"13 of 15 genes absent from the ranking". Both were **false**, caused by
judging absence from a field that holds only the top ten rows. Corrected, from
the full-ranking fields:

| | n of 15 |
|---|---:|
| cis pair proposed as a pair **somewhere in the full ranking**, at ranks **16–24** | **11** |
| cis pair **inside the top ten** | **0** |
| gene at rank 1 as a **single-allele** hypothesis (MDH2, DYNC2I2) | 2 |
| **absent from the full ranking** (COL17A1, DCDC2) | 2 |

So the ×0.05 *cis* penalty **demotes** the pair out of the reportable range; it
does not suppress it. "No cis pair in the top ten" is supported. "Never
proposed" is not.

**Hypothesis ranking versus submission rows.** The benchmark scores the
pipeline's *hypothesis ranking*. The real submission CSV is a different object:
it carries the leading pair plus **duplicate single-allele rows** for each of
its two alleles (a backstop, because the evaluator computes a partial match
only when the truth has exactly two variants), plus one secondary finding.
Benchmark rank points therefore do **not** translate directly into the
submission's score.

---

## 7. The sentence to use

> On 45 constructed cases with a complete causal truth, the method returned the
> exact variant configuration at rank 1 in 37 (82%, Wilson 95% CI 69–91%).
> This measures recovery of planted answers in cases built to be findable; it
> is **not** a probability that the prediction for this child is correct.
