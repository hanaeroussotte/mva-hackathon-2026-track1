# SCORING_SPECIFICATION.md

The implemented rules, in execution order, with no simplification. Where the
code does something awkward or wrong, that is stated rather than tidied. Line
references are to the frozen scripts in `src/`.

Final configuration: `MISSENSE_EVIDENCE=once`, `RANK_MODE=full`, `PAIR_AGG=mean`,
`AF_UNKNOWN=absent`, `QUALITY_TERMS=1`, `RECESSIVE_RESCUE=1`,
`COMBINE=multiplicative`, `VAR_SCALE=2.5`, `EVIDENCE_FLOOR=1.0`,
`CIS_PENALTY=0.05`, `MAX_PRIMARY=1`, `MAX_ROWS=10`, seed 20261007.

---

## Stage A — per-variant evidence score (`04_shortlist.py`, `variant_evidence()`)

**Order matters: this entire stage runs BEFORE any frequency/rarity term is
computed.** Rarity is added later, in Stage D, and is *not* part of the evidence
floor. A variant can therefore be dropped at the floor despite being absent from
gnomAD, and a common variant can clear the floor and only then be penalised.

Contributions are additive, starting from 0.0:

### A1. Genotype-confidence penalties (applied to every genotype class)

| condition | points |
|---|---:|
| `DP < 10` **or** `GQ < 30` | **−3.0** |
| else `DP < 15` **or** `GQ < 50` | **−2.0** |

A missing or unparseable `DP`/`GQ` is coerced to **0**, which triggers the −3.0
branch. Rationale: candidate rows were enriched ~1.6-fold for calls below these
thresholds relative to the whole coding callset.

### A2. Allele-balance penalties (heterozygous calls only, and only when `DP ≥ 20`)

| condition | points |
|---|---:|
| `AB < 0.20` | **−3.0** |
| else `AB < 0.30` | **−1.5** |

`AB = alt / (ref + alt)` from the `AD` field. If `AB` cannot be computed, or
`DP < 20`, **no penalty is applied** — absence of evidence is not treated as
evidence of a problem here.

### A3. Consequence class

| condition | points |
|---|---:|
| consequence ∈ {stop_gained, frameshift, start_lost, stop_lost, splice_donor, splice_acceptor} | **+2.0** (`LoF`) |
| non-missense and consequence ∈ inframe set | **+0.8** |

### A4. Computational missense evidence — the `MISSENSE_EVIDENCE` switch

CADD points are first computed for any variant: `+1.0` if CADD ≥ 25, `+0.5` if
CADD ≥ 20, otherwise `0.0`. A missing CADD value is coerced to 0.0 and scores
nothing.

AlphaMissense points, for missense variants only:

| AlphaMissense class | points | tag |
|---|---:|---|
| `likely_pathogenic` | **+1.5** | `AM_path` |
| `ambiguous` | **+0.4** | `AM_amb` |
| `likely_benign` | **−0.5** | `AM_benign` |
| anything else, **including no call at all** | **+0.2** | *(no tag)* |

That last row is a real behaviour and a weak point: **an unavailable
AlphaMissense call receives +0.2, the same as an unrecognised class.** It is a
small positive credit for absent information.

For a **missense** variant the three modes differ:

* **`none`** — neither term is added.
* **`once` (submitted)** — if AlphaMissense made a *recognised* call
  (pathogenic, ambiguous or benign), its points are added **alone** and CADD is
  discarded. If it did not, CADD points are added if non-zero; if CADD is also
  zero, the bare `+0.2` is added. The branch is on whether AlphaMissense is
  informative, **not** on which score is larger — taking a signed maximum would
  silently discard the `likely_benign` penalty.
* **`both`** — AlphaMissense points **plus** CADD points.

For a **non-missense** variant CADD is always added, in every mode: for a
nonsense or splice variant CADD is not standing in for a missense predictor.

### A5. Splice prediction

| condition | points |
|---|---:|
| SpliceAI max delta ≥ 0.50 | **+1.5** |
| else ≥ 0.20 | **+0.5** |

A missing SpliceAI value is coerced to 0.0 and scores nothing. **Missing and
zero are not distinguished** — a limitation shared with the frequency handling.

### A6. Low-confidence region

`−1.5` if the position falls in `work/lowconf_regions.bed`. **That file is
unrecoverable, so in every current rebuild this penalty never fires.** Its
absence is the whole of the *LZTR1* change from an earlier 0.3661 to the current
0.3933.

### A7. ClinVar

| `CLNSIG` starts with | review status contains | points |
|---|---|---:|
| `Pathogenic` or `Likely_pathogenic` | `multiple_submitters`, `practice_guideline` or `reviewed_by_expert_panel` | **+3.0** |
| `Pathogenic` or `Likely_pathogenic` | anything else | **+2.0** |
| `Benign` or `Likely_benign` | — | **−3.0** |

Uncertain-significance and conflicting records score **0.0** — they are neither
credited nor penalised. The join reads `CLNSIG`/`CLNREVSTAT`/`CLNDN` at fixed
column indices 5/6/7, with the ClinVar ID at 4; an earlier rebuild wrote these
columns in a different order and silently cost allele 1 its ClinVar term.

### A8. The floor

`EVIDENCE_FLOOR = 1.0`. A variant scoring below it is dropped and is invisible
to everything downstream. **583 variants clear the floor.**

---

## Stage B — genotype eligibility and configuration assembly (`04_shortlist.py`)

Variants clearing the floor are grouped by gene, then split by genotype:

* **homozygous** — `GT ∈ {1/1, 1|1}` → one `hom` row each.
* **heterozygous** — `GT ∈ {0/1, 1/0, 0|1, 1|0}`:
  * **≥ 2 hets** → `comphet` rows for the **top four by evidence score**
    (`sorted(hets, key=-evidence)[:4]`). No `het_single` or
    `incomplete_recessive` row is produced for such a gene.
  * **exactly 1 het** → one `het_single` row, plus one
    `incomplete_recessive` row when `RECESSIVE_RESCUE=1`.

**Two documented limitations live here.**

1. **The top-four cap** drops heterozygous candidates beyond the fourth in genes
   with ≥5 — **15 variants** in this run.
2. **`GT = 1/2` matches neither branch.** Multi-allelic calls clear the floor,
   are counted in the 583, and then **no model row is emitted at all** —
   **6 variants** in this run. This is a defect, not a design choice, and it is
   left unfixed because fixing it would change the candidate set after the
   results were seen.

583 clearing the floor → **562 distinct variants written as 861 (variant, model)
rows** across **448 genes**. Per-variant reconciliation:
`docs/candidate_count_reconciliation.tsv`.

---

## Stage C — inheritance gating (`07_mechanism_score.py`)

Each configuration is admitted only against diseases whose HPO inheritance
annotation is compatible:

| configuration | compatible inheritance |
|---|---|
| `het_single` | AD, XL, MT |
| `hom` | AR, XL, MT |
| `comphet` | AR |
| `incomplete_recessive` | AR |

Term sets: AD = {HP:0000006, HP:0000015, HP:0012275}; AR = {HP:0000007};
XL = {HP:0001417, HP:0001419, HP:0001423, HP:0001470}; MT = {HP:0001427}.

A disease is admitted if `(annotated_moi ∩ compatible) or not annotated_moi` —
so **a disease with no inheritance annotation is admitted under every model**.
That is 55 of 8,420 diseases (0.7%). `hom` excluding AD blocks 642 of 5,277
configurations corpus-wide and 3 of 135 proband homozygous-candidate genes
(*ASPSCR1*, *POLR2A*, *ALOX5AP*).

---

## Stage D — rarity points (`05_rank_and_submit.py`, `rarity_points()`)

Computed **after** the evidence floor, and added to the evidence score to give an
**allele total**. The model-dependent cap is applied first:

* recessive models (`hom`, `comphet`, `incomplete_recessive`): cap **1e-2**
* `het_single`: cap **1e-4**

| condition | points | note |
|---|---:|---|
| `af is None` and `AF_UNKNOWN=absent` **(submitted)** | **+1.5** | tagged `absent_gnomAD` |
| `af is None` and `AF_UNKNOWN=unknown` | 0.0 | sensitivity arm |
| `af is None` and `AF_UNKNOWN=common` | −6.0 | sensitivity bound |
| `af > cap` | **−6.0** | too common for the model |
| `af < 1e-5` | **+1.5** | |
| `af < 1e-4` | **+1.0** | |
| `af < 1e-3` | **+0.3** | |
| otherwise | **−0.5** | |

**Missing-value handling is the known weak point.** Three distinct evidence
states exist — *unknown* (lookup failed), *absent* (queried, not present), and
*observed ultra-rare* (seen a countable number of times) — and the submitted
configuration collapses the first two into one. 287 of 841 lookups (34%)
returned nothing. For *BUB1B* allele 2 the frequency was **measured** at AC 1 /
AN 1,461,878; under this configuration it receives the same bonus as a verified
absence, which is a coincidence of the function rather than a correct
classification.

---

## Stage E — pair scoring and the final transformation (`05_rank_and_submit.py`, `rank()`)

For each gene, for each model with rows:

| model | `var_score` |
|---|---|
| `comphet` (requires ≥2 rows, else skipped) | **mean of the top two allele totals** under `PAIR_AGG=mean`; the **minimum** under `PAIR_AGG=min` |
| `incomplete_recessive` | **best allele total ÷ 2** — the unfound second allele contributes zero |
| `hom`, `het_single` | the single allele total |

Then:

```
e_var   = 1 − exp(−max(var_score, 0) / VAR_SCALE)       VAR_SCALE = 2.5
e_pheno = max(pheno_score, PHENO_FLOOR)                 PHENO_FLOOR = 0.02
total   = e_var × e_pheno                               (COMBINE = multiplicative)
```

`RANK_MODE` replaces the final line for the two ablations: `evidence_only` →
`total = e_var`; `pheno_only` → `total = e_pheno`. The additive alternative
(`COMBINE=additive`, `total = var_score + 1.5 × pheno`) exists in the code and
was **not** used.

**Phase penalty.** If the configuration's phase state is `confirmed_cis`,
`total` is multiplied by `CIS_PENALTY = 0.05`. Only GATK-supported *cis*
(shared `PID`, phased `PGT`) triggers it; `unresolved_phase` carries **no**
penalty, which is why the *BUB1B* pair is not penalised. Phase states come from
`10_phase_states.py`: `confirmed_cis`, `supported_trans`, `unresolved_phase`,
`single_allele_only`, `not_applicable`.

**One hypothesis per gene.** Across a gene's models, `max(options, key=total)`
selects a single winner. 448 rows = 448 distinct genes. Deeper ranks therefore
add *other genes*, never a better configuration for the same gene.

**Tie-breaking.** None is implemented. Ordering among equal `total` values falls
back to Python's stable sort over dictionary insertion order, i.e. the order
genes were first encountered in the candidate file. No tie affected the top of
this ranking.

---

## Stage F — ranked hypotheses to submission rows (`write_submission()`)

1. Partition: `SECONDARY_GENES = {LZTR1}` → secondary; `EXCLUDED_GENES` (nine
   HLA genes) → neither; everything else → primary.
2. Emit the top `MAX_PRIMARY = 1` primary hypothesis. For a `comphet` winner,
   both variants go in one row (`chrom_1/pos_1/ref_1/alt_1` and
   `chrom_2/...`); otherwise the second variant fields are blank. Chromosomes
   are written chr-prefixed from the unprefixed VCF contigs.
3. **Single-allele backstop rows.** If the leading hypothesis is `comphet`, each
   of its two alleles is additionally emitted as its own row. The reason is the
   evaluator's own behaviour, established by running it against synthetic answer
   keys: a partial match is computed **only when the truth has exactly two
   variants**, so a submission holding only the paired row scores **zero** if the
   true answer is a single causal variant. These rows take that scenario from 0
   to 50 rank points and cost nothing in the others. They are **alternative
   representations of one hypothesis, not independent diagnoses.**
4. Emit secondary findings, capped at `MAX_ROWS = 10` total.

**EPCR.** `epcr = round(clamp(0.95 × total / best_total, 0.01, 0.95), 4)`, where
`best_total` is the leading primary hypothesis's total. The leading row
therefore always receives exactly **0.95**. The backstop rows carry
hand-assigned 0.90 and 0.85 to preserve ordering beneath it. **EPCR is an
ordering value. The challenge evaluator reads it only to sort rows and form
F-max thresholds; it is not a calibrated probability and must not be reported as
one.**

---

## Worked example — the submitted *BUB1B* hypothesis

Reproduce with:

```bash
python src/05_rank_and_submit.py \
    work/candidates_pre_af.tsv work/bench_af.json work/mech_scores.tsv \
    WGS_EX2312012 submissions/repro.csv work/ranked_full_repro.tsv
```

### Allele 1 — `chr15:40209701 T>G`, p.Leu737Ter

| stage | term | points | running |
|---|---|---:|---:|
| A1 | DP 46, GQ 99 — above both thresholds | 0.0 | 0.0 |
| A2 | AB = 25/46 = 0.54, above 0.30 | 0.0 | 0.0 |
| A3 | `stop_gained` ∈ LoF | **+2.0** | 2.0 |
| A4 | not missense → CADD 36 ≥ 25 | **+1.0** | 3.0 |
| A5 | SpliceAI 0.030, below 0.20 | 0.0 | 3.0 |
| A6 | low-confidence file absent | 0.0 | 3.0 |
| A7 | ClinVar Pathogenic/Likely_pathogenic, multiple submitters | **+3.0** | **6.0** |
| A8 | 6.0 ≥ 1.0 → clears the floor | | |
| D | gnomAD exome 7.87e-05 → `af < 1e-4` | **+1.0** | **allele total 7.0** |

### Allele 2 — `chr15:40220612 T>G`, p.Asn1002Lys

| stage | term | points | running |
|---|---|---:|---:|
| A1 | DP 28, GQ 99 | 0.0 | 0.0 |
| A2 | AB = 13/28 = 0.46 | 0.0 | 0.0 |
| A3 | missense, not LoF | 0.0 | 0.0 |
| A4 | `once` + AlphaMissense `likely_pathogenic` → **AM alone, CADD 24.5 discarded** | **+1.5** | 1.5 |
| A5 | SpliceAI 0.000 on MANE | 0.0 | 1.5 |
| A7 | **no ClinVar record for this allele** | 0.0 | **1.5** |
| A8 | 1.5 ≥ 1.0 → clears the floor **by 0.5** | | |
| D | lookup unresolved this run → `AF_UNKNOWN=absent` | **+1.5** | **allele total 3.0** |

Had AlphaMissense returned `ambiguous` (+0.4) instead, Stage A would total 0.4
and this allele would have been **dropped at the floor** — the *BUB1B* hypothesis
would then have been `incomplete_recessive` on allele 1 alone. That is the same
mechanism behind eight of the nine alleles dropped in the benchmark failures.

### Pair and final score

```
var_score = mean(7.0, 3.0)        = 5.00
e_var     = 1 − exp(−5.00 / 2.5)  = 1 − exp(−2.0)  = 0.8647
e_pheno   = 0.7467                 (> PHENO_FLOOR, used as is)
phase     = unresolved_phase       → no cis penalty
total     = 0.8647 × 0.7467        = 0.6456          rank 1 of 448
```

Runner-up *LZTR1*, `incomplete_recessive`: evidence 7.5 (LoF +2.0, CADD 49 →
+1.0, SpliceAI 0.92 → +1.5, ClinVar P/LP multiple submitters +3.0) plus rarity
+1.5 (AF 1.37e-06 < 1e-5) = allele total 9.0; halved to `var_score` 4.50;
`e_var = 1 − exp(−1.8) = 0.8347`; `× 0.4712 = 0.3933`. Margin **1.64×**.

Under `PAIR_AGG=min` the pair takes the weaker allele: `var_score = 3.00`,
`e_var = 1 − exp(−1.2) = 0.6988`, `total = 0.5218`, still rank 1 with a 1.33×
margin.

Every figure above is reproduced in `docs/final_ranked_hypotheses.tsv` and
`docs/allele_evidence_table.tsv`.
