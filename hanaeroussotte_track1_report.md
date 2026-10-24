# Track 1 — causal variant prediction for proband WGS_EX2312012

**Rare Disease, Real Kid: The MVA Hackathon 2026**
Submitted by: hanaeroussotte · Approach name: `mechanism-gated`
Prediction file: `hanaeroussotte_mechanism-gated_FINAL.csv` (md5 `d1b74e84be505744149e1cab1877d27a`)

> **AI assistance disclosure** (per the challenge update of 28 August 2026):
> Anthropic, Claude Opus 5, accessed through Claude Science. *The plan/tier and
> the applicable data-handling setting must be confirmed by the submitter before
> this report is filed — see `docs/CONTRIBUTIONS.md`.*

---

## 1. Summary and objective

**Objective.** Given one child's singleton whole-genome VCF and eight HPO
phenotype terms, nominate the variant configuration most likely to explain the
presentation, and state honestly what the evidence does and does not support.

**Prediction.** Two variants in *BUB1B*, proposed as a compound-heterozygous
pair:

| | coordinate (GRCh38) | HGVS | consequence |
|---|---|---|---|
| allele 1 | `chr15:40209701 T>G` | `NM_001211.6:c.2210T>G` | p.Leu737Ter, stop_gained, exon 17/23 |
| allele 2 | `chr15:40220612 T>G` | `NM_001211.6:c.3006T>G` | p.Asn1002Lys, missense, exon 23/23 |

Implied disease: **Mosaic variegated aneuploidy syndrome 1** (OMIM:257300), the
recessive chromosome-instability disorder caused by biallelic *BUB1B* variants
[1]. Rank 1 of 448 ranked hypotheses, score 0.6456, 1.64× the runner-up.

**Two requirements of this hypothesis are unresolved, and they are stated here
rather than in a closing caveat:**

1. **Phase.** Whether the two variants lie on opposite copies of *BUB1B* is
   unresolved in the analysed reads and marker set. A recessive pair requires
   *trans*; a *cis* configuration would undermine this proposed two-allele
   explanation.
2. **The second allele.** No allele-specific functional experiment for
   p.Asn1002Lys was identified in the reviewed evidence. Its support is one
   computational call plus extreme rarity.

One supplied clinical feature — nephrocalcinosis — is **not accounted for** by
this diagnosis.

This is a prioritisation of candidates from a singleton genome, not a diagnosis,
and not a clinical report.

---

## 2. Input data and quality assessment

One GATK HaplotypeCaller VCF (GRCh38, hard filters, no VQSR) [16], md5
`a04ea354141cfb032872d67a11c8d2d8`, 315,153,971 bytes. No parental samples.
Reference contigs are unprefixed (`1`…`22`,`X`,`Y`,`M`); submission coordinates
are chr-prefixed, and the three submitted REF/ALT pairs were re-read from the
source VCF to confirm that mapping.

| metric | value |
|---|---|
| records / PASS | 5,012,204 / 4,740,790 |
| mean depth | 45.1× |
| Ti/Tv | 1.94 |
| chrM calls | 13 |
| runs of homozygosity >1 Mb | 3 segments, 3.47 Mb total, none on chr15 |

Genotype quality enters the scoring directly as a penalty rather than a filter
(§4 of `docs/SCORING_SPECIFICATION.md`): candidate rows were enriched ~1.6-fold
for calls below DP 15 / GQ 50 relative to the whole coding callset, which is the
signature of low-depth calls generating spurious rare-and-damaging annotations.

---

## 3. Annotation, filtering and candidate generation

| step | remaining | output |
|---|---:|---|
| PASS records | 4,740,790 | — |
| exonic ±10 bp | 180,564 | — |
| coding, after `bcftools norm -m -any` [19] | **29,907** | `work/cds_vep.tsv` |
| clearing the evidence floor | **583 variants** | — |
| written to the candidate file | **562 distinct variants as 861 (variant, model) rows** | `work/candidates_pre_af.tsv` |
| genes with ≥1 candidate | **448** | — |

The 583 → 562 difference is not a discrepancy between runs; both figures come
from the same run and measure different things. 15 of the 21 unwritten variants
are heterozygous candidates beyond a top-four cap in genes with ≥5 heterozygous
candidates; **6 carry a `1/2` multi-allelic genotype that the shortlist code
matches to neither the homozygous nor the heterozygous branch, so no row is
emitted at all.** That second case is an implementation defect, disclosed in §7
and reconciled variant-by-variant in `docs/candidate_count_reconciliation.tsv`.

Annotation is Ensembl VEP [9]: consequence class, CADD [10], AlphaMissense [11],
SpliceAI [12], gnomAD r4 frequency [14], and a ClinVar join [13] (release
2026-10-19). Each candidate receives an additive evidence score — loss of
function +2.0, AlphaMissense likely_pathogenic +1.5, CADD ≥25 +1.0, SpliceAI
≥0.5 +1.5, ClinVar Pathogenic/Likely-pathogenic with multiple submitters +3.0,
with quality and low-confidence-region penalties — and the floor is 1.0. Every
contribution, its sign, and its order relative to the rarity term is specified
in `docs/SCORING_SPECIFICATION.md`.

For missense variants the pipeline counts **one** computational line of
evidence, not two: AlphaMissense's call stands and CADD is discarded, because
the two predictors share training features and calibration sets and adding them
double-counts one signal. This correction has a measured cost (§6).

**The consequence class, the predictor scores and the ClinVar assertion are not
independent of one another.** A truncating change drives all three. The additive
score is a ranking device, not a likelihood.

---

## 4. Inheritance-compatible configuration ranking and phenotype matching

This is the central methodological choice. The pipeline ranks **allele
configurations**, not genes. No systematic comparison against other
prioritisation tools was run, so no claim is made about what other tools do; the
motivation is internal to the task.

A gene is not a hypothesis. *BUB1B*-related MVA1 is recessive, so one damaged
copy is a carrier state and a single-allele finding in that gene is not a
candidate explanation for a recessive disease. Each gene's candidates are
therefore grouped into configurations — `hom`, `comphet` (≥2 heterozygous
alleles), `het_single`, and `incomplete_recessive` (one damaging allele in a
recessive gene, partner not found, scored at **half** the allele total) — and
each configuration is admitted only against diseases whose inheritance
annotation is compatible with it. Diseases with no inheritance annotation (55 of
8,420, 0.7%) are admitted under every model.

**What this implements is disease-level phenotype matching gated by
compatibility with the observed allele configuration.** It is not a mechanistic
model of disease biology.

The phenotype term is an information-content-weighted match between the
proband's eight HPO terms [15] and each admissible disease's annotation. Final
score `total = e_var × pheno`, with `e_var = 1 − exp(−var_score / 2.5)`; for a
pair, `var_score` is the mean of the two allele totals. **One hypothesis is
retained per gene** — the best-scoring configuration — giving 448 hypotheses
from 448 genes.

| rank | gene | configuration | var_score | pheno | total |
|---:|---|---|---:|---:|---:|
| 1 | *BUB1B* | comphet | 5.00 | 0.7467 | **0.6456** |
| 2 | *LZTR1* | incomplete_recessive | 4.50 | 0.4712 | 0.3933 |
| 3 | *HLA-DRB1* | comphet | 4.25 | 0.4632 | 0.3786 |

### What the phenotype score actually represents

Of the eight supplied features, **three are specific matches** and one of those
carries most of the weight:

| feature | match | information content |
|---|---|---:|
| Rhabdomyosarcoma | → embryonal rhabdomyosarcoma (disease term more specific) | **6.42** |
| Small for gestational age | exact | 3.65 |
| Short stature | exact | 1.75 |
| Premature birth | generic ancestor only | 2.37 |
| Nephrocalcinosis | generic ancestor only | 2.16 |
| Failure to thrive | generic ancestor only | 1.93 |
| Skeletal muscle atrophy | generic ancestor only | 0.86 |
| Recurrent miscarriage (**parental**) | ontology root — no real overlap | **0.00** |

Corpus maximum is 9.46. The summary statistic is a **similarity of 0.75**; it is
not a count of features biologically explained, and the pipeline cannot
distinguish a proband feature from a parental one.

**Nephrocalcinosis is unexplained.** A PubMed search for *BUB1B* with
nephrocalcinosis returns no records, and the renal phenotype reported in MVA is
cystic rather than calcific. It occurs in a substantial minority of infants born
before 32 weeks [6]; this child was born at 32 weeks at about 1 kg, so
prematurity is a possible explanation requiring clinical context, not a
demonstrated one. Features common in *BUB1B*-positive MVA that do not appear in
the supplied list — microcephaly, Dandy-Walker complex, cataract [1,5] — are
**not reported**, which is not the same as clinically absent.

### Worked example of the configuration gate

*LZTR1* carries one truncating allele, p.Tyr748Ter. Two distinct things follow,
and the submission keeps them apart. The **ranking hypothesis** is
`incomplete_recessive` toward autosomal-recessive Noonan syndrome: one damaging
allele found, partner not found, allele total 9.0 halved to `var_score` 4.50,
giving rank 2 at 0.3933. The **reporting rationale** is different — a single
truncating *LZTR1* allele is separately relevant to dominant *LZTR1*-related
schwannomatosis, which is why it appears in the submission at all, as
`finding_type=secondary` with an EPCR of 0.02. Neither is proposed as this
child's diagnosis, and neither is a management recommendation.

---

## 5. Leading prediction and allele-specific investigation

The two alleles are **not** equally supported, and the submission says so. The
complete field-by-field table, each row labelled OBSERVED / DATABASE / PREDICTED
/ HYPOTHESIS, is `docs/allele_evidence_table.tsv`.

| | p.Leu737Ter | p.Asn1002Lys |
|---|---|---|
| consequence | stop_gained, exon 17/23 | missense, exon 23/23 (last) |
| transcript fate | premature stop 747 nt before the last exon–exon junction → nonsense-mediated decay **predicted from transcript geometry**, not measured | **no premature stop is introduced**; a single amino-acid substitution, so no decay rule applies |
| ClinVar | **VCV000533901.9** Pathogenic/Likely pathogenic, criteria provided, multiple submitters, no conflicts (SCV000762865, SCV004031025), annotated to MVA syndrome 1, last evaluated 2024-10-09 | **no record for this exact allele** |
| gnomAD r4 exomes | AC 115 / AN 1,461,846 = 7.87e-05, 0 homozygotes | AC 1 / AN 1,461,878 = 6.84e-07 — **observed once** |
| CADD | 36 | 24.5 |
| AlphaMissense | n/a | likely_pathogenic 0.9229 — a **model output for research prioritisation**, not an ACMG/AMP classification [8] |
| SpliceAI | 0.030 | 0.000 on the MANE transcript; 0.020 is the maximum across all 24 scored transcripts |
| evidence / rarity / allele total | 6.0 / +1.0 / **7.0** | 1.5 / +1.5 / **3.0** |

Allele 1 supplies about 70% of the pair score. **Allele 2 clears the evidence
floor only because AlphaMissense calls it likely_pathogenic**; an `ambiguous`
call would have dropped it below the floor and the *BUB1B* hypothesis would have
been `incomplete_recessive` on the nonsense allele alone. That is the same
mechanism behind most benchmark failures (§6).

### Is p.Asn1002Lys destabilising? A hypothesis, not a result

Suijkerbuijk *et al.* [2] (full text read) report **two distinct classes** of MVA
substitution, and the five they tested must not be described collectively as
destabilising:

| substitution | protein abundance | function |
|---|---|---|
| I909T | reduced; turnover ~2× wild type | checkpoint response restored by forced overexpression |
| R727C | reduced; turnover ~2× wild type | — |
| L1012P | reduced | checkpoint response restored by forced overexpression |
| **Y155C** | **not reduced** | **separation-of-function**: restored chromosome alignment to BUBR1-depleted cells but failed to reconstitute mitotic checkpoint activity |
| **R550Q** | **not reduced** | rescued both BUBR1 functions |

Residue 1002 was not tested. Y155C is why the two questions must stay apart:
**normal abundance would weaken a destabilisation hypothesis for p.Asn1002Lys
but would not exclude a functional defect.** The paper's generalisation is
positional — destabilising changes "consistently occur in or near the BUBR1
kinase domain" — but the boundary is itself ambiguous (InterPro IPR011009:
791–984, excluding residue 1002; UniProt O60566: 766–1050, including it), and
p.Leu1012Pro ten residues away is a ClinVar variant of **uncertain
significance** (VCV000006765.53). Adjacency is not evidence.

The same protein substitution from a **different DNA allele**, `c.3006T>A`, is
ClinVar **VCV004600147.1, Uncertain significance**, single submitter, last
evaluated 2025-09-19. That record was retrieved and confirmed here. Its
classification is **not transferred**: different nucleotide change, different
submission, and an uncertain-significance assertion carries no directional
weight in either direction.

### Read support and phase

Three checks: the GATK genotypes (0/1, AD 21,25 and 15,13, GQ 99, PASS); exact
31-mer counting on raw FASTQ, involving no aligner; and a targeted re-alignment
with minimap2 [18]. **The latter two read the same sequencing material**, so
they are independent of the original alignment, not independent observations.
Read evidence supports the variants' presence; it does not establish causality.
A cautionary case from the same callset: *TMEM216* `c.432-1G>C` had allele
balance 0.49, DP 55, GQ 99 and PASS, yet returned **zero** alternate 31-mers.

`whatshap 2.8` [17] over the three heterozygous markers at the locus returned
**0 reads covering any two markers and 0 phase blocks**. Marker gaps are 6,769 bp
and 4,142 bp; the fragment distribution at the locus (n = 1,085) has median
436 bp and **maximum 798 bp**. More depth on this library would not help — the
limit is fragment length against gap size — though a longer-insert library or
informative markers inside the gaps could resolve it. No general claim is made
about short-read data.

### What was searched for an alternative second allele, and what was not

*Searched.* Called small variants (SNVs and short indels) in the 60 kb upstream
of *BUB1B*: 85 PASS records, 52 heterozygous, all annotated, **no CADD filter
applied**, none both rare and heterozygous. A regulatory small variant of the
class reported 44 kb upstream of a *BUB1B* transcription start site [4] was
therefore **within the search space in that interval**, and none was found. This
does not exclude regulatory mechanisms generally: the interval is finite, only
called variants were examined, and a regulatory effect carried by a structural
change would not appear in this scan.

*Not assessable.* Structural variants and mobile-element insertions at the
locus — soft-clip clusters reach 7–75× local depth and the insertion-anchor
background is 0.578, so no call is supportable either way; a 3 kb intronic
mobile-element insertion of the kind reported in another case would not have
been detected. Single-exon copy-number changes were underpowered (power 0.349 at
one exon, 0.973 at 1 kb, under an idealised Poisson model).

---

## 6. Synthetic benchmark, baselines and sensitivity analyses

### Why two endpoints were measured

The original validation measured whether the **correct gene** reached rank 1.
Gene prioritisation is a reasonable intermediate endpoint, but it is
insufficient on its own for a task scored on variant identity: a correct gene
with the wrong allele pair earns partial credit at best. Gene recovery and
exact-configuration recovery were therefore measured separately. On
complete-truth cases the measured difference is **39/45 correct gene versus
37/45 exact configuration**; on the `comphet_trans` stratum, 26/30 versus 24/30.

### Construction

80 cases, seed 20261007. Each plants one or two variants of known identity into
this proband's own callset; the truth label is exactly those variants. Truth
variants are ClinVar Pathogenic/Likely-pathogenic with a coding molecular
consequence and frequency <1e-3; genes must have a recessive disease with ≥8
annotated terms and must not be an MVA gene, a submitted gene, or a gene already
in this proband's candidate set. Strata: `comphet_trans` 30 (≥1,000 bp apart),
`comphet_cis` 15 (≤60 bp, shared phase tag), `single_allele` 20, `hom` 15.
Genotype quality is resampled from the proband's own calls. Phenotypes are
sampled from the planted disease's annotation. All arms run with the planted
variants **removed from the ClinVar evidence file**.

This is a **new benchmark, not a reproduction** of an earlier one: the ClinVar
release changed, the earlier runner never recorded its injected alleles, and the
current cases share only 4 of 80 genes with the earlier set. Earlier
gene-level figures are superseded and are not quoted here. The compensating
advantage is that these 80 cases were used to tune nothing.

### Results — primary configuration

| stratum | n | correct gene @1 | **exact config @1** | @3 | one allele only | absent from ranking |
|---|---:|---:|---:|---:|---:|---:|
| comphet_trans | 30 | 26 | **24** | 24 | 5 | 1 |
| hom | 15 | 13 | **13** | 13 | 0 | 2 |
| single_allele | 20 | 17 | 17 | **19** | 0 | 1 |
| comphet_cis | 15 | 2 | **0** | 0 | 2 | 2 |

**Complete-truth: 37 of 45 (82%), Wilson 95% CI 69–91%.** This is recovery of
planted answers in cases built to be findable. **It is not a probability that
the prediction for this child is correct.** Applying the challenge's own
`evaluation.py` to the hypothesis ranking gives mean rank points 85.17 and mean
F-max 0.883 on `comphet_trans`, 86.67 and 0.867 on `hom`.

Depth helps only where the correct configuration was already formed but ranked
below first: exact recovery improves with depth in the `single_allele` stratum
(17 → 19 at rank 3) and **not** in the complete-truth strata, because the ranking
holds one hypothesis per gene.

*Cis* controls are **demoted, not suppressed**: 11 of 15 are still proposed as
pairs, at ranks 16–24; **0** reach the top ten; 2 rank first as single-allele
hypotheses; 2 are absent. "No *cis* pair in the top ten" is supported; "never
proposed" is not.

### Baselines

Three arms over identical cases, pool, annotations, inheritance gate and truth
definitions. They are **internal ablations of the same pipeline**, not
independent external tools.

| arm | rule | exact @1 (45) | gene @1 (45) |
|---|---|---:|---:|
| evidence + rarity only | `total = e_var` | **0/45** | 0/45 |
| phenotype only | `total = e_pheno` | **37/45** | **42/45** |
| full (submitted) | `total = e_var × e_pheno` | **37/45** | 39/45 |

Paired exact McNemar, full versus phenotype-only: exact configuration **37 vs
37, zero discordant pairs, p = 1.000**; correct gene **39 vs 42, p = 0.250**,
all three discordant cases favouring the baseline.

**Within the shared filtered pipeline, the phenotype term alone matched exact
recovery.** Two things this does not establish: it says nothing about the value
of the filtering and annotation stages, because every arm is built on the same
shortlist and none varies it; and non-significance on 45 cases does not prove
equivalence. The evidence-only arm's 0/45 shows only that variant evidence and
rarity alone cannot order these candidates.

### Sensitivity analyses, kept separate from the primary result

| analysis | primary | alternative | effect |
|---|---|---|---|
| missense evidence | `once` 37/45 | `both` 39/45 · `none` 26/45 | `both` recovers exactly HMGCR and PTRHD1, the two the failure trace predicts; p = 0.50, so `both` is **not** claimed better. `none` being far worse confirms missense evidence matters somewhere |
| pair aggregation | `mean` | `min` (weaker allele governs) | identical on every reported endpoint, but **not** identical rankings: 25 of 80 scores differ, 2 ranks change. `mean` retained because nothing measured favours either |
| unresolved frequency | `absent` 37/45 | `unknown` **36/45** · `common` 17/45 | the one change under `unknown` is **COQ7**, rank 1 → 3. `common` is an implausible bound, not a scenario |

*BUB1B* is rank 1 under the submitted configuration, under `PAIR_AGG=min`
(0.5218), under `pheno_only` (0.7467) and under `AF_UNKNOWN=unknown` (0.6103);
it is **rank 3 under `evidence_only`**, where *LZTR1* leads at 0.9727, and
**rank 2 under `AF_UNKNOWN=common`**. The same *BUB1B* pair remains the selected
configuration across these analyses.

### Failure analysis

Of the 8 complete-truth cases that did not return the exact configuration at
rank 1: **5** lost one allele during filtering, so `comphet` was unavailable and
`incomplete_recessive` was selected; **3** never entered the candidate set; **0**
formed a complete pair but ranked it below first; **0** selected a wrong partner.

**One mechanism explains eight of the nine dropped alleles.** They are missense
variants AlphaMissense calls *ambiguous* or *likely_benign*; under the `once`
policy CADD is not added for missense, so a CADD of 21–28 contributes nothing
and the allele falls below the 1.0 floor. These are ClinVar
Pathogenic/Likely-pathogenic variants by construction, so the policy discards
the predictor that agrees with ClinVar in favour of one that disagrees.

The ninth dropped allele exposed a **defect in this pipeline's transcript
selection** (§7).

**Zero wrong partners carries no information about partner selection.** In 80 of
80 cases the injected gene contained zero competing candidate alleles, a direct
consequence of excluding genes already in the proband's candidate set. The
method was never offered a wrong partner.

---

## 7. Limitations, reproducibility and conclusion

### Implementation defects found during audit and deliberately left unfixed

Fixing either would change the candidate set after the results had been seen.
Both are documented instead.

1. **Transcript selection.** The annotation step reduces VEP's transcript list
   to one consequence by taking the **first MANE-tagged transcript**. Where two
   genes overlap, both carry a MANE transcript and the neighbour's can come
   first, so a coding consequence in the intended gene is discarded for a
   MODIFIER consequence in the neighbour. Confirmed on four benchmark alleles
   across three cases by re-querying all transcripts: the intended gene's MANE
   transcript carried `stop_gained` or `missense_variant` in every case, and
   VEP's own `most_severe_consequence` was the coding term. **275 of 29,907
   proband coding calls (0.92%)** carry a selected upstream/downstream
   consequence and are candidate instances
   (`docs/transcript_selection_suspects.tsv`); whether each is a genuine
   misassignment was not individually checked. **The three submitted alleles are
   not affected** — both *BUB1B* alleles are annotated on ENST00000287598
   (`NM_001211.6`) and *LZTR1* on ENST00000646124 (`NM_006767.4`).
2. **Unhandled `1/2` genotypes.** Six multi-allelic calls clear the evidence
   floor but match neither genotype branch of the shortlist, so no model row is
   emitted. None is in a gene relevant to the prediction.

### Genuinely unresolved science

* **Phase** — unresolved in the analysed reads and marker set. *Cis* would
  undermine the proposed pair.
* **p.Asn1002Lys** — no allele-specific functional experiment identified;
  residue 1002 untested; destabilisation is a hypothesis.
* **Nephrocalcinosis** — no genetic candidate.
* **An alternative second allele** — structural and mobile-element classes not
  assessable, and these are where a missing partner would most plausibly sit.

### Benchmark limitations

No competing alleles existed in any injected gene (80/80). Phenotypes derive
from the planted disease's own annotation. One shared genomic background.
Coding SNVs with an existing ClinVar pathogenic assertion only — the real second
allele is absent from ClinVar, a class the benchmark contains none of. Phase is
never tested. The `MISSENSE_EVIDENCE` default was selected using an **earlier**
benchmark, which makes that choice exploratory, not pre-specified.

### Reproducibility

Three different things are called reproduction here, and only the first runs
from the submitted repository alone.

| level | status |
|---|---|
| **ranking rerun** | **executed in this pass** from the packaged tree: exit 0, *BUB1B* 0.6456 at rank 1, *LZTR1* 0.3933, 448 hypotheses; emitted rows match the submitted CSV on variants, ordering, EPCR and finding type |
| **cached-input reproduction** | previously executed in a clean directory (`docs/REPRO_LOG.txt`); needs two withheld genotype-level files |
| **raw-data reproduction** | **not performed.** Nine of eleven scripts not run, each missing input named. `08_kmer_validate.py` cannot be run — the FASTQ is deleted under the dataset terms |

A fresh `conda env create` from `environment.yml` was **not** tested.
`work/lowconf_regions.bed` is unrecoverable, so its −1.5 penalty never fires and
*LZTR1* scores 0.3933 rather than an earlier 0.3661. Three of the eleven scripts
are reconstructions written after a workspace cleanup; their equivalence to the
originals was never tested and is not claimed.

### Conclusion

A compound-heterozygous *BUB1B* pair is the best-supported account of this
child's tumour and growth phenotype, and it is prioritised first under every
scoring configuration tested except the two named in §6. The nonsense allele
carries independent database support; the missense allele does not, and the
phase that a recessive explanation requires could not be resolved. The two
measurements that would move this furthest are parental genotyping or
long-read sequencing for phase, and allele-specific BUBR1 protein
quantification in the child's cells for the second allele. Until then this is a
prioritised hypothesis for follow-up.

---

## 8. References

Two kinds of checking are reported and must not be conflated. **Bibliographic
verification**: every DOI was resolved through Crossref in this work
(`docs/reference_verification.tsv`, 20 of 20) — the reference exists and is
correctly cited. **Substantive verification**: whether the source supports the
claim attached to it was checked only where the article text was read, which
applies to **[2] and [21]**. Entries marked *record only* are bibliographically
verified but substantively not; the claim rests on the citing literature.

1. Hanks S, *et al.* Constitutional aneuploidy and cancer predisposition caused by biallelic mutations in *BUB1B*. *Nat Genet* 2004;36:1159–1161. doi:10.1038/ng1449 — record only.
2. Suijkerbuijk SJE, *et al.* Molecular causes for BUBR1 dysfunction in the human cancer predisposition syndrome mosaic variegated aneuploidy. *Cancer Res* 2010;70:4891–4900. doi:10.1158/0008-5472.CAN-09-4319 — **full text read.**
3. Matsuura S, *et al.* Monoallelic *BUB1B* mutations and defective mitotic-spindle checkpoint in seven families with premature chromatid separation syndrome. *Am J Med Genet A* 2006;140A:358–367. doi:10.1002/ajmg.a.31069 — record only.
4. Ochiai H, *et al.* TALEN-mediated single-base-pair editing identification of an intergenic mutation upstream of *BUB1B* as causative of PCS (MVA) syndrome. *PNAS* 2014;111:1461–1466. doi:10.1073/pnas.1317008111 — record only.
5. García-Castillo H, *et al.* Clinical and genetic heterogeneity in patients with mosaic variegated aneuploidy. *Am J Med Genet A* 2008;146A:1687–1695. doi:10.1002/ajmg.a.32315 — record only.
6. Narendra A, *et al.* Nephrocalcinosis in preterm babies. *Arch Dis Child Fetal Neonatal Ed* 2001;85:F207–F213. doi:10.1136/fn.85.3.F207 — record only.
7. Yost S, *et al.* Biallelic *TRIP13* mutations predispose to Wilms tumor and chromosome missegregation. *Nat Genet* 2017;49:1148–1151. doi:10.1038/ng.3883 — record only.
8. Richards S, *et al.* Standards and guidelines for the interpretation of sequence variants (ACMG/AMP). *Genet Med* 2015;17:405–424. doi:10.1038/gim.2015.30 — record only. Cited to state what this work did **not** do: no ACMG criteria were applied.
9. McLaren W, *et al.* The Ensembl Variant Effect Predictor. *Genome Biol* 2016;17:122. doi:10.1186/s13059-016-0974-4 — record only.
10. Rentzsch P, *et al.* CADD. *Nucleic Acids Res* 2019;47:D886–D894. doi:10.1093/nar/gky1016 — record only.
11. Cheng J, *et al.* Accurate proteome-wide missense variant effect prediction with AlphaMissense. *Science* 2023;381:eadg7492. doi:10.1126/science.adg7492 — record only.
12. Jaganathan K, *et al.* Predicting splicing from primary sequence with deep learning. *Cell* 2019;176:535–548.e24. doi:10.1016/j.cell.2018.12.015 — record only.
13. Landrum MJ, *et al.* ClinVar. *Nucleic Acids Res* 2018;46:D1062–D1067. doi:10.1093/nar/gkx1153 — record only. Individual records were retrieved directly: VCV000533901.9, VCV004600147.1, VCV000006765.53, VCV004842697.
14. Chen S, *et al.* A genomic mutational constraint map using variation in 76,156 human genomes. *Nature* 2024;625:92–100. doi:10.1038/s41586-023-06045-0 — record only.
15. Gargano MA, *et al.* The Human Phenotype Ontology in 2024. *Nucleic Acids Res* 2024;52:D1333–D1346. doi:10.1093/nar/gkad1005 — record only.
16. Van der Auwera GA, O'Connor BD. *Genomics in the Cloud*. O'Reilly, 2020 — GATK, the caller behind the supplied VCF. Not re-run here.
17. Martin M, *et al.* WhatsHap: fast and accurate read-based phasing. *bioRxiv* 2016. doi:10.1101/085050 — record only; `whatshap 2.8` was run.
18. Li H. Minimap2. *Bioinformatics* 2018;34:3094–3100. doi:10.1093/bioinformatics/bty191 — record only.
19. Danecek P, *et al.* Twelve years of SAMtools and BCFtools. *GigaScience* 2021;10:giab008. doi:10.1093/gigascience/giab008 — record only.
20. Martin AR, *et al.* PanelApp. *Nat Genet* 2019;51:1560–1565. doi:10.1038/s41588-019-0528-2 — record only. No tool comparison was performed.
21. Stenton SL, *et al.* Critical assessment of variant prioritization methods for rare disease diagnosis within the Rare Genomes Project. *Hum Genomics* 2024;18:44. doi:10.1186/s40246-024-00604-w — **full text read.**

Software and database **versions** are tabulated separately in
`docs/METHODS_FINAL.md` §6: a citation names a method, a version names what was
actually run.

---

*Research output for a hackathon. Not a clinical report and not a diagnosis.
Any finding here would require confirmation on an accredited platform and
interpretation by the treating clinical team before it could bear on care.*
