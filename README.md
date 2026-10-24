# Track 1 — causal variant prediction for proband WGS_EX2312012

Rare Disease, Real Kid: The MVA Hackathon 2026 · submitted by `hanaeroussotte`

**Prediction:** two variants in *BUB1B*, proposed as a compound-heterozygous pair
with **unresolved phase**.

```
chr15:40209701 T>G   NM_001211.6:c.2210T>G   p.Leu737Ter   stop_gained, exon 17/23
chr15:40220612 T>G                c.3006T>G   p.Asn1002Lys  missense,   exon 23/23
```

Implied disease: Mosaic variegated aneuploidy syndrome 1 (OMIM:257300). Rank 1
of 448 hypotheses, score 0.6456, 1.64× the runner-up.

**This is a prioritised hypothesis for follow-up, not a diagnosis.** Two
requirements are unresolved: whether the variants lie on opposite copies of the
gene, and whether p.Asn1002Lys damages the protein. One supplied clinical
feature, nephrocalcinosis, is unexplained.

**Not submitted to the leaderboard. Nothing published.**

---

## Start here

| | |
|---|---|
| **The prediction** | [`submissions/hanaeroussotte_mechanism-gated_FINAL.csv`](submissions/hanaeroussotte_mechanism-gated_FINAL.csv) — md5 `d1b74e84be505744149e1cab1877d27a` |
| **The methods** | [`hanaeroussotte_track1_report.md`](hanaeroussotte_track1_report.md) — **this is the submitted report** |
| **How scoring works, exactly** | [`docs/SCORING_SPECIFICATION.md`](docs/SCORING_SPECIFICATION.md) — every term, penalty and transformation, with a worked *BUB1B* example |
| **What the validation shows** | [`docs/BENCHMARK_EXPLAINED.md`](docs/BENCHMARK_EXPLAINED.md) — read before quoting any benchmark figure |
| **How to reproduce** | §Reproducing below |
| **Submission requirements** | [`docs/SUBMISSION_CHECKLIST.md`](docs/SUBMISSION_CHECKLIST.md) |
| **What is and is not in this repository** | [`docs/REPO_SCOPE.md`](docs/REPO_SCOPE.md) |

---

## Method in one paragraph

The pipeline ranks **allele configurations**, not genes. A gene is not a
hypothesis: *BUB1B*-related MVA1 is recessive, so one damaged copy is a carrier
state. Each gene's candidate variants are grouped into `hom`, `comphet`,
`het_single` or `incomplete_recessive` configurations, each admitted only
against diseases whose inheritance annotation is compatible with it, and then
scored against the child's eight HPO terms weighted by information content. What
this implements is disease-level phenotype matching gated by compatibility with
the observed allele configuration — not a mechanistic model of disease biology.

```
5,012,204 VCF records -> 4,740,790 PASS -> 180,564 exonic ±10bp -> 29,907 coding
  -> 583 clear the evidence floor -> 562 distinct variants (861 variant/model rows)
  -> 448 genes -> 448 ranked hypotheses (one per gene)
```

## Validation

80 constructed cases (seed 20261007) planting variants of known identity into
this proband's own callset, with the planted variants removed from the ClinVar
evidence file.

**Complete-truth: exact configuration at rank 1 in 37 of 45 (82%, Wilson 95% CI
69–91%); correct gene in 39 of 45.** These are different endpoints and are not
interchangeable. **Neither is a probability that this child's prediction is
correct.**

Three findings a reviewer should have before quoting that figure: a
**phenotype-only ablation of the same pipeline matches it** (37/45, zero
discordant pairs, p = 1.000); **80 of 80** injected genes contained zero
competing alleles, so partner selection was never tested; and phenotypes were
sampled from the planted disease's own annotation.

## Reproducing

Three different things are called reproduction, and only the first runs from
this repository.

**1. Ranking rerun** — needs `work/candidates_pre_af.tsv`, which is held back
pending the clarification in `docs/REPO_SCOPE.md`:

```bash
python src/05_rank_and_submit.py \
    work/candidates_pre_af.tsv work/bench_af.json work/mech_scores.tsv \
    WGS_EX2312012 submissions/repro.csv work/ranked_full_repro.tsv
# -> genes_ranked 448 / BUB1B comphet 0.6456 rank 1 / LZTR1 0.3933 rank 2
```

**2. Cached-input reproduction** — additionally needs `work/cds_vep.tsv` and
`work/cds_genotypes.tsv`, both withheld as individual-level genotype data.

**3. Raw-data reproduction** — **not performed.** Nine of eleven scripts were not
run; `08_kmer_validate.py` cannot be run, as the FASTQ is deleted under the
dataset terms.

`environment.yml` pins the versions used (Python 3.11.16, bcftools 1.24). **A
fresh `conda env create` from it was not tested.**

## Known defects, disclosed and unfixed

Fixing either would change the candidate set after the results were seen.

1. **Transcript selection** takes the first MANE-tagged transcript, so at
   overlapping loci a neighbouring gene's transcript can displace a coding
   consequence in the intended gene. 275 of 29,907 coding calls (0.92%) are
   candidate instances. **The three submitted alleles are verified unaffected.**
2. **`GT = 1/2` genotypes** clear the evidence floor but produce no model row.
   Six variants; none in a gene relevant to the prediction.

## Data

No raw sequencing data in this repository, and none will be added. See
`docs/REPO_SCOPE.md` for what is held back and the exact question outstanding
with the organisers.

---

*Research output for a hackathon. Not a clinical report and not a diagnosis.*
