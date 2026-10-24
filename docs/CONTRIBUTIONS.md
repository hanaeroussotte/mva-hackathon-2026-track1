# CONTRIBUTIONS.md — draft for Hanae Roussotte to review and edit

This is a **draft**. The AI-assistance section is factual and can be filed as
written once you confirm the two bracketed items. **The human-contribution
section is deliberately left for you to complete** — it is not mine to write,
and I have not guessed at it.

---

## AI assistance disclosure

The challenge requires this (submission tab, update of 28 August 2026): *"If you
used an LLM or AI assistant, please record the provider, the plan or tier, and
the relevant data-handling setting in your methods description."*

**Provider and model:** Anthropic, Claude Opus 5, accessed through Claude
Science.

**Plan or tier:** `[TO CONFIRM — I cannot see your subscription tier. Fill in the
plan under which you accessed the assistant.]`

**Data-handling setting:** `[TO CONFIRM — state the applicable setting for your
plan, e.g. whether customer content is excluded from training. I cannot verify
this from inside the session and will not assert it.]`

A one-line form, once completed, satisfies the requirement — for example:
*"Anthropic API, Claude Opus 5, <plan>, <data-handling setting>."*

### What the AI assistance actually did

Stated plainly so a reviewer can calibrate it:

* Wrote and ran the analysis code — the eleven pipeline scripts, the benchmark
  generator and runner, the recovery classification, and the baseline and
  sensitivity arms.
* Queried external resources: Ensembl VEP, ClinVar (local release and NCBI
  esummary), gnomAD via VEP, Crossref for reference metadata, PubMed, and the
  challenge's own evaluator.
* Read one primary paper in full (Suijkerbuijk et al. 2010) and verified
  bibliographic metadata for the rest.
* Drafted the methods, README, specification and supporting documents.
* Conducted the audit that found and documented the two implementation defects,
  and retracted the claims listed in `docs/CHANGELOG.md`.

**Three sub-agents** were used for bounded tasks and their outputs are recorded
with their own declared deviations: a literature review, a Track 2 candidate
generation pass, and a simulation (the latter two belong to Track 2).

### Errors the AI made in this work, corrected and recorded

Listed because a reviewer should see them, and because their presence is more
informative than their absence:

* Misread the MVA acronym at the outset (mevalonic aciduria rather than mosaic
  variegated aneuploidy) and then ran a confirmation-seeking literature search
  on the wrong meaning.
* Drafted an ablation table before running it.
* Introduced a ClinVar column-order bug during a rebuild that silently cost
  allele 1 its ClinVar evidence term.
* Judged "absent from the ranking" from a field holding only the top ten rows,
  producing two false claims about the *cis* controls.
* Presented "zero wrong partners" as evidence about partner selection when no
  wrong partner had ever been available.
* Attributed the F5/COL17A1/RNF31 annotation mismatches to the benchmark before
  demonstrating they were a defect in its own transcript-selection code.
* Attributed a candidate-count difference to a missing file without checking
  that the explanation was even directionally possible.

Full list with corrections: `docs/CHANGELOG.md` (27 numbered items).

---

## Human contribution

`[FOR YOU TO COMPLETE. Suggested headings only — delete, rewrite or ignore.]`

* Scientific direction and the questions pursued:
* Domain judgement applied to intermediate results:
* Review of the methods, audit briefs and corrections directed:
* Decisions on what to submit and what to exclude:
* Writing:

What I can state factually, if it helps you write this: the project proceeded
through **eight successive review briefs** you authored, each of which changed
the analysis or the documents. Several of the corrections above were found
because you challenged a specific claim — the *cis*-control interpretation, the
"zero wrong partners" inference, the annotation-mismatch attribution, the
reference-verification language, and the data-sharing statement were each
reopened at your instruction, not spontaneously. The standing instruction never
to write a number that had not been computed and read back from actual output
came from you.

---

## Statement of what this work is not

No wet-lab validation. No independent cohort. No ACMG/AMP classification —
Richards et al. 2015 is cited only to state that the framework was **not**
applied. No clinical validation. No comparison against other prioritisation
tools. Nothing here is a diagnosis or a basis for clinical decisions.

---

## Data handling

Raw sequencing data are not redistributed. The FASTQ has been deleted under the
dataset terms. Per the gated agreement, Sage Bionetworks must be notified at
**cnb@sagebase.org** once all data copies are deleted — **the VCF is still
present locally and that notification has not been sent.** The Hugging Face
access token used to download the dataset should be revoked after submission.
