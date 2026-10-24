# SUBMISSION_CHECKLIST.md

Each official Track 1 requirement mapped to the exact file. Source for the
requirements: the challenge Space submission tab
(`SageBio/rare-disease-real-kid-mva-hackathon-2026`, `tabs/submit_track1.py`,
retrieved 2026-10-24) and the gated data-use agreement on
`SageBio/mva-hackathon-2026-data`.

**Nothing has been submitted and nothing has been published.**

## Required by the submission form

| # | Requirement | Satisfied by | Status |
|---|---|---|---|
| 1 | Predictions CSV, schema `proband_id, chrom_1, pos_1, ref_1, alt_1, chrom_2, pos_2, ref_2, alt_2, epcr, finding_type, notes` | `submissions/hanaeroussotte_mechanism-gated_FINAL.csv` | **READY** — accepted by the challenge's own `load_submission()`; md5 `d1b74e84be505744149e1cab1877d27a` |
| 2 | Filename includes username and a short approach name | `hanaeroussotte_mechanism-gated_FINAL.csv` | **READY** |
| 3 | `epcr` in (0, 1] | 0.95 / 0.90 / 0.85 / 0.02 | **READY** — range-checked by the loader |
| 4 | `finding_type` is `primary` or `secondary` | 3 primary, 1 secondary | **READY** |
| 5 | Report in PDF or Markdown | `hanaeroussotte_track1_report.md` | **READY** — Markdown is accepted |
| 6 | Report filename includes username | `hanaeroussotte_track1_report.md` | **READY** |
| 7 | AI-assistance line: provider, plan/tier, data-handling setting | report header + `docs/CONTRIBUTIONS.md` | **INCOMPLETE** — provider and model stated; **plan/tier and data-handling setting must be supplied by the submitter** |
| 8 | GitHub URL starting `https://github.com/` | `https://github.com/hanaeroussotte/mva-hackathon-2026-track1` | **UNVERIFIED IN THIS PASS** — see below |
| 9 | Display name | `hanaeroussotte` | **READY** |

## Required of the repository

| # | Requirement | Satisfied by | Status |
|---|---|---|---|
| 10 | Documented, reproducible code — all scripts and configuration needed to reproduce the results | `src/*.py` (11 scripts), `config/`, `environment.yml` | **READY** in the scoped package |
| 11 | Submission files *(optional)* | `submissions/` | included |
| 12 | Methods description *(optional but recommended)* | `hanaeroussotte_track1_report.md`, `docs/METHODS_FINAL.md` | included |
| 13 | May remain private during the hackathon; **must be made public once it ends** so the judging panel can review code and methods | — | **ACTION REQUIRED AFTER THE HACKATHON** |

## Required by the data-use agreement

| # | Requirement | Status |
|---|---|---|
| 14 | Do not release or grant data access to anyone | **MET** — no raw data in any package |
| 15 | Delete all data after completing the hackathon | **OUTSTANDING** — the VCF is still on local disk |
| 16 | Notify Sage Bionetworks at `cnb@sagebase.org` once deleted | **OUTSTANDING** |
| 17 | Do not attempt to recontact the data subject, family or MVA Society contacts | **MET** |
| 18 | Submissions are published CC-BY 4.0 with your name | acknowledged |

## Which document to submit

**Submit `hanaeroussotte_track1_report.md`.** It is the reader-facing report, in
the eight sections requested, and it carries the AI-disclosure line the rules
require.

`docs/METHODS_FINAL.md` is retained as the **technical supporting document** —
longer, with the full version table, the complete reference annotations and the
audit detail. It is referenced from the report but is not the submission.

## Genuinely unresolved

1. **AI-disclosure plan/tier and data-handling setting** (#7). I cannot see your
   subscription terms and will not assert them. One line from you closes this.
2. **Repository existence and visibility** (#8). `github.com` and
   `api.github.com` were **unreachable from this environment** in this pass —
   the local proxy accepts the connection and then times out after 25–28 s,
   while other allowlisted domains return 200. The repositories were created and
   pushed earlier in this project and were recorded then as **private**, but
   that is a **previously recorded** state, not a check executed now. Open the
   URL yourself, or push the scoped package from your own machine, before
   entering the link on the form.
3. **Which individual-level derived files may be published** (#13 vs #14). The
   rules require a public repository with code and methods; the data agreement
   prohibits releasing the Data and does not define whether genotype-scale
   derived outputs count. The scoped package resolves this conservatively by
   containing **code, configuration and methods only**. The affected categories
   and the exact question to ask the organisers are in
   `docs/REPO_SCOPE.md`.
4. **Data deletion and token revocation** (#15, #16). Both outstanding; both are
   yours to perform.
