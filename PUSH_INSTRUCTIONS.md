# PUSH_INSTRUCTIONS.md — getting this repository onto GitHub

I could not reach `github.com` from my environment (three hosts, 25 s timeouts
each, while other allowlisted domains returned HTTP 200), so I could not check
your account or push for you. Everything below runs on **your** machine.

Two files are provided and they contain the same repository:

| file | use it if |
|---|---|
| `mva-hackathon-2026-track1.bundle` | you want the git history too — clone from it in one command |
| `mva-hackathon-2026-track1-repo.zip` | you prefer plain files, or want to drag-and-drop into the GitHub web UI |

Verified before shipping: the bundle records a complete history, clones cleanly
to 33 files at commit `c2fe1c9` on branch `main`, all 11 scripts parse, and the
prediction CSV inside the clone has md5 `d1b74e84be505744149e1cab1877d27a`.

---

## Route A — from the bundle (recommended)

```bash
# 1. clone the bundle into a working directory
git clone mva-hackathon-2026-track1.bundle mva-hackathon-2026-track1
cd mva-hackathon-2026-track1

# 2. drop the bundle as the origin and point at GitHub
git remote remove origin
git remote add origin https://github.com/hanaeroussotte/mva-hackathon-2026-track1.git

# 3. push
git push -u origin main
```

**If step 3 says the repository does not exist**, create it first — either on
github.com → New repository (name it `mva-hackathon-2026-track1`, **Private**,
do *not* add a README or .gitignore), or from the command line if you have the
`gh` CLI:

```bash
gh repo create hanaeroussotte/mva-hackathon-2026-track1 --private --source=. --push
```

**If step 3 says a repository already exists with different history**, the
earlier push did land and you are looking at a naming or visibility problem
rather than a missing repo. Check
`https://github.com/hanaeroussotte?tab=repositories` while logged in — a private
repository returns 404 to anyone not signed in, which would explain a dead link.
In that case either push with `--force` to replace it with this clean version,
or create a new repository under a different name and use that URL on the form.

## Route B — from the ZIP

```bash
unzip mva-hackathon-2026-track1-repo.zip
cd mva-hackathon-2026-track1
git init -b main
git add -A
git commit -m "Track 1: causal variant prediction for proband WGS_EX2312012"
git remote add origin https://github.com/hanaeroussotte/mva-hackathon-2026-track1.git
git push -u origin main
```

## Authentication

GitHub stopped accepting account passwords for git operations. Use a **personal
access token** as the password when prompted (Settings → Developer settings →
Personal access tokens → Fine-grained or Classic, scope `repo`), or sign in
through the GitHub CLI with `gh auth login`.

## Commit authorship

The commit is attributed to `hanaeroussotte <hanaeroussotte@users.noreply.github.com>`.
That is GitHub's own no-reply address form; I did not use a personal email
address. If you want your own identity on it:

```bash
git -c user.name="Your Name" -c user.email="you@example.com" commit --amend --reset-author --no-edit
```

## Before you make it public

The hackathon rules require the repository to be **public once the hackathon
ends** so the judging panel can review code and methods. It is currently set up
to be pushed private, which the rules permit during the event.

This repository contains **no proband sequencing data** and no genotype-scale
derived tables — see `docs/REPO_SCOPE.md` for exactly what is held back and the
one question that remains open with the organisers. A `.gitignore` blocks
`data/`, `*.vcf*`, `*.fastq*`, `*.bam`, `*.cram` and `*.docx` so proband data
cannot be committed by accident.

## After pushing

1. Open the repository URL in a **private/incognito window**. Logged out, a
   private repo shows 404 and a public repo renders — that is the quickest way
   to confirm what a judge will actually see.
2. Paste that exact URL into the submission form's GitHub field. The form
   requires it to start with `https://github.com/`.
3. Submit `hanaeroussotte_track1_report.md` as the report and
   `hanaeroussotte_mechanism-gated_FINAL.csv` as the predictions.
