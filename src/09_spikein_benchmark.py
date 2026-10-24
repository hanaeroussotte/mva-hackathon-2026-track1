"""Spike-in benchmark: where does the pipeline rank a known causal genotype?

The pipeline was developed on one case whose answer was already believed, so
its ranking is not evidence that the scoring generalises. This benchmark
supplies the missing control.

Construction of one synthetic case:
  * pick a gene G with an autosomal-recessive disease D annotated with >=8 HPO
    terms, and two distinct ClinVar Pathogenic/Likely-pathogenic SNVs in G
  * synthetic phenotype = 8 terms sampled from D's annotations
  * synthetic genotype = the two variants injected as heterozygous, with depth
    and quality fields sampled from the proband's own heterozygous CDS calls
  * background distractors = the proband's real candidate variants, unchanged

The pipeline then runs end to end, blind to which gene was spiked, and the rank
of G is recorded.

Two arms, because the scoring awards up to +3.0 for a ClinVar P/LP annotation
and the spiked variants were selected for having one:
  visible  ClinVar annotation present  (optimistic: the answer is pre-labelled)
  masked   ClinVar annotation withheld (the pipeline must reach it from
           consequence, CADD, AlphaMissense, SpliceAI and rarity alone)

Usage:
    python src/09_spikein_benchmark.py handoff/spikein_cases.json work/spikein_results.tsv
"""
import json, os, subprocess, sys, tempfile
from importlib import import_module

sys.path.insert(0, "src")
_r = import_module("05_rank_and_submit")

REAL_VEP = "work/cds_vep.tsv"
REAL_GT = "work/cds_genotypes.tsv"
REAL_CLINVAR = "work/clinvar_hits.tsv"
SPIKE_VEP = os.environ.get("SPIKE_VEP", "work/spike_vep.tsv")


def load_spike_vep():
    rows = {}
    with open(SPIKE_VEP) as fh:
        header = fh.readline()
        for line in fh:
            f = line.rstrip("\n").split("\t")
            rows[(f[0], f[1], f[2], f[3])] = line
    return header, rows


def run_case(case, spike_header, spike_rows, clinvar_extra, af_path, tmp, arm):
    keys = [(v["chrom"], str(v["pos"]), v["ref"], v["alt"]) for v in case["variants"]]
    if not all(k in spike_rows for k in keys):
        return None

    vep_p = os.path.join(tmp, "vep.tsv")
    with open(vep_p, "w") as out:
        out.write(open(REAL_VEP).read())
        for k in keys:
            out.write(spike_rows[k])

    gt_p = os.path.join(tmp, "gt.tsv")
    with open(gt_p, "w") as out:
        out.write(open(REAL_GT).read())
        for v in case["variants"]:
            out.write(f"{v['chrom']}\t{v['pos']}\t{v['ref']}\t{v['alt']}\t500.0\t"
                      f"{v['gt']}\t{v['ad']}\t{v['dp']}\t{v['gq']}\n")

    cv_p = os.path.join(tmp, "clinvar.tsv")
    with open(cv_p, "w") as out:
        out.write(open(REAL_CLINVAR).read())
        if arm == "visible":
            for k in keys:
                if k in clinvar_extra:
                    out.write(clinvar_extra[k])

    hpo_p = os.path.join(tmp, "hpo.txt")
    with open(hpo_p, "w") as out:
        out.write("\n".join(case["hpo"]) + "\n")

    cand_p = os.path.join(tmp, "cand.tsv")
    mech_p = os.path.join(tmp, "mech.tsv")
    env = dict(os.environ)
    subprocess.run([sys.executable, "src/04_shortlist.py", vep_p, gt_p, cv_p,
                    cand_p, os.path.join(tmp, "genes.txt")],
                   check=True, capture_output=True, env=env)
    subprocess.run([sys.executable, "src/07_mechanism_score.py", hpo_p, mech_p],
                   check=True, capture_output=True, env=env)
    # Per-case phase states. GATK PGT/PID for the real background calls comes
    # from the real callset; the spiked variants contribute their own tags via
    # work/bench_phase_tags.tsv, written by the case builder.
    phase_p = os.path.join(tmp, "phase.tsv")
    tag_p = os.path.join(tmp, "tags.tsv")
    with open(tag_p, "w") as out:
        out.write(open("work/cds_phase.tsv").read())
        for v in case["variants"]:
            out.write(f"{v['chrom']}\t{v['pos']}\t{v['ref']}\t{v['alt']}\t{v['gt']}\t"
                      f"{v.get('pgt','.')}\t{v.get('pid','.')}\n")
    if os.environ.get("USE_PHASE", "1") == "1":
        subprocess.run([sys.executable, "src/10_phase_states.py", cand_p, tag_p, phase_p,
                        "/dev/null"],
                       check=True, capture_output=True, env=env)
        env["PHASE_FILE"] = phase_p
        _r.PHASE = _r.load_phase(phase_p)
    else:
        _r.PHASE = {}

    sweep = os.environ.get("VAR_SCALE_SWEEP")
    if sweep:
        out = []
        for vs in (float(x) for x in sweep.split(",")):
            _r.VAR_SCALE = vs
            sc = _r.rank(cand_p, af_path, mech_p)
            rk = next((i for i, s in enumerate(sc, 1) if s["gene"] == case["gene"]), None)
            out.append((vs, rk))
        _r.VAR_SCALE = 2.5
        return {"truth_variants": "", "n_truth": "", "hit_variants": "", "pred_top10": "",
"gene": case["gene"], "stratum": case.get("stratum",""),
                "disease": case["disease"], "arm": arm,
                "n_ranked": len(sc), "rank": ";".join(f"{v}:{r}" for v, r in out),
                "top_gene": "", "score": "", "top_score": "", "model": "",
                "coverage": "", "var_score": "", "implied_disease": ""}

    scored = _r.rank(cand_p, af_path, mech_p)
    rank = next((i for i, s in enumerate(scored, 1) if s["gene"] == case["gene"]), None)
    hit = next((s for s in scored if s["gene"] == case["gene"]), None)

    # --- variant-level capture -------------------------------------------
    # This benchmark previously recorded only whether the right GENE ranked
    # first. The challenge scores VARIANTS, so the injected truth and the exact
    # predicted alleles of the top rows are retained, letting exact-pair
    # recovery be measured instead of inferred from gene rank.
    def _vk(r):
        return f"{r['chrom']}:{r['pos']}:{r['ref']}:{r['alt']}"
    _truth = "+".join(sorted(f"{v['chrom']}:{v['pos']}:{v['ref']}:{v['alt']}"
                             for v in case["variants"]))
    _pred = ";".join("{}|{}|{}".format(x["gene"], x["model"],
                     "+".join(sorted(_vk(v) for v in x["variants"])))
                     for x in scored[:10])
    # ---------------------------------------------------------------------

    return {
        "truth_variants": _truth, "n_truth": len(case["variants"]),
        "hit_variants": ("+".join(sorted(_vk(v) for v in hit["variants"])) if hit else ""),
        "pred_top10": _pred,
        "gene": case["gene"], "stratum": case.get("stratum",""),
        "disease": case["disease"], "arm": arm,
        "n_ranked": len(scored),
        "rank": rank if rank else "", "top_gene": scored[0]["gene"] if scored else "",
        "score": f"{hit['total']:.4f}" if hit else "",
        "top_score": f"{scored[0]['total']:.4f}" if scored else "",
        "model": hit["model"] if hit else "",
        "coverage": f"{hit['pheno']:.4f}" if hit else "",
        "var_score": f"{hit['var_score']:.2f}" if hit else "",
        "implied_disease": hit["disease"] if hit else "",
    }


def main(cases_path, out_path, af_path="work/spikein_af.json"):
    cases = json.load(open(cases_path))
    spike_header, spike_rows = load_spike_vep()
    clinvar_extra = json.load(open(os.environ.get("SPIKE_CLINVAR", "handoff/spikein_clinvar.json")))
    clinvar_extra = {tuple(k.split("|")): v for k, v in clinvar_extra.items()}

    cols = ["gene","stratum","disease","arm","n_ranked","rank","top_gene","score","top_score",
            "model","coverage","var_score","implied_disease",
            "truth_variants","n_truth","hit_variants","pred_top10"]
    with open(out_path, "w") as out:
        out.write("\t".join(cols) + "\n")
        for i, case in enumerate(cases, 1):
            with tempfile.TemporaryDirectory() as tmp:
                for arm in os.environ.get("ARMS", "visible,masked").split(","):
                    r = run_case(case, spike_header, spike_rows, clinvar_extra,
                                 af_path, tmp, arm)
                    if r:
                        out.write("\t".join(str(r[c]) for c in cols) + "\n")
                        out.flush()
            if i % 5 == 0:
                print(f"progress\t{i}/{len(cases)}", flush=True)


if __name__ == "__main__":
    main(*sys.argv[1:4])
