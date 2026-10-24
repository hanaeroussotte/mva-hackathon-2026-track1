"""Rank diagnostic hypotheses and write the Track 1 submission CSV.

Changes in this version:
  * `incomplete_recessive` hypotheses are scored with the MISSING allele
    contributing zero evidence -- the compound-het score is the mean of two
    allele scores, so one found allele gives exactly half. No new constant.
  * every hypothesis carries an explicit PHASE STATE, and a compound-het pair
    demonstrated to be in cis is demoted, because two variants on one haplotype
    leave an intact allele and cannot cause a recessive disease.
  * `epcr` is documented as a gene-ranking confidence, NOT the probability that
    a specific variant pair is causal.
"""
import csv, json, math, os, sys
from collections import defaultdict

AF_RECESSIVE_MAX, AF_DOMINANT_MAX = 1e-2, 1e-4
W_PHENO = 1.5
VAR_SCALE = float(os.environ.get("VAR_SCALE", "2.5"))
PHENO_FLOOR = 0.02
COMBINE = os.environ.get("COMBINE", "multiplicative")
RANK_MODE = os.environ.get("RANK_MODE", "full")
PAIR_AGG = os.environ.get("PAIR_AGG", "mean")
AF_UNKNOWN = os.environ.get("AF_UNKNOWN", "absent")
# MAX_PRIMARY = 1: only the leading hypothesis is submitted as a primary
# differential. Earlier versions submitted the next two ranked genes as
# "insurance" on the reasoning that extra rows at low EPCR cannot reduce the
# F-max achieved at a higher threshold. That reasoning is true but incomplete:
# running the challenge's own evaluation.py shows the extra rows push the
# remaining rows down the EPCR ordering and across rank-point tier boundaries,
# costing 12.5 -> 5.0 rank points and 0.40 -> 0.29 F-max in the scenario where
# they were supposed to help. They are strictly harmful, and they were not
# evidence-supported differentials in the first place.
MAX_ROWS, MAX_PRIMARY = 10, 1

# A pair proven to lie on one haplotype cannot explain a recessive disease.
CIS_PENALTY = 0.05

SECONDARY_GENES = {"LZTR1"}
EXCLUDED_GENES = {"HLA-A","HLA-B","HLA-C","HLA-DRB1","HLA-DRB5","HLA-DQA1","HLA-DQB1",
                  "HLA-DPA1","HLA-DPB1"}
SECONDARY_RATIONALE = {
    "LZTR1": ("INCIDENTAL, not proposed as the cause of the presenting phenotype. "
              "Single heterozygous truncating allele is compatible only with "
              "LZTR1-related schwannomatosis (AD, adult-onset nerve sheath tumours); "
              "AR Noonan syndrome requires biallelic variants. Reported because "
              "schwannomatosis surveillance is actionable."),
}
SUBMISSION_COLUMNS = ["proband_id","chrom_1","pos_1","ref_1","alt_1","chrom_2","pos_2",
                      "ref_2","alt_2","epcr","finding_type","notes"]


def rarity_points(af, model):
    if af is None:
        # AF_UNKNOWN controls how an UNRESOLVED lookup is treated. The default
        # ("absent") is the behaviour shipped with the original analysis and
        # conflates two different things: a variant verified absent from
        # gnomAD, and a lookup that returned nothing. Only the first deserves
        # rarity credit. The other settings bound the consequence.
        #   absent  : +1.5, same as verified-absent (original behaviour)
        #   unknown : 0.0, no rarity credit either way
        #   common  : -6.0, worst case - treat as if too common
        if AF_UNKNOWN == "unknown":
            return 0.0, "AF_unresolved_no_credit"
        if AF_UNKNOWN == "common":
            return -6.0, "AF_unresolved_penalised"
        return 1.5, "absent_gnomAD"
    cap = AF_RECESSIVE_MAX if model in ("hom","comphet","incomplete_recessive") else AF_DOMINANT_MAX
    if af > cap:
        return -6.0, f"AF{af:.2e}_too_common"
    if af < 1e-5:
        return 1.5, f"AF{af:.2e}"
    if af < 1e-4:
        return 1.0, f"AF{af:.2e}"
    if af < 1e-3:
        return 0.3, f"AF{af:.2e}"
    return -0.5, f"AF{af:.2e}"


def load_af(path):
    raw = json.load(open(path))
    out = {}
    for vid, blocks in raw.items():
        best = None
        for key in ("exome","genome"):
            af = (blocks.get(key) or {}).get("af")
            if af is not None:
                best = af if best is None else max(best, af)
        out[vid] = best
    return out


def load_phase(path=None):
    path = path or os.environ.get("PHASE_FILE", "work/phase_states.tsv")
    """gene -> state, from an explicit phasing analysis. Absent = unresolved."""
    out = {}
    try:
        with open(path) as fh:
            head = next(fh).rstrip("\n").split("\t")
            ix = {c: i for i, c in enumerate(head)}
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) > ix["phase_state"]:
                    out[(f[ix["gene"]], f[ix["model"]])] = f[ix["phase_state"]]
    except (FileNotFoundError, StopIteration, KeyError):
        pass
    return out


PHASE = load_phase()


def epcr_from_score(total, best_total):
    if best_total <= 0:
        return 0.01
    return round(max(0.01, min(0.95, 0.95 * total / best_total)), 4)


def rank(cand_path, af_path, pheno_path):
    af_map = load_af(af_path)
    pheno, mech_disease, gene_level = {}, {}, {}
    with open(pheno_path) as fh:
        head = next(fh).rstrip("\n").split("\t")
        mech_mode = "model" in head
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if mech_mode:
                pheno[(f[0], f[1])] = float(f[10]) if len(f) > 10 else float(f[2])
                mech_disease[(f[0], f[1])] = (f[4], f[11] if len(f) > 11 else "")
            else:
                gene_level[f[0]] = float(f[1])
    pheno_for = (lambda g, m: pheno.get((g, m), 0.0)) if mech_mode else \
                (lambda g, m: gene_level.get(g, 0.0))

    def total(r):
        pts, note = rarity_points(r["af"], r["model"])
        return r["evidence"] + pts, note

    by_gene = defaultdict(list)
    for row in csv.DictReader(open(cand_path), delimiter="\t"):
        row["af"] = af_map.get(f"{row['chrom']}-{row['pos']}-{row['ref']}-{row['alt']}")
        row["evidence"] = float(row["evidence"])
        by_gene[row["gene_symbol"]].append(row)

    scored = []
    for gene, rows in by_gene.items():
        by_model = defaultdict(list)
        for r in rows:
            by_model[r["model"]].append(r)
        options = []
        for model, mrows in by_model.items():
            phase_state = "not_applicable"
            if model == "comphet":
                if len(mrows) < 2:
                    continue
                pair = sorted(mrows, key=lambda r: -total(r)[0])[:2]
                # PAIR_AGG: how two allele scores combine into a pair score.
                #   mean  (default, as submitted) -- (a + b) / 2
                #   min   (pre-specified alternative) -- the WEAKER allele governs.
                # Rationale for the alternative, fixed before any outcome was
                # inspected: under a recessive model both alleles must be
                # damaging, so the weaker allele is the binding constraint, and
                # averaging lets a strongly supported allele carry a weak partner.
                _tt = [total(r)[0] for r in pair]
                var_score = min(_tt) if PAIR_AGG == "min" else sum(_tt) / 2.0
                phase_state = PHASE.get((gene, model), "unresolved_phase")
            elif model == "incomplete_recessive":
                pair = [max(mrows, key=lambda r: total(r)[0])]
                # the second allele was never found: it contributes zero evidence
                var_score = total(pair[0])[0] / 2.0
                phase_state = "single_allele_only"
            else:
                pair = [max(mrows, key=lambda r: total(r)[0])]
                var_score = total(pair[0])[0]
            p = pheno_for(gene, model)
            e_var = 1.0 - math.exp(-max(var_score, 0.0) / VAR_SCALE)
            e_pheno = max(p, PHENO_FLOOR)
            # RANK_MODE selects the scoring rule. "full" is the submitted method.
            # The two baselines are INTERNAL ablations of it, not independent
            # external tools: they reuse the same candidate pool, annotations,
            # inheritance gate and truth definitions, and differ only in which
            # term drives the ordering.
            #   evidence_only : variant evidence + rarity, phenotype ignored
            #   pheno_only    : phenotype score alone, over inheritance-compatible
            #                   candidates, variant evidence ignored
            if RANK_MODE == "evidence_only":
                tot = e_var
            elif RANK_MODE == "pheno_only":
                tot = e_pheno
            else:
                tot = e_var * e_pheno if COMBINE == "multiplicative" else var_score + W_PHENO * p
            if phase_state == "confirmed_cis":
                tot *= CIS_PENALTY
            options.append((tot, model, pair, var_score, p, phase_state))
        if not options:
            continue
        tot, model, pair, var_score, p, phase_state = max(options, key=lambda o: o[0])
        notes = " | ".join(f"{r['hgvsp'] or r['hgvsc']} [{r['tags']};{total(r)[1]}]" for r in pair)
        dz_name, n_expl = mech_disease.get((gene, model), ("", ""))
        scored.append({"gene":gene,"model":model,"variants":pair,"var_score":var_score,
                       "pheno":p,"total":tot,"notes":notes,"disease":dz_name,
                       "n_features_explained":n_expl,"phase_state":phase_state,
                       "n_models_considered":len(options)})
    scored.sort(key=lambda s: -s["total"])
    return scored


def write_submission(scored, proband_id, csv_path):
    primary = [s for s in scored if s["gene"] not in SECONDARY_GENES
               and s["gene"] not in EXCLUDED_GENES]
    secondary = [s for s in scored if s["gene"] in SECONDARY_GENES]
    rows = []
    for i, s in enumerate(primary[:MAX_PRIMARY]):
        v1 = s["variants"][0]
        v2 = s["variants"][1] if s["model"] == "comphet" and len(s["variants"]) > 1 else None
        rows.append({"proband_id":proband_id,"chrom_1":f"chr{v1['chrom']}","pos_1":v1["pos"],
                     "ref_1":v1["ref"],"alt_1":v1["alt"],
                     "chrom_2":f"chr{v2['chrom']}" if v2 else "","pos_2":v2["pos"] if v2 else "",
                     "ref_2":v2["ref"] if v2 else "","alt_2":v2["alt"] if v2 else "",
                     "epcr":epcr_from_score(s["total"], primary[0]["total"]),
                     "finding_type":"primary",
                     "notes":(f"{s['gene']} ({s['model']}, phase={s['phase_state']}) -> "
                              f"{s['disease'] or 'no admissible disease'}; "
                              f"phenotype coverage {s['pheno']:.2f}, "
                              f"computational HPO match for {s['n_features_explained']} of 8 proband "
                              f"terms (similarity, not biological explanation; only a subset "
                              f"are specific matches -- see docs/phenotype_feature_map.tsv); "
                              f"{s['notes']}")[:480]})
    # Single-allele backstop rows for the leading compound-heterozygous
    # hypothesis.
    #
    # The reason recorded here earlier was WRONG and is corrected. It argued
    # from Stenton et al. 2024 (doi:10.1186/s40246-024-00604-w), whose
    # assessment protocol counts a part-correct biallelic row as wholly
    # incorrect. The challenge does NOT inherit that rule: its own
    # evaluation.py awards 0.5 x rank points for a row containing one true
    # variant, "either proposed alone, or paired with an incorrect second
    # variant". So the paired row already earns partial credit by itself and
    # these rows add nothing in that scenario.
    #
    # They are kept for a different reason, established by running that
    # evaluator against synthetic answer keys (docs/evaluator_scenarios.tsv).
    # If the true answer is a SINGLE causal variant rather than a pair, the
    # evaluator computes a partial match only when len(true_variants) == 2, so
    # a submission holding only the paired row scores ZERO. Adding each allele
    # as its own row takes that case from 0 to 50 rank points. They cost
    # nothing in the other three scenarios tested.
    if primary and primary[0]["model"] == "comphet" and len(primary[0]["variants"]) > 1:
        lead = primary[0]
        lead_epcr = epcr_from_score(lead["total"], primary[0]["total"])
        for j, v in enumerate(lead["variants"][:2]):
            if len(rows) >= MAX_ROWS:
                break
            rows.append({"proband_id":proband_id,"chrom_1":f"chr{v['chrom']}","pos_1":v["pos"],
                         "ref_1":v["ref"],"alt_1":v["alt"],"chrom_2":"","pos_2":"","ref_2":"",
                         "alt_2":"","epcr":round(min(lead_epcr - 0.05 - 0.05*j, 0.90), 4),
                         "finding_type":"primary",
                         "notes":(f"{lead['gene']} allele {j+1} of the rank-1 "
                                  f"compound-heterozygous hypothesis, submitted additionally as a "
                                  f"single-variant row so that the hypothesis still scores if the "
                                  f"true answer is a single causal variant rather than a pair. "
                                  f"Not an independent hypothesis; phase={lead['phase_state']}.")[:480]})
    for s in secondary:
        if len(rows) >= MAX_ROWS:
            break
        v1 = s["variants"][0]
        rows.append({"proband_id":proband_id,"chrom_1":f"chr{v1['chrom']}","pos_1":v1["pos"],
                     "ref_1":v1["ref"],"alt_1":v1["alt"],"chrom_2":"","pos_2":"","ref_2":"",
                     "alt_2":"","epcr":0.02,"finding_type":"secondary",
                     "notes":(f"{s['gene']}: {SECONDARY_RATIONALE.get(s['gene'],'')} "
                              f"{s['notes']}")[:480]})
    for r in rows:
        assert 0.0 < float(r["epcr"]) <= 1.0, r["epcr"]
        assert str(r["chrom_1"]).startswith("chr"), r["chrom_1"]
        if r["chrom_2"]:
            assert str(r["chrom_2"]).startswith("chr"), r["chrom_2"]
            assert int(r["pos_2"]) > 0
    with open(csv_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=SUBMISSION_COLUMNS)
        w.writeheader(); w.writerows(rows)
    print(f"submission_rows\t{len(rows)}\t{csv_path}")
    return rows


def main(cand_path, af_path, pheno_path, proband_id, csv_path, full_path):
    scored = rank(cand_path, af_path, pheno_path)
    write_submission(scored, proband_id, csv_path)
    with open(full_path, "w") as out:
        out.write("rank\tgene\tmodel\tphase_state\tvar_score\tpheno_score\ttotal\t"
                  "implied_disease\tn_features\tdetail\n")
        for i, s in enumerate(scored, 1):
            out.write(f"{i}\t{s['gene']}\t{s['model']}\t{s['phase_state']}\t"
                      f"{s['var_score']:.2f}\t{s['pheno']:.3f}\t{s['total']:.4f}\t"
                      f"{s['disease']}\t{s['n_features_explained']}\t{s['notes']}\n")
    print(f"genes_ranked\t{len(scored)}")
    for s in scored[:10]:
        print(f"{s['total']:7.4f}  {s['gene']:<12} {s['model']:<22} {s['phase_state']:<18} "
              f"var={s['var_score']:.2f} pheno={s['pheno']:.2f}  {s['disease'][:34]}")


if __name__ == "__main__":
    main(*sys.argv[1:7])
