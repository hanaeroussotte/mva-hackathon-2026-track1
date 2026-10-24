"""Apply inheritance models, quality terms, read-support exclusions, and a
recessive rescue route.

Change from the previous version (motivated by the failure mode the spike-in
benchmark exposed): a gene with exactly ONE heterozygous variant clearing the
evidence floor previously produced only a `het_single` hypothesis, which the
inheritance gate then correctly ruled inadmissible for an autosomal-recessive
disease -- so a true recessive gene whose second allele was missed upstream
collapsed to the phenotype floor instead of ranking. Such genes now ALSO emit an
`incomplete_recessive` hypothesis, which 07 scores against AR diseases and 05
scores with the missing allele contributing zero evidence. The route makes the
gene visible for follow-up without treating an incomplete hypothesis as a
confirmed biallelic one.
"""
import os, sys
from collections import defaultdict

QUALITY_TERMS = os.environ.get("QUALITY_TERMS", "1") != "0"
# "once" is the default after audit: AlphaMissense and CADD are correlated and
# must not both be counted for the same missense variant.
MISSENSE_EVIDENCE = os.environ.get("MISSENSE_EVIDENCE", "once")
# Evidence floor. Lower it only to build a FIXED candidate pool for scoring
# sensitivity analysis, where the shortlist must not change between arms.
EVIDENCE_FLOOR = float(os.environ.get("EVIDENCE_FLOOR", "1.0"))
RESCUE = os.environ.get("RECESSIVE_RESCUE", "1") != "0"

LOF = {"stop_gained","frameshift_variant","splice_acceptor_variant","splice_donor_variant",
       "start_lost","stop_lost","transcript_ablation"}
DAMAGING_MODERATE = {"missense_variant","inframe_deletion","inframe_insertion",
                     "protein_altering_variant"}


def load_bed(path):
    regions = defaultdict(list)
    try:
        for line in open(path):
            f = line.split("\t")
            regions[f[0]].append((int(f[1]), int(f[2])))
    except FileNotFoundError:
        pass
    return regions


def load_rejects(path="work/kmer_excluded.tsv"):
    out = set()
    try:
        for line in open(path):
            if line.startswith("#") or not line.strip():
                continue
            f = line.split("\t")
            out.add((f[0], f[1], f[2], f[3].strip()))
    except FileNotFoundError:
        pass
    return out


LOWCONF = load_bed("work/lowconf_regions.bed")
READ_REJECTS = load_rejects()
in_lowconf = lambda c, p: any(s <= p < e for s, e in LOWCONF.get(c, ()))


def float_or(x, d=0.0):
    try:
        return float(x)
    except (TypeError, ValueError):
        return d


def allele_balance(row):
    try:
        r, a = (int(x) for x in row["ad"].split(",")[:2])
    except (ValueError, KeyError):
        return None
    return a / (r + a) if (r + a) else None


def variant_evidence(row):
    cons = set(row["consequence"].split(","))
    score, tags = 0.0, []
    if QUALITY_TERMS:
        # Genotype confidence, applied to every genotype class. GATK's own
        # guidance for a confident call is DP >= 10 / GQ >= 20; clinical
        # reporting practice is stricter. Candidate rows are enriched ~1.6-fold
        # for calls below these thresholds relative to the whole CDS callset,
        # which is the signature of low-depth calls generating spurious
        # rare-and-damaging annotations (most often a low-depth het miscalled
        # as homozygous-alt).
        dp_i, gq_i = int(row.get("dp") or 0), int(row.get("gq") or 0)
        if dp_i < 10 or gq_i < 30:
            score -= 3.0; tags.append(f"DP{dp_i}_GQ{gq_i}_unreliable")
        elif dp_i < 15 or gq_i < 50:
            score -= 2.0; tags.append(f"DP{dp_i}_GQ{gq_i}_low")
    if QUALITY_TERMS and row["gt"] in ("0/1","1/0","0|1","1|0"):
        ab = allele_balance(row)
        if ab is not None and int(row.get("dp") or 0) >= 20:
            if ab < 0.20:
                score -= 3.0; tags.append(f"AB{ab:.2f}_artifact")
            elif ab < 0.30:
                score -= 1.5; tags.append(f"AB{ab:.2f}_low")
    if cons & LOF:
        score += 2.0; tags.append("LoF")
    # Computational missense evidence. AlphaMissense and CADD are NOT
    # independent -- CADD's training features overlap heavily with the
    # conservation and structural signals AlphaMissense uses, and both are
    # calibrated on overlapping variant sets. Adding them was double-counting.
    # MISSENSE_EVIDENCE controls how they are combined for missense rows:
    #   both  - legacy behaviour, AM + CADD summed (retained for the audit only)
    #   once  - the STRONGER of the two is counted, the other discarded
    #   none  - no computational missense evidence at all
    # Non-missense rows keep their CADD term: for a nonsense or splice variant
    # CADD is not standing in for a missense predictor.
    am = row.get("alphamissense_class", "")
    cadd = float_or(row.get("cadd_phred"))
    cadd_pts = 1.0 if cadd >= 25 else (0.5 if cadd >= 20 else 0.0)
    is_missense = "missense_variant" in cons
    if is_missense:
        if am == "likely_pathogenic":
            am_pts = 1.5; am_tag = "AM_path"
        elif am == "ambiguous":
            am_pts = 0.4; am_tag = "AM_amb"
        elif am == "likely_benign":
            am_pts = -0.5; am_tag = "AM_benign"
        else:
            am_pts = 0.2; am_tag = ""
        if MISSENSE_EVIDENCE == "none":
            pass
        elif MISSENSE_EVIDENCE == "once":
            # One line of evidence, chosen by specificity rather than by
            # magnitude: AlphaMissense is a missense-specific classifier, so
            # when it makes a call (pathogenic or benign) that call stands
            # alone and CADD is discarded as redundant. Taking a max over
            # signed scores would silently drop the benign penalty, so the
            # branch is on whether AM is informative, not on which is larger.
            if am_tag:
                score += am_pts; tags.append(am_tag + "(sole)")
            elif cadd_pts:
                score += cadd_pts; tags.append(f"CADD{cadd:.0f}(sole)")
            else:
                score += am_pts
        else:
            score += am_pts
            if am_tag: tags.append(am_tag)
            if cadd_pts: score += cadd_pts; tags.append(f"CADD{cadd:.0f}")
    else:
        if cons & DAMAGING_MODERATE:
            score += 0.8; tags.append("inframe")
        if cadd_pts:
            score += cadd_pts; tags.append(f"CADD{cadd:.0f}")
    sai = float_or(row.get("spliceai_max"))
    if sai >= 0.5:
        score += 1.5; tags.append(f"SpliceAI{sai:.2f}")
    elif sai >= 0.2:
        score += 0.5; tags.append(f"SpliceAI{sai:.2f}")
    if QUALITY_TERMS and in_lowconf(row["chrom"], int(row["pos"])):
        score -= 1.5; tags.append("lowconf_region")
    cln = row.get("clnsig", "")
    if cln.startswith(("Pathogenic","Likely_pathogenic")):
        rev = row.get("clnrevstat", "")
        score += 3.0 if any(k in rev for k in ("multiple_submitters","practice_guideline",
                                               "reviewed_by_expert_panel")) else 2.0
        tags.append(f"ClinVar:{cln}")
    elif cln.startswith(("Benign","Likely_benign")):
        score -= 3.0; tags.append("ClinVar:benign")
    return score, ";".join(tags)


OUT_COLS = ["gene_symbol","model","chrom","pos","ref","alt","gt","ad","dp","gq","consequence",
            "impact","hgvsc","hgvsp","cadd_phred","alphamissense_class","spliceai_max",
            "clnsig","clnrevstat","evidence","tags"]


def main(vep_path, gt_path, clinvar_path, out_path, genes_path):
    gts = {}
    for line in open(gt_path):
        f = line.rstrip("\n").split("\t")
        gts[(f[0],f[1],f[2],f[3])] = {"qual":f[4],"gt":f[5],"ad":f[6],"dp":f[7],"gq":f[8]}
    clinvar = {}
    for line in open(clinvar_path):
        f = line.rstrip("\n").split("\t")
        clinvar[(f[0],f[1],f[2],f[3])] = {"clnsig":f[5],"clnrevstat":f[6],"clndn":f[7]}

    rows = []
    with open(vep_path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) != len(header):
                continue
            row = dict(zip(header, f))
            key = (row["chrom"], row["pos"], row["ref"], row["alt"])
            row.update(clinvar.get(key, {"clnsig":"","clnrevstat":"","clndn":""}))
            g = gts.get(key)
            if not g or not row["gene_symbol"] or key in READ_REJECTS:
                continue
            row.update(g)
            sc, tags = variant_evidence(row)
            if sc < EVIDENCE_FLOOR:
                continue
            row["evidence"], row["tags"] = sc, tags
            rows.append(row)

    by_gene = defaultdict(list)
    for r in rows:
        by_gene[r["gene_symbol"]].append(r)

    kept, n_rescue = set(), 0
    with open(out_path, "w") as out:
        out.write("\t".join(OUT_COLS) + "\n")
        for gene, vs in by_gene.items():
            homs = [v for v in vs if v["gt"] in ("1/1","1|1")]
            hets = [v for v in vs if v["gt"] in ("0/1","1/0","0|1","1|0")]
            sel = [(v, "hom") for v in homs]
            if len(hets) >= 2:
                sel += [(v, "comphet") for v in sorted(hets, key=lambda v: -v["evidence"])[:4]]
            elif len(hets) == 1:
                sel += [(hets[0], "het_single")]
                if RESCUE:
                    sel += [(hets[0], "incomplete_recessive")]
                    n_rescue += 1
            for v, model in sel:
                kept.add(gene)
                w = dict(v); w["model"] = model
                out.write("\t".join(str(w.get(c, "")) for c in OUT_COLS) + "\n")

    with open(genes_path, "w") as out:
        for g in sorted(kept):
            out.write(g + "\n")
    print(f"variants_passing_evidence\t{len(rows)}")
    print(f"genes_with_candidates\t{len(kept)}")
    print(f"genes_with_rescue_route\t{n_rescue}")


if __name__ == "__main__":
    main(*sys.argv[1:6])
