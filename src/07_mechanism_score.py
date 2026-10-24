"""Disease-level, inheritance-gated phenotype scoring with coverage."""
import math, sys
from collections import defaultdict
from importlib import import_module
sys.path.insert(0, "src")
_p3 = import_module("03_phenotype_score")
parse_obo, ancestors_of = _p3.parse_obo, _p3.ancestors_of

AD = {"HP:0000006","HP:0000015","HP:0012275"}
AR = {"HP:0000007"}
XL = {"HP:0001417","HP:0001419","HP:0001423","HP:0001470"}
MT = {"HP:0001427"}
MODEL_COMPATIBLE = {"het_single": AD | XL | MT, "hom": AR | XL | MT, "comphet": AR,
                    "incomplete_recessive": AR}

def main(hpo_path, out_path, obo="ref/hp.obo", hpoa="ref/phenotype.hpoa",
         g2p="ref/genes_to_phenotype.txt"):
    parents, names, alt, obsolete = parse_obo(obo)
    cache = {}; norm = lambda t: alt.get(t, t)
    dzp, dzm, dzn = defaultdict(set), defaultdict(set), {}
    col = None
    for line in open(hpoa):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if col is None:
            col = {c: i for i, c in enumerate(f)}; continue
        dz, term, asp = f[col["database_id"]], f[col["hpo_id"]], f[col["aspect"]]
        dzn.setdefault(dz, f[col["disease_name"]])
        if asp == "I":
            dzm[dz].add(norm(term))
        elif asp == "P":
            dzp[dz].add(norm(term))
    counts = defaultdict(int)
    for terms in dzp.values():
        cl = set()
        for t in terms:
            cl |= ancestors_of(t, parents, cache)
        for t in cl:
            counts[t] += 1
    n_dz = len(dzp)
    ic = {t: -math.log(c/n_dz) for t, c in counts.items() if c > 0}
    proband = []
    for line in open(hpo_path):
        line = line.strip()
        if line and not line.startswith("#"):
            t = norm(line.split()[0])
            if t not in obsolete:
                proband.append(t)
    pa = {t: ancestors_of(t, parents, cache) for t in proband}
    self_ic = sum(ic.get(t, 0.0) for t in proband) or 1.0
    sim = lambda a, b: max((ic.get(c, 0.0) for c in (a & b)), default=0.0)
    sc, mt_, cov, nexp = {}, {}, {}, {}
    for dz, terms in dzp.items():
        ga = {t: ancestors_of(t, parents, cache) for t in terms}
        if not ga:
            continue
        fwd, pairs = [], []
        for a in proband:
            best, arg = 0.0, ""
            for b, bs in ga.items():
                s = sim(pa[a], bs)
                if s > best:
                    best, arg = s, b
            fwd.append(best); pairs.append((a, arg, best))
        rev = [max((sim(pa[a], bs) for a in proband), default=0.0) for bs in ga.values()]
        sc[dz] = 0.5*(sum(fwd)/len(fwd)+sum(rev)/len(rev))
        cov[dz] = min(1.0, sum(fwd)/self_ic)
        nexp[dz] = sum(1 for s in fwd if s > 0)
        mt_[dz] = ";".join(f"{names.get(a,a)}->{names.get(b,b)}:{s:.2f}"
                           for a, b, s in sorted(pairs, key=lambda x: -x[2])[:3] if s > 0)
    gdz = defaultdict(set)
    with open(g2p) as fh:
        next(fh)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 6 and f[5]:
                gdz[f[1]].add(f[5])
    with open(out_path, "w") as out:
        out.write("gene_symbol\tmodel\tmech_score\tdisease_id\tdisease_name\tdisease_moi\t"
                  "n_admissible\tn_total_diseases\tbest_incompatible_score\ttop_matches\t"
                  "pheno_coverage\tn_features_explained\n")
        for gene, dzs in sorted(gdz.items()):
            for model, ok in MODEL_COMPATIBLE.items():
                adm, inc = [], []
                for dz in dzs:
                    if dz not in sc:
                        continue
                    moi = dzm.get(dz, set())
                    (adm if (moi & ok) or not moi else inc).append(dz)
                best = max(adm, key=lambda d: sc[d], default=None)
                miss = max(inc, key=lambda d: sc[d], default=None)
                fl = [gene, model]
                fl += ([f"{sc[best]:.4f}", best, dzn.get(best,""),
                        ",".join(sorted(dzm.get(best, set())))] if best else [f"{0.0:.4f}","","",""])
                fl += [str(len(adm)), str(len(dzs))]
                fl.append(f"{sc[miss]:.4f}" if miss is not None else f"{0.0:.4f}")
                fl.append(mt_.get(best,"") if best else "")
                fl += ([f"{cov[best]:.4f}", str(nexp[best])] if best else [f"{0.0:.4f}","0"])
                out.write("\t".join(fl) + "\n")
    print(f"diseases_scored\t{len(sc)}\ngenes\t{len(gdz)}")

if __name__ == "__main__":
    main(*sys.argv[1:3])
