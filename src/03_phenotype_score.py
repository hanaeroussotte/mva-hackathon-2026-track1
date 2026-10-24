"""Gene-level HPO similarity (pre-gate baseline; provides ontology helpers for 07)."""
import math, sys
from collections import defaultdict

def parse_obo(path):
    parents, names, alt, obsolete = defaultdict(set), {}, {}, set()
    term = None
    for line in open(path):
        line = line.rstrip("\n")
        if line == "[Term]":
            term = None
        elif line.startswith("id: HP:"):
            term = line[4:]
        elif term and line.startswith("name: "):
            names[term] = line[6:]
        elif term and line.startswith("is_a: HP:"):
            parents[term].add(line[6:].split(" ")[0])
        elif term and line.startswith("alt_id: HP:"):
            alt[line[8:]] = term
        elif term and line.startswith("is_obsolete: true"):
            obsolete.add(term)
    return parents, names, alt, obsolete

def ancestors_of(term, parents, cache):
    if term in cache:
        return cache[term]
    out = {term}; stack = list(parents.get(term, ()))
    while stack:
        p = stack.pop()
        if p not in out:
            out.add(p); stack.extend(parents.get(p, ()))
    cache[term] = out
    return out

def main(hpo_path, out_path, obo="ref/hp.obo", hpoa="ref/phenotype.hpoa",
         g2p="ref/genes_to_phenotype.txt"):
    parents, names, alt, obsolete = parse_obo(obo)
    cache = {}; norm = lambda t: alt.get(t, t)
    dzt = defaultdict(set)
    for line in open(hpoa):
        if line.startswith("#") or line.startswith("database_id"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) >= 4 and f[3].startswith("HP:"):
            dzt[f[0]].add(norm(f[3]))
    counts = defaultdict(int)
    for terms in dzt.values():
        cl = set()
        for t in terms:
            cl |= ancestors_of(t, parents, cache)
        for t in cl:
            counts[t] += 1
    n_dz = len(dzt)
    ic = {t: -math.log(c/n_dz) for t, c in counts.items() if c > 0}
    gt = defaultdict(set)
    with open(g2p) as fh:
        next(fh)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 4:
                gt[f[1]].add(norm(f[2]))
    proband = []
    for line in open(hpo_path):
        line = line.strip()
        if line and not line.startswith("#"):
            t = norm(line.split()[0])
            if t not in obsolete:
                proband.append(t)
    pa = {t: ancestors_of(t, parents, cache) for t in proband}
    sim = lambda a, b: max((ic.get(c, 0.0) for c in (a & b)), default=0.0)
    rows = []
    for gene, terms in gt.items():
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
        rows.append((gene, 0.5*(sum(fwd)/len(fwd)+sum(rev)/len(rev)), sum(fwd)/len(fwd), len(terms),
                     ";".join(f"{names.get(a,a)}->{names.get(b,b)}:{s:.2f}"
                              for a, b, s in sorted(pairs, key=lambda x: -x[2])[:4] if s > 0)))
    rows.sort(key=lambda r: -r[1])
    with open(out_path, "w") as out:
        out.write("gene_symbol\tpheno_score\tfwd_score\tn_gene_terms\ttop_matches\n")
        for g, s, f, n, m in rows:
            out.write(f"{g}\t{s:.4f}\t{f:.4f}\t{n}\t{m}\n")
    print(f"proband_terms\t{len(proband)}\ndiseases_for_IC\t{n_dz}\ngenes_scored\t{len(rows)}")

if __name__ == "__main__":
    main(*sys.argv[1:3])
