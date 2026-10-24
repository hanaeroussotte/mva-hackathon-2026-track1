"""Assign an explicit phase state to every compound-heterozygous hypothesis.

Three independent sources of phase information, in decreasing strength:

  1. GATK HaplotypeCaller physical phasing (PGT/PID). Two calls sharing a PID
     were assembled on the same local haplotype graph; identical PGT means the
     alternate alleles are on the SAME haplotype (cis), opposite PGT means
     opposite haplotypes (trans). This is read-derived and is the only
     experimental phase available genome-wide in this callset.
  2. A dedicated read-based phasing run (whatshap) over a targeted alignment,
     recorded in work/phase_manual.tsv for loci where it was attempted.
  3. Nothing. State is then `unresolved_phase` -- NOT an assumption of trans.

A `confirmed_cis` pair cannot be a compound heterozygote and is penalised
downstream; `supported_trans` is the only state that experimentally supports
the recessive hypothesis.
"""
import sys
from collections import defaultdict


def load_phase_tags(path):
    tags = {}
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if len(f) < 7:
            continue
        pgt, pid = f[5], f[6]
        if pid and pid != ".":
            tags[(f[0], f[1], f[2], f[3])] = (pgt, pid)
    return tags


def main(cand_path, phase_tag_path, out_path, manual_path="work/phase_manual.tsv"):
    tags = load_phase_tags(phase_tag_path)
    manual = {}
    try:
        with open(manual_path) as fh:
            next(fh)
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) >= 4:
                    manual[(f[0], f[1])] = (f[2], f[3])
    except (FileNotFoundError, StopIteration):
        pass

    per = defaultdict(list)
    with open(cand_path) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        ix = {c: i for i, c in enumerate(head)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            per[(f[ix["gene_symbol"]], f[ix["model"]])].append(
                (f[ix["chrom"]], f[ix["pos"]], f[ix["ref"]], f[ix["alt"]]))

    counts = defaultdict(int)
    with open(out_path, "w") as out:
        out.write("gene\tmodel\tpos_1\tpos_2\tphase_state\tevidence\n")
        for (gene, model), vs in sorted(per.items()):
            if model != "comphet" or len(vs) < 2:
                continue
            a, b = vs[0], vs[1]
            if (gene, model) in manual:
                state, ev = manual[(gene, model)]
            elif a in tags and b in tags and tags[a][1] == tags[b][1]:
                pgt_a, pgt_b = tags[a][0], tags[b][0]
                if pgt_a == pgt_b:
                    state = "confirmed_cis"
                    ev = f"GATK_PID={tags[a][1]};PGT={pgt_a}={pgt_b}"
                else:
                    state = "supported_trans"
                    ev = f"GATK_PID={tags[a][1]};PGT={pgt_a}vs{pgt_b}"
            else:
                state, ev = "unresolved_phase", "no_shared_PID;no_targeted_phasing"
            counts[state] += 1
            out.write(f"{gene}\t{model}\t{a[1]}\t{b[1]}\t{state}\t{ev}\n")
    for k in sorted(counts):
        print(f"{k}\t{counts[k]}")


if __name__ == "__main__":
    main(*sys.argv[1:5])
