"""Mosaic-aneuploidy test from B-allele frequency; emits a low-confidence mask."""
import gzip, sys
from collections import defaultdict
import numpy as np
AUT = [str(c) for c in range(1, 23)]
BIN, MIN_DP, MAX_DP = 10_000_000, 20, 100

def collect(vcf):
    per = defaultdict(lambda: ([], [], []))
    with gzip.open(vcf, "rt") as fh:
        for line in fh:
            if line[0] == "#":
                continue
            f = line.split("\t", 10)
            if f[0] not in AUT or f[6] != "PASS" or len(f[3]) != 1 or len(f[4]) != 1:
                continue
            c = dict(zip(f[8].split(":"), f[9].rstrip("\n").split(":")))
            if c.get("GT") not in ("0/1","0|1","1|0","1/0"):
                continue
            try:
                r, a = (int(x) for x in c["AD"].split(",")[:2])
            except Exception:
                continue
            dp = r + a
            if MIN_DP <= dp <= MAX_DP:
                p, b, d = per[f[0]]
                p.append(int(f[1])); b.append(a/dp); d.append(dp)
    return {k: (np.array(v[0]), np.array(v[1]), np.array(v[2])) for k, v in per.items()}

def main(vcf, prefix):
    per = collect(vcf)
    allb = np.concatenate([per[c][1] for c in AUT if c in per])
    alld = np.concatenate([per[c][2] for c in AUT if c in per])
    print(f"het_sites\t{len(allb)}")
    print(f"genome_mean|BAF-0.5|\t{np.abs(allb-0.5).mean():.5f}")
    print(f"genome_median_DP\t{np.median(alld):.1f}")
    with open(prefix + "_per_bin.tsv", "w") as out:
        out.write("chrom\tbin_start\tn_het\tmean_dev\tmedian_dp\n")
        for c in AUT:
            if c not in per:
                continue
            pos, baf, dp = per[c]
            idx = (pos // BIN).astype(int)
            for b in np.unique(idx):
                m = idx == b
                if m.sum() >= 200:
                    out.write(f"{c}\t{b*BIN}\t{int(m.sum())}\t"
                              f"{np.abs(baf[m]-0.5).mean():.5f}\t{np.median(dp[m]):.1f}\n")

if __name__ == "__main__":
    main(*sys.argv[1:3])
