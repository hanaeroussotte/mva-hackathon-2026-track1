"""Build genic and CDS interval BEDs from a GENCODE GTF."""
import gzip, re, sys
from collections import defaultdict
ATTR = re.compile(r'(\S+) "([^"]*)"')
PRIMARY = {str(c) for c in range(1, 23)} | {"X", "Y", "M"}
KEEP = {"protein_coding","lncRNA","miRNA","snRNA","snoRNA","rRNA","Mt_tRNA","Mt_rRNA"}
key = lambda c: (0, int(c)) if c.isdigit() else (1, c)

def merge_write(per, path):
    n = tot = 0
    with open(path, "w") as out:
        for c in sorted(per, key=key):
            m = []
            for s, e in sorted(per[c]):
                if m and s <= m[-1][1]:
                    m[-1][1] = max(m[-1][1], e)
                else:
                    m.append([s, e])
            for s, e in m:
                out.write(f"{c}\t{s}\t{e}\n"); n += 1; tot += e - s
    return n, tot

def main(gtf, genic_bed, spans_path, cds_bed):
    genic, cds, span = defaultdict(list), defaultdict(list), {}
    for line in gzip.open(gtf, "rt"):
        if line[0] == "#":
            continue
        f = line.rstrip("\n").split("\t")
        c = f[0][3:] if f[0].startswith("chr") else f[0]
        c = "M" if c in ("MT", "M") else c
        if c not in PRIMARY:
            continue
        if f[2] == "exon":
            a = dict(ATTR.findall(f[8]))
            if a.get("gene_type","") not in KEEP:
                continue
            genic[c].append((max(0, int(f[3])-1-10), int(f[4])+10))
            nm = a.get("gene_name", a.get("gene_id","."))
            r = span.get(nm)
            if r is None:
                span[nm] = [c, int(f[3]), int(f[4]), a.get("gene_id","."), a.get("gene_type","")]
            elif r[0] == c:
                r[1] = min(r[1], int(f[3])); r[2] = max(r[2], int(f[4]))
        elif f[2] in ("CDS","start_codon","stop_codon"):
            cds[c].append((max(0, int(f[3])-1-12), int(f[4])+12))
    n1, t1 = merge_write(genic, genic_bed)
    n2, t2 = merge_write(cds, cds_bed)
    with open(spans_path, "w") as out:
        out.write("gene_name\tcontig\tstart\tend\tgene_id\tgene_type\n")
        for nm, (c, s, e, gid, bt) in sorted(span.items(), key=lambda kv: (key(kv[1][0]), kv[1][1])):
            out.write(f"{nm}\t{c}\t{s}\t{e}\t{gid}\t{bt}\n")
    print(f"genic_intervals\t{n1}\tgenic_Mb\t{t1/1e6:.1f}")
    print(f"cds_intervals\t{n2}\tcds_Mb\t{t2/1e6:.1f}")
    print(f"genes\t{len(span)}")

if __name__ == "__main__":
    main(*sys.argv[1:5])
