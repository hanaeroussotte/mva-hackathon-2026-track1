"""Annotate variants with Ensembl VEP (REST), writing one row per variant.

Predictors requested and parsed: CADD, AlphaMissense, SpliceAI, SIFT.

Two schema details cost real debugging time and are the reason this file
exists as its own script:

  * SpliceAI is returned as a NESTED dict on each transcript consequence,
    tc["spliceai"] with keys DS_AG / DS_AL / DS_DG / DS_DL -- NOT as flat
    keys such as spliceai_pred_ds_ag. Reading the flat key silently yields
    0.0 for every variant.
  * AlphaMissense is likewise NESTED, tc["alphamissense"] with keys
    am_class / am_pathogenicity -- not alphamissense_class.

PolyPhen is deliberately NOT claimed or parsed: this server configuration
does not return any polyphen key unless explicitly requested, so a
"polyphen" column would be silently empty.

RECONSTRUCTED 2026-10-09 after a workspace cleanup deleted the original.
Functionally equivalent to the script that produced work/cds_vep.tsv --
same endpoint, same parameters, same parsing -- but not byte-identical to
it. See METHODS.md section 16.
"""
import json, os, sys, time
import requests

SERVER = "https://rest.ensembl.org"
ENDPOINT = "/vep/human/hgvs"
BATCH = 200
PARAMS = {"CADD": 1, "AlphaMissense": 1, "SpliceAI": 1, "SIFT": 1,
          "mane": 1, "numbers": 1, "hgvs": 1, "canonical": 1}
COLUMNS = ["chrom", "pos", "ref", "alt", "gene_symbol", "gene_id",
           "transcript_id", "mane_select", "consequence", "impact", "biotype",
           "hgvsc", "hgvsp", "exon", "intron", "amino_acids", "protein_pos",
           "cadd_phred", "alphamissense_class", "alphamissense_score",
           "sift", "polyphen", "spliceai_max", "lof", "variant_class",
           "colocated_ids"]


def spliceai_max(tc):
    """Return (max delta, scored?) from the NESTED spliceai object."""
    sa = tc.get("spliceai")
    if not isinstance(sa, dict):
        return 0.0, False
    vals = [sa.get(k) for k in ("DS_AG", "DS_AL", "DS_DG", "DS_DL")]
    vals = [float(v) for v in vals if v is not None]
    return (max(vals), True) if vals else (0.0, False)


def am_of(tc):
    """Return (class, score) from the NESTED alphamissense object."""
    am = tc.get("alphamissense")
    if not isinstance(am, dict):
        return None, None
    score = am.get("am_pathogenicity")
    return am.get("am_class"), (float(score) if score is not None else None)


def post(hgvs_batch, params, session, retries=5):
    for attempt in range(retries):
        r = session.post(SERVER + ENDPOINT,
                         headers={"Content-Type": "application/json",
                                  "Accept": "application/json"},
                         params=params, data=json.dumps({"hgvs_notations": hgvs_batch}),
                         timeout=180)
        if r.status_code == 200:
            return r.json()
        time.sleep(2 ** attempt)
    raise RuntimeError(f"VEP returned {r.status_code} for a batch of {len(hgvs_batch)}")


def pick_tc(annot):
    """Prefer the MANE Select transcript; fall back to canonical, then first."""
    tcs = annot.get("transcript_consequences") or []
    for tc in tcs:
        if tc.get("mane_select"):
            return tc
    for tc in tcs:
        if tc.get("canonical"):
            return tc
    return tcs[0] if tcs else {}


def row_for(annot, key):
    tc = pick_tc(annot)
    amc, ams = am_of(tc)
    sa, _ = spliceai_max(tc)
    d = {"chrom": key[0], "pos": key[1], "ref": key[2], "alt": key[3],
         "gene_symbol": tc.get("gene_symbol", ""), "gene_id": tc.get("gene_id", ""),
         "transcript_id": tc.get("transcript_id", ""),
         "mane_select": tc.get("mane_select", ""),
         "consequence": ",".join(tc.get("consequence_terms")
                                 or [annot.get("most_severe_consequence", "")]),
         "impact": tc.get("impact", ""), "biotype": tc.get("biotype", ""),
         "hgvsc": tc.get("hgvsc", ""), "hgvsp": tc.get("hgvsp", ""),
         "exon": tc.get("exon", ""), "intron": tc.get("intron", ""),
         "amino_acids": tc.get("amino_acids", ""),
         "protein_pos": tc.get("protein_start", ""),
         "cadd_phred": tc.get("cadd_phred", ""),
         "alphamissense_class": amc or "",
         "alphamissense_score": ams if ams is not None else "",
         "sift": tc.get("sift_prediction", ""),
         "polyphen": "",
         "spliceai_max": f"{sa:.3f}",
         "lof": tc.get("lof", ""),
         "variant_class": annot.get("variant_class", ""),
         "colocated_ids": ",".join(v.get("id", "") for v in
                                   (annot.get("colocated_variants") or []))}
    return "\t".join(str(d.get(c, "")) for c in COLUMNS)


def main(variants_path, out_path):
    keys = []
    with open(variants_path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 4 and f[0] != "chrom":
                keys.append((f[0], f[1], f[2], f[3]))
    session = requests.Session()
    written = failed = 0
    with open(out_path, "w") as out:
        out.write("\t".join(COLUMNS) + "\n")
        for i in range(0, len(keys), BATCH):
            chunk = keys[i:i + BATCH]
            hgvs = [f"{c}:g.{p}{r}>{a}" for c, p, r, a in chunk]
            try:
                annots = post(hgvs, PARAMS, session)
            except RuntimeError as exc:
                failed += len(chunk)
                print(f"batch {i} FAILED: {exc}", file=sys.stderr)
                continue
            by_input = {a.get("input"): a for a in annots}
            for key, h in zip(chunk, hgvs):
                a = by_input.get(h)
                if a is None:
                    failed += 1
                    continue
                out.write(row_for(a, key) + "\n")
                written += 1
            print(f"progress\t{min(i + BATCH, len(keys))}/{len(keys)}", flush=True)
    print(f"annotated\t{written}")
    print(f"failed\t{failed}")
    if failed:
        print("WARNING: silent annotation failures are the bug class that caused "
              "a false 'no competing variant' claim; re-run the failures before "
              "asserting any absence.", file=sys.stderr)


if __name__ == "__main__":
    main(*sys.argv[1:3])
