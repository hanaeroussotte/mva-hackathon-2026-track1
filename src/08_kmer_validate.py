"""Aligner-independent allele confirmation by exact k-mer counting on FASTQ.

For each candidate variant, build the 31-mer centred on the reference allele
and the 31-mer centred on the alternate allele, then count exact occurrences
(and reverse complements) directly in the raw reads. No aligner is involved,
so a genotype that depends on how reads were placed will not reproduce here.

Scope limits that must travel with any result from this script:
  * a 31-mer must fall entirely within a 150 bp read, and only lane 1 was
    used, so effective coverage is well below the sequencing depth -- counts
    of 6-15 per site are expected;
  * the assay is QUALITATIVE. It establishes presence or absence of allele
    support, not an allele fraction;
  * absence of the alternate k-mer is consistent with mismapping but does not
    demonstrate it. Demonstrating a paralogue origin needs the k-mer screened
    against the paralogue and against the whole genome, at several offsets.

RECONSTRUCTED 2026-10-09 after a workspace cleanup deleted the original.
Functionally equivalent, not byte-identical. See METHODS.md section 16.

Performance note: pyahocorasick over 64 MB chunks sustains ~143 MB/s, so a
~100 GB FASTQ pass takes ~11 minutes. BSD grep is unusable at this scale.
"""
import gzip, sys
import ahocorasick
import pysam

K = 31
FLANK = K // 2
COMP = str.maketrans("ACGTN", "TGCAN")


def revcomp(s):
    return s.translate(COMP)[::-1]


def kmers_for(fasta, chrom, pos, ref, alt):
    """Return (ref_kmer, alt_kmer) centred on the variant position."""
    start = pos - 1 - FLANK
    left = fasta.fetch(chrom, start, pos - 1).upper()
    right = fasta.fetch(chrom, pos - 1 + len(ref), pos - 1 + len(ref) + FLANK).upper()
    return (left + ref.upper() + right)[:K], (left + alt.upper() + right)[:K]


def build(variants, fasta_path):
    fa = pysam.FastaFile(fasta_path)
    auto = ahocorasick.Automaton()
    meta = {}
    for vid, (chrom, pos, ref, alt) in variants.items():
        r, a = kmers_for(fa, chrom, int(pos), ref, alt)
        for label, seq in (("ref", r), ("alt", a)):
            for strand_seq in {seq, revcomp(seq)}:
                auto.add_word(strand_seq, (vid, label))
        meta[vid] = {"ref_kmer": r, "alt_kmer": a, "ref": 0, "alt": 0}
    auto.make_automaton()
    return auto, meta


def scan(auto, meta, fastq_path, chunk_bytes=64 * 1024 * 1024):
    with gzip.open(fastq_path, "rt") as fh:
        carry = ""
        while True:
            block = fh.read(chunk_bytes)
            if not block:
                break
            buf = carry + block
            for _, (vid, label) in auto.iter(buf):
                meta[vid][label] += 1
            carry = buf[-(K - 1):]
    return meta


def main(variants_path, fasta_path, fastq_path, out_path):
    variants = {}
    with open(variants_path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 4 and f[0] != "chrom":
                variants[f"{f[0]}:{f[1]}{f[2]}>{f[3]}"] = (f[0], f[1], f[2], f[3])
    auto, meta = build(variants, fasta_path)
    meta = scan(auto, meta, fastq_path)
    with open(out_path, "w") as out:
        out.write("variant\tref_kmer\talt_kmer\tref_count\talt_count\tverdict\n")
        for vid, m in meta.items():
            verdict = ("no alternate support" if m["alt"] == 0 and m["ref"] > 0
                       else "both alleles present" if m["alt"] and m["ref"]
                       else "alternate only" if m["alt"] else "uninformative")
            out.write(f"{vid}\t{m['ref_kmer']}\t{m['alt_kmer']}\t"
                      f"{m['ref']}\t{m['alt']}\t{verdict}\n")
    print(f"variants\t{len(meta)}")


if __name__ == "__main__":
    main(*sys.argv[1:5])
