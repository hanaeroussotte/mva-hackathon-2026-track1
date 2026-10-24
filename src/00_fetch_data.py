"""Fetch the challenge inputs from the gated Hugging Face dataset.

RECONSTRUCTED 2026-10-09 after a workspace cleanup deleted the original.
Functionally equivalent to the script that fetched the inputs, not
byte-identical to it. See METHODS.md section 16.

Requires a Hugging Face token with access to the gated dataset, read from
the HUGGINGFACE environment variable. The token is never printed.

Note: on this machine the huggingface_hub python client fails TLS through
the sandbox proxy (OSStatus -26276); curl works. Hence the curl path.
"""
import os, subprocess, sys

REPO = "SageBio/mva-hackathon-2026-data"
BASE = f"https://huggingface.co/datasets/{REPO}/resolve/main"
FILES = ["WGS_EX2312012_HGWCNDSX7.vcf.gz",
         "WGS_EX2312012_HGWCNDSX7.vcf.gz.tbi",
         "WGS_EX2312012_HGWCNDSX7_R1.fastq.gz",
         "WGS_EX2312012_HGWCNDSX7_R2.fastq.gz",
         "Challenge_Clinical_Phenotype_1.docx"]


def fetch(name, outdir="data"):
    token = os.environ.get("HUGGINGFACE") or os.environ.get("HF_TOKEN")
    if not token:
        sys.exit("No Hugging Face token in HUGGINGFACE or HF_TOKEN")
    os.makedirs(outdir, exist_ok=True)
    dst = os.path.join(outdir, name)
    cmd = ["curl", "-sSL", "--max-time", "7200",
           "-H", f"Authorization: Bearer {token}",
           "-o", dst, f"{BASE}/{name}"]
    subprocess.run(cmd, check=True)
    return dst, os.path.getsize(dst)


if __name__ == "__main__":
    wanted = sys.argv[1:] or FILES
    for name in wanted:
        path, size = fetch(name)
        print(f"{name}\t{size}")
