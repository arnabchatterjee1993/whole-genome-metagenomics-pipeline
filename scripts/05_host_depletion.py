"""Stage 5: Host (human) DNA contamination removal via Bowtie2 + Samtools.

Aligns reads against the human reference genome and keeps only the pairs
that DON'T map -- those are the microbial reads we actually care about
for the downstream taxonomy/resistome/assembly steps.
"""

import os
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import config
from utils import detect_sample, require_paths

config.HOST_DEPLETED_DIR.mkdir(parents=True, exist_ok=True)
require_paths(config.BOWTIE2_BIN, config.SAMTOOLS_BIN, config.HUMAN_INDEX_DIR)
print("Bowtie2 and Samtools binaries found")

# --- 1. Point Bowtie2 at the pre-built human genome index ---
if not any(config.HUMAN_INDEX_DIR.glob("*.bt2")):
    sys.exit(f"ERROR: No .bt2 files found in {config.HUMAN_INDEX_DIR}")
bt2_prefix = str(config.HUMAN_INDEX_DIR / config.HUMAN_INDEX_PREFIX)
os.environ["BT2_INDEX"] = bt2_prefix
print(f"Human genome index found. Prefix: {bt2_prefix}")

# --- 2. Locate stage 04's complexity-filtered FASTQs ---
SAMPLE = detect_sample(config.MASTER_DIR)
f1 = config.FILTERED_DIR / f"{SAMPLE}_R1_final.fastq"
f2 = config.FILTERED_DIR / f"{SAMPLE}_R2_final.fastq"
require_paths(f1, f2)
print("Filtered FASTQs verified")

bam_file = config.HOST_DEPLETED_DIR / f"{SAMPLE}_unmapped.bam"
r1_depleted = config.HOST_DEPLETED_DIR / f"{SAMPLE}_R1_depleted.fastq"
r2_depleted = config.HOST_DEPLETED_DIR / f"{SAMPLE}_R2_depleted.fastq"

# --- 3. Align to the human genome, keep only unmapped pairs (samtools flag -f 12) ---
#    Piped directly from Bowtie2's stdout into `samtools view` -- no huge
#    intermediate SAM file ever touches disk.
print("Running host depletion with Bowtie2...")
bowtie2_cmd = [
    str(config.BOWTIE2_BIN), "-x", bt2_prefix,
    "-1", str(f1), "-2", str(f2),
    "--fast-local", "-p", str(config.THREADS),
]
with open(bam_file, "wb") as bam_out:
    bowtie_proc = subprocess.Popen(bowtie2_cmd, stdout=subprocess.PIPE)
    try:
        subprocess.run(
            [str(config.SAMTOOLS_BIN), "view", "-b", "-f", "12"],
            stdin=bowtie_proc.stdout, stdout=bam_out, check=True,
        )
    finally:
        if bowtie_proc.stdout is not None:
            bowtie_proc.stdout.close()

    bowtie_returncode = bowtie_proc.wait()
    if bowtie_returncode != 0:
        sys.exit(f"ERROR: Bowtie2 failed during host depletion (exit code {bowtie_returncode})")

if not bam_file.exists() or bam_file.stat().st_size < 1000:
    sys.exit("ERROR: BAM file empty or missing - check alignment")
print("Host-depleted BAM created")

# --- 4. Convert the unmapped-reads BAM back to FASTQ for downstream tools ---
print("Converting BAM to host-depleted FASTQ...")
subprocess.run([
    str(config.SAMTOOLS_BIN), "fastq",
    "-1", str(r1_depleted), "-2", str(r2_depleted),
    "-0", "/dev/null", "-s", "/dev/null", "-n",
    str(bam_file),
], check=True)

for fq in [r1_depleted, r2_depleted]:
    if not fq.exists():
        sys.exit(f"ERROR: Host-depleted FASTQ missing: {fq}")
print("Host removal complete.")
