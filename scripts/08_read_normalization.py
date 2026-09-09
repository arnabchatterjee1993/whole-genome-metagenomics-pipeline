"""Stage 8: Read normalization (down-sampling) to a target depth.

The normalized dataset is the canonical input for the final validation and
downstream analysis stages. If the dataset is already below the target
depth, the paired reads are copied unchanged so the pipeline remains
continuous and deterministic.
"""

import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import config
from utils import count_reads_gz, detect_sample, require_paths, require_tools

config.NORM_DIR.mkdir(parents=True, exist_ok=True)
require_tools("seqtk", "pigz")

# --- 1. Locate stage 07's deduplicated reads ---
SAMPLE = detect_sample(config.MASTER_DIR)
final_gz_r1 = config.FINAL_CLEAN_DIR / f"{SAMPLE}_clean_R1.fastq.gz"
final_gz_r2 = config.FINAL_CLEAN_DIR / f"{SAMPLE}_clean_R2.fastq.gz"
require_paths(final_gz_r1, final_gz_r2)
print(f"Found paired FASTQs in: {config.FINAL_CLEAN_DIR}")

total_reads = count_reads_gz(final_gz_r1)
print(f"Current total microbial read pairs: {total_reads}")

norm_r1 = config.NORM_DIR / f"{SAMPLE}_R1_norm.fastq.gz"
norm_r2 = config.NORM_DIR / f"{SAMPLE}_R2_norm.fastq.gz"

# --- 2. Only down-sample if we're actually above the target depth ---
if total_reads > config.TARGET_READ_DEPTH:
    print(f"Normalizing paired-end reads to {config.TARGET_READ_DEPTH} read pairs...")
    r1_tmp = config.NORM_DIR / f"{SAMPLE}_R1_norm.fastq"
    r2_tmp = config.NORM_DIR / f"{SAMPLE}_R2_norm.fastq"

    # seed -s100 keeps R1/R2 sub-sampling in sync so pairs stay matched
    subprocess.run(["seqtk", "sample", "-s100", str(final_gz_r1), str(config.TARGET_READ_DEPTH)],
                    stdout=open(r1_tmp, "w"), check=True)
    subprocess.run(["seqtk", "sample", "-s100", str(final_gz_r2), str(config.TARGET_READ_DEPTH)],
                    stdout=open(r2_tmp, "w"), check=True)

    subprocess.run(["pigz", "-f", str(r1_tmp)], check=True)
    subprocess.run(["pigz", "-f", str(r2_tmp)], check=True)
    r1_tmp.with_suffix(".fastq.gz").rename(norm_r1)
    r2_tmp.with_suffix(".fastq.gz").rename(norm_r2)
    print("Normalization complete")
elif total_reads > 0:
    print("Total reads already below target depth. Keeping all reads.")
    import shutil
    shutil.copy2(final_gz_r1, norm_r1)
    shutil.copy2(final_gz_r2, norm_r2)
    print("Normalization not required; canonical normalized dataset created.")
else:
    sys.exit("ERROR: No reads detected; cannot create normalized dataset.")

for fq in [norm_r1, norm_r2]:
    if not fq.exists() or fq.stat().st_size == 0:
        sys.exit(f"ERROR: Normalized FASTQ missing or empty: {fq}")

print(f"Normalized dataset ready: {norm_r1.name}, {norm_r2.name}")
