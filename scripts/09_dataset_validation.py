"""Stage 9: Final pre-processed dataset validation and gold-standard packaging."""

import gzip
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import config
from utils import count_reads_gz, detect_sample, require_paths

config.GOLD_DIR.mkdir(parents=True, exist_ok=True)

# --- 1. Locate stage 08's canonical normalized dataset ---
SAMPLE = detect_sample(config.MASTER_DIR)
final_r1 = config.NORM_DIR / f"{SAMPLE}_R1_norm.fastq.gz"
final_r2 = config.NORM_DIR / f"{SAMPLE}_R2_norm.fastq.gz"
raw_r1 = config.RAW_DIR / f"{SAMPLE}_R1.fastq"
require_paths(final_r1, final_r2)
print("Canonical normalized paired FASTQs detected.")

final_count_r1 = count_reads_gz(final_r1)
final_count_r2 = count_reads_gz(final_r2)

# raw_count is used below to report what fraction of reads survived the whole pipeline
raw_count = None
if raw_r1.exists():
    with open(raw_r1) as f:
        raw_count = sum(1 for _ in f) // 4

print("\n" + "=" * 40)
print("FINAL DATASET VALIDATION REPORT")
print("=" * 40)

# --- 2. R1 and R2 must have identical read counts, or something upstream desynced them ---
if final_count_r1 == final_count_r2:
    print(f"Paired-end integrity: PASSED ({final_count_r1} read pairs)")
else:
    print(f"Paired-end integrity: FAILED (R1={final_count_r1}, R2={final_count_r2})")

if raw_count:
    retention = (final_count_r1 / raw_count) * 100
    print(f"Data retention: {retention:.2f}% of raw reads retained")

# --- 3. Quick sanity metric: mean read length after all processing ---
print("Computing average read length...")
total_bases = 0
sequence_count = 0
with gzip.open(final_r1, "rt") as fh:
    for line_number, line in enumerate(fh):
        if line_number % 4 == 1:
            total_bases += len(line.strip())
            sequence_count += 1

if sequence_count == 0:
    sys.exit("ERROR: No sequences found in final R1 FASTQ")
mean_length = total_bases / sequence_count
print(f"Average read length: {mean_length:.2f} bp")

# --- 4. Copy the validated reads into GOLD_DIR -- this is what stage 10 analyzes ---
print("\nCreating final gold-standard copies...")
gold_r1 = config.GOLD_DIR / f"{SAMPLE}_clean_R1.fastq.gz"
gold_r2 = config.GOLD_DIR / f"{SAMPLE}_clean_R2.fastq.gz"
shutil.copy2(final_r1, gold_r1)
shutil.copy2(final_r2, gold_r2)
print(f"Gold-standard dataset saved to: {config.GOLD_DIR}")

for f in config.GOLD_DIR.iterdir():
    print(f"{f.name} - {f.stat().st_size / 1e6:.2f} MB")

print("\nFINAL VALIDATION COMPLETE.")
