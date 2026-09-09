"""Stage 9: Final pre-processed dataset validation and gold-standard packaging."""

import shutil
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import config
from utils import count_reads_gz, detect_sample, require_paths

config.GOLD_DIR.mkdir(parents=True, exist_ok=True)

# --- 1. Locate stage 07's final clean reads (see stage 08's note on normalization) ---
SAMPLE = detect_sample(config.MASTER_DIR)
final_r1 = config.FINAL_CLEAN_DIR / f"{SAMPLE}_clean_R1.fastq.gz"
final_r2 = config.FINAL_CLEAN_DIR / f"{SAMPLE}_clean_R2.fastq.gz"
raw_r1 = config.RAW_DIR / f"{SAMPLE}_R1.fastq"
require_paths(final_r1, final_r2)
print("Final paired FASTQs detected.")

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

# --- 3. Quick sanity metric: mean read length after all the trimming/filtering ---
print("Computing average read length...")
result = subprocess.run(
    f"zcat {final_r1} | awk 'NR%4==2 {{sum+=length($0); n++}} END {{if(n>0) print sum/n}}'",
    shell=True, stdout=subprocess.PIPE, check=True,
)
print(f"Average read length: {result.stdout.decode().strip()} bp")

# --- 4. Copy the validated reads into GOLD_DIR -- this is what stage 10 analyzes ---
print("\nCreating final gold-standard copies...")
shutil.copy(final_r1, config.GOLD_DIR / final_r1.name)
shutil.copy(final_r2, config.GOLD_DIR / final_r2.name)
print(f"Gold-standard dataset saved to: {config.GOLD_DIR}")

for f in config.GOLD_DIR.iterdir():
    print(f"{f.name} - {f.stat().st_size / 1e6:.2f} MB")

print("\nFINAL VALIDATION COMPLETE.")
