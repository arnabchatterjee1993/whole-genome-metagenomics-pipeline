"""Stage 7: Repair paired-end reads and remove duplicates with Clumpify.

PCR/optical duplicates inflate abundance estimates downstream, so this
stage collapses exact-duplicate read pairs to one representative each.
"""

import shutil
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import config
from utils import count_reads_gz, detect_sample, require_paths

config.DEDUP_DIR.mkdir(parents=True, exist_ok=True)
require_paths(config.REPAIR_SH, config.CLUMPIFY)

# --- 1. Locate stage 06's compressed clean reads ---
SAMPLE = detect_sample(config.MASTER_DIR)
depleted_r1 = config.CLEAN_BACKUP_DIR / f"{SAMPLE}_clean_R1.fastq.gz"
depleted_r2 = config.CLEAN_BACKUP_DIR / f"{SAMPLE}_clean_R2.fastq.gz"
require_paths(depleted_r1, depleted_r2)
print("Input FASTQs found. Starting repair step...")

repaired_r1 = config.DEDUP_DIR / f"{SAMPLE}_repaired_R1.fastq.gz"
repaired_r2 = config.DEDUP_DIR / f"{SAMPLE}_repaired_R2.fastq.gz"
singletons = config.DEDUP_DIR / f"{SAMPLE}_singletons.fastq.gz"

# --- 2. Re-sync R1/R2 pairing (Clumpify requires strictly paired input) ---
#    Any read whose mate got dropped upstream goes into `singletons` instead.
print("Running repair.sh...")
subprocess.run([
    str(config.REPAIR_SH),
    f"in1={depleted_r1}", f"in2={depleted_r2}",
    f"out1={repaired_r1}", f"out2={repaired_r2}",
    f"outs={singletons}", "repair",
], check=True)
for fq in [repaired_r1, repaired_r2]:
    if not fq.exists():
        sys.exit(f"ERROR: Repair step failed: {fq}")
print("Repair completed successfully")

unique_r1 = config.DEDUP_DIR / f"{SAMPLE}_clean_R1.fastq.gz"
unique_r2 = config.DEDUP_DIR / f"{SAMPLE}_clean_R2.fastq.gz"

# --- 3. Deduplicate: subs=0 means only exact-sequence duplicates are collapsed ---
print("Running Clumpify deduplication...")
subprocess.run([
    str(config.CLUMPIFY),
    f"in1={repaired_r1}", f"in2={repaired_r2}",
    f"out1={unique_r1}", f"out2={unique_r2}",
    "dedupe=t", "optical=t", "subs=0",
], check=True)
for fq in [unique_r1, unique_r2]:
    if not fq.exists():
        sys.exit(f"ERROR: Clumpify failed to produce output: {fq}")
print("Deduplication completed")

# --- 4. Back up the deduplicated reads (this becomes the "final clean" dataset) ---
config.FINAL_CLEAN_DIR.mkdir(parents=True, exist_ok=True)
shutil.copy2(unique_r1, config.FINAL_CLEAN_DIR / unique_r1.name)
shutil.copy2(unique_r2, config.FINAL_CLEAN_DIR / unique_r2.name)
shutil.copy2(singletons, config.FINAL_CLEAN_DIR / f"{SAMPLE}_orphans.fastq.gz")
for name in [unique_r1.name, unique_r2.name, f"{SAMPLE}_orphans.fastq.gz"]:
    if not (config.FINAL_CLEAN_DIR / name).exists():
        sys.exit(f"ERROR: Backup failed for {name}")
print(f"Final cleaned data backed up to: {config.FINAL_CLEAN_DIR}")

# --- 5. The dedup working files are now redundant with the backup; remove them ---
for f in [repaired_r1, repaired_r2, unique_r1, unique_r2, singletons]:
    if f.exists():
        f.unlink()
print("Intermediate files cleaned")

# --- 6. Sanity-check final read counts ---
print("Final Read Count Verification:")
for fq in [config.FINAL_CLEAN_DIR / f"{SAMPLE}_clean_R1.fastq.gz",
           config.FINAL_CLEAN_DIR / f"{SAMPLE}_clean_R2.fastq.gz"]:
    print(f"{fq.name}: {count_reads_gz(fq)} reads")

print(f"SUCCESS! Deduplicated data saved in {config.FINAL_CLEAN_DIR}")
