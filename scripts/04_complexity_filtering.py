"""Stage 4: Low-complexity read filtering with BBDuk.

Reads that are mostly repeats (e.g. poly-A runs, low-diversity artifacts)
carry little biological signal and can skew downstream taxonomic/assembly
results, so this stage drops them using a sequence-entropy cutoff.
"""

import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import config
from utils import detect_sample, require_paths

config.FILTERED_DIR.mkdir(parents=True, exist_ok=True)
require_paths(config.BBDUK)
print("BBTools is available")

# --- 1. Locate stage 03's trimmed FASTQs ---
SAMPLE = detect_sample(config.MASTER_DIR)
t1 = config.TRIMMED_DIR / f"{SAMPLE}_R1_trimmed.fastq"
t2 = config.TRIMMED_DIR / f"{SAMPLE}_R2_trimmed.fastq"
require_paths(t1, t2)
print("Trimmed FASTQs verified")

f1 = config.FILTERED_DIR / f"{SAMPLE}_R1_final.fastq"
f2 = config.FILTERED_DIR / f"{SAMPLE}_R2_final.fastq"
stats_file = config.FILTERED_DIR / f"{SAMPLE}_complexity_stats.txt"

# --- 2. Drop reads below the entropy threshold (entropy=0.7 over a 50 bp window) ---
print("Running BBDuk complexity filtering...")
subprocess.run([
    str(config.BBDUK),
    f"in1={t1}", f"in2={t2}",
    f"out1={f1}", f"out2={f2}",
    "entropy=0.7", "entropywindow=50", "entropymask=f", "maxns=0",
    f"stats={stats_file}",
], check=True)

for out in [f1, f2, stats_file]:
    if not out.exists():
        sys.exit(f"ERROR: BBDuk output missing: {out}")
print("Complexity filtering completed successfully")

print("--- Complexity Filtering Summary ---")
print(stats_file.read_text())
