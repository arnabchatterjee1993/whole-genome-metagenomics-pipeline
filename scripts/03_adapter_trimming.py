"""Stage 3: Adapter trimming and quality filtering with fastp."""

import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import config
from utils import detect_sample, require_tools

require_tools("fastp")
config.TRIMMED_DIR.mkdir(parents=True, exist_ok=True)

# --- 1. Locate stage 01's raw FASTQs ---
SAMPLE = detect_sample(config.MASTER_DIR)
r1 = config.RAW_DIR / f"{SAMPLE}_R1.fastq"
r2 = config.RAW_DIR / f"{SAMPLE}_R2.fastq"
for fq in [r1, r2]:
    if not fq.exists():
        sys.exit(f"ERROR: Input FASTQ not found: {fq}")
print("Input FASTQ files verified")

r1_trimmed = config.TRIMMED_DIR / f"{SAMPLE}_R1_trimmed.fastq"
r2_trimmed = config.TRIMMED_DIR / f"{SAMPLE}_R2_trimmed.fastq"
html_report = config.TRIMMED_DIR / f"{SAMPLE}_fastp_report.html"
json_report = config.TRIMMED_DIR / f"{SAMPLE}_fastp_report.json"

# --- 2. Auto-detect and trim adapters, then quality-trim both read ends ---
#    --cut_front/--cut_tail: sliding-window quality trim from both ends
#    --length_required 50:   drop reads shorter than 50 bp after trimming
print("Running fastp trimming...")
subprocess.run([
    "fastp",
    "-i", str(r1), "-I", str(r2),
    "-o", str(r1_trimmed), "-O", str(r2_trimmed),
    "--detect_adapter_for_pe",
    "--cut_front", "--cut_tail",
    "--cut_window_size", "4",
    "--cut_mean_quality", "20",
    "--length_required", "50",
    "--html", str(html_report),
    "--json", str(json_report),
], check=True)

# --- 3. Confirm fastp actually produced everything it was supposed to ---
for out in [r1_trimmed, r2_trimmed, html_report, json_report]:
    if not out.exists():
        sys.exit(f"ERROR: fastp output missing: {out}")
print("fastp trimming completed successfully")

print("Trimmed reads and reports:")
for f in config.TRIMMED_DIR.iterdir():
    print(f" - {f.name}")
