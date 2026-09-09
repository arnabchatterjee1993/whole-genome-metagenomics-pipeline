"""Stage 2: Initial quality assessment of raw reads (FastQC + MultiQC)."""

import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import config
from utils import detect_sample, require_tools

require_tools("fastqc", "multiqc")
config.FASTQC_DIR.mkdir(parents=True, exist_ok=True)
config.MULTIQC_DIR.mkdir(parents=True, exist_ok=True)

# --- 1. Locate stage 01's output ---
SAMPLE = detect_sample(config.MASTER_DIR)
r1 = config.RAW_DIR / f"{SAMPLE}_R1.fastq"
r2 = config.RAW_DIR / f"{SAMPLE}_R2.fastq"
for fq in [r1, r2]:
    if not fq.exists():
        sys.exit(f"ERROR: FASTQ file not found: {fq}")
print("FASTQ files verified")

# --- 2. Per-read QC metrics (per-base quality, adapter content, GC%, etc.) ---
print("Running FastQC...")
subprocess.run(["fastqc", str(r1), str(r2), "-o", str(config.FASTQC_DIR)], check=True)
if not list(config.FASTQC_DIR.glob("*_fastqc.html")):
    sys.exit("ERROR: FastQC did not generate output files")
print("FastQC completed successfully")

# --- 3. Roll the R1/R2 FastQC reports into one summary HTML page ---
print("Aggregating reports with MultiQC...")
subprocess.run(["multiqc", str(config.FASTQC_DIR), "-o", str(config.MULTIQC_DIR)], check=True)
if not (config.MULTIQC_DIR / "multiqc_report.html").exists():
    sys.exit("ERROR: MultiQC report not generated")
print("MultiQC report generated successfully")

# --- 4. FastQC's .zip files just duplicate the .html report; drop them ---
print("Cleaning up FastQC zip files...")
for f in config.FASTQC_DIR.iterdir():
    if f.suffix == ".zip":
        f.unlink()
        print(f" Removed: {f.name}")
