"""Stage 6: Post-host-removal quality assessment, compression, and backup."""

import shutil
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import config
from utils import detect_sample, require_paths

for d in [config.POST_QC_DIR, config.FINAL_MULTIQC_DIR, config.CLEAN_BACKUP_DIR]:
    d.mkdir(parents=True, exist_ok=True)

# --- 1. Locate stage 05's host-depleted FASTQs ---
SAMPLE = detect_sample(config.MASTER_DIR)
r1 = config.HOST_DEPLETED_DIR / f"{SAMPLE}_R1_depleted.fastq"
r2 = config.HOST_DEPLETED_DIR / f"{SAMPLE}_R2_depleted.fastq"
require_paths(r1, r2)
print("Host-depleted FASTQs verified")

# --- 2. Re-run FastQC now that host reads are gone, to see the "clean" quality profile ---
print("Running FastQC...")
subprocess.run(["fastqc", str(r1), str(r2), "-o", str(config.POST_QC_DIR)], check=True)
if not any(config.POST_QC_DIR.glob("*_fastqc.html")):
    sys.exit("ERROR: FastQC reports missing")
print("Post-removal FastQC completed")

# --- 3. Combine the raw-read QC (stage 02) and post-host-removal QC into one report ---
print("Running MultiQC...")
qc_inputs = [str(d) for d in [config.FASTQC_DIR, config.POST_QC_DIR] if d.exists()]
subprocess.run(["multiqc", *qc_inputs, "-o", str(config.FINAL_MULTIQC_DIR)], check=True)
final_report = config.FINAL_MULTIQC_DIR / "multiqc_report.html"
if not final_report.exists():
    sys.exit("ERROR: Final MultiQC report not generated")
print(f"MultiQC report generated: {final_report}")

# --- 4. Compress the clean FASTQs for storage/backup ---
r1_gz = config.CLEAN_BACKUP_DIR / f"{SAMPLE}_clean_R1.fastq.gz"
r2_gz = config.CLEAN_BACKUP_DIR / f"{SAMPLE}_clean_R2.fastq.gz"
print("Compressing host-depleted FASTQs with pigz...")
for fq, out in zip([r1, r2], [r1_gz, r2_gz]):
    with open(out, "wb") as fh:
        subprocess.run(["pigz", "-c", str(fq)], stdout=fh, check=True)
for gz in [r1_gz, r2_gz]:
    if not gz.exists():
        sys.exit(f"ERROR: Compressed FASTQ missing: {gz}")
print("Compression successful")

# --- 5. Keep a copy of the final QC report alongside the clean reads ---
shutil.copy2(final_report, config.CLEAN_BACKUP_DIR / "Final_QC_Report.html")
print("MultiQC report backed up successfully")

# --- 6. FastQC's .zip files just duplicate the .html report; drop them ---
for f in config.POST_QC_DIR.iterdir():
    if f.suffix == ".zip":
        f.unlink()

print("SUCCESS! Cleaned, compressed data and final QC report are safely stored.")
