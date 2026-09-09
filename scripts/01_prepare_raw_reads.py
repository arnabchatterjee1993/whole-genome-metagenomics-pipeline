"""
Stage 1: Locate the input FASTQ pair, decompress, checksum, and stage it.

Input:  <sample>_R1.fastq.gz and <sample>_R2.fastq.gz placed directly in
        config.MASTER_DIR (the pipeline's root working directory).
Output: decompressed + verified copies in config.RAW_DIR.
"""

import shutil
import subprocess
import sys
from pathlib import Path

# Make the repo root importable so `import config` / `from utils import ...`
# work no matter where this script is invoked from.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import config
from utils import detect_sample, require_tools

config.RAW_DIR.mkdir(parents=True, exist_ok=True)
require_tools("pigz", "md5sum")

# --- 1. Auto-detect the sample from the FASTQ pair in MASTER_DIR ---
SAMPLE = detect_sample(config.MASTER_DIR)
print(f"Detected sample: {SAMPLE}")

r1_gz = config.MASTER_DIR / f"{SAMPLE}_R1.fastq.gz"
r2_gz = config.MASTER_DIR / f"{SAMPLE}_R2.fastq.gz"

r1_src = config.MASTER_DIR / f"{SAMPLE}_R1.fastq"
r2_src = config.MASTER_DIR / f"{SAMPLE}_R2.fastq"

# --- 2. Decompress to plain FASTQ (skip if already done) ---
for gz, out in [(r1_gz, r1_src), (r2_gz, r2_src)]:
    if not out.exists():
        print(f"Decompressing {gz.name} -> {out.name}")
        subprocess.run(["pigz", "-dc", str(gz)], stdout=open(out, "wb"), check=True)

missing = [str(f) for f in [r1_src, r2_src] if not f.exists()]
if missing:
    sys.exit("ERROR: Missing FASTQ files after decompression:\n" + "\n".join(missing))
print("Input FASTQ files found")


# --- 3. Verify integrity with MD5 checksums (printed for the run log / provenance) ---
def md5sum(file_path):
    result = subprocess.run(["md5sum", str(file_path)], capture_output=True, text=True, check=True)
    return result.stdout.strip()


print("MD5 checksums:")
r1_md5 = md5sum(r1_src)
r2_md5 = md5sum(r2_src)
print(r1_md5)
print(r2_md5)

# Save checksums alongside the staged raw reads for reproducibility.
checksum_file = config.RAW_DIR / "checksums.md5"
checksum_file.write_text(
    f"{r1_md5.split()[0]}  {r1_src.name}\n"
    f"{r2_md5.split()[0]}  {r2_src.name}\n"
)
print(f"Checksums saved to: {checksum_file}")

# --- 4. Copy the decompressed FASTQs into RAW_DIR, the input for stage 02 ---
r1_dest = config.RAW_DIR / r1_src.name
r2_dest = config.RAW_DIR / r2_src.name
for src, dest in [(r1_src, r1_dest), (r2_src, r2_dest)]:
    if not dest.exists():
        shutil.copy2(src, dest)

print(f"Raw FASTQs ready in: {config.RAW_DIR}")
for f in config.RAW_DIR.iterdir():
    print(f" - {f.name}")
