"""
Stage 10: Integrated downstream analyses - Kraken2, MetaPhlAn, DeepARG, MEGAHIT.

NOTE: MetaPhlAn and DeepARG are run on R1 only, even though the pipeline is
paired-end throughout. That was true of the original script too. If you want
paired-mode profiling, pass a comma-separated "r1,r2" string to MetaPhlAn's
positional input and adjust the DeepARG call accordingly.

Each of the four tools runs in its own try/except so that one tool failing
(e.g. a missing database) doesn't stop the others from running.
"""

import os
import shutil
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import config
from utils import detect_sample, require_paths

config.ANALYSIS_DIR.mkdir(parents=True, exist_ok=True)
status = {}

# --- 1. Locate stage 09's gold-standard reads ---
SAMPLE = detect_sample(config.MASTER_DIR)
clean_r1 = config.GOLD_DIR / f"{SAMPLE}_clean_R1.fastq.gz"
clean_r2 = config.GOLD_DIR / f"{SAMPLE}_clean_R2.fastq.gz"
require_paths(clean_r1, clean_r2)

# Stale *.bowtie2out.txt files from a previous MetaPhlAn run will make the
# next run crash (MetaPhlAn treats them as cached, incompatible output).
for f in config.GOLD_DIR.glob("*.bowtie2out.txt"):
    f.unlink()

# --- 2. Kraken2: fast k-mer based taxonomic classification of every read ---
print("\n[1/4] Running Kraken2...")
kraken_report = config.ANALYSIS_DIR / f"{SAMPLE}.kraken.report"
kraken_out = config.ANALYSIS_DIR / f"{SAMPLE}.kraken.out"
try:
    subprocess.run([
        config.KRAKEN_EXE, "--db", config.KRAKEN_DB, "--threads", str(config.THREADS), "--paired",
        "--report", str(kraken_report), "--output", str(kraken_out),
        str(clean_r1), str(clean_r2),
    ], check=True)
    print("Kraken2 complete.")
    status["Kraken2"] = "PASS"
except (OSError, subprocess.CalledProcessError) as e:
    print(f"Kraken2 failed: {e}")
    status["Kraken2"] = "FAIL"

# --- 3. MetaPhlAn: marker-gene based species-level relative abundance ---
print("\n[2/4] Running MetaPhlAn...")
metaphlan_out = config.ANALYSIS_DIR / f"{SAMPLE}_metaphlan_profile.txt"
mp_env = os.environ.copy()
mp_env["PATH"] = f"{config.METAPHLAN_DIR}:" + mp_env["PATH"]
try:
    subprocess.run([
        f"{config.METAPHLAN_DIR}/metaphlan", str(clean_r1),
        "--input_type", "fastq",
        "--bowtie2db", config.METAPHLAN_BT2DB,
        "--index", config.METAPHLAN_INDEX,
        "--offline", "--nproc", str(config.THREADS),
        "-o", str(metaphlan_out),
    ], check=True, env=mp_env)
    print("MetaPhlAn complete.")
    status["MetaPhlAn"] = "PASS"
except (OSError, subprocess.CalledProcessError):
    print("MetaPhlAn failed. Ensure .bowtie2out.txt files were deleted.")

# --- 4. DeepARG: deep-learning based antibiotic resistance gene prediction ---
print("\n[3/4] Running DeepARG...")
arg_env = os.environ.copy()
arg_env["PATH"] = f"{config.DEEPARG_DIR}:" + arg_env["PATH"]
try:
    subprocess.run([
        f"{config.DEEPARG_DIR}/deeparg", "predict",
        "--input", str(clean_r1),
        "--output", str(config.ANALYSIS_DIR / f"{SAMPLE}_resistome"),
        "--type", "reads",
        "--data-path", config.DEEPARG_DB,
        "--model", "v2",
    ], check=True, env=arg_env)
    print("DeepARG complete.")
    status["DeepARG"] = "PASS"
except (OSError, subprocess.CalledProcessError) as e:
    print(f"DeepARG failed: {e}")
    status["DeepARG"] = "FAIL"

# --- 5. MEGAHIT: de novo assembly of reads into contigs ---
print("\n[4/4] Running MEGAHIT (Conservative Mode)...")
megahit_dir = config.ANALYSIS_DIR / f"{SAMPLE}_megahit_assembly"
if megahit_dir.exists():
    shutil.rmtree(megahit_dir)  # MEGAHIT refuses to write into an existing output dir
try:
    subprocess.run([
        "megahit", "-1", str(clean_r1), "-2", str(clean_r2),
        "-o", str(megahit_dir), "-t", str(config.THREADS),
        "--mem-flag", "1", "--no-mercy",
    ], check=True)
    print("MEGAHIT complete.")
except (OSError, subprocess.CalledProcessError) as e:
    print(f"MEGAHIT failed: {e}")
    status["MEGAHIT"] = "FAIL"

print("\n" + "=" * 40)
print("DOWNSTREAM ANALYSIS SUMMARY")
print("=" * 40)
for tool, result in status.items():
    print(f"{tool:<12}: {result}")
failed = [tool for tool, result in status.items() if result == "FAIL"]
print(f"Overall      : {'FAILED' if failed else 'SUCCESS'}")
print("=" * 40)

# Keep running the remaining tools when one fails, but signal failure to
# an orchestrator/CI job after all independent analyses have been attempted.
if failed:
    sys.exit(1)

