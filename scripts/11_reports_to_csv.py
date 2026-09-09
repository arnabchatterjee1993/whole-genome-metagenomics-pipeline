"""Stage 11: Convert raw analysis outputs (Kraken2, MetaPhlAn, DeepARG, MEGAHIT) to CSV."""

import shutil
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import config
from utils import detect_sample, read_metaphlan_profile

config.RAW_REPORTS_DIR.mkdir(parents=True, exist_ok=True)
SAMPLE = detect_sample(config.MASTER_DIR)

print("Converting analysis reports to CSV format...")

# --- Kraken2 -> CSV: report has no header row, so column names are supplied here ---
kraken_report = config.ANALYSIS_DIR / f"{SAMPLE}.kraken.report"
if kraken_report.exists():
    cols = ["percentage", "reads_at_taxon", "reads_directly_at_taxon", "rank", "taxid", "name"]
    df = pd.read_csv(kraken_report, sep="\t", header=None, names=cols)
    df["name"] = df["name"].str.strip()
    out_csv = config.ANALYSIS_DIR / f"{SAMPLE}_taxonomic_classification.csv"
    df.to_csv(out_csv, index=False)
    print(f"Created: {out_csv.name}")
else:
    print("Kraken2 report not found.")

# --- MetaPhlAn -> CSV ---
# Uses read_metaphlan_profile() (utils.py) rather than a plain pd.read_csv(comment='#')
# call, which would silently corrupt the output -- see that function's docstring.
metaphlan_file = config.ANALYSIS_DIR / f"{SAMPLE}_metaphlan_profile.txt"
if metaphlan_file.exists():
    df_meta = read_metaphlan_profile(metaphlan_file)
    out_csv = config.ANALYSIS_DIR / f"{SAMPLE}_species_abundance.csv"
    df_meta.to_csv(out_csv, index=False)
    print(f"Created: {out_csv.name}")
else:
    print("MetaPhlAn profile not found.")

# --- DeepARG -> CSV ---
deeparg_file = config.ANALYSIS_DIR / f"{SAMPLE}_resistome.mapping.ARG"
if deeparg_file.exists():
    df_arg = pd.read_csv(deeparg_file, sep="\t")
    out_csv = config.ANALYSIS_DIR / f"{SAMPLE}_antibiotic_resistance.csv"
    df_arg.to_csv(out_csv, index=False)
    print(f"Created: {out_csv.name}")
else:
    print("DeepARG ARG mapping not found.")

# --- MEGAHIT: just confirm the assembly exists; the .fa is used as-is, not converted ---
contigs_file = config.ANALYSIS_DIR / f"{SAMPLE}_megahit_assembly" / "final.contigs.fa"
if contigs_file.exists():
    print(f"Assembly found: {contigs_file.name}")
else:
    print("MEGAHIT assembly file not found.")

# --- Organize: move the raw/intermediate text reports out of the way, into raw_reports/ ---
print("\nOrganizing raw files...")
files_to_move = [
    f"{SAMPLE}.kraken.report",
    f"{SAMPLE}.kraken.out",
    f"{SAMPLE}_metaphlan_profile.txt",
    f"{SAMPLE}_resistome.mapping.ARG",
    f"{SAMPLE}_resistome.align.daa.tsv",
    f"{SAMPLE}_resistome.potential.ARG",
]
for filename in files_to_move:
    f = config.ANALYSIS_DIR / filename
    if f.exists():
        shutil.move(str(f), str(config.RAW_REPORTS_DIR / filename))

print("\nFINAL CSV DATASETS READY:")
csv_files = list(config.ANALYSIS_DIR.glob(f"{SAMPLE}_*.csv"))
if csv_files:
    for f in csv_files:
        print(f"{f.name} ({f.stat().st_size / 1024:.2f} KB)")
else:
    print("No CSV files found.")
