"""Stage 12: Final analysis report formatting.

Turns the CSV outputs from stage 11 into three tidy, analysis-ready tables:
  1. MetaPhlAn taxonomy, split into Domain/Phylum/.../Species columns
  2. A wide Kraken2 species-by-sample abundance matrix
  3. DeepARG's resistome hits, reformatted into a fixed column layout

Note: this reuses utils.read_metaphlan_profile(), which works around a real
bug in the original script -- see that function's docstring, or the README's
"Bug fix vs. the original script" section, for details.
"""

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import config
from utils import detect_sample, read_metaphlan_profile

SAMPLE = detect_sample(config.MASTER_DIR)

raw_reports_dir = config.RAW_REPORTS_DIR
raw_reports_dir.mkdir(exist_ok=True)

final_reports_dir = config.ANALYSIS_DIR / f"{SAMPLE}_Final_Reports_For_Analysis"
final_reports_dir.mkdir(parents=True, exist_ok=True)
print(f"Final reports directory: {final_reports_dir}")

# Inputs were moved here by stage 11's "Organize" step
metaphlan_input = raw_reports_dir / f"{SAMPLE}_metaphlan_profile.txt"
tax_path = config.ANALYSIS_DIR / f"{SAMPLE}_taxonomic_classification.csv"
resistome_input = raw_reports_dir / f"{SAMPLE}_resistome.align.daa.tsv"

metaphlan_output = final_reports_dir / f"{SAMPLE}_metaphlan_profile_parsed.csv"
matrix_output = final_reports_dir / f"{SAMPLE}_taxonomy_final_report.csv"
resistome_output_base = final_reports_dir / f"{SAMPLE}_resistome_final_report"


def parse_metaphlan_profile(input_file, output_file):
    """Split MetaPhlAn's pipe-delimited lineage string into taxonomy columns.

    e.g. "k__Bacteria|p__Firmicutes|...|s__Lactobacillus_reuteri"
      -> Domain=d__Bacteria, Phylum=p__Firmicutes, ..., Species=s__Lactobacillus_reuteri
    Only rows with a species-level assignment (containing "s__") are kept.
    """
    print("\n[Step 1] Parsing MetaPhlAn species taxonomy...")
    if not input_file.exists():
        print(f"Skipping: {input_file} not found.")
        return

    df = read_metaphlan_profile(input_file)
    tax_col = "clade_name" if "clade_name" in df.columns else df.columns[0]

    rows = []
    for clade in df[tax_col].astype(str):
        if "s__" not in clade:
            continue
        tax = {"Domain": "", "Phylum": "", "Class": "", "Order": "", "Family": "", "Genus": "", "Species": ""}
        for lvl in clade.split("|"):
            if lvl.startswith("k__"): tax["Domain"] = "d__" + lvl.replace("k__", "")
            elif lvl.startswith("p__"): tax["Phylum"] = lvl
            elif lvl.startswith("c__"): tax["Class"] = lvl
            elif lvl.startswith("o__"): tax["Order"] = lvl
            elif lvl.startswith("f__"): tax["Family"] = lvl
            elif lvl.startswith("g__"): tax["Genus"] = lvl
            elif lvl.startswith("s__"): tax["Species"] = lvl
        rows.append(tax)

    df_parsed = pd.DataFrame(rows)
    df_parsed.to_csv(output_file, index=False)
    print(f"Saved MetaPhlAn parsed taxonomy: {output_file.name} ({len(df_parsed)} species)")


def generate_wide_abundance_matrix(tax_file, output_path):
    """Pivot Kraken2's species rows into one wide row: species become columns,
    read counts become values -- the standard "abundance matrix" shape most
    downstream ecology/stats tools (e.g. vegan, scikit-bio) expect."""
    print("\n[Step 2] Generating Wide Abundance Matrix (Kraken2)...")
    if not tax_file.exists():
        print(f"Skipping: {tax_file} not found.")
        return

    tax_df = pd.read_csv(tax_file)
    # rank 'S' = species, 'S1' = sub-species/strain -- both count as species-level here
    species_tax = tax_df[tax_df["rank"].isin(["S", "S1"])].copy()
    species_tax["matrix_header"] = "s__" + species_tax["name"].astype(str).str.strip()

    matrix_df = species_tax.pivot_table(
        index=None, columns="matrix_header", values="reads_at_taxon", aggfunc="sum"
    )
    matrix_df.index = [SAMPLE]
    matrix_df.fillna(0).to_csv(output_path)
    print(f"Saved taxonomy matrix with {matrix_df.shape[1]} species")


def process_resistome_to_custom_format(input_path, output_path_base):
    """Reformat DeepARG's raw alignment output into a fixed, analysis-friendly
    column layout (strand, gene, coverage, resistance class, etc.), saved as
    both CSV and TSV."""
    print("\n[Step 3] Processing Resistome data...")
    if not input_path.exists():
        print(f"Skipping: {input_path} not found.")
        return

    headers = ["SEQUENCE", "Subject", "%IDENTITY", "Length", "Mismatches", "GAPS",
               "Q_start", "Q_end", "START", "END", "E_value", "Bit_score"]
    df = pd.read_csv(input_path, sep="\t", names=headers)

    df["#FILE"] = input_path.name
    df["STRAND"] = df.apply(lambda x: "+" if x["START"] < x["END"] else "-", axis=1)

    # DeepARG packs gene/resistance-class/product into the pipe-delimited Subject field
    split_sub = df["Subject"].str.split("|", expand=True)
    df["ACCESSION"] = split_sub[0].str.split(":").str[-2]
    df["GENE"] = split_sub[2]
    df["RESISTANCE"] = split_sub[3].str.upper()
    df["PRODUCT"] = split_sub[4]
    df["DATABASE"] = "ncbi"

    df["COVERAGE"] = "1-" + df["Length"].astype(str) + "/" + df["Length"].astype(str)
    df["COVERAGE_MAP"] = "==============="
    df["%COVERAGE"] = 100.00

    final_cols = ["#FILE", "SEQUENCE", "START", "END", "STRAND", "GENE",
                  "COVERAGE", "COVERAGE_MAP", "GAPS", "%COVERAGE",
                  "%IDENTITY", "DATABASE", "ACCESSION", "PRODUCT", "RESISTANCE"]
    final_df = df[final_cols]

    final_df.to_csv(f"{output_path_base}.csv", index=False)
    final_df.to_csv(f"{output_path_base}.tsv", index=False, sep="\t")
    print(f"Saved resistome formatted reports: {output_path_base.name}.csv / .tsv")


if __name__ == "__main__":
    parse_metaphlan_profile(metaphlan_input, metaphlan_output)
    generate_wide_abundance_matrix(tax_path, matrix_output)
    process_resistome_to_custom_format(resistome_input, resistome_output_base)
    print("\nSUCCESS: STAGE 12 FINAL REPORTS READY FOR ANALYSIS")
