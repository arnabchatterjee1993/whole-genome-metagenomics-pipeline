"""
Central configuration for the metagenomics pipeline.

Everything that used to be a hardcoded, machine-specific absolute path
now lives here and can be overridden with environment variables, so the
pipeline runs on anyone's machine without editing the scripts themselves.

Set these in a .env file (see .env.example) or export them in your shell
before running any script in scripts/.
"""

import os
from pathlib import Path


def _env_path(name: str, default: str) -> Path:
    return Path(os.environ.get(name, default)).expanduser().resolve()


# -------------------------------------------------------------------
# Top-level working directories
# -------------------------------------------------------------------
# MASTER_DIR: where raw input FASTQs live and all pipeline outputs are written
MASTER_DIR = _env_path("PIPELINE_MASTER_DIR", "./pipeline_run")

# TOOLS_DIR: where third-party binaries (BBTools, Bowtie2, Samtools, etc.) live
TOOLS_DIR = _env_path("PIPELINE_TOOLS_DIR", "./tools")

# -------------------------------------------------------------------
# Stage directories (all relative to MASTER_DIR)
# -------------------------------------------------------------------
RAW_DIR = MASTER_DIR / "01_Raw_Reads"
FASTQC_DIR = MASTER_DIR / "02_FastQC_Reports"
MULTIQC_DIR = MASTER_DIR / "03_MultiQC_Reports"
TRIMMED_DIR = MASTER_DIR / "04_Trimmed_Reads"
FILTERED_DIR = MASTER_DIR / "05_Filtered_Reads"
HOST_DEPLETED_DIR = MASTER_DIR / "06_Host_Depleted"
POST_QC_DIR = MASTER_DIR / "07_Post_Host_QC"
FINAL_MULTIQC_DIR = MASTER_DIR / "08_Final_MultiQC"
CLEAN_BACKUP_DIR = MASTER_DIR / "09_Clean_Reads"          # renamed from "01_Clean_Reads" (was confusing)
DEDUP_DIR = MASTER_DIR / "10_Deduplicated_Data"
FINAL_CLEAN_DIR = MASTER_DIR / "11_Final_Clean_Microbial_Data"
NORM_DIR = MASTER_DIR / "12_Normalized_Data"
GOLD_DIR = MASTER_DIR / "13_Gold_Standard_Dataset"
ANALYSIS_DIR = MASTER_DIR / "14_Analysis_Results"
RAW_REPORTS_DIR = ANALYSIS_DIR / "raw_reports"

# NOTE: In the original script, stage 12 (read normalization / down-sampling)
# produced output that stages 13-14 never consumed -- they read straight from
# the pre-normalization "final clean" directory instead. That gap is preserved
# here (not silently "fixed") because it's not clear which was intended.
# See README.md "Known gaps" section before you rely on normalized output.

# -------------------------------------------------------------------
# Third-party binaries / scripts
# -------------------------------------------------------------------
BBMAP_DIR = TOOLS_DIR / "bbmap"
BBDUK = BBMAP_DIR / "bbduk.sh"
REPAIR_SH = BBMAP_DIR / "repair.sh"
CLUMPIFY = BBMAP_DIR / "clumpify.sh"

BOWTIE2_BIN = _env_path("PIPELINE_BOWTIE2_BIN", str(TOOLS_DIR / "bowtie2" / "bowtie2"))
SAMTOOLS_BIN = _env_path("PIPELINE_SAMTOOLS_BIN", str(TOOLS_DIR / "samtools" / "samtools"))
HUMAN_INDEX_DIR = _env_path("PIPELINE_HUMAN_INDEX_DIR", str(TOOLS_DIR / "bowtie2" / "GRCh38_noalt_as"))
HUMAN_INDEX_PREFIX = "GRCh38_noalt_as"

KRAKEN_EXE = os.environ.get("PIPELINE_KRAKEN_EXE", str(TOOLS_DIR / "kraken2" / "kraken2"))
KRAKEN_DB = os.environ.get("PIPELINE_KRAKEN_DB", str(TOOLS_DIR / "kraken2" / "kraken2_db"))

METAPHLAN_DIR = os.environ.get("PIPELINE_METAPHLAN_DIR", str(TOOLS_DIR / "metaphlan" / "bin"))
METAPHLAN_BT2DB = os.environ.get("PIPELINE_METAPHLAN_BT2DB", str(TOOLS_DIR / "metaphlan" / "metaphlan_bt2db"))
METAPHLAN_INDEX = os.environ.get("PIPELINE_METAPHLAN_INDEX", "mpa_vJun23_CHOCOPhlAnSGB_202307")

DEEPARG_DIR = os.environ.get("PIPELINE_DEEPARG_DIR", str(TOOLS_DIR / "deeparg" / "bin"))
DEEPARG_DB = os.environ.get("PIPELINE_DEEPARG_DB", str(TOOLS_DIR / "deeparg_db"))

# -------------------------------------------------------------------
# Runtime parameters
# -------------------------------------------------------------------
THREADS = int(os.environ.get("PIPELINE_THREADS", "4"))
TARGET_READ_DEPTH = int(os.environ.get("PIPELINE_TARGET_DEPTH", "1000000"))
