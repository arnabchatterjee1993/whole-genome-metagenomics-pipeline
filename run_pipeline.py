#!/usr/bin/env python3
"""Run the complete whole-genome metagenomics pipeline in stage order."""

import argparse
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPTS = ROOT / "scripts"

STAGES = [
    SCRIPTS / "01_prepare_raw_reads.py",
    SCRIPTS / "02_initial_qc.py",
    SCRIPTS / "03_adapter_trimming.py",
    SCRIPTS / "04_complexity_filtering.py",
    SCRIPTS / "05_host_depletion.py",
    SCRIPTS / "06_post_host_qc.py",
    SCRIPTS / "07_deduplication.py",
    SCRIPTS / "08_read_normalization.py",
    SCRIPTS / "09_dataset_validation.py",
    SCRIPTS / "10_downstream_analysis.py",
    SCRIPTS / "11_reports_to_csv.py",
    SCRIPTS / "12_final_report_formatting.py",
]


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Run the whole-genome metagenomics pipeline."
    )
    parser.add_argument(
        "--from-stage",
        type=int,
        choices=range(1, 13),
        default=1,
        help="Start at this stage number (default: 1).",
    )
    args = parser.parse_args()

    for stage in STAGES[args.from_stage - 1:]:
        print("\n" + "=" * 72)
        print(f"RUNNING {stage.name}")
        print("=" * 72)
        result = subprocess.run([sys.executable, str(stage)], cwd=ROOT)
        if result.returncode != 0:
            print(f"\nPIPELINE STOPPED: {stage.name} exited with code {result.returncode}.")
            return result.returncode

    print("\n" + "=" * 72)
    print("PIPELINE COMPLETED SUCCESSFULLY")
    print("=" * 72)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
