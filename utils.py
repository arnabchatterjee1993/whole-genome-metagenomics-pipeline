"""
Shared helpers used across the pipeline stage scripts.

The original single-file version relied on a `SAMPLE` variable set once at
the top and shared via Python's global scope across all 12 "parts" because
they were literally one script. Now that each stage is its own script/process,
`detect_sample()` re-derives it the same way every time it's needed.
"""

import shutil
import subprocess
import sys
from pathlib import Path


def detect_sample(raw_dir: Path) -> str:
    """Find the single paired-end FASTQ set in raw_dir and return its sample name.

    Expects exactly one `<sample>_R1.fastq.gz` and one `<sample>_R2.fastq.gz`.
    Exits with an error message if that's not the case.
    """
    r1_candidates = sorted(raw_dir.glob("*_R1.fastq.gz"))
    r2_candidates = sorted(raw_dir.glob("*_R2.fastq.gz"))

    if len(r1_candidates) != 1 or len(r2_candidates) != 1:
        sys.exit(
            f"ERROR: Exactly ONE paired FASTQ.GZ expected in {raw_dir}:\n"
            "  <sample>_R1.fastq.gz and <sample>_R2.fastq.gz"
        )

    r1_gz, r2_gz = r1_candidates[0], r2_candidates[0]
    sample = r1_gz.name.replace("_R1.fastq.gz", "")

    if r2_gz.name != f"{sample}_R2.fastq.gz":
        sys.exit("ERROR: R1 and R2 sample names do not match")

    return sample


def require_tools(*tool_names: str) -> None:
    """Exit with a clear message if any of the given tools aren't on PATH."""
    missing = [t for t in tool_names if shutil.which(t) is None]
    if missing:
        sys.exit(f"ERROR: required tool(s) not found on PATH: {', '.join(missing)}")


def require_paths(*paths: Path) -> None:
    """Exit with a clear message if any of the given paths don't exist."""
    missing = [str(p) for p in paths if not p.exists()]
    if missing:
        sys.exit("ERROR: required file(s)/binary(ies) not found:\n" + "\n".join(missing))


def count_reads_gz(fq_path: Path) -> int:
    """Count reads in a gzipped FASTQ file (4 lines per read)."""
    result = subprocess.run(["zcat", str(fq_path)], stdout=subprocess.PIPE, check=True)
    return len(result.stdout.decode().splitlines()) // 4


def read_metaphlan_profile(path: Path):
    """Read a MetaPhlAn profile TSV into a DataFrame, with correct headers.

    Bug this works around: MetaPhlAn profiles start with a metadata line
    AND a header line that both begin with '#' (e.g. `#clade_name\\t...`).
    Naively reading with `pd.read_csv(path, sep='\\t', comment='#')` strips
    both, so pandas is left using the first *species* row as the column
    header -- silently dropping that species and mislabeling every column.
    This reads the header explicitly instead of comment-stripping it away.
    """
    import pandas as pd
    from io import StringIO

    with open(path) as fh:
        lines = fh.readlines()

    header_idx = next((i for i, l in enumerate(lines) if l.startswith("#clade_name")), None)
    if header_idx is None:
        raise ValueError(f"No '#clade_name' header found in {path}")

    header_line = lines[header_idx].lstrip("#")
    data_lines = [l for l in lines[header_idx + 1:] if not l.startswith("#")]
    return pd.read_csv(StringIO(header_line + "".join(data_lines)), sep="\t")
