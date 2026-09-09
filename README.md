# Whole-genome Metagenomics Preprocessing & Analysis Pipeline

A 12-stage command-line pipeline that takes raw paired-end FASTQ reads from
a Whole-genome metagenomic sample through quality control, trimming, host
(human) DNA depletion, deduplication, and into taxonomic classification,
species abundance profiling, antibiotic-resistance gene detection, and
metagenomic assembly — with tidy CSV reports at the end.

```
Raw FASTQ
   │
   ▼
01  Prepare raw reads (decompress, checksum)
02  Initial QC (FastQC + MultiQC)
03  Adapter trimming (fastp)
04  Complexity filtering (BBDuk)
05  Host DNA depletion (Bowtie2 + Samtools)
06  Post-depletion QC & compression
07  Deduplication (BBTools repair + Clumpify)
08  Read normalization / down-sampling (seqtk)      *see "Known gaps"
09  Dataset validation & gold-standard packaging
10  Downstream analysis: Kraken2, MetaPhlAn, DeepARG, MEGAHIT
11  Raw outputs → CSV
12  Final report formatting (taxonomy, abundance matrix, resistome)
```

Note: the numbered output folders inside `pipeline_run/` (`01_Raw_Reads`
... `14_Analysis_Results`) don't map 1:1 to the 12 script numbers above —
a few stages produce more than one folder (e.g. stage 06 produces three:
`07_Post_Host_QC`, `08_Final_MultiQC`, `09_Clean_Reads`). The folder
numbers just reflect the order data flows through the pipeline; see
`config.py` for the full stage→folder mapping.

## Requirements

Python packages:

```
pip install -r requirements.txt
```

External bioinformatics tools (install separately, not bundled here).
"Invoked in stage" means that stage runs the tool's binary; stages 11–12
only parse the text/report files those tools already produced, using
pandas — you don't need Kraken2/MetaPhlAn/DeepARG installed just to run
stages 11–12 on existing output.

| Tool | Invoked in stage |
|---|---|
| [pigz](https://zlib.net/pigz/) | 01, 06, 08 |
| [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) | 02, 06 |
| [MultiQC](https://multiqc.info/) | 02, 06 |
| [fastp](https://github.com/OpenGene/fastp) | 03 |
| [BBTools](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/) (bbduk, repair, clumpify) | 04, 07 |
| [Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/) | 05 |
| [Samtools](https://www.htslib.org/) | 05 |
| [seqtk](https://github.com/lh3/seqtk) | 08 |
| [Kraken2](https://ccb.jhu.edu/software/kraken2/) + a Kraken2 database | 10 |
| [MetaPhlAn](https://github.com/biobakery/MetaPhlAn) + its marker database | 10 |
| [DeepARG](https://github.com/gaarangoa/deeparg) + its database | 10 |
| [MEGAHIT](https://github.com/voutcn/megahit) | 10 |

You'll also need a Bowtie2 index for the host genome (e.g. GRCh38) for
stage 05. These tools and databases are large (many GB) and are
intentionally not part of this repo.

## Setup

1. Clone the repo and install Python deps:
   ```bash
   git clone https://github.com/<your-username>/<repo-name>.git
   cd <repo-name>
   pip install -r requirements.txt
   ```
2. Install the external tools above and note their install paths.
3. Copy `.env.example` to `.env`, edit the paths for your machine, and
   export them (or use a tool like `python-dotenv` / `direnv`):
   ```bash
   cp .env.example .env
   export $(cat .env | xargs)
   ```
   All paths default to sensible relative locations (`./tools/...`) if you
   don't set anything.
4. Place your input reads directly in `PIPELINE_MASTER_DIR` (default
   `./pipeline_run/`), named `<sample>_R1.fastq.gz` and
   `<sample>_R2.fastq.gz`. Exactly one pair per run — the sample name is
   auto-detected from the filenames.

## Usage

Run the complete pipeline with the launcher:

```bash
python run_pipeline.py
```

To resume from a particular stage:

```bash
python run_pipeline.py --from-stage 5
```

You can also run stages individually while inspecting intermediate output:

```bash
python scripts/01_prepare_raw_reads.py
python scripts/02_initial_qc.py
# ... etc.
```

Each script re-detects the sample name from the FASTQ pair in
`PIPELINE_MASTER_DIR`, so they can be run independently as long as the
previous stage's output already exists.

## Design notes

- The pipeline currently processes exactly one paired-end sample per run.
- Stage 08 is part of the canonical data flow. It creates the normalized
  dataset consumed by stages 09–12. If the read count is already below the
  target depth, the reads are copied unchanged so the downstream stages
  always have a defined input.
- MetaPhlAn and DeepARG (stage 10) are currently run on R1 only. Confirm this
  choice against the requirements of your study and the versions of those
  tools that you validate.
- Stage 10 runs its four analyses independently, reports PASS/FAIL status,
  and exits non-zero if one or more analyses fail.
- The repository intentionally excludes FASTQ/BAM/SAM files, third-party
  binaries, and large reference databases. Configure their locations in
  `.env` rather than committing machine-specific paths.

See `docs/` for installation, workflow, software-version tracking, and
troubleshooting notes.

## Bug fix vs. the original script

Both stage 11 and stage 12 parsed MetaPhlAn's profile with
`pd.read_csv(file, sep='\t', comment='#')`. MetaPhlAn's real header row
(`#clade_name\t...`) also starts with `#`, so pandas stripped it along with
the metadata line above it and silently promoted the first species row into
the column header — permanently dropping that species from every report
(confirmed by running it against a synthetic MetaPhlAn-format file: 2
species in, 1 species out). Fixed once in `utils.read_metaphlan_profile()`,
which locates the header line explicitly instead of comment-stripping it
away, and both stages now call that shared function.

## License

MIT — see [LICENSE](LICENSE).
