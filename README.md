# Shotgun Metagenomics Preprocessing & Analysis Pipeline

A 12-stage command-line pipeline that takes raw paired-end FASTQ reads from
a shotgun metagenomic sample through quality control, trimming, host
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

Run the stages in order:

```bash
for s in scripts/*.py; do python "$s" || break; done
```

or run them one at a time while you inspect intermediate output:

```bash
python scripts/01_prepare_raw_reads.py
python scripts/02_initial_qc.py
# ... etc.
```

Each script re-detects the sample name from the FASTQ pair in
`PIPELINE_MASTER_DIR`, so they can be run independently as long as the
previous stage's output already exists.

## Known gaps

Carried over from the original pipeline and not silently changed:

- **Stage 08's normalized/down-sampled reads aren't consumed downstream.**
  Stages 09–12 read from the pre-normalization "final clean" directory, not
  from stage 08's output. If you want normalized reads to feed the
  downstream analysis, point stage 09's `final_r1/final_r2` at
  `config.NORM_DIR` instead.
- **MetaPhlAn and DeepARG (stage 10) are run on R1 only**, despite the
  pipeline being paired-end throughout. This may be intentional (a common
  simplification for these tools) but is worth confirming against your
  use case.

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
