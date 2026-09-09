# Installation

## 1. Python dependencies

```bash
python -m pip install -r requirements.txt
```

## 2. External tools

Install the bioinformatics tools listed in the main README and make their
paths available through `.env` or environment variables.

## 3. Configure paths

Copy the example configuration:

```bash
cp .env.example .env
```

Edit the paths for the local machine. Do not commit `.env`.

## 4. Input

Place exactly one paired-end sample in `PIPELINE_MASTER_DIR`:

```text
<sample>_R1.fastq.gz
<sample>_R2.fastq.gz
```

The sample name is detected automatically.

## 5. Run

```bash
python run_pipeline.py
```

To resume from a later stage:

```bash
python run_pipeline.py --from-stage 5
```
