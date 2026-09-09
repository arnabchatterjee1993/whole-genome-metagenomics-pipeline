# Troubleshooting

## Required tool not found

Check the paths in `.env` and verify that the relevant executable exists.

## Exactly one paired FASTQ expected

The pipeline is currently designed for one sample per run. The input
directory must contain exactly:

```text
<sample>_R1.fastq.gz
<sample>_R2.fastq.gz
```

## Stage 10 fails

Stage 10 runs Kraken2, MetaPhlAn, DeepARG and MEGAHIT independently. Check
the tool's database path and environment. A summary is printed at the end;
the stage exits non-zero if one or more analyses fail.

## Resuming

After fixing an issue, resume from the affected stage:

```bash
python run_pipeline.py --from-stage 10
```
