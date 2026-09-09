# Pipeline workflow

The pipeline processes one paired-end whole-genome metagenomic sample per run.

```text
Raw FASTQ
  -> Initial QC
  -> Adapter trimming
  -> Complexity filtering
  -> Human-host depletion
  -> Post-host QC
  -> Pair repair + deduplication
  -> Read normalization
  -> Final validation
  -> Kraken2 / MetaPhlAn / DeepARG / MEGAHIT
  -> CSV conversion
  -> Final analysis reports
```

Stage 08 is now part of the canonical data flow. If the dataset is already
below `PIPELINE_TARGET_DEPTH`, the reads are copied unchanged into the
normalized-data directory.

The pipeline intentionally keeps external tools and large reference
databases outside Git. Their locations are configured through environment
variables in `.env`.
