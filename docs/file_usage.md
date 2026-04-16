# File Usage Reference

This document explains how to use key pipeline files and artifacts.

## Input files

- `data/.../*.fastq.gz`
  - Supported forward-read patterns:
    - `*_1.fastq.gz` (ENA-style)
    - `*_R1_001.fastq.gz` (demux-export style)
  - Reverse reads are ignored by this pipeline.

## Deblur outputs (`run_deblur.py`)

- `work/<run_name>/workflow/all.biom`
  - Deblur feature table in BIOM format.
  - Main handoff artifact to `build_table.py`.
- `work/<run_name>/workflow/all.seqs.fa`
  - Representative feature sequences from Deblur.

## Table outputs (`build_table.py`)

- `results/<run_name>/feature_table.tsv`
  - Sample x feature matrix exported from `all.biom` via Python `biom.load_table`.
  - Demux sample IDs are normalized by removing lane/run suffixes:
    - `L5S222_17_L001_R1_001` -> `L5S222`

## Diversity / ordination outputs (`diversity.py`)

- `results/<run_name>/distance_matrix_braycurtis.tsv`
  - Sample x sample Bray-Curtis distance matrix.
- `results/<run_name>/pcoa_coordinates.tsv`
  - Ordination coordinates for each sample (PC axes).
- `results/<run_name>/pcoa_plot.png`
  - Quick static plot of PC1 vs PC2.

## Comparing to QIIME moving picture artifacts

- Compare this pipeline's `feature_table.tsv` to:
  - `data/moving_picture/table.qza` (post-Deblur table)
- Do not compare directly to:
  - `data/moving_picture/diversity-core-metrics-phylogenetic/rarefied_table.qza`
    when expecting unrarefied parity, because that table is depth-normalized and has fewer samples.

## Metadata usage

- Biological labels (body site, month, subject, intervention) come from metadata, not from feature tables.
- Use sample ID as join key:
  - moving picture metadata file: `data/moving_picture/sample-metadata.tsv`
