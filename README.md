# PCOA prototype

A prototype project to process ENA 16S forward reads, denoise with Deblur, and generate beta-diversity + PCoA outputs.

## Tooling

- `fastq-dl` for downloading FASTQ data from ENA
- `seqkit` for read utilities
- `deblur` for QC / denoising
- `scikit-bio` for distance metrics and PCoA analysis

## Environment setup

The repository's `environment.yml` is exported from the solved Conda environment, not maintained by hand.

Create the environment from the repository root:

```bash
conda env create -f environment.yml
conda activate pcoa-prototype
```

If you update the environment, regenerate the file from the installed env with:

```bash
conda env export -n pcoa-prototype --from-history | sed '/^prefix: /d' > environment.yml
```

Check that the main tools resolve correctly:

```bash
fastq-dl --help
seqkit version
deblur --help
python -c "import skbio; print(skbio.__version__)"
```

## First end-to-end forward-only pipeline

Run the forward-only PRJNA825639 ordination pipeline:

```bash
python src/run_deblur_pcoa.py
```

### Inputs used

The pipeline intentionally uses only these per-run forward reads under `data/`:

- `SRR27336825_1.fastq.gz`
- `SRR27336826_1.fastq.gz`
- `SRR27336827_1.fastq.gz`
- `SRR27336828_1.fastq.gz`
- `SRR27336829_1.fastq.gz`
- `SRR27336830_1.fastq.gz`
- `SRR27336831_1.fastq.gz`
- `SRR27336832_1.fastq.gz`

It ignores the merged `SAMN27531837.fastq.gz` file and does not do paired-end merging or subsampling.

### What the script does

1. Validates expected input files in `data/`.
2. Converts each run FASTQ to trimmed (250 bp) FASTA.
3. Runs `deblur dereplicate` and `deblur deblur-seqs` per run.
4. Builds one merged sample-by-feature table from Deblur-cleaned outputs.
5. Computes beta diversity (default: `braycurtis`).
6. Runs PCoA.
7. Writes tabular results and a PNG plot.

### Outputs

- Working outputs: `work/deblur/`
  - per-run trimmed FASTA, dereplicated FASTA, and Deblur-clean FASTA files
- Final outputs: `results/forward_only/`
  - `feature_table.tsv`
  - `distance_matrix_braycurtis.tsv`
  - `pcoa_coordinates.tsv`
  - `pcoa_plot.png`

## Notes

This is the first-pass reproducible implementation and is intentionally minimal so future work can branch into paired-end workflows, subsampling comparisons, richer metadata integration, and alternative distance metrics.
