# PCOA prototype

A prototype project to process ENA 16S forward reads, denoise with Deblur, and generate beta-diversity + PCoA outputs.

## Tooling

- `fastq-dl` for downloading FASTQ data from ENA
- `deblur` for QC / denoising
- `biom-format` Python package for reading BIOM feature tables
- `scikit-bio` for distance metrics and PCoA analysis

## Environment setup

The `pcoa-prototype` conda environment covers Deblur, BIOM, scikit-bio, and plotting.
The UniFrac step additionally requires a [QIIME2 amplicon environment](https://docs.qiime2.org/2024.10/install/) with the [q2-greengenes2 plugin](https://github.com/biocore/q2-greengenes2) installed — any recent QIIME2 amplicon distribution will work.

Create the pcoa-prototype environment from the repository root:

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
deblur --help
biom --help
python -c "import skbio; print(skbio.__version__)"
```

The pipeline code also expects Python packages that are now listed in `environment.yml`, including `pandas`, `matplotlib`, and `pytest`.

## Forward-only pipeline

Run the pipeline in three steps:

```bash
python src/run_deblur.py --data-dir data
python src/build_table.py
python src/diversity.py
```

Use `--work-dir` and `--results-dir` when running multiple experiments side-by-side:

```bash
python src/run_deblur.py --data-dir data/demux_moving_picture --work-dir work/deblur_moving_picture --trim-length 120
python src/build_table.py --work-dir work/deblur_moving_picture --results-dir results/moving_picture
python src/diversity.py --results-dir results/moving_picture --metric braycurtis
```


### Inputs used

The pipeline recursively scans the supplied `--data-dir` and uses forward-read sample FASTQs matching either `*_1.fastq.gz` (ENA-style) or `*_R1_001.fastq.gz` (demux-export style). Reverse reads like `*_2.fastq.gz` are ignored. This works whether the FASTQs are directly under `data/` or nested under a study subdirectory.

It does not do paired-end merging or subsampling.

### What the scripts do

1. `run_deblur.py` recursively discovers forward reads matching `*_1.fastq.gz` or `*_R1_001.fastq.gz` under the input directory.
2. `run_deblur.py` calls `deblur workflow` directly on that directory.
3. The primary Deblur artifact is `work/deblur/workflow/all.biom`.
4. `build_table.py` reads `all.biom` via `biom.load_table` and writes a sample x feature TSV (`results/forward_only/feature_table.tsv`).
5. `build_table.py` normalizes demux-style sample IDs by removing lane/run suffixes (for example `L5S222_17_L001_R1_001` to `L5S222`).
6. `diversity.py` computes beta diversity (default: `braycurtis`) and runs PCoA from the feature table.
7. `diversity.py` writes tabular outputs and a PNG plot.

### Outputs

- Working outputs: `work/deblur/`
  - `workflow/` Deblur outputs including `all.biom`, `all.seqs.fa`, and workflow metadata/logs
- Final outputs: `results/forward_only/`
  - `feature_table.tsv`
  - `distance_matrix_braycurtis.tsv`
  - `pcoa_coordinates.tsv`
  - `pcoa_plot.png`

### File usage quick guide

- `work/.../workflow/all.biom`
  - Primary Deblur output (feature table in BIOM format).
  - Use this when you want to re-export tables or compare Deblur runs.
- `results/.../feature_table.tsv`
  - Sample x feature matrix used by `diversity.py`.
  - This is the direct input for scikit-bio beta-diversity/PCoA in this repo.
- `results/.../distance_matrix_braycurtis.tsv`
  - Pairwise sample dissimilarity matrix used to derive PCoA coordinates.
- `results/.../pcoa_coordinates.tsv`
  - Coordinates for each sample (PC1, PC2, ...); use this for custom plotting or metadata joins.
- `results/.../pcoa_plot.png`
  - Static quick-look visualization of the ordination.

Metadata note:

- Feature tables do not contain disease/group labels.
- Join `pcoa_coordinates.tsv` or `feature_table.tsv` to a metadata TSV by sample ID.

## Notes

This is the first-pass reproducible implementation and is intentionally minimal so future work can branch into paired-end workflows, subsampling comparisons, richer metadata integration, and alternative distance metrics.

Deblur defaults used by `run_deblur.py` are standardized as:

- `--trim-length 150`
- `--min-reads 0` to disable the final cross-sample abundance filter and keep per-sample denoising independent

## PRJEB44533 Subsampling

Use `src/subsample_fastq.py` to subsample forward reads (`*_1.fastq.gz` or `*_R1_001.fastq.gz`) with `seqkit sample2` in two-pass mode.

Expected PRJEB44533 layout:

```text
data/fastq_data/PRJEB44533/
  PRJEB44533-run-info.tsv
  full/
  subsample_50/
  subsample_25/
```

Example commands:

```bash
python src/subsample_fastq.py \
  --input-dir data/fastq_data/PRJEB44533/full \
  --output-dir data/fastq_data/PRJEB44533/subsample_50 \
  --percent 50 \
  --seed 11

python src/subsample_fastq.py \
  --input-dir data/fastq_data/PRJEB44533/full \
  --output-dir data/fastq_data/PRJEB44533/subsample_25 \
  --percent 25 \
  --seed 11
```

Implementation detail:

- `subsample_fastq.py` always calls `seqkit sample2` with `-2` (`--two-pass`) for stable fraction sampling on large FASTQ files.
