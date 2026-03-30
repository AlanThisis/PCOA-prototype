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
fastp --help
seqkit version
deblur --help
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

If you are not already inside an activated Conda shell, run it explicitly with the environment interpreter:

```bash
/opt/conda/envs/pcoa-prototype/bin/python src/run_deblur.py --data-dir data
/opt/conda/envs/pcoa-prototype/bin/python src/build_table.py
/opt/conda/envs/pcoa-prototype/bin/python src/diversity.py
```

### Inputs used

The pipeline recursively scans the supplied `--data-dir` and uses every file matching `*_1.fastq.gz` as a forward-read sample input. Reverse reads like `*_2.fastq.gz` are ignored. This works whether the FASTQs are directly under `data/` or nested under a study subdirectory such as `data/PRJNA825639/`.

It does not do paired-end merging or subsampling.

### What the scripts do

1. `run_deblur.py` recursively discovers forward reads matching `*_1.fastq.gz` under the input directory.
2. `run_deblur.py` stages each discovered run through `fastp` into `work/deblur/fastp/`.
3. `run_deblur.py` calls `deblur workflow` once on that staged directory. The primary Deblur artifact is `work/deblur/workflow/all.biom`.
4. `build_table.py` is currently incompatible with the refactored `run_deblur.py`, because it still expects per-sample `.clean` files rather than Deblur workflow BIOM output.
5. `diversity.py` is therefore also currently incompatible with the refactored `run_deblur.py`, because it depends on `build_table.py` producing `results/forward_only/feature_table.tsv`.
6. Once `build_table.py` is updated or replaced, `diversity.py` will still compute beta diversity (default: `braycurtis`) and run PCoA from the resulting feature table.
7. `diversity.py` writes tabular outputs and a PNG plot.

### Outputs

- Working outputs: `work/deblur/`
  - `fastp/` staged forward-read files
  - `workflow/` Deblur outputs including `all.biom`, `all.seqs.fa`, and workflow metadata/logs
- Final outputs: `results/forward_only/`
  - `feature_table.tsv`
  - `distance_matrix_braycurtis.tsv`
  - `pcoa_coordinates.tsv`
  - `pcoa_plot.png`

## Notes

This is the first-pass reproducible implementation and is intentionally minimal so future work can branch into paired-end workflows, subsampling comparisons, richer metadata integration, and alternative distance metrics.
