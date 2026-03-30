# PCOA prototype

A prototype project to download 16S data from ENA, subsample reads, run QC, and generate PCoA plots quickly.

## Tooling

- `fastq-dl` for downloading FASTQ data from ENA
- `seqkit` for subsampling reads
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
