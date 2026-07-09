# PCOA prototype

A prototype project to process ENA 16S forward reads, denoise with Deblur, and generate UniFrac PCoA outputs.

## Environment setup

This pipeline runs inside a single [QIIME2 Rachis environment](https://library.qiime2.org/quickstart/qiime2). QIIME2 supplies Deblur, BIOM, scikit-bio, seqkit, matplotlib, pandas, requests, vsearch, and UniFrac support. This repo adds `fastq-dl` for ENA downloads and `q2-greengenes2` for Greengenes2-backed UniFrac.

The environment is assembled in three layers:

1. Official QIIME2 Rachis conda environment. Treat this as the base because QIIME2 pins a large stack of Python, compiled, R, and plugin dependencies.
2. Repo extras from `environment.yml`. This adds small helper dependencies, including `joblib`, `msgpack-python`, `nltk`, and a pinned `setuptools` compatible with the current Greengenes2 plugin.
3. Pip-installed tools. `fastq-dl==4.0.1` is installed with pip to avoid expanding the conda solve for the QIIME2 environment. `q2-greengenes2` and `redbiom` are installed with `pip --no-deps` so pip does not replace QIIME2's conda-managed packages such as `scikit-bio`.

This is intentionally not a fully hand-written conda YAML for every package. A custom all-in-one YAML is more fragile because QIIME2's distribution packages need to stay internally consistent.

### One-shot setup

Run from the repo root:

```bash
bash scripts/setup_qiime2_gg2_env.sh
```

By default this creates `rachis-qiime2-2026.4`, matching the current QIIME2 Library quickstart. Override the version or env name if needed:

```bash
QIIME2_VERSION=2026.1 QIIME2_ENV_NAME=pcoa-qiime2-gg2 bash scripts/setup_qiime2_gg2_env.sh
```

Activate before running any pipeline script:

```bash
conda activate rachis-qiime2-2026.4
```

Check the main tools resolve:

```bash
fastq-dl --help
deblur --help
qiime info
python -c "import skbio; print(skbio.__version__)"
```

### Existing QIIME2 env

If QIIME2 is already installed, update that env with the repo extras. Install the GG2 plugin with `--no-deps` so pip does not replace QIIME2's conda-managed packages.

```bash
conda env update -n <your-qiime2-env-name> -f environment.yml
conda install -n <your-qiime2-env-name> -c conda-forge "setuptools<81"
conda run -n <your-qiime2-env-name> python -m pip install --force-reinstall --no-deps setuptools==80.9.0
conda run -n <your-qiime2-env-name> python -m pip install fastq-dl==4.0.1
conda run -n <your-qiime2-env-name> python -m pip install --no-deps redbiom==0.3.9 q2-greengenes2==2024.1
```

`fastq-dl --check` may report `sracha` missing. That only affects SRA downloads; the documented PRJEB44533 command uses `--provider ena`, which needs `wget`.

### Greengenes2 reference files

Download once into `data/gg2/` before running UniFrac (~175 MB total):

```bash
mkdir -p data/gg2
wget -P data/gg2 https://ftp.microbio.me/greengenes_release/current/2024.09.backbone.full-length.fna.qza
wget -P data/gg2 https://ftp.microbio.me/greengenes_release/current/2024.09.phylogeny.id.nwk.qza
```

---

## Server runbook: PRJEB44533

PRJEB44533 is a 16S gut microbiome study (204 samples). This example runs the pipeline from raw ENA FASTQs to a Greengenes2-backed unweighted UniFrac PCoA plot using the 10% subsample.

Start from a fresh server clone:

```bash
git clone <repo-url>
cd PCOA-prototype
bash scripts/setup_qiime2_gg2_env.sh
conda activate rachis-qiime2-2026.4
```

Download the Greengenes2 artifacts once:

```bash
mkdir -p data/gg2
wget -P data/gg2 https://ftp.microbio.me/greengenes_release/current/2024.09.backbone.full-length.fna.qza
wget -P data/gg2 https://ftp.microbio.me/greengenes_release/current/2024.09.phylogeny.id.nwk.qza
```

### 1. Download

```bash
fastq-dl --accession PRJEB44533 --provider ena --outdir data/fastq_data/PRJEB44533/full
```

### 2. Fetch Metadata

```bash
mkdir -p data/PRJEB44533

python src/get_ENA_metadata.py PRJEB44533 \
  -o data/PRJEB44533/metadata.csv
```

For PRJEB44533, use the `description` column for disease/control-style labels when plotting metadata-colored PCoA output.

### 3. Subsample

```bash
python src/subsample_fastq.py \
  --input-dir data/fastq_data/PRJEB44533/full \
  --output-dir data/fastq_data/PRJEB44533/subsample_10 \
  --percent 10 \
  --seed 11
```

`subsample_fastq.py` uses `seqkit sample2` with two-pass mode (`-2`) for stable fraction sampling. Run it again with `--percent 25` / `--percent 50` for additional subsampling levels.

### 4. Deblur

```bash
python src/run_deblur.py \
  --data-dir data/fastq_data/PRJEB44533/subsample_10 \
  --work-dir work/PRJEB44533_sub10 \
  --trim-length 120 \
  --min-reads 0 \
  --jobs-to-start 10
```

Outputs land in `work/PRJEB44533_sub10/workflow/`, including `all.biom` and `all.seqs.fa`.

Use `--trim-length 120` for PRJEB44533 because the reads are about 122 bp. `--min-reads 0` disables Deblur's cross-sample feature-count filter, so samples are treated independently. Adjust `--jobs-to-start` to match available server cores.

### 5. UniFrac PCoA

```bash
python src/unifrac.py \
  --deblur-dir work/PRJEB44533_sub10/workflow \
  --results-dir results/PRJEB44533_sub10 \
  --threads 10
```

Outputs written to `results/PRJEB44533_sub10/`:
- `distance_matrix_unweighted_unifrac.tsv`
- `pcoa_coordinates_unweighted_unifrac.txt`
- `pcoa_plot_unweighted_unifrac.png`

`unifrac.py` uses GG2's `non-v4-16s` closed-reference action (vsearch at 99%) to map Deblur ASVs onto the GG2 backbone, then computes UniFrac against the GG2 ID phylogeny. Rarefaction depth defaults to the minimum sample depth after backbone mapping; override with `--sampling-depth`.

For a smaller local test, use `--jobs-to-start 4` and `--threads 4`.

### 6. Plot PCoA With Metadata Labels

After UniFrac finishes, use `plot_pcoa.py` to color the ordination by the PRJEB44533 metadata `description` column:

```bash
python src/plot_pcoa.py \
  --pcoa results/PRJEB44533_sub10/pcoa_coordinates_unweighted_unifrac.txt \
  --metadata data/PRJEB44533/metadata.csv \
  --color-by description \
  --out results/PRJEB44533_sub10/pcoa_description.png \
  --title "PRJEB44533 - description"
```

`plot_pcoa.py` matches PCoA sample IDs to the metadata CSV through the `run_accessions` column written by `get_ENA_metadata.py`. Samples not found in metadata are labeled `Unknown`.

---

## Scripts

| Script | Env | Description |
|--------|-----|-------------|
| `src/subsample_fastq.py` | QIIME2 + repo extras | Subsample FASTQs to a given percent with seqkit |
| `src/run_deblur.py` | QIIME2 + repo extras | Run Deblur on a directory of forward FASTQs |
| `src/build_table.py` | QIIME2 + repo extras | Export Deblur BIOM to TSV feature table |
| `src/diversity.py` | QIIME2 + repo extras | Bray-Curtis beta diversity + PCoA from feature table |
| `src/merge_biom.py` | QIIME2 + repo extras | Merge BIOM tables across studies |
| `src/unifrac.py` | QIIME2 + repo extras | UniFrac PCoA via GG2 and QIIME2 |
| `src/get_ENA_metadata.py` | QIIME2 + repo extras | Fetch sample metadata CSV for an ENA project |

### Inputs

Scripts recursively scan for forward reads matching `*_1.fastq.gz` (ENA-style) or `*_R1_001.fastq.gz` (demux-export style). Reverse reads are ignored.

### Metadata

Feature tables do not contain disease/group labels. Join `pcoa_coordinates.tsv` or `feature_table.tsv` to a metadata TSV by sample ID.
