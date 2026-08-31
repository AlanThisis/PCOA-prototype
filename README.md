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
python src/download_ena_amplicon.py PRJEB44533 \
  --output-dir data/fastq_data/PRJEB44533/full \
  --workers 8
```

The downloader queries ENA first and selects only `AMPLICON` runs. For paired
runs it downloads `_1.fastq.gz`. When ENA documents a paired-layout run with
only an unsuffixed archive FASTQ, that file contains unpaired reads and is
retained. Single-end and unpaired FASTQs are normalized to the `_1.fastq.gz`
naming convention.

Downloads are concurrent, resumable through `.part` files, and checked against
ENA's reported byte count and MD5. Failed files receive multiple recovery
passes with backoff. The command exits successfully only after every selected
FASTQ passes final validation, then writes `.ena_download_complete.json` in the
output directory. The manifest, per-run status, and project summary TSVs are
updated on every invocation. `ena_amplicon_forward_history.tsv` retains the
timing and outcome of every invocation, including the initial download and any
later repair or validation passes. Rerunning the same command is safe: verified
files are skipped and only incomplete files are repaired.

For an unattended job that should keep retrying within a wall-time budget:

```bash
python src/download_ena_amplicon.py PRJEB44533 \
  --output-dir data/fastq_data/PRJEB44533/full \
  --workers 4 \
  --retry-until-complete \
  --max-runtime 82800
```

`--max-runtime` is a soft limit in seconds; `82800` leaves one hour of margin in
a 24-hour SLURM allocation. If the command exits incomplete, run the identical
command again to resume. Use `--dry-run` to inspect the selected files and total
size without downloading or invalidating an existing completion marker.

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
  --seed 11 \
  --timings-tsv data/fastq_data/PRJEB44533/subsample_10_timings.tsv
```

`subsample_fastq.py` uses `seqkit sample2` with two-pass mode (`-2`) for stable fraction sampling. Run it again with `--percent 25` / `--percent 50` for additional subsampling levels.

### 4. Run the Prepared FASTQs

`run_pipeline.py` is the canonical entry point once the forward-read FASTQs are
downloaded and, if desired, subsampled. It runs Deblur, GG2-backed unweighted
UniFrac, PCoA, and requested metadata-colored plots in one self-contained run
directory.

```bash
python src/run_pipeline.py \
  --study PRJEB44533=data/fastq_data/PRJEB44533/subsample_10 \
  --metadata data/PRJEB44533/metadata.csv \
  --color-by description \
  --run-dir runs/prjeb44533/sub10-001 \
  --trim-length 120 \
  --min-reads 0 \
  --threads 10
```

Use a unique `--run-dir` for every server job, normally including the SLURM job
ID. The directory contains immutable provenance, resume state, per-attempt
timings, retained handoff artifacts, and final results:

```text
runs/prjeb44533/sub10-001/
├── run_manifest.json
├── run_state.json
├── timings/attempt-001/
├── work/
│   ├── deblur/PRJEB44533/workflow/
│   └── qiime2/
└── results/
    ├── distance_matrix_unweighted_unifrac.tsv
    ├── pcoa_coordinates_unweighted_unifrac.txt
    ├── pcoa_plot_unweighted_unifrac.png
    └── pcoa_description.png
```

Use `--trim-length 120` for PRJEB44533 because the reads are about 122 bp.
`--min-reads 0` disables Deblur's cross-sample feature-count filter. Deblur's
bulky internal temporary files are discarded after successful processing by
default; add `--keep-deblur-tmp-files` only when diagnosing Deblur.

### Cross-study Runs

Repeat `--study NAME=FASTQ_DIR` to process multiple studies. Each study is
Deblurred separately, then their `all.biom` and `all.seqs.fa` artifacts are
merged before one shared GG2/UniFrac analysis:

```bash
python src/run_pipeline.py \
  --study ERP005534=data/fastq_data/CRC_cross_study/ERP005534/sub10 \
  --study PRJEB46665=data/fastq_data/CRC_cross_study/PRJEB46665/sub10 \
  --metadata data/fastq_data/CRC_cross_study/crc_cross_study_metadata.tsv \
  --color-by disease_status \
  --color-by study \
  --run-dir runs/crc-cross-study/sub10-450671 \
  --trim-length 120 \
  --min-reads 0 \
  --threads 16
```

Study names and sample identifiers must be unique. For a single study, the
merge stage is skipped and that study's Deblur workflow is passed directly to
UniFrac.

`unifrac.py` uses GG2's `non-v4-16s` closed-reference action (vsearch at 99%) to map Deblur ASVs onto the GG2 backbone, then computes UniFrac against the GG2 ID phylogeny. Rarefaction depth defaults to 1,000 reads after backbone mapping; override with `--sampling-depth`.

### Resume a Failed Run

Existing runs are never resumed implicitly. Rerun the same command with
`--resume`; inputs, metadata, GG2 path, Git commit, and scientific parameters
must match the saved manifest. FASTQs, metadata, and GG2 artifacts are compared
using their resolved path, size, and modification time. `--threads` may differ
between attempts.

```bash
python src/run_pipeline.py \
  --study PRJEB44533=data/fastq_data/PRJEB44533/subsample_10 \
  --metadata data/PRJEB44533/metadata.csv \
  --color-by description \
  --run-dir runs/prjeb44533/sub10-001 \
  --trim-length 120 \
  --min-reads 0 \
  --threads 4 \
  --resume
```

Completed stages are skipped only when all expected outputs remain non-empty.
If a stage must rerun, every downstream stage reruns as well. Each invocation
gets a new `timings/attempt-NNN/` directory, so failed-attempt timing records
are preserved. When Deblur or merge reruns, the orchestrator refreshes QIIME2's
imported table, representative sequences, and GG2 mapping instead of reusing
stale upstream artifacts.

### Generic SLURM Wrapper

SLURM should only allocate resources, activate the environment, choose a unique
run path, and invoke the Python orchestrator. Keep project-specific job scripts
local rather than adding one script per dataset to the repository.

```bash
#!/usr/bin/env bash
#SBATCH --partition=short
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=1-00:00:00
#SBATCH --job-name=pcoa
#SBATCH --output=slurm-%x-%j.out

set -euo pipefail
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate qiime2-amplicon-2024.10
cd /path/to/PCOA-prototype

python src/run_pipeline.py \
  --study PRJEB44533=data/fastq_data/PRJEB44533/subsample_10 \
  --metadata data/PRJEB44533/metadata.csv \
  --color-by description \
  --run-dir "runs/prjeb44533/sub10-${SLURM_JOB_ID}" \
  --trim-length 120 \
  --min-reads 0 \
  --threads "${SLURM_CPUS_PER_TASK}"
```

Adjust the partition, environment name, memory, and wall time for the server.

---

## Scripts

| Script | Env | Description |
|--------|-----|-------------|
| `src/run_pipeline.py` | QIIME2 + repo extras | Canonical prepared-FASTQ pipeline orchestrator |
| `src/subsample_fastq.py` | QIIME2 + repo extras | Subsample FASTQs to a given percent with seqkit |
| `src/run_deblur.py` | QIIME2 + repo extras | Run Deblur on a directory of forward FASTQs |
| `src/merge_biom.py` | QIIME2 + repo extras | Merge BIOM tables across studies |
| `src/unifrac.py` | QIIME2 + repo extras | UniFrac PCoA via GG2 and QIIME2 |
| `src/plot_pcoa.py` | QIIME2 + repo extras | Plot PCoA coordinates colored by metadata |
| `src/get_ENA_metadata.py` | QIIME2 + repo extras | Fetch sample metadata CSV for an ENA project |

### Runtime Profiling

The orchestrator always writes an attempt-level `pipeline.tsv` plus detailed
component timing TSVs under `<run-dir>/timings/attempt-NNN/`. These separate
Deblur, merge, Greengenes2 mapping, rarefaction, UniFrac, PCoA, exports, and
metadata plotting rather than labeling the entire QIIME stage as UniFrac.

When running a component CLI directly for debugging, enable the same detailed
timing format with its optional `--timings-tsv` argument:

```bash
python src/unifrac.py \
  --deblur-dir work/PRJEB44533_sub10/workflow \
  --results-dir results/PRJEB44533_sub10 \
  --threads 10 \
  --timings-tsv results/PRJEB44533_sub10/unifrac_timings.tsv
```

The common timing schema records the component, operation, optional input item,
UTC start/end timestamps, elapsed seconds, status, exit code, external command,
and error or skip reason. Records are written when an operation starts and
an additional terminal row is appended when it finishes. This append-only event
log keeps timing I/O linear and leaves an interrupted operation's latest event
marked `running`. Use `completed`, `failed`, and `skipped` rows for elapsed-time
analysis; cached operations are marked `skipped` instead of being reported as
zero-second work.

For `unifrac.py`, the detailed operations distinguish:

- QIIME2 feature-table and representative-sequence imports
- Greengenes2 `non-v4-16s` backbone mapping
- mapped-table export and rarefaction-depth handling
- feature-table rarefaction
- `diversity beta-phylogenetic`, the actual unweighted UniFrac calculation
- `diversity pcoa`
- artifact exports and local plotting

For a complete benchmark, use a fresh run directory. A validated resume can
reuse QIIME2 imports and GG2 mapping artifacts, so its component timings do not
represent a clean full-pipeline runtime.

### Component CLIs for Debugging

The orchestrator invokes component scripts using the active Python executable.
They remain available for inspecting or rerunning individual boundaries:

```bash
python src/run_deblur.py --help
python src/merge_biom.py --help
python src/unifrac.py --help
python src/plot_pcoa.py --help
```

### Inputs

Scripts recursively scan for forward reads matching `*_1.fastq.gz` (ENA-style) or `*_R1_001.fastq.gz` (demux-export style). Reverse reads are ignored.

### Metadata

Metadata labels are not stored in the feature table. `plot_pcoa.py` joins the
UniFrac PCoA sample IDs to either the metadata `sample-id` or `run_accessions`
column.
