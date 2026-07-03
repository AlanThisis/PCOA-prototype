# PCOA prototype

A prototype project to process ENA 16S forward reads, denoise with Deblur, and generate UniFrac PCoA outputs.

## Environment setup

### pcoa-prototype

Covers Deblur, BIOM, scikit-bio, seqkit, and plotting. Required for all steps except UniFrac.

```bash
conda env create -f environment.yml
conda activate pcoa-prototype
```

Check the main tools resolve:

```bash
fastq-dl --help
deblur --help
biom --help
python -c "import skbio; print(skbio.__version__)"
```

### QIIME2 (UniFrac only)

The UniFrac step requires a separate [QIIME2 amplicon environment](https://docs.qiime2.org/2024.10/install/) with the [q2-greengenes2 plugin](https://github.com/biocore/q2-greengenes2) installed. Any recent QIIME2 amplicon distribution works.

### Greengenes2 reference files

Download once into `data/gg2/` before running UniFrac (~175 MB total):

```bash
mkdir -p data/gg2
wget -P data/gg2 https://ftp.microbio.me/greengenes_release/current/2024.09.backbone.full-length.fna.qza
wget -P data/gg2 https://ftp.microbio.me/greengenes_release/current/2024.09.phylogeny.id.nwk.qza
```

If you update the pcoa-prototype environment, regenerate `environment.yml` with:

```bash
conda env export -n pcoa-prototype --from-history | sed '/^prefix: /d' > environment.yml
```

---

## End-to-end example: PRJEB44533

PRJEB44533 is a 16S gut microbiome study (204 samples). This example runs the full pipeline from download to UniFrac PCoA at 10% read subsampling.

### 1. Download

```bash
fastq-dl --accession PRJEB44533 --provider ena --outdir data/fastq_data/PRJEB44533/full
```

### 2. Subsample

```bash
python src/subsample_fastq.py \
  --input-dir data/fastq_data/PRJEB44533/full \
  --output-dir data/fastq_data/PRJEB44533/subsample_10 \
  --percent 10 \
  --seed 11
```

`subsample_fastq.py` uses `seqkit sample2` with two-pass mode (`-2`) for stable fraction sampling. Run it again with `--percent 25` / `--percent 50` for additional subsampling levels.

### 3. Deblur

```bash
python src/run_deblur.py \
  --data-dir data/fastq_data/PRJEB44533/subsample_10 \
  --work-dir work/PRJEB44533_sub10
```

Outputs land in `work/PRJEB44533_sub10/workflow/`, including `all.biom` and `all.seqs.fa`.

Deblur defaults: `--trim-length 150`, `--min-reads 0` (per-sample denoising, no cross-sample filter).

### 4. UniFrac PCoA

Activate the QIIME2 amplicon environment, then:

```bash
python src/unifrac.py \
  --deblur-dir work/PRJEB44533_sub10/workflow \
  --results-dir results/PRJEB44533_sub10
```

Outputs written to `results/PRJEB44533_sub10/`:
- `unweighted_unifrac_distance_matrix.tsv`
- `unweighted_unifrac_pcoa_results.tsv`
- `pcoa_plot.png`

`unifrac.py` uses GG2's `non-v4-16s` closed-reference action (vsearch at 99%) to map Deblur ASVs onto the GG2 backbone, then computes UniFrac against the GG2 ID phylogeny. Rarefaction depth defaults to the minimum sample depth after backbone mapping; override with `--sampling-depth`.

---

## Scripts

| Script | Env | Description |
|--------|-----|-------------|
| `src/subsample_fastq.py` | pcoa-prototype | Subsample FASTQs to a given percent with seqkit |
| `src/run_deblur.py` | pcoa-prototype | Run Deblur on a directory of forward FASTQs |
| `src/build_table.py` | pcoa-prototype | Export Deblur BIOM to TSV feature table |
| `src/diversity.py` | pcoa-prototype | Bray-Curtis beta diversity + PCoA from feature table |
| `src/merge_biom.py` | pcoa-prototype | Merge BIOM tables across studies |
| `src/unifrac.py` | QIIME2 amplicon | UniFrac PCoA via GG2 and QIIME2 |
| `src/get_ENA_metadata.py` | pcoa-prototype | Fetch sample metadata CSV for an ENA project |

### Inputs

Scripts recursively scan for forward reads matching `*_1.fastq.gz` (ENA-style) or `*_R1_001.fastq.gz` (demux-export style). Reverse reads are ignored.

### Metadata

Feature tables do not contain disease/group labels. Join `pcoa_coordinates.tsv` or `feature_table.tsv` to a metadata TSV by sample ID.
