# File Usage Reference

This document describes the prepared inputs and self-contained output layout
used by `src/run_pipeline.py`.

## Prepared FASTQ Inputs

- Forward reads may match `*_1.fastq.gz` (ENA style) or
  `*_R1_001.fastq.gz` (demultiplexed export style).
- Input discovery is recursive. Reverse reads are ignored.
- Every file must be non-empty, and sample identifiers must be unique across
  every `--study NAME=FASTQ_DIR` supplied to one run.
- Downloading and subsampling happen before orchestration. They intentionally
  remain separate so server jobs can reuse the same prepared FASTQs.

## Canonical Run Layout

Every invocation requires a unique `--run-dir`. `runs/` is ignored by Git, but
the run directory may also point to server scratch or project storage outside
the clone.

```text
<run-dir>/
├── run_manifest.json
├── run_state.json
├── timings/
│   ├── attempt-001/
│   │   ├── pipeline.tsv
│   │   ├── deblur-<study>.tsv
│   │   ├── merge.tsv                  # cross-study only
│   │   ├── unifrac.tsv
│   │   └── plot-<metadata-column>.tsv
│   └── attempt-002/                   # created by --resume
├── work/
│   ├── deblur/<study>/workflow/
│   ├── merged/                        # cross-study only
│   └── qiime2/
└── results/
    ├── distance_matrix_unweighted_unifrac.tsv
    ├── pcoa_coordinates_unweighted_unifrac.txt
    ├── pcoa_plot_unweighted_unifrac.png
    └── pcoa_<metadata-column>.png
```

## Provenance and Resume Files

- `run_manifest.json` records resolved FASTQ inputs, scientific parameters,
  metadata and GG2 paths, input sizes and modification times, Git commit,
  executable paths, Conda environment, hostname, SLURM job ID, and creation
  time.
- `run_state.json` records each stage as `pending`, `running`, `completed`, or
  `failed`, plus every execution attempt and its thread count.
- `--resume` is explicit. It requires a matching manifest and skips completed
  stages only when all expected outputs remain non-empty.
- Thread count may change between attempts. It is recorded per attempt rather
  than treated as a scientific compatibility parameter.
- If Deblur or merge reruns, input-derived QIIME2 imports and GG2 mapping are
  invalidated before UniFrac resumes. A failed UniFrac stage can still reuse
  its partial QIIME2 work when no upstream stage changed.

## Deblur Handoff Artifacts

- `work/deblur/<study>/workflow/all.biom` is the per-study feature table.
- `work/deblur/<study>/workflow/all.seqs.fa` contains representative feature
  sequences.
- Cross-study runs merge these into `work/merged/all.biom` and
  `work/merged/all.seqs.fa` before GG2 mapping.
- Deblur internals are not retained by default. Use
  `--keep-deblur-tmp-files` only for diagnosis because these files can consume
  several gigabytes.

## QIIME2 Work Artifacts

`work/qiime2/` contains imported tables and sequences, GG2 backbone-mapped
artifacts, the rarefied table, UniFrac distance artifact, PCoA artifact, and
exports. These are retained to diagnose a failure and may be reused only by a
manifest-compatible resume of the same run.

## Final Results

- `distance_matrix_unweighted_unifrac.tsv` is the sample-by-sample unweighted
  UniFrac distance matrix.
- `pcoa_coordinates_unweighted_unifrac.txt` is the QIIME2 PCoA ordination
  export and is the input to metadata plotting.
- `pcoa_plot_unweighted_unifrac.png` is an unlabeled quick-look plot.
- `pcoa_<metadata-column>.png` is generated for each repeated `--color-by`.

Metadata must include either a `sample-id` column or a `run_accessions` column.
The latter may contain semicolon-separated ENA run accessions. Samples not
found in metadata are plotted as `Unknown`.
