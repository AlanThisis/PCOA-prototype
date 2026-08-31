# MMC1 human 16S metadata cleanup

The source metadata is preserved outside this repository at
`../ENA-metadata-harmonization/data/mmc1_16s_amplicon_filtered.csv`. Generated
pipeline inputs are written under `data/mmc1_16s_cleaning/`, which is ignored
by Git because it contains data-derived artifacts.

Run the reproducible cleanup from the repository root:

```bash
python src/prepare_mmc1_16s_metadata.py \
  ../ENA-metadata-harmonization/data/mmc1_16s_amplicon_filtered.csv \
  --output-dir data/mmc1_16s_cleaning
```

## Selection policy

1. Retain only studies labeled `study_hmr=Human` and specimens labeled
   `host_species_resolved=Human`. Mixed human/animal studies are excluded.
2. Expand semicolon-delimited accessions to one record per ENA run. The run
   accession is the pipeline `sample-id`, including for umbrella BioSamples
   whose submitter attached multiple distinct specimens to one BioSample.
3. Exclude explicit laboratory blanks and negative controls based on run- and
   experiment-level aliases; the word `control` alone is not excluded because
   healthy biological controls are valid samples.
4. Exclude explicitly identified ITS, 18S, and fungal experiments. ENA's
   `AMPLICON` strategy is not by itself evidence that a run is bacterial 16S.
5. For `PRJNA597168`, retain only the separately deposited forward-direction
   run. The reverse-direction run is not a second biological specimen.
6. Exclude `PRJNA339931`, `PRJNA564314`, and `PRJNA985881` because their
   BioSample structure or sample-type labeling cannot support this automated
   cross-study analysis.
7. Exclude `PRJNA1212986` because it provides PacBio subread archives that are
   incompatible with the fixed-length Deblur workflow. Exclude `PRJNA381268`
   because its FASTQs are protected by dbGaP and cannot be downloaded from ENA.
8. Exclude `SRR10377565` from `PRJNA580506`: ENA advertises the archive file,
   but repeated HTTPS requests return an HTML error document and FTP reports
   that the file is absent. The other 47 runs from that project are retained.
9. Exclude `PRJEB75547`. ENA identifies its runs as Oxford Nanopore GridION
   full-length amplicons. Their quality-score range is incompatible with the
   current short-read Deblur workflow; this is a platform mismatch rather than
   a corrupt or incomplete download.

The exact project-level reasons are emitted in `excluded_projects.tsv`; all
run-level decisions are recorded in `excluded_runs.tsv`.

## Generated audit artifacts

- `cleaned_metadata.csv`: run-resolved metadata accepted by `run_pipeline.py`.
- `approved_runs.tsv`: exact project/run allowlist consumed by the downloader.
- `projects.txt`: sorted BioProject accessions used by SLURM arrays.
- `excluded_projects.tsv` and `excluded_runs.tsv`: exclusion decisions and
  reasons.
- `cleanup_summary.json`: selected/excluded totals used by preflight checks.

The downloader validates ENA byte counts and MD5 checksums before recording a
project as complete. It is invoked with `--run-manifest approved_runs.tsv`, so
ENA records not present in the reviewed allowlist are not downloaded. A listed
run that ENA does not return as a downloadable AMPLICON forward read makes the
project fail rather than silently shrinking the cohort. Any pre-existing
top-level FASTQ outside the allowlist is moved under `unapproved_fastqs/`, where
the pipeline will not discover it. Greengenes2 mapping remains a downstream
compatibility check rather than the primary method for removing known non-16S
inputs.

The current generated set contains 15,438 unique approved runs in 108 projects.
The exact counts are authoritative in `cleanup_summary.json` and are checked
again by the Barnacle preflight before submission.
