# SLURM scripts

These scripts are the versioned batch wrappers used for the repository's CRC,
MMC20, and MMC1 analyses. They are intentionally plain Bash files: edit and
test them locally, copy the required file to Barnacle, and submit it from the
repository root.

## Scripts

| Script | Purpose | Invocation |
| --- | --- | --- |
| `crc_download.sh` | Download the two CRC studies | `sbatch` |
| `crc_subsample.sh` | Build the CRC 50%, 25%, and 10% FASTQs | `sbatch` |
| `crc_pcoa.sh` | Run one CRC level selected by `LEVEL` | `sbatch` |
| `mmc20_download.sh` | Download the 20-project MMC20 cohort as an array | `sbatch` |
| `mmc20_analysis.sh` | Submit the MMC20 full/subsampling/PCoA dependency graph | `bash` |
| `mmc1_human.sh` | Submit the MMC1 download/subsampling/PCoA dependency graph | `bash` |

## Copy to Barnacle

Create the destination once, then copy only the script being changed:

```bash
ssh barnacle 'mkdir -p ~/PCOA-prototype/scripts/slurm'
scp scripts/slurm/mmc1_human.sh \
  barnacle:~/PCOA-prototype/scripts/slurm/
```

The repository source files and metadata expected by a script must already be
present on the server. Copying a wrapper does not synchronize `src/`, metadata,
FASTQs, or GG2 artifacts.

## Submit

Run commands from `~/PCOA-prototype` so relative input and output paths resolve
consistently.

```bash
ssh barnacle
cd ~/PCOA-prototype
```

Single-job and array wrappers are submitted directly:

```bash
sbatch scripts/slurm/mmc20_download.sh
sbatch scripts/slurm/crc_download.sh
sbatch scripts/slurm/crc_subsample.sh
sbatch --job-name=crc-pcoa-sub10 \
  --export=ALL,LEVEL=sub10 \
  scripts/slurm/crc_pcoa.sh
```

The MMC20 and MMC1 orchestrators submit their own arrays and dependencies:

```bash
bash scripts/slurm/mmc20_analysis.sh
bash scripts/slurm/mmc1_human.sh
```

Configuration remains overridable through environment variables. For example:

```bash
PCOA_CONDA_ENV=qiime2-amplicon-2024.10 \
SAMPLING_DEPTH=1000 \
bash scripts/slurm/mmc1_human.sh
```

## Monitor and inspect completed jobs

```bash
squeue -u "$USER"
seff JOB_ID
sacct -j JOB_ID --format=JobID,State,Elapsed,TotalCPU,MaxRSS,ExitCode
```

SLURM `*.out` and `*.err` files, run directories, downloaded data, and analysis
outputs remain ignored by Git. Keep usernames, SSH endpoints beyond the local
`barnacle` alias, allocation/account names, notification emails, credentials,
and absolute home or lab filesystem paths out of these public scripts.
