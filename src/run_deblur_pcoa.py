#!/usr/bin/env python3
from __future__ import annotations
import argparse
import gzip
import re
import shutil
import subprocess
from collections import Counter
from pathlib import Path


EXPECTED_RUN_FILES = [
    "SRR27336825_1.fastq.gz",
    "SRR27336826_1.fastq.gz",
    "SRR27336827_1.fastq.gz",
    "SRR27336828_1.fastq.gz",
    "SRR27336829_1.fastq.gz",
    "SRR27336830_1.fastq.gz",
    "SRR27336831_1.fastq.gz",
    "SRR27336832_1.fastq.gz",
]
SIZE_PATTERN = re.compile(r";size=(\d+);?$")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run a forward-only Deblur -> beta diversity -> PCoA workflow."
    )
    parser.add_argument("--data-dir", type=Path, default=Path("data"))
    parser.add_argument("--work-dir", type=Path, default=Path("work/deblur"))
    parser.add_argument("--results-dir", type=Path, default=Path("results/forward_only"))
    parser.add_argument("--trim-length", type=int, default=250)
    parser.add_argument("--metric", default="braycurtis")
    parser.add_argument(
        "--error-dist",
        default="1,0.06,0.02,0.02,0.01,0.005,0.005,0.005,0.001,0.001,0.001,0.0005",
    )
    return parser.parse_args()


def ensure_executable(name: str) -> None:
    if shutil.which(name) is None:
        raise RuntimeError(f"Required executable not found in PATH: {name}")


def discover_inputs(data_dir: Path) -> list[Path]:
    merged_file = data_dir / "SAMN27531837.fastq.gz"
    if merged_file.exists():
        print(f"Ignoring merged sample file: {merged_file}")

    missing = [f for f in EXPECTED_RUN_FILES if not (data_dir / f).exists()]
    if missing:
        raise FileNotFoundError(
            "Missing expected per-run forward FASTQs: " + ", ".join(missing)
        )

    run_files = [data_dir / f for f in EXPECTED_RUN_FILES]
    empty = [str(fp) for fp in run_files if fp.stat().st_size == 0]
    if empty:
        raise RuntimeError("Input FASTQ(s) are empty: " + ", ".join(empty))

    return run_files


def fastq_to_trimmed_fasta(fastq_fp: Path, fasta_fp: Path, trim_length: int) -> int:
    total = 0
    kept = 0
    with gzip.open(fastq_fp, "rt") as fin, fasta_fp.open("w") as fout:
        while True:
            header = fin.readline().strip()
            if not header:
                break
            seq = fin.readline().strip()
            fin.readline()
            fin.readline()
            total += 1
            if len(seq) < trim_length:
                continue
            kept += 1
            seq_id = header[1:].split()[0]
            fout.write(f">{seq_id}\n{seq[:trim_length]}\n")

    if kept == 0:
        raise RuntimeError(f"No reads >= {trim_length} bp in {fastq_fp}")
    return kept


def run_deblur_for_sample(fastq_fp: Path, work_dir: Path, trim_length: int, error_dist: str) -> Path:
    sample_id = fastq_fp.name.replace("_1.fastq.gz", "")
    sample_out = work_dir / sample_id
    sample_out.mkdir(parents=True, exist_ok=True)

    fasta_fp = sample_out / f"{sample_id}.trim{trim_length}.fasta"
    kept_reads = fastq_to_trimmed_fasta(fastq_fp, fasta_fp, trim_length)
    print(f"  {sample_id}: converted {kept_reads} reads to trimmed FASTA", flush=True)

    derep_fp = sample_out / f"{fasta_fp.name}.derep"
    subprocess.run(
        [
            "deblur",
            "dereplicate",
            str(fasta_fp),
            str(derep_fp),
            "--min-size",
            "2",
        ],
        check=True,
    )

    subprocess.run(
        [
            "deblur",
            "deblur-seqs",
            str(derep_fp),
            "--error-dist",
            error_dist,
        ],
        check=True,
    )

    clean_fp = Path(str(derep_fp) + ".clean")
    if not clean_fp.exists():
        raise FileNotFoundError(f"Deblur clean output missing: {clean_fp}")
    return clean_fp


def parse_deblur_clean_fasta(clean_fp: Path) -> Counter:
    counts = Counter()
    with clean_fp.open() as fin:
        seq = None
        current_count = 0
        for line in fin:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"): 
                m = SIZE_PATTERN.search(line)
                current_count = int(m.group(1)) if m else 1
                seq = None
            else:
                seq = line
                counts[seq] += current_count
    return counts


def build_feature_table(clean_fastas: dict[str, Path]) -> pd.DataFrame:
    import pandas as pd

    by_sample = {}
    for sample_id, clean_fp in clean_fastas.items():
        by_sample[sample_id] = parse_deblur_clean_fasta(clean_fp)

    all_features = sorted({seq for counts in by_sample.values() for seq in counts})
    table = pd.DataFrame(0, index=sorted(by_sample.keys()), columns=all_features, dtype=int)
    for sample_id, counts in by_sample.items():
        for seq, n in counts.items():
            table.at[sample_id, seq] = n

    table.index.name = "sample_id"
    table.columns.name = "feature_sequence"
    return table


def compute_pcoa(feature_table: pd.DataFrame, metric: str):
    from skbio.diversity import beta_diversity
    from skbio.stats.ordination import pcoa

    ids = feature_table.index.tolist()
    matrix = feature_table.to_numpy(dtype=float)
    dm = beta_diversity(metric=metric, counts=matrix, ids=ids)
    ord_res = pcoa(dm)
    return dm, ord_res


def write_outputs(feature_table: pd.DataFrame, dm, ord_res, results_dir: Path, metric: str) -> None:
    import matplotlib.pyplot as plt
    import pandas as pd

    results_dir.mkdir(parents=True, exist_ok=True)

    feature_fp = results_dir / "feature_table.tsv"
    dm_fp = results_dir / f"distance_matrix_{metric}.tsv"
    coords_fp = results_dir / "pcoa_coordinates.tsv"
    plot_fp = results_dir / "pcoa_plot.png"

    feature_table.to_csv(feature_fp, sep="\t")

    dm_df = pd.DataFrame(dm.data, index=dm.ids, columns=dm.ids)
    dm_df.index.name = "sample_id"
    dm_df.to_csv(dm_fp, sep="\t")

    coords = ord_res.samples.copy()
    coords.index.name = "sample_id"
    coords.to_csv(coords_fp, sep="\t")

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.scatter(coords["PC1"], coords["PC2"])
    for sample_id, xv, yv in zip(coords.index, coords["PC1"], coords["PC2"]):
        ax.annotate(sample_id, (xv, yv), fontsize=8)
    ax.set_xlabel(f"PC1 ({ord_res.proportion_explained['PC1'] * 100:.2f}%)")
    ax.set_ylabel(f"PC2 ({ord_res.proportion_explained['PC2'] * 100:.2f}%)")
    ax.set_title("PCoA (forward-only Deblur table)")
    fig.tight_layout()
    fig.savefig(plot_fp, dpi=150)


def main() -> int:
    args = parse_args()
    ensure_executable("deblur")
    ensure_executable("vsearch")

    inputs = discover_inputs(args.data_dir)
    args.work_dir.mkdir(parents=True, exist_ok=True)

    clean_fastas = {}
    for fp in inputs:
        print(f"Running Deblur for: {fp}", flush=True)
        sample_id = fp.name.replace("_1.fastq.gz", "")
        clean_fastas[sample_id] = run_deblur_for_sample(
            fp, args.work_dir, args.trim_length, args.error_dist
        )

    feature_table = build_feature_table(clean_fastas)
    print(f"Feature table shape: {feature_table.shape}")
    print("Samples in feature table:", ", ".join(feature_table.index.tolist()))

    dm, ord_res = compute_pcoa(feature_table, metric=args.metric)
    write_outputs(feature_table, dm, ord_res, args.results_dir, args.metric)

    print(f"Finished. Outputs written under: {args.results_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
