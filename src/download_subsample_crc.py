#!/usr/bin/env python3
"""Download and 10%-subsample R1 FASTQs for CRC cross-study.

Selects 25 samples per study (balanced cancer/control where possible):
  - ERP005534: labels from ENA sample XML <TAG>diagnosis</TAG>
  - PRJEB46665: labels from sample_title prefix (cancer_ / nocancer_)

Downloads R1 only, subsamples 10% in-place, deletes originals.
Writes crc_cross_study_metadata.tsv.
"""
from __future__ import annotations

import csv
import subprocess
import sys
import time
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path

SEED = 42
FRACTION = 0.1
N_PER_STUDY = 25          # total per study; split as evenly as possible
LARGE_FILE_BYTES = 8 * 1024 * 1024  # >8 MB = not yet subsampled

BASE = Path(__file__).resolve().parent.parent
CRC_DIR = BASE / "data" / "fastq_data" / "CRC_cross_study"
ERP_TSV = CRC_DIR / "ERP005534" / "ERP005534-run-info.tsv"
PRJEB_TSV = CRC_DIR / "PRJEB46665" / "PRJEB46665-run-info.tsv"
METADATA_OUT = CRC_DIR / "crc_cross_study_metadata.tsv"


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def ftp_to_https(ftp_url: str) -> str:
    return ftp_url.replace("ftp.sra.ebi.ac.uk", "ftp.sra.ebi.ac.uk", 1)


def r1_url(fastq_ftp: str) -> str:
    """Pick the R1 URL from a semicolon-separated fastq_ftp field."""
    urls = fastq_ftp.split(";")
    for u in urls:
        if u.endswith("_1.fastq.gz") or u.endswith("_R1_001.fastq.gz"):
            return "https://" + u if not u.startswith("http") else u
    # fallback: first URL
    u = urls[0]
    return "https://" + u if not u.startswith("http") else u


def seqkit_path() -> str:
    env_bin = Path(sys.executable).resolve().parent / "seqkit"
    if env_bin.exists():
        return str(env_bin)
    import shutil
    resolved = shutil.which("seqkit")
    if resolved is None:
        raise RuntimeError("seqkit not found in PATH")
    return resolved


def subsample_inplace(fastq: Path, seqkit: str) -> None:
    orig = fastq.with_name(fastq.name.replace(".fastq.gz", "_orig.fastq.gz"))
    fastq.rename(orig)
    subprocess.run(
        [seqkit, "sample2", "-p", str(FRACTION), "-s", str(SEED), "-2",
         str(orig), "-o", str(fastq)],
        check=True,
    )
    orig.unlink()


def wget_file(url: str, dest: Path) -> None:
    print(f"  downloading {dest.name} ...", flush=True)
    dest.parent.mkdir(parents=True, exist_ok=True)
    urllib.request.urlretrieve(url, dest)


def fetch_diagnosis(sample_acc: str, retries: int = 3) -> str | None:
    url = f"https://www.ebi.ac.uk/ena/browser/api/xml/{sample_acc}"
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=20) as resp:
                tree = ET.fromstring(resp.read())
            for attr in tree.iter("SAMPLE_ATTRIBUTE"):
                tag = attr.findtext("TAG") or ""
                if tag.strip().lower() == "diagnosis":
                    return (attr.findtext("VALUE") or "").strip()
            return None
        except Exception as e:
            if attempt == retries - 1:
                print(f"  WARNING: could not fetch XML for {sample_acc}: {e}", flush=True)
                return None
            time.sleep(2)
    return None


# ---------------------------------------------------------------------------
# ERP005534 — XML-based labels
# ---------------------------------------------------------------------------

def select_erp005534(n: int) -> list[dict]:
    """Return up to n rows (balanced Cancer/Normal) with diagnosis from XML."""
    rows = read_tsv(ERP_TSV)

    # one run per sample_accession (first occurrence)
    seen: dict[str, dict] = {}
    for row in rows:
        sa = row["sample_accession"]
        if sa not in seen and row.get("fastq_ftp"):
            seen[sa] = row

    samples = list(seen.values())
    print(f"ERP005534: {len(samples)} unique samples — fetching XML diagnoses...", flush=True)

    labelled: dict[str, list[dict]] = {"Cancer": [], "Normal": []}
    for i, row in enumerate(samples):
        sa = row["sample_accession"]
        diag = fetch_diagnosis(sa)
        if diag in labelled:
            labelled[diag].append(row)
        if (i + 1) % 20 == 0:
            print(f"  ...{i+1}/{len(samples)} XMLs fetched "
                  f"(Cancer={len(labelled['Cancer'])}, Normal={len(labelled['Normal'])})",
                  flush=True)
        # early-stop once we have enough of both
        if len(labelled["Cancer"]) >= n // 2 + 2 and len(labelled["Normal"]) >= n // 2 + 2:
            print(f"  early stop after {i+1} XMLs", flush=True)
            break

    n_cancer = min(n // 2, len(labelled["Cancer"]))
    n_normal = min(n - n_cancer, len(labelled["Normal"]))
    selected = labelled["Cancer"][:n_cancer] + labelled["Normal"][:n_normal]
    print(f"ERP005534 selected: {n_cancer} Cancer + {n_normal} Normal", flush=True)

    result = []
    for row in selected:
        diag = "Cancer" if row in labelled["Cancer"][:n_cancer] else "Normal"
        result.append({
            "run_accession": row["run_accession"],
            "study": "ERP005534",
            "disease_status": "cancer" if diag == "Cancer" else "control",
            "url": r1_url(row["fastq_ftp"]),
            "fastq_dir": CRC_DIR / "ERP005534" / "fastq",
        })
    return result


# ---------------------------------------------------------------------------
# PRJEB46665 — sample_title labels
# ---------------------------------------------------------------------------

def select_prjeb46665(n: int) -> list[dict]:
    rows = read_tsv(PRJEB_TSV)
    cancer, control = [], []
    for row in sorted(rows, key=lambda r: r["run_accession"]):
        title = row.get("sample_title", "")
        if title.startswith("cancer_"):
            cancer.append(row)
        elif title.startswith("no_cancer_"):
            control.append(row)

    n_cancer = min(n // 2, len(cancer))
    n_control = min(n - n_cancer, len(control))
    selected_cancer = cancer[:n_cancer]
    selected_control = control[:n_control]
    print(f"PRJEB46665 selected: {n_cancer} cancer + {n_control} control", flush=True)

    result = []
    for row in selected_cancer + selected_control:
        label = "cancer" if row in selected_cancer else "control"
        result.append({
            "run_accession": row["run_accession"],
            "study": "PRJEB46665",
            "disease_status": label,
            "url": r1_url(row["fastq_ftp"]),
            "fastq_dir": CRC_DIR / "PRJEB46665" / "fastq",
        })
    return result


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def main() -> int:
    seqkit = seqkit_path()

    erp_samples = select_erp005534(N_PER_STUDY)
    prjeb_samples = select_prjeb46665(N_PER_STUDY)
    all_samples = erp_samples + prjeb_samples

    print(f"\nTotal selected: {len(all_samples)} samples. Starting download + subsample...\n",
          flush=True)

    for s in all_samples:
        run_acc = s["run_accession"]
        fastq_dir: Path = s["fastq_dir"]
        dest = fastq_dir / f"{run_acc}_1.fastq.gz"

        if dest.exists() and dest.stat().st_size <= LARGE_FILE_BYTES:
            print(f"  {run_acc}: already subsampled, skipping", flush=True)
            continue

        if dest.exists() and dest.stat().st_size > LARGE_FILE_BYTES:
            print(f"  {run_acc}: large file exists, subsampling in-place...", flush=True)
            subsample_inplace(dest, seqkit)
            continue

        wget_file(s["url"], dest)
        print(f"  {run_acc}: subsampling...", flush=True)
        subsample_inplace(dest, seqkit)

    # write metadata
    with METADATA_OUT.open("w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["sample-id", "study", "disease_status"])
        for s in all_samples:
            writer.writerow([s["run_accession"], s["study"], s["disease_status"]])

    print(f"\nMetadata written to {METADATA_OUT}")
    print(f"Done. {len(all_samples)} samples ready.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
