from __future__ import annotations

import gzip
from pathlib import Path

from pipeline_lib import fastq_to_trimmed_fasta, parse_deblur_clean_fasta


def test_fastq_to_trimmed_fasta_keeps_only_reads_at_trim_length(tmp_path: Path) -> None:
    fastq_fp = tmp_path / "sample_1.fastq.gz"
    fasta_fp = tmp_path / "trimmed.fasta"
    trim_length = 5

    with gzip.open(fastq_fp, "wt") as fout:
        fout.write(
            "@read_short\n"
            "ACG\n"
            "+\n"
            "!!!\n"
            "@read_exact\n"
            "ACGTA\n"
            "+\n"
            "!!!!!\n"
            "@read_long\n"
            "ACGTACGT\n"
            "+\n"
            "!!!!!!!!\n"
        )

    kept = fastq_to_trimmed_fasta(fastq_fp, fasta_fp, trim_length)

    assert kept == 2
    assert fasta_fp.read_text() == ">read_exact\nACGTA\n>read_long\nACGTA\n"


def test_parse_deblur_clean_fasta_parses_sizes_and_counts_last_sequence(tmp_path: Path) -> None:
    clean_fp = tmp_path / "mock.clean"
    clean_fp.write_text(
        ">seqA;size=4;\n"
        "ACGT\n"
        ">seqB;size=2;\n"
        "TGCA\n"
        ">seqC;size=7;\n"
        "ACGT\n"
    )

    counts = parse_deblur_clean_fasta(clean_fp)

    assert counts["ACGT"] == 11
    assert counts["TGCA"] == 2
