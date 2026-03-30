from pathlib import Path

from src.run_deblur_pcoa import parse_deblur_clean_fasta


def test_parse_deblur_clean_fasta_reads_size_with_trailing_semicolon(tmp_path: Path) -> None:
    clean_fp = tmp_path / "sample.derep.clean"
    clean_fp.write_text(
        ">seq1;size=7;\n"
        "ACGT\n"
        ">seq2;size=3;\n"
        "TGCA\n"
    )

    counts = parse_deblur_clean_fasta(clean_fp)

    assert counts["ACGT"] == 7
    assert counts["TGCA"] == 3


def test_parse_deblur_clean_fasta_still_supports_no_trailing_semicolon(tmp_path: Path) -> None:
    clean_fp = tmp_path / "sample.derep.clean"
    clean_fp.write_text(
        ">seq1;size=5\n"
        "ACGT\n"
    )

    counts = parse_deblur_clean_fasta(clean_fp)

    assert counts["ACGT"] == 5
