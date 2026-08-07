from get_ENA_metadata import ENAClient, PROJECT_RE, resolve_sample_accessions, study_query


def test_project_accession_validation_accepts_ena_study_prefixes() -> None:
    assert PROJECT_RE.fullmatch("ERP005534")
    assert PROJECT_RE.fullmatch("PRJEB46665")


def test_project_accession_validation_rejects_malformed_values() -> None:
    assert PROJECT_RE.fullmatch("ER005534") is None
    assert PROJECT_RE.fullmatch("ERP") is None
    assert PROJECT_RE.fullmatch("ERP-005534") is None


def test_run_linked_samples_are_fallback_for_empty_study_sample_search() -> None:
    sample_to_runs = {"SAME1": ["ERR1"], "SAME2": ["ERR2"]}

    assert resolve_sample_accessions([], sample_to_runs) == ["SAME1", "SAME2"]


def test_study_sample_search_takes_precedence_when_available() -> None:
    sample_to_runs = {"SAME1": ["ERR1"]}

    assert resolve_sample_accessions(["SAME2", "SAME2"], sample_to_runs) == ["SAME2"]


def test_study_query_uses_primary_and_secondary_accession_fields() -> None:
    assert study_query("PRJEB46665") == 'study_accession="PRJEB46665"'
    assert study_query("ERP005534") == 'secondary_study_accession="ERP005534"'


def test_secondary_accession_samples_come_from_run_mapping() -> None:
    client = ENAClient()

    assert client.fetch_sample_accessions("ERP005534") == []
