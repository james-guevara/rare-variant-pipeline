from pathlib import Path


REPO = Path(__file__).resolve().parents[1]


def test_cohort_analysis_subworkflow_connects_rare_and_analysis_boundaries():
    source = (REPO / "workflows" / "cohort_analysis.nf").read_text()

    assert "COHORT_RARE_BURDEN_WORKFLOW" in source
    assert "ANALYSIS_DATASET_WORKFLOW" in source
    assert "tuple(true, burdens)" in source
    assert "analysis_dataset =" in source
    assert "rare_burdens =" in source


def test_cohort_analysis_entrypoint_has_portable_optional_inputs():
    source = (REPO / "cohort_analysis.nf").read_text()

    assert "params.pgs_dataset" in source
    assert "params.cnv_dataset" in source
    assert "params.missing_rare_policy = 'error'" in source
    assert "empty-analysis-pgs-dataset.tsv" in source
    assert "empty-analysis-cnv-dataset.tsv" in source
    assert "--pgs_dataset and --pgs_dictionary must be supplied together" in source
    assert "--cnv_dataset and --cnv_dictionary must be supplied together" in source
