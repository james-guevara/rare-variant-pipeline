from pathlib import Path


REPO = Path(__file__).resolve().parents[1]


def test_reusable_analysis_workflow_has_optional_component_contract():
    workflow = (REPO / "workflows" / "analysis_dataset.nf").read_text()
    assert "workflow ANALYSIS_DATASET_WORKFLOW" in workflow
    assert "tuple val(pgs_enabled), path(pgs_dataset), path(pgs_dictionary)" in workflow
    assert "tuple val(rare_enabled), path(rare_burdens)" in workflow
    assert "tuple val(cnv_enabled), path(cnv_dataset), path(cnv_dictionary)" in workflow
    assert "--participant-manifest" in workflow
    assert "--missing-pgs-policy" in workflow
    assert "--missing-rare-policy" in workflow
    assert "--missing-cnv-policy" in workflow
    assert "analysis_dataset.tsv" in workflow
    assert "analysis_dictionary.tsv" in workflow
    assert "analysis_qc.json" in workflow
    assert "analysis_exclusions.tsv" in workflow


def test_standalone_analysis_wrapper_validates_optional_pairs():
    wrapper = (REPO / "analysis.nf").read_text()
    assert "include { ANALYSIS_DATASET_WORKFLOW }" in wrapper
    assert "--pgs_dataset and --pgs_dictionary must be supplied together" in wrapper
    assert "--cnv_dataset and --cnv_dictionary must be supplied together" in wrapper
    assert "At least one PGS, rare-burden, or CNV dataset is required" in wrapper
    assert "empty-analysis-pgs-dataset.tsv" in wrapper
    assert "empty-analysis-cnv-dictionary.tsv" in wrapper
