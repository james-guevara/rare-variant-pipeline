import csv
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).parents[1]


def rows(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_hc_output_retains_plofs_outside_genebayes_tiers(tmp_path):
    source = tmp_path / "loftee.tsv"
    source.write_text(
        "record_id\tGene\tSYMBOL\tLoF\n"
        "v1\tENSG1\tONE\tHC\n"
        "v2\tENSG2\tTWO\tHC\n"
        "v3\tENSG3\tTHREE\tLC\n"
    )
    genebayes = tmp_path / "genebayes.tsv"
    genebayes.write_text(
        "ensg\tobs_lof\texp_lof\tprior_mean\tpost_mean\tpost_lower_95\tpost_upper_95\n"
        "ENSG1\t0\t1\t0.2\t0.20\t0.1\t0.3\n"
        "ENSG2\t0\t1\t0.01\t0.01\t0\t0.02\n"
    )
    all_output = tmp_path / "all.tsv"
    tiered = tmp_path / "tiered.tsv"
    hc = tmp_path / "hc.tsv"
    subprocess.run([
        sys.executable, str(ROOT / "scripts" / "join_genebayes_lof_tiers.py"),
        "--input", str(source), "--genebayes", str(genebayes),
        "--output", str(all_output), "--qualifying-output", str(tiered),
        "--hc-output", str(hc),
    ], check=True)

    tiered_rows = rows(tiered)
    assert [row["record_id"] for row in tiered_rows] == ["v1"]
    assert tiered_rows[0]["lof_tier"] == "lof_t1"
    assert [row["record_id"] for row in rows(hc)] == ["v1", "v2"]
    assert [row["record_id"] for row in rows(all_output)] == ["v1", "v2", "v3"]
