import csv
import subprocess
import sys
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "select_synonymous_tiered.py"


def test_selects_synonymous_in_gene_tiers(tmp_path):
    picked = tmp_path / "picked.tsv"
    picked.write_text(
        "record_id\tGene\tSYMBOL\tConsequence\n"
        "v1\tENSG1.2\tONE\tsynonymous_variant\n"
        "v2\tENSG2\tTWO\tmissense_variant\n"
        "v3\tENSG3\tTHREE\tsynonymous_variant&splice_region_variant\n"
    )
    gb = tmp_path / "gb.tsv"
    gb.write_text(
        "ensg\tobs_lof\texp_lof\tprior_mean\tpost_mean\tpost_lower_95\tpost_upper_95\n"
        "ENSG1\t1\t2\t0.2\t0.18\t0.1\t0.3\n"
        "ENSG2\t1\t2\t0.1\t0.10\t0.1\t0.2\n"
        "ENSG3\t1\t2\t0.01\t0.02\t0\t0.1\n"
    )
    output = tmp_path / "out.tsv"
    subprocess.run([sys.executable, str(SCRIPT), "--picked", str(picked), "--genebayes", str(gb), "--output", str(output)], check=True)
    rows = list(csv.DictReader(output.open(), delimiter="\t"))
    assert [(row["record_id"], row["lof_tier"]) for row in rows] == [("v1", "lof_t1")]
