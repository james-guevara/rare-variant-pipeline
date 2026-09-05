import csv
import subprocess
import sys
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq


ROOT = Path(__file__).parents[1]
SCRIPT = ROOT / "scripts" / "build_gene_set_burdens.py"


def test_gene_set_burdens_are_independent_of_lof_tiers(tmp_path):
    samples = tmp_path / "samples.tsv"
    samples.write_text("FID\tIID\nF1\tS1\nF2\tS2\n")
    memberships = tmp_path / "sets.tsv"
    memberships.write_text(
        "gene_set_id\tgene_symbol\tensembl_gene_id\n"
        "set_a\tGENE1\tENSG1\n"
        "set_a\tGENE2\t.\n"
        "set_b\tGENE1\tENSG1\n"
    )
    plof = tmp_path / "plof.parquet"
    pq.write_table(pa.Table.from_pylist([
        {"sample_id": "S1", "record_id": "v1", "Gene": "ENSG1.2", "SYMBOL": "GENE1", "LoF": "HC", "lof_tier": None},
        {"sample_id": "S1", "record_id": "v2", "Gene": "OTHER", "SYMBOL": "GENE2", "LoF": "HC", "lof_tier": None},
        {"sample_id": "S1", "record_id": "v3", "Gene": "ENSG1", "SYMBOL": "GENE1", "LoF": "LC", "lof_tier": None},
    ]), plof)
    missense = tmp_path / "miss.parquet"
    pq.write_table(pa.Table.from_pylist([
        {"sample_id": "S1", "record_id": "m1", "Gene": "ENSG1", "SYMBOL": "GENE1", "miss_tier": "miss_t1"},
        {"sample_id": "S1", "record_id": "m1", "Gene": "ENSG1", "SYMBOL": "GENE1", "miss_tier": "miss_t1"},
        {"sample_id": "S2", "record_id": "m2", "Gene": "ENSG1", "SYMBOL": "GENE1", "miss_tier": "miss_t4"},
    ]), missense)
    output = tmp_path / "burdens.tsv"

    subprocess.run([
        sys.executable, str(SCRIPT), "--samples", str(samples), "--gene-sets", str(memberships),
        "--plof", str(plof), "--missense", str(missense), "--output", str(output),
    ], check=True)
    with output.open() as handle:
        rows = {(r["IID"], r["gene_set_id"]): r for r in csv.DictReader(handle, delimiter="\t")}

    assert rows[("S1", "set_a")]["plof"] == "2"
    assert rows[("S1", "set_a")]["miss_t1"] == "1"
    assert rows[("S1", "set_b")]["plof"] == "1"
    assert rows[("S2", "set_a")]["plof"] == "0"
    assert rows[("S2", "set_a")]["miss_t4"] == "1"
    assert len(rows) == 4
