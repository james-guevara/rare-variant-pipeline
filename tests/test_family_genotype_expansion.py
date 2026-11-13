import subprocess
import sys
from pathlib import Path

import numpy as np
import pyarrow.parquet as pq
import zarr


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "extract_zarr_allele_genotypes.py"


def test_family_expansion_uses_fid_membership(tmp_path):
    store = zarr.open_group(str(tmp_path / "calls.zarr"), mode="w")
    store.create_array("sample_id", data=np.array(["CHILD", "MOM", "SIB"], dtype="U5"))
    store.create_array("filter_id", data=np.array(["PASS"], dtype="U4"))
    genotypes = np.array([[[0, 1], [0, 0], [-1, -1]]], dtype=np.int8)
    store.create_array("call_genotype", data=genotypes, chunks=(1, 3, 2))
    store.create_array("call_genotype_mask", data=genotypes < 0, chunks=(1, 3, 2))

    alleles = tmp_path / "alleles.tsv"
    alleles.write_text(
        "record_id\tCHROM\tPOS\tREF\tALT\tGene\tSYMBOL\tConsequence\tlof_tier\n"
        "z0_a1\t22\t100\tA\tG\tENSG1\tONE\tsynonymous_variant\tlof_t1\n"
    )
    manifest = tmp_path / "samples.tsv"
    manifest.write_text(
        "#FID\tIID\tSEX\n"
        "F1\tCHILD\tF\n"
        "F1\tMOM\tF\n"
        "F1\tSIB\tM\n"
    )
    carriers = tmp_path / "carriers.parquet"
    summary = tmp_path / "summary.parquet"
    family = tmp_path / "family.parquet"
    subprocess.run([
        sys.executable, str(SCRIPT), "--zarr", str(tmp_path / "calls.zarr"),
        "--alleles", str(alleles), "--carriers-output", str(carriers),
        "--summary-output", str(summary), "--sample-manifest", str(manifest),
        "--family-output", str(family),
    ], check=True)

    rows = {row["sample_id"]: row for row in pq.read_table(family).to_pylist()}
    assert set(rows) == {"CHILD", "MOM", "SIB"}
    assert rows["CHILD"]["family_call_state"] == "carrier"
    assert rows["CHILD"]["is_index_carrier"] is True
    assert rows["MOM"]["family_call_state"] == "reference_or_other_alt"
    assert rows["MOM"]["is_index_carrier"] is False
    assert rows["SIB"]["family_call_state"] == "missing_genotype"
