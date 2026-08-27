import subprocess
import sys
from pathlib import Path

import duckdb


REPO = Path(__file__).resolve().parents[1]
SCRIPT = REPO / "scripts/select_missense_candidates.py"


def test_requires_picked_missense_and_matching_candidate_gene(tmp_path: Path):
    picked = tmp_path / "picked.tsv"
    picked.write_text(
        "record_id\tCHROM\tPOS\tREF\tALT\tGene\tConsequence\n"
        "v1\t22\t100\tA\tG\tENSG1.5\tmissense_variant\n"
        "v2\t22\t200\tC\tT\tENSG2\tdownstream_gene_variant\n"
        "v3\t22\t300\tG\tA\tENSG_OTHER\tmissense_variant\n"
    )
    candidates = tmp_path / "candidates.parquet"
    output = tmp_path / "selected.parquet"
    con = duckdb.connect()
    con.execute(
        """
        COPY (
          SELECT * FROM (VALUES
            ('chr22', 100::BIGINT, 'A', 'G', 'ENSG1', 4, 'miss_t1', 1.0, 1.0, 1.0, 1.0),
            ('chr22', 200::BIGINT, 'C', 'T', 'ENSG2', 2, 'miss_t3', 1.0, 1.0, NULL, NULL),
            ('chr22', 300::BIGINT, 'G', 'A', 'ENSG3', 1, 'miss_t4', 1.0, NULL, NULL, NULL)
          ) t(chrom,pos,ref,alt,candidate_genes,miss_n_flag,miss_tier,
              ClinPred_rankscore,AlphaMissense_rankscore,
              popEVE_converted_rankscore,MPC_rankscore)
        ) TO ? (FORMAT PARQUET)
        """,
        [str(candidates)],
    )
    result = subprocess.run(
        [sys.executable, str(SCRIPT), "--picked", str(picked),
         "--candidates", str(candidates), "--output", str(output)],
        capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stderr
    assert "selected_missense=1" in result.stdout
    assert con.execute(
        "SELECT record_id, Gene, miss_tier FROM read_parquet(?)", [str(output)]
    ).fetchall() == [("v1", "ENSG1.5", "miss_t1")]
