import subprocess
import sys
from pathlib import Path

import duckdb


REPO = Path(__file__).resolve().parents[1]
SCRIPT = REPO / "scripts/apply_population_af_filter.py"


def test_population_filter_keeps_rare_and_missing_without_mutating_input(tmp_path: Path):
    source = tmp_path / "annotated.parquet"
    output = tmp_path / "eligible.parquet"
    con = duckdb.connect()
    con.execute(
        """
        COPY (
          SELECT * FROM (VALUES
            ('rare', '0.0009'), ('boundary', '0.001'),
            ('common', '0.2'), ('missing', NULL), ('dot_missing', '.')
          ) t(label, "gnomAD4.1_joint_AF")
        ) TO ? (FORMAT PARQUET)
        """,
        [str(source)],
    )
    subprocess.run(
        [sys.executable, str(SCRIPT), "--input", str(source), "--output", str(output),
         "--max-af", "0.001"],
        check=True,
    )
    assert con.execute(
        "SELECT label FROM read_parquet(?) ORDER BY label", [str(output)]
    ).fetchall() == [("dot_missing",), ("missing",), ("rare",)]
    assert con.execute(
        "SELECT COUNT(*) FROM read_parquet(?)", [str(source)]
    ).fetchone()[0] == 5
