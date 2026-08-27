import json
import subprocess
import sys
from pathlib import Path

import duckdb


REPO = Path(__file__).resolve().parents[1]


def _write_population_inputs(tmp_path: Path) -> tuple[Path, Path, Path]:
    input_path = tmp_path / "input.parquet"
    dbnsfp_dir = tmp_path / "dbnsfp"
    dbnsfp_dir.mkdir()
    dbnsfp_path = dbnsfp_dir / "chr22.parquet"

    con = duckdb.connect()
    con.execute(
        """
        COPY (
          SELECT * FROM (VALUES
            ('chr22', 100::BIGINT, 'A', 'G', 'common'),
            ('chr22', 200::BIGINT, 'C', 'T', 'rare'),
            ('chr22', 300::BIGINT, 'G', 'A', 'absent')
          ) AS t("#CHROM", POS, REF, ALT, label)
        ) TO ? (FORMAT PARQUET)
        """,
        [str(input_path)],
    )

    af_columns = [
        "gnomAD4.1_joint_AF",
        "gnomAD4.1_joint_POPMAX_AF",
        "gnomAD4.1_joint_nhomalt",
        "gnomAD4.1_joint_flag",
        "AllofUs_ALL_AF",
        "AllofUs_POPMAX_AF",
        "1000Gp3_AF",
        "ALFA_Total_AF",
        "RegeneronME_ALL_AF",
        "TOPMed_frz8_AC",
        "TOPMed_frz8_AF",
        "TOPMed_frz8_AN",
        "dbNSFP_POPMAX_AC",
        "dbNSFP_POPMAX_AF",
    ]
    select_af = ", ".join(
        [
            'af AS "gnomAD4.1_joint_AF"',
            *[f'NULL::VARCHAR AS "{column}"' for column in af_columns[1:]],
        ]
    )
    con.execute(
        f"""
        COPY (
          SELECT chrom AS "#chr", pos AS "pos(1-based)", ref, alt, {select_af}
          FROM (VALUES
            ('22', 100::BIGINT, 'A', 'G', '0.20'),
            ('22', 200::BIGINT, 'C', 'T', '0.0001')
          ) AS t(chrom, pos, ref, alt, af)
        ) TO ? (FORMAT PARQUET)
        """,
        [str(dbnsfp_path)],
    )
    return input_path, dbnsfp_dir, tmp_path / "annotated.parquet"


def test_population_af_join_annotates_without_filtering(tmp_path: Path):
    input_path, dbnsfp_dir, output_path = _write_population_inputs(tmp_path)
    config_path = tmp_path / "resources.json"
    config_path.write_text(
        json.dumps(
            {
                "cohorts": {"g2mh": {}},
                "output_base": str(tmp_path),
                "dbnsfp_af_dir": str(dbnsfp_dir),
            }
        )
    )

    subprocess.run(
        [
            sys.executable,
            str(REPO / "scripts/postprocess/join_pop_af.py"),
            "--cohort",
            "g2mh",
            "--chrom",
            "chr22",
            "--resources",
            str(config_path),
            "--input",
            str(input_path),
            "--output",
            str(output_path),
        ],
        check=True,
    )

    rows = duckdb.connect().execute(
        'SELECT label, "gnomAD4.1_joint_AF" FROM read_parquet(?) ORDER BY POS',
        [str(output_path)],
    ).fetchall()
    assert rows == [
        ("common", "0.20"),
        ("rare", "0.0001"),
        ("absent", None),
    ]


def test_cohort_af_filter_keeps_complete_annotated_output(tmp_path: Path):
    carriers = tmp_path / "carriers.parquet"
    summary = tmp_path / "summary.parquet"
    annotated = tmp_path / "annotated.parquet"
    eligible = tmp_path / "eligible.parquet"
    con = duckdb.connect()
    con.execute(
        "COPY (SELECT * FROM (VALUES (1, 'S1'), (2, 'S2')) t(record_id, sample_id)) "
        "TO ? (FORMAT PARQUET)",
        [str(carriers)],
    )
    con.execute(
        """
        COPY (
          SELECT * FROM (VALUES
            (1, 1, 200, 0.005::DOUBLE, 1, 0),
            (2, 20, 200, 0.10::DOUBLE, 20, 0)
          ) t(record_id, genotype_ac, genotype_an, genotype_af, carrier_count, hom_alt_count)
        ) TO ? (FORMAT PARQUET)
        """,
        [str(summary)],
    )

    subprocess.run(
        [
            sys.executable,
            str(REPO / "scripts/apply_cohort_af_filter.py"),
            "--input",
            str(carriers),
            "--allele-summary",
            str(summary),
            "--max-af",
            "0.01",
            "--annotated-output",
            str(annotated),
            "--eligible-output",
            str(eligible),
        ],
        check=True,
    )

    assert con.execute("SELECT COUNT(*) FROM read_parquet(?)", [str(annotated)]).fetchone()[0] == 2
    assert con.execute("SELECT COUNT(*) FROM read_parquet(?)", [str(eligible)]).fetchone()[0] == 1
    assert con.execute(
        "SELECT record_id FROM read_parquet(?)", [str(eligible)]
    ).fetchone()[0] == 1
