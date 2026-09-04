import subprocess
import sys
from pathlib import Path

import duckdb


SCRIPT = Path(__file__).parents[1] / "scripts" / "build_missense_candidate_alleles.py"


def test_observed_vcf_join_uses_normalized_allele(tmp_path):
    dbnsfp = tmp_path / "dbnsfp.parquet"
    output = tmp_path / "observed.parquet"
    vcf = tmp_path / "normalized.vcf"
    duckdb.connect().execute(
        """
        COPY (SELECT '1' AS "#chr", '100' AS "pos(1-based)",
                     'A' AS ref, 'G' AS alt, 'ENSG1' AS Ensembl_geneid,
                     '0.5' AS ClinPred_rankscore,
                     '0.1' AS AlphaMissense_rankscore,
                     '0.1' AS popEVE_converted_rankscore,
                     '0.1' AS MPC_rankscore)
        TO ? (FORMAT PARQUET)
        """,
        [str(dbnsfp)],
    )
    vcf.write_text(
        "##fileformat=VCFv4.3\n"
        "##source=test\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "1\t100\tz1_a2\tA\tG\t.\tPASS\t.\n"
    )

    subprocess.run(
        [
            sys.executable, str(SCRIPT), "--dbnsfp", str(dbnsfp),
            "--all-genes", "--chrom", "chr1", "--observed-vcf", str(vcf),
            "--observed-output", str(output),
        ],
        check=True,
    )

    row = duckdb.connect().execute(
        "SELECT chrom, pos, ref, alt, miss_tier FROM read_parquet(?)",
        [str(output)],
    ).fetchone()
    assert row == ("1", 100, "A", "G", "miss_t4")


def test_score_thresholds_are_configurable(tmp_path):
    dbnsfp = tmp_path / "dbnsfp.parquet"
    output = tmp_path / "candidates.parquet"
    duckdb.connect().execute(
        '''
        COPY (SELECT '1' AS "#chr", '100' AS "pos(1-based)",
                     'A' AS ref, 'G' AS alt, 'ENSG1' AS Ensembl_geneid,
                     '0.5' AS ClinPred_rankscore,
                     '0.5' AS AlphaMissense_rankscore,
                     '0.5' AS popEVE_converted_rankscore,
                     '0.5' AS MPC_rankscore)
        TO ? (FORMAT PARQUET)
        ''',
        [str(dbnsfp)],
    )
    subprocess.run(
        [
            sys.executable, str(SCRIPT), "--dbnsfp", str(dbnsfp),
            "--all-genes", "--chrom", "chr1", "--output", str(output),
            "--clinpred-rankscore-min", "0.4",
            "--alphamissense-rankscore-min", "0.4",
            "--popeve-converted-rankscore-min", "0.4",
            "--mpc-rankscore-min", "0.4",
        ],
        check=True,
    )
    tier = duckdb.connect().execute(
        "SELECT miss_tier FROM read_parquet(?)", [str(output)]
    ).fetchone()[0]
    assert tier == "miss_t1"
