#!/usr/bin/env python3
"""Build participant-level burdens for arbitrary gene sets.

The pLoF count is independent of GeneBayes LoF tiers: every eligible HC-pLoF
carrier row in a member gene contributes. Missense counts retain the existing
four classifier tiers. A variant is counted once per participant and gene set,
while a variant may count in multiple overlapping gene sets.
"""

import argparse
import csv
from pathlib import Path

import duckdb


def sql_paths(paths: list[Path]) -> str:
    values = ", ".join("'{}'".format(str(path).replace("'", "''")) for path in paths)
    return f"[{values}]"


def read_samples(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        first = handle.readline()
        handle.seek(0)
        delimiter = "\t" if "\t" in first else None
        reader = csv.DictReader(handle, delimiter=delimiter) if delimiter else csv.DictReader(
            handle, delimiter=" ", skipinitialspace=True
        )
        fields = {name.lstrip("#"): name for name in (reader.fieldnames or [])}
        if "IID" not in fields:
            raise ValueError("sample file must contain IID (or #IID)")
        rows = []
        seen = set()
        for row in reader:
            iid = row[fields["IID"]]
            if not iid or iid in seen:
                raise ValueError(f"missing or duplicate IID: {iid!r}")
            seen.add(iid)
            rows.append({"FID": row.get(fields.get("FID", ""), ""), "IID": iid})
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--samples", required=True, type=Path)
    parser.add_argument("--gene-sets", required=True, type=Path)
    parser.add_argument("--plof", required=True, type=Path, nargs="+")
    parser.add_argument("--missense", required=True, type=Path, nargs="+")
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    samples = read_samples(args.samples)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    con = duckdb.connect()
    con.execute("CREATE TABLE samples(FID VARCHAR, IID VARCHAR)")
    con.executemany("INSERT INTO samples VALUES (?, ?)", [(r["FID"], r["IID"]) for r in samples])
    gene_sets = str(args.gene_sets).replace("'", "''")
    con.execute(
        "CREATE VIEW memberships AS SELECT DISTINCT gene_set_id, "
        "NULLIF(ensembl_gene_id, '.') AS ensembl_gene_id, "
        "upper(NULLIF(gene_symbol, '.')) AS gene_symbol "
        f"FROM read_csv('{gene_sets}', delim='\\t', header=true, all_varchar=true)",
    )

    # Prefer stable Ensembl IDs. Symbol matching is a fallback for catalog rows
    # without an Ensembl ID, preventing aliases from double-counting a gene.
    match = """(
        (m.ensembl_gene_id IS NOT NULL AND regexp_replace(c.Gene, '\\.[0-9]+$', '') = m.ensembl_gene_id)
        OR (m.ensembl_gene_id IS NULL AND m.gene_symbol = upper(c.SYMBOL))
    )"""
    con.execute(f"""
        CREATE TABLE events AS
        SELECT DISTINCT c.sample_id AS IID, m.gene_set_id, c.record_id, 'plof' AS variant_class
        FROM read_parquet({sql_paths(args.plof)}, union_by_name=true) c
        JOIN memberships m ON {match}
        WHERE c.LoF = 'HC'
        UNION
        SELECT DISTINCT c.sample_id AS IID, m.gene_set_id, c.record_id, c.miss_tier AS variant_class
        FROM read_parquet({sql_paths(args.missense)}, union_by_name=true) c
        JOIN memberships m ON {match}
        WHERE c.miss_tier IN ('miss_t1', 'miss_t2', 'miss_t3', 'miss_t4')
    """)

    unknown = con.execute(
        "SELECT count(DISTINCT IID) FROM events e WHERE NOT EXISTS "
        "(SELECT 1 FROM samples s WHERE s.IID=e.IID)"
    ).fetchone()[0]
    if unknown:
        raise ValueError(f"carrier inputs contain {unknown} sample IDs absent from sample file")

    output = str(args.output).replace("'", "''")
    con.execute(f"""
        COPY (
          WITH gene_sets AS (SELECT DISTINCT gene_set_id FROM memberships),
          counts AS (
            SELECT IID, gene_set_id, variant_class, count(*) AS burden_count
            FROM events GROUP BY ALL
          )
          SELECT s.FID, s.IID, g.gene_set_id,
                 coalesce(max(c.burden_count) FILTER (WHERE c.variant_class='plof'), 0)::BIGINT AS plof,
                 coalesce(max(c.burden_count) FILTER (WHERE c.variant_class='miss_t1'), 0)::BIGINT AS miss_t1,
                 coalesce(max(c.burden_count) FILTER (WHERE c.variant_class='miss_t2'), 0)::BIGINT AS miss_t2,
                 coalesce(max(c.burden_count) FILTER (WHERE c.variant_class='miss_t3'), 0)::BIGINT AS miss_t3,
                 coalesce(max(c.burden_count) FILTER (WHERE c.variant_class='miss_t4'), 0)::BIGINT AS miss_t4
          FROM samples s CROSS JOIN gene_sets g
          LEFT JOIN counts c ON c.IID=s.IID AND c.gene_set_id=g.gene_set_id
          GROUP BY s.FID, s.IID, g.gene_set_id
          ORDER BY s.IID, g.gene_set_id
        ) TO '{output}' (FORMAT CSV, DELIMITER '\t', HEADER)
    """)
    n_sets = con.execute("SELECT count(DISTINCT gene_set_id) FROM memberships").fetchone()[0]
    n_events = con.execute("SELECT count(*) FROM events").fetchone()[0]
    print(f"samples={len(samples):,}")
    print(f"gene_sets={n_sets:,}")
    print(f"distinct_participant_gene_set_variant_events={n_events:,}")
    print(f"output_rows={len(samples) * n_sets:,}")


if __name__ == "__main__":
    main()
