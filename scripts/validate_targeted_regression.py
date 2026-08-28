#!/usr/bin/env python3
"""Validate current LoF/missense outputs against pinned counts and hashes."""

import argparse
import hashlib
import json
from pathlib import Path

import duckdb


def file_sha256(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def normalized(value):
    if value is None:
        return None
    if isinstance(value, float):
        return format(value, ".17g")
    return str(value)


def row_digest(rows) -> str:
    payload = json.dumps(
        sorted(tuple(normalized(value) for value in row) for row in rows),
        separators=(",", ":"),
    ).encode()
    return hashlib.sha256(payload).hexdigest()


def quote(value: str) -> str:
    return '"{}"'.format(value.replace('"', '""'))


def columns(connection, path: Path) -> set[str]:
    return {
        row[0]
        for row in connection.execute(
            "DESCRIBE SELECT * FROM read_parquet(?)", [str(path)]
        ).fetchall()
    }


def canonical_hashes(connection, path: Path, tier_column: str) -> dict[str, str]:
    fixed = [
        "#CHROM", "POS", "REF", "ALT", "sample_id", "GT", "GQ", "DP",
        "AD", "FILTER", "Gene", "SYMBOL", "Consequence", "LoF",
        tier_column, "genebayes_post_mean",
    ]
    rows = connection.execute(
        "SELECT {} FROM read_parquet(?)".format(",".join(map(quote, fixed))),
        [str(path)],
    ).fetchall()
    variants = {
        tuple(row[index] for index in (0, 1, 2, 3, 10, 11, 12, 13, 14, 15))
        for row in rows
    }
    coordinates = {tuple(row[index] for index in (1, 2, 3)) for row in rows}
    variant_core = {
        tuple(row[index] for index in (1, 2, 3, 10, 11, 12, 13, 14))
        for row in rows
    }
    carrier_core = {
        tuple(row[index] for index in (1, 2, 3, 4, 5, 14)) for row in rows
    }
    return {
        "carrier_sha256": row_digest(rows),
        "variant_sha256": row_digest(variants),
        "coordinate_sha256": row_digest(coordinates),
        "variant_core_sha256": row_digest(variant_core),
        "carrier_core_sha256": row_digest(carrier_core),
    }


def count_rows(connection, path: Path, predicate: str = "TRUE") -> int:
    return connection.execute(
        f"SELECT count(*) FROM read_parquet(?) WHERE {predicate}", [str(path)]
    ).fetchone()[0]


def eligible_predicates(connection, path: Path) -> tuple[str, str]:
    available = columns(connection, path)
    if "primary_analysis_eligible" not in available:
        return "TRUE", "FALSE"
    primary = "COALESCE(primary_analysis_eligible, FALSE)"
    sensitivity = f"NOT {primary}"
    if "burden_count_available" in available:
        sensitivity += " AND COALESCE(burden_count_available, FALSE)"
    return primary, sensitivity


def add_legacy_missense_aliases(observed: dict[str, int]) -> None:
    """Expose names used by the existing pinned regression documents."""
    for suffix in (
        "burden_eligible_rows", "burden_eligible_alleles",
        "burden_eligible_samples", "primary_burden_rows", "sensitivity_only_rows",
    ):
        observed[suffix] = observed[f"missense_{suffix}"]
    observed["sensitivity_only_burden_rows"] = observed[
        "missense_sensitivity_only_rows"
    ]
    for tier in ("miss_t1", "miss_t2", "miss_t3", "miss_t4"):
        observed[f"primary_{tier}_rows"] = observed[f"missense_primary_{tier}_rows"]


def observe(args, expected_document: dict) -> tuple[dict[str, int], dict[str, dict]]:
    bindings = json.loads(args.bindings.read_text())["resources"]
    run = args.run_root
    if not (run / "_SUCCESS").is_file():
        raise ValueError(f"run is not complete: {run}")
    for name in ("target_bed", "missense_candidates"):
        if f"{name}_sha256" not in expected_document:
            continue
        observed_digest = file_sha256(Path(bindings[name]))
        expected_digest = expected_document[f"{name}_sha256"]
        if observed_digest != expected_digest:
            raise ValueError(
                f"{name} checksum differs: expected={expected_digest} observed={observed_digest}"
            )

    connection = duckdb.connect()
    pq = lambda name: run / name
    scalar = lambda query, values=(): connection.execute(query, values).fetchone()[0]
    target_bed = Path(bindings["target_bed"])
    candidates = Path(bindings["missense_candidates"])
    observed = {
        "target_intervals": sum(
            bool(line and not line.startswith("#"))
            for line in target_bed.read_text().splitlines()
        ),
        "targeted_alleles": count_rows(connection, pq("01.target-alleles.parquet")),
        "plof_rows": scalar(
            "SELECT count(*) FROM read_csv(?, delim='\\t', header=true)",
            [str(pq("06.plof-genebayes.tsv"))],
        ),
        "qualifying_lof": scalar(
            "SELECT count(*) FROM read_csv(?, delim='\\t', header=true)",
            [str(pq("06.plof-tiered.tsv"))],
        ),
        "observed_missense_candidates": count_rows(connection, candidates),
        "selected_missense": count_rows(connection, pq("06.missense-tiered.parquet")),
    }
    for tier in ("lof_t1", "lof_t2"):
        observed[tier] = scalar(
            "SELECT count(*) FROM read_csv(?, delim='\\t', header=true) WHERE lof_tier=?",
            [str(pq("06.plof-tiered.tsv")), tier],
        )
    for tier in ("miss_t1", "miss_t2", "miss_t3", "miss_t4"):
        observed[f"selected_{tier}"] = count_rows(
            connection, pq("06.missense-tiered.parquet"), f"miss_tier='{tier}'"
        )

    hashes = {}
    for branch, stem, tier_column in (
        ("lof", "plof", "lof_tier"),
        ("missense", "missense", "miss_tier"),
    ):
        raw = pq(f"07.{stem}-tiered.carriers.parquet")
        region = pq(f"08.{stem}-region-filtered.parquet")
        genotype = pq(f"09.{stem}-genotype-qc.parquet")
        population = pq(f"10.{stem}-population-af-eligible.parquet")
        burden = pq(f"11.{stem}-burden-eligible.parquet")
        prefix = "lof_" if branch == "lof" else "missense_"
        observed[f"{prefix}raw_carrier_rows"] = count_rows(connection, raw)
        observed[f"{prefix}region_filtered_rows"] = count_rows(connection, region)
        observed[f"{prefix}genotype_qc_rows"] = count_rows(connection, genotype)
        population_key = (
            "lof_population_af_eligible_rows"
            if branch == "lof" else "population_af_eligible_rows"
        )
        observed[population_key] = count_rows(connection, population)
        observed[f"{prefix}burden_eligible_rows"] = count_rows(connection, burden)
        observed[f"{prefix}burden_eligible_alleles"] = scalar(
            "SELECT count(DISTINCT record_id) FROM read_parquet(?)", [str(burden)]
        )
        observed[f"{prefix}burden_eligible_samples"] = scalar(
            "SELECT count(DISTINCT sample_id) FROM read_parquet(?)", [str(burden)]
        )
        primary, sensitivity = eligible_predicates(connection, burden)
        observed[f"{prefix}primary_burden_rows"] = count_rows(connection, burden, primary)
        observed[f"{prefix}sensitivity_only_rows"] = count_rows(
            connection, burden, sensitivity
        )
        tiers = (
            ("lof_t1", "lof_t2")
            if branch == "lof"
            else ("miss_t1", "miss_t2", "miss_t3", "miss_t4")
        )
        for tier in tiers:
            primary_key = (
                f"lof_primary_{tier.removeprefix('lof_')}_rows"
                if branch == "lof"
                else f"missense_primary_{tier}_rows"
            )
            observed[primary_key] = count_rows(
                connection, burden, f"({primary}) AND {tier_column}='{tier}'"
            )
        hashes[branch] = canonical_hashes(connection, burden, tier_column)

    # Older expectation files use unprefixed names for the missense final stage.
    add_legacy_missense_aliases(observed)
    return observed, hashes


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--regression", required=True, type=Path)
    parser.add_argument("--bindings", required=True, type=Path)
    parser.add_argument("--run-root", required=True, type=Path)
    args = parser.parse_args()
    expected_document = json.loads(args.regression.read_text())
    observed, hashes = observe(args, expected_document)
    expected = expected_document["expected_counts"]
    differences = {
        key: {"expected": value, "observed": observed.get(key)}
        for key, value in expected.items()
        if observed.get(key) != value
    }
    for branch in ("lof", "missense"):
        expected_hashes = expected_document.get(f"canonical_{branch}_hashes", {})
        for key, value in expected_hashes.items():
            if hashes[branch].get(key) != value:
                differences[f"canonical_{branch}_hashes.{key}"] = {
                    "expected": value,
                    "observed": hashes[branch].get(key),
                }
    if differences:
        raise ValueError(f"regression differs: {json.dumps(differences, sort_keys=True)}")
    hash_count = sum(
        len(expected_document.get(f"canonical_{branch}_hashes", {}))
        for branch in ("lof", "missense")
    )
    (args.run_root / "_REGRESSION_SUCCESS").write_text(
        json.dumps(
            {
                "regression": str(args.regression),
                "count_checks": len(expected),
                "hash_checks": hash_count,
            },
            sort_keys=True,
        ) + "\n"
    )
    print(
        f"regression=PASS chromosome={args.regression.stem} "
        f"counts={len(expected)} hashes={hash_count}"
    )


if __name__ == "__main__":
    main()
