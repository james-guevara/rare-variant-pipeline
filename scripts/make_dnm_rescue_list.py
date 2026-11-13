#!/usr/bin/env python3
"""Build the DNM rescue list for post-rarity QC from one or more DNM callset Parquets (e.g. dnm-callsets-revised-20261021).
Writes sample,CHROM,POS,REF,ALT with both the raw and the normalized allele of each call (deduplicated, sorted).
Usage: python make_dnm_rescue_list.py --out rescue.tsv.gz calls.parquet [calls.parquet ...]"""
import argparse
import gzip
import duckdb

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--out', required=True); p.add_argument('calls', nargs='+')
a = p.parse_args()
src = "[" + ",".join("'" + c.replace("'", "''") + "'" for c in a.calls) + "]"
rows = duckdb.sql(f"""SELECT DISTINCT sample, CHROM, POS, REF, ALT FROM (
        SELECT child_id AS sample, chromosome AS CHROM, CAST(position AS BIGINT) AS POS, ref AS REF, alt AS ALT FROM read_parquet({src})
        UNION ALL SELECT child_id, chromosome, CAST(norm_pos AS BIGINT), norm_ref, norm_alt FROM read_parquet({src}) WHERE norm_pos IS NOT NULL)
    WHERE sample IS NOT NULL AND REF IS NOT NULL AND ALT IS NOT NULL ORDER BY CHROM, POS, REF, ALT, sample""").fetchall()
with gzip.GzipFile(a.out, 'wb', mtime=0) as raw:
    raw.write(b"sample\tCHROM\tPOS\tREF\tALT\n")
    for r in rows:
        raw.write(("\t".join(str(x) for x in r) + "\n").encode())
print(f"{len(rows)} rescue keys -> {a.out}")
