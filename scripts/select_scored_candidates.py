#!/usr/bin/env python3
"""Candidate-only DuckDB joins; no genotype reads and no full scored annotation export."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import time

import duckdb


def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for data in iter(lambda: f.read(8*1024*1024), b''):
            h.update(data)
    return h.hexdigest()


def lit(value):
    return "'"+str(value).replace("'", "''")+"'"


def ident(value):
    return '"'+value.replace('"','""')+'"'


def check(item):
    stat = Path(item['path']).stat()
    if (stat.st_size, stat.st_mtime_ns) != (item['bytes'], item['mtime_ns']):
        raise ValueError('Shared resource changed since checksum verification')


def run(a):
    sys.path.insert(0, str(Path(a.postprocess_dir).resolve()))
    from tier_variants import T_STARS, LOF_T1, LOF_T2
    from join_scores import SCORES, extract_expr
    meta = json.loads(Path(a.metadata).read_text())
    out = Path(a.outdir); out.mkdir(parents=True, exist_ok=True)
    receipt = dict(status='failed', unit_id=meta['unit_id'], chromosome=meta['chromosome'],
                   sources=meta, duckdb_version=duckdb.__version__,
                   thresholds=T_STARS, lof_thresholds=dict(lof_t1=LOF_T1, lof_t2=LOF_T2),
                   selection='missense consequence token and n_flag >= 1; all HC LoF rows retained',
                   allele_policy='sequence and spanning_deletion labelled separately; no biological star-burden decision',
                   genebayes_join='Gene == ensg, exact; no symbol fallback or gene-ID rewriting',
                   dbnsfp_policy='unfiltered parquet_expanded; no MANE inclusion filter; exact allele join; no picked-transcript or gene match',
                   score_policy='reuse join_scores.SCORES/extract_expr; ranks scalar; raw list MAX except popEVE MIN; duplicate exact keys aggregate as existing code')
    start = time.perf_counter()
    con = None
    try:
        if meta.get('dbnsfp_representation') != 'parquet_expanded':
            raise ValueError('Candidate scoring requires a rebuilt unfiltered parquet_expanded lock')
        for resource in [meta['dbnsfp'], meta['genebayes']]: check(resource)
        con = duckdb.connect()
        con.execute(f'SET threads={int(a.threads)}')
        con.execute('SET memory_limit='+lit(a.memory))
        con.execute('SET temp_directory='+lit(out/'.duckdb-spill'))
        def count(sql): return con.execute(sql).fetchone()[0]
        def tsv(path):
            return f"read_csv({lit(path)}, delim='\t', header=true, all_varchar=true, nullstr=['','.'])"
        for name, path in [('picked', a.picked), ('loftee', a.loftee)]:
            con.execute(f'CREATE VIEW {name}_input AS SELECT * FROM {tsv(path)}')
            fields = {r[0] for r in con.execute(f'DESCRIBE {name}_input').fetchall()}
            required = {'CHROM','POS','REF','ALT','Gene','Feature','Consequence'}
            if name == 'loftee': required.add('LoF')
            if not required <= fields:
                raise ValueError('Required annotation fields missing')
        normalize = "regexp_replace(CHROM, '^chr', '')"
        allele_class = "CASE WHEN ALT='*' THEN 'spanning_deletion' ELSE 'sequence' END"
        # Only these reduced relations are materialized. Input views never produce a full scored copy.
        con.execute(f"CREATE TEMP TABLE miss AS SELECT * REPLACE(TRY_CAST(POS AS BIGINT) AS POS), {normalize} AS chrom_key, {allele_class} AS allele_class FROM picked_input WHERE list_contains(string_split(Consequence,'&'),'missense_variant')")
        con.execute(f"CREATE TEMP TABLE hc AS SELECT * REPLACE(TRY_CAST(POS AS BIGINT) AS POS), {normalize} AS chrom_key, {allele_class} AS allele_class FROM loftee_input WHERE LoF='HC'")
        for table in ['miss', 'hc']:
            if count(f"SELECT count(*) FROM {table} WHERE chrom_key IS DISTINCT FROM {lit(meta['chromosome'].removeprefix('chr'))} OR POS IS NULL OR POS < 1 OR REF IS NULL OR NOT regexp_full_match(REF,'[ACGTN]+') OR ALT IS NULL OR NOT (ALT='*' OR regexp_full_match(ALT,'[ACGTN]+')) OR REF=ALT"):
                raise ValueError('Invalid candidate chromosome/position/biallelic allele')
            if count(f'SELECT count(*) FROM (SELECT chrom_key,POS,REF,ALT FROM {table} GROUP BY ALL HAVING count(*)>1)'):
                raise ValueError('Duplicate candidate exact allele key')
        if count("SELECT count(*) FROM hc h ANTI JOIN picked_input p ON h.chrom_key=regexp_replace(p.CHROM,'^chr','') AND h.POS=TRY_CAST(p.POS AS BIGINT) AND h.REF=p.REF AND h.ALT=p.ALT AND h.Gene IS NOT DISTINCT FROM p.Gene AND h.Feature IS NOT DISTINCT FROM p.Feature"):
            raise ValueError('HC annotations do not match picked allele/gene/transcript inputs')
        resource = lit(meta['dbnsfp']['path'])
        per_scores = ', '.join(f'{extract_expr(src, encoding, agg)} AS {ident(target)}' for src,target,encoding,agg in SCORES)
        con.execute(f'''CREATE TEMP TABLE matched_scores AS
            SELECT regexp_replace(d."#chr",'^chr','') AS chrom_key, TRY_CAST(d."pos(1-based)" AS BIGINT) AS POS,
                   d.ref AS REF, d.alt AS ALT, {per_scores}
            FROM read_parquet({resource}) d
            SEMI JOIN miss m ON regexp_replace(d."#chr",'^chr','')=m.chrom_key
               AND TRY_CAST(d."pos(1-based)" AS BIGINT)=m.POS AND d.ref=m.REF AND d.alt=m.ALT''')
        aggregate = ', '.join(f'{agg}({ident(target)}) AS {ident(target)}' for src,target,encoding,agg in SCORES)
        con.execute(f'CREATE TEMP TABLE scores AS SELECT chrom_key,POS,REF,ALT,count(*) AS dbnsfp_source_rows,{aggregate} FROM matched_scores GROUP BY chrom_key,POS,REF,ALT')
        columns = ', '.join('s.'+ident(target) for src,target,encoding,agg in SCORES)
        flags = ' + '.join(f'CAST(COALESCE({ident(col)}>={threshold},FALSE) AS INTEGER)' for col,threshold in T_STARS.items())
        available = ' + '.join(f'CAST({ident(col)} IS NOT NULL AS INTEGER)' for col in T_STARS)
        con.execute(f'''CREATE TEMP TABLE miss_scored AS
            WITH joined AS (SELECT m.*,s.dbnsfp_source_rows IS NOT NULL AS dbnsfp_matched,
                           COALESCE(s.dbnsfp_source_rows,0) AS dbnsfp_source_rows,{columns}
                           FROM miss m LEFT JOIN scores s USING(chrom_key,POS,REF,ALT)),
                 flagged AS (SELECT *, {flags} AS n_flag, {available} AS n_scored FROM joined)
            SELECT *, n_flag AS miss_n_flag,
                   CASE n_flag WHEN 4 THEN 'miss_t1' WHEN 3 THEN 'miss_t2' WHEN 2 THEN 'miss_t3' WHEN 1 THEN 'miss_t4' END AS tier
            FROM flagged''')
        gb = ['obs_lof','exp_lof','prior_mean','post_mean','post_lower_95','post_upper_95']
        projection = ', '.join(f'TRY_CAST({ident(c)} AS DOUBLE) AS genebayes_{c}' for c in gb)
        con.execute(f'CREATE TEMP TABLE gb AS SELECT ensg,{projection} FROM {tsv(meta["genebayes"]["path"])} WHERE ensg IS NOT NULL')
        if count('SELECT count(*) FROM (SELECT ensg FROM gb GROUP BY ensg HAVING count(*)>1)'):
            raise ValueError('Duplicate GeneBayes gene IDs would multiply candidates')
        con.execute(f'''CREATE TEMP TABLE hc_scored AS
            SELECT h.*, b.ensg IS NOT NULL AS genebayes_matched,
                   {', '.join('b.genebayes_'+c for c in gb)},
                   CASE WHEN b.genebayes_post_mean >= {LOF_T1} THEN 'lof_t1'
                        WHEN b.genebayes_post_mean >= {LOF_T2} AND b.genebayes_post_mean < {LOF_T1} THEN 'lof_t2' END AS tier
            FROM hc h LEFT JOIN gb b ON h.Gene=b.ensg''')
        outputs = [('missense.parquet','miss_scored','n_flag>=1'),('lof_hc.parquet','hc_scored','true')]
        for name, table, predicate in outputs:
            con.execute(f'COPY (SELECT * EXCLUDE(chrom_key) FROM {table} WHERE {predicate} ORDER BY chrom_key,POS,REF,ALT) TO {lit(out/name)} (FORMAT PARQUET, COMPRESSION ZSTD)')
        by_class = {}
        for kind in ['sequence','spanning_deletion']:
            where = 'allele_class='+lit(kind)
            by_class[kind] = dict(
                missense_input_rows=count(f'SELECT count(*) FROM miss WHERE {where}'),
                dbnsfp_matches=count(f'SELECT count(*) FROM miss_scored WHERE {where} AND dbnsfp_matched'),
                dbnsfp_nonmatches=count(f'SELECT count(*) FROM miss_scored WHERE {where} AND NOT dbnsfp_matched'),
                n_flag_counts={str(n):count(f'SELECT count(*) FROM miss_scored WHERE {where} AND n_flag={n}') for n in range(5)},
                n_scored_counts={str(n):count(f'SELECT count(*) FROM miss_scored WHERE {where} AND n_scored={n}') for n in range(5)},
                rankscore_missing_counts={col:count(f'SELECT count(*) FROM miss_scored WHERE {where} AND {ident(col)} IS NULL') for col in T_STARS},
                matched_without_rankscores=count(f'SELECT count(*) FROM miss_scored WHERE {where} AND dbnsfp_matched AND n_scored=0'),
                fully_scored_below_thresholds=count(f'SELECT count(*) FROM miss_scored WHERE {where} AND n_scored=4 AND n_flag=0'),
                missense_selected_rows=count(f'SELECT count(*) FROM miss_scored WHERE {where} AND n_flag>=1'),
                hc_rows=count(f'SELECT count(*) FROM hc_scored WHERE {where}'),
                genebayes_matches=count(f'SELECT count(*) FROM hc_scored WHERE {where} AND genebayes_matched'),
                genebayes_nonmatches=count(f'SELECT count(*) FROM hc_scored WHERE {where} AND NOT genebayes_matched'),
                lof_tier_counts={tier:count(f'SELECT count(*) FROM hc_scored WHERE {where} AND '+('tier IS NULL' if tier=='untiered' else 'tier='+lit(tier))) for tier in ['lof_t1','lof_t2','untiered']})
        for resource in [meta['dbnsfp'], meta['genebayes']]: check(resource)
        receipt.update(status='passed', picked_input_rows=count('SELECT count(*) FROM picked_input'),
                       loftee_input_rows=count('SELECT count(*) FROM loftee_input'),
                       missense_input_rows=count('SELECT count(*) FROM miss'), hc_rows=count('SELECT count(*) FROM hc'),
                       matched_dbnsfp_resource_rows=count('SELECT count(*) FROM matched_scores'),
                       overlapping_missense_hc_keys=count('SELECT count(*) FROM miss SEMI JOIN hc USING(chrom_key,POS,REF,ALT)'),
                       duplicate_dbnsfp_keys=count('SELECT count(*) FROM scores WHERE dbnsfp_source_rows>1'),
                       by_allele_class=by_class,
                       input_sha256=dict(picked=sha(a.picked), loftee=sha(a.loftee)),
                       reused_code_sha256={name:sha(Path(a.postprocess_dir)/name) for name in ['tier_variants.py','join_scores.py']},
                       outputs={name:dict(bytes=(out/name).stat().st_size,sha256=sha(out/name),rows=count(f'SELECT count(*) FROM {table} WHERE {predicate}')) for name,table,predicate in outputs})
    except Exception as exc:
        receipt['error_type'] = type(exc).__name__
        for name in ['missense.parquet','lof_hc.parquet']:
            (out/name).unlink(missing_ok=True)
        raise
    finally:
        if con is not None: con.close()
        receipt['wall_seconds'] = time.perf_counter()-start
        (out/'receipt.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n')
    print(json.dumps({k:receipt[k] for k in ['status','unit_id','missense_input_rows','hc_rows','by_allele_class']}))


if __name__ == '__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['picked','loftee','metadata','postprocess-dir','outdir']: p.add_argument('--'+name,required=True)
    p.add_argument('--threads',type=int,default=2)
    p.add_argument('--memory',default='3GB')
    try: run(p.parse_args())
    except Exception: sys.exit('Candidate selection failed; inspect the task receipt. No individual records are printed.')
