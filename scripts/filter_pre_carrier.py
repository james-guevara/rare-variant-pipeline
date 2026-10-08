#!/usr/bin/env python3
"""Filter existing candidate Parquets using sites-only AC/AN, POPmax, and BEDs."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import re
import time
import duckdb

TRACKS = ('genomicSuperDups', 'simpleRepeat', 'rmsk')
AF = 'gnomAD4.1_joint_POPMAX_AF'


def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(8*1024*1024), b''): h.update(chunk)
    return h.hexdigest()


def lit(x): return "'" + str(x).replace("'", "''") + "'"
def chrom(x): return str(x).removeprefix('chr')


def unchanged(item):
    s = Path(item['path']).stat()
    if (s.st_size, s.st_mtime_ns) != (item['bytes'], item['mtime_ns']):
        raise ValueError('Shared resource changed after locking')


def read_sites(path, keys):
    """Reject genotyped input at its header, before reading any data records."""
    found = {}
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt') as f:
        header = False
        for line in f:
            if line.startswith('##'): continue
            if line.startswith('#CHROM'):
                if line.rstrip('\n\r').split('\t') != ['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO']:
                    raise ValueError('Expected genotype-free eight-column sites VCF')
                header = True
                break
            raise ValueError('Invalid sites VCF header')
        if not header: raise ValueError('Missing sites VCF header')
        for line in f:
            fields = line.rstrip('\n\r').split('\t')
            if len(fields) != 8: raise ValueError('Expected eight sites-only columns')
            key = (chrom(fields[0]), int(fields[1]), fields[3], fields[4])
            if key not in keys: continue
            if key in found: raise ValueError('Duplicate exact allele in sites VCF')
            info = {}
            for token in fields[7].split(';'):
                name, sep, value = token.partition('=')
                if name in ('AC','AN','AF'):
                    if name in info: raise ValueError('Duplicate frequency INFO field')
                    info[name] = value if sep else '.'
            ac, an = info.get('AC'), info.get('AN')
            valid = bool(ac and an and re.fullmatch(r'[0-9]+', ac) and re.fullmatch(r'[0-9]+', an))
            ac, an = (int(ac), int(an)) if valid else (None, None)
            valid = valid and an > 0 and ac <= an
            found[key] = (ac, an, ac/an if valid else None, info.get('AC'), info.get('AN'), info.get('AF'))
    return found


def read_regions(path, track, chromosome):
    rows = []
    with open(path) as f:
        for line in f:
            if not line.strip() or line.startswith(('#','track ','browser ')): continue
            fields = line.rstrip('\n\r').split('\t')
            if len(fields) < (7 if track == 'rmsk' else 3):
                raise ValueError('BED schema mismatch; rmsk repClass must be column 7')
            start, end = int(fields[1]), int(fields[2])
            if start < 0 or end < start: raise ValueError('Invalid BED interval')
            if chrom(fields[0]) != chromosome: continue
            if track == 'rmsk' and fields[6] not in ('Simple_repeat','Low_complexity'): continue
            rows.append((start, end))
    return rows


def run(a):
    meta = json.loads(Path(a.metadata).read_text())
    out = Path(a.outdir).resolve(); out.mkdir(parents=True, exist_ok=True)
    inputs = {'missense': Path(a.missense).resolve(), 'lof_hc': Path(a.lof_hc).resolve(), 'sites': Path(a.sites).resolve()}
    resources = [meta['popmax'], *meta['regions'].values()]
    output_names = ['missense.filtered.parquet','lof_hc.filtered.parquet','filter_audit.parquet','receipt.json']
    protected = set(inputs.values()) | {Path(r['path']).resolve() for r in resources}
    if any(out/name in protected for name in output_names): raise ValueError('Output would overwrite an input')
    receipt = dict(status='failed', unit_id=meta['unit_id'], chromosome=meta['chromosome'],
                   policy=dict(cohort='uncorrected INFO/AC / INFO/AN < 0.005; preliminary eligibility only; missing/invalid fails',
                               gnomad=AF+' missing OR < 0.001', regions='BED start < POS <= end; POS only',
                               rmsk_classes=['Simple_repeat','Low_complexity']), sources=meta)
    start = time.perf_counter(); con = None
    try:
        for item in resources: unchanged(item)
        con = duckdb.connect(); con.execute(f'SET threads={int(a.threads)}')
        con.execute('SET memory_limit='+lit(a.memory)); con.execute('SET temp_directory='+lit(out/'.duckdb-spill'))
        def count(sql): return con.execute(sql).fetchone()[0]
        for kind in ('missense','lof_hc'):
            con.execute(f'CREATE TABLE {kind} AS SELECT * FROM read_parquet({lit(inputs[kind])})')
            fields = {r[0] for r in con.execute(f'DESCRIBE {kind}').fetchall()}
            if not {'CHROM','POS','REF','ALT','allele_class'} <= fields or any(c.startswith('pcf_') for c in fields):
                raise ValueError('Candidate schema missing keys or contains reserved pcf_ fields')
            if count(f"SELECT count(*) FROM {kind} WHERE regexp_replace(CHROM,'^chr','') IS DISTINCT FROM {lit(chrom(meta['chromosome']))} OR TRY_CAST(POS AS BIGINT) IS NULL OR TRY_CAST(POS AS BIGINT)<1 OR REF IS NULL OR NOT regexp_full_match(REF,'[ACGTN]+') OR ALT IS NULL OR NOT (ALT='*' OR regexp_full_match(ALT,'[ACGTN]+')) OR allele_class IS DISTINCT FROM CASE WHEN ALT='*' THEN 'spanning_deletion' ELSE 'sequence' END"):
                raise ValueError('Invalid candidate identity')
            if count(f"SELECT count(*) FROM (SELECT regexp_replace(CHROM,'^chr',''),TRY_CAST(POS AS BIGINT),REF,ALT FROM {kind} GROUP BY ALL HAVING count(*)>1)"):
                raise ValueError('Duplicate candidate exact allele')
        con.execute("CREATE TABLE keys AS SELECT regexp_replace(CHROM,'^chr','') AS c,CAST(POS AS BIGINT) AS p,REF AS r,ALT AS a FROM missense UNION SELECT regexp_replace(CHROM,'^chr',''),CAST(POS AS BIGINT),REF,ALT FROM lof_hc")
        keys = set(con.execute('SELECT * FROM keys').fetchall())
        sites = read_sites(inputs['sites'], keys)
        hashes = {k: sha(v) for k,v in inputs.items()}
        con.execute('CREATE TABLE sites(c VARCHAR,p BIGINT,r VARCHAR,a VARCHAR,ac BIGINT,an BIGINT,af DOUBLE,source_ac VARCHAR,source_an VARCHAR,source_af VARCHAR)')
        if sites: con.executemany('INSERT INTO sites VALUES (?,?,?,?,?,?,?,?,?,?)', [(*key,*value) for key,value in sites.items()])
        con.execute(f'''CREATE TABLE pop_raw AS SELECT regexp_replace(d."#chr",'^chr','') AS c,
                       TRY_CAST(d."pos(1-based)" AS BIGINT) AS p,d.ref AS r,d.alt AS a,
                       NULLIF(NULLIF(trim(CAST(d."{AF}" AS VARCHAR)),'.'),'') AS raw
                       FROM read_parquet({lit(meta['popmax']['path'])}) d SEMI JOIN keys k
                       ON regexp_replace(d."#chr",'^chr','')=k.c AND TRY_CAST(d."pos(1-based)" AS BIGINT)=k.p AND d.ref=k.r AND d.alt=k.a''')
        if count("SELECT count(*) FROM pop_raw WHERE raw IS NOT NULL AND (TRY_CAST(raw AS DOUBLE) IS NULL OR NOT isfinite(TRY_CAST(raw AS DOUBLE)) OR TRY_CAST(raw AS DOUBLE)<0 OR TRY_CAST(raw AS DOUBLE)>1)"):
            raise ValueError('Malformed nonmissing POPmax value')
        if count('SELECT count(*) FROM (SELECT c,p,r,a FROM pop_raw GROUP BY ALL HAVING count(DISTINCT TRY_CAST(raw AS DOUBLE))>1)'):
            raise ValueError('Conflicting POPmax values at exact allele')
        con.execute('CREATE TABLE pop AS SELECT c,p,r,a,max(TRY_CAST(raw AS DOUBLE)) AS af,count(*) AS source_rows FROM pop_raw GROUP BY c,p,r,a')
        for track in TRACKS:
            con.execute(f'CREATE TABLE {track}(start_ BIGINT,end_ BIGINT)')
            rows = read_regions(meta['regions'][track]['path'],track,chrom(meta['chromosome']))
            if rows: con.executemany(f'INSERT INTO {track} VALUES (?,?)',rows)
        receipt['counts'] = {}; receipt['outputs'] = {}
        for kind in ('missense','lof_hc'):
            overlaps = ','.join(f'EXISTS(SELECT 1 FROM {t} b WHERE CAST(v.POS AS BIGINT)>b.start_ AND CAST(v.POS AS BIGINT)<=b.end_) AS pcf_overlap_{t}' for t in TRACKS)
            union = ' OR '.join('pcf_overlap_'+t for t in TRACKS)
            con.execute(f'''CREATE TABLE {kind}_audit AS WITH joined AS (
                SELECT v.*,s.c IS NOT NULL AS pcf_site_matched,s.ac AS pcf_cohort_ac,s.an AS pcf_cohort_an,s.af AS pcf_cohort_af,
                       s.source_ac AS pcf_source_info_ac,s.source_an AS pcf_source_info_an,s.source_af AS pcf_source_info_af,
                       g.c IS NOT NULL AS pcf_gnomad_matched,g.af AS pcf_gnomad_popmax_af,{overlaps}
                FROM {kind} v LEFT JOIN sites s ON regexp_replace(v.CHROM,'^chr','')=s.c AND CAST(v.POS AS BIGINT)=s.p AND v.REF=s.r AND v.ALT=s.a
                LEFT JOIN pop g ON regexp_replace(v.CHROM,'^chr','')=g.c AND CAST(v.POS AS BIGINT)=g.p AND v.REF=g.r AND v.ALT=g.a),
                flags AS (SELECT *,COALESCE(pcf_cohort_af<0.005,FALSE) AS pcf_cohort_pass,
                          pcf_gnomad_popmax_af IS NULL OR pcf_gnomad_popmax_af<0.001 AS pcf_gnomad_pass,
                          NOT ({union}) AS pcf_region_pass FROM joined)
                SELECT *,pcf_cohort_pass AND pcf_gnomad_pass AND pcf_region_pass AS pcf_retained FROM flags''')
            def n(where='true'): return count(f'SELECT count(*) FROM {kind}_audit WHERE {where}')
            receipt['counts'][kind] = dict(input_candidates=n(),
                independent={name:dict(pass_count=n('pcf_'+name+'_pass'),fail_count=n('NOT pcf_'+name+'_pass')) for name in ['cohort','gnomad','region']},
                sequential=dict(after_cohort=n('pcf_cohort_pass'),after_gnomad=n('pcf_cohort_pass AND pcf_gnomad_pass'),final_retained=n('pcf_retained')),
                final_retained=n('pcf_retained'),missing_site_matches=n('NOT pcf_site_matched'),
                missing_cohort_af=n('pcf_cohort_af IS NULL'),matched_missing_or_invalid_cohort_af=n('pcf_site_matched AND pcf_cohort_af IS NULL'),
                missing_gnomad_popmax=n('pcf_gnomad_popmax_af IS NULL'),missing_gnomad_matches=n('NOT pcf_gnomad_matched'),
                overlaps={t:n('pcf_overlap_'+t) for t in TRACKS},overlap_union=n('NOT pcf_region_pass'),
                allele_classes={c:dict(input_candidates=n('allele_class='+lit(c)),retained=n('allele_class='+lit(c)+' AND pcf_retained')) for c in ['sequence','spanning_deletion']})
            name=kind+'.filtered.parquet'
            con.execute(f'COPY (SELECT * FROM {kind}_audit WHERE pcf_retained ORDER BY CHROM,POS,REF,ALT) TO {lit(out/name)} (FORMAT PARQUET,COMPRESSION ZSTD)')
        columns=[r[0] for r in con.execute('DESCRIBE missense_audit').fetchall() if r[0].startswith('pcf_')]
        projection=','.join(['CHROM','POS','REF','ALT','allele_class',*columns])
        con.execute(f"COPY (SELECT 'missense' AS candidate_type,{projection} FROM missense_audit UNION ALL SELECT 'lof_hc',{projection} FROM lof_hc_audit) TO {lit(out/'filter_audit.parquet')} (FORMAT PARQUET,COMPRESSION ZSTD)")
        for item in resources: unchanged(item)
        if hashes != {k:sha(v) for k,v in inputs.items()}: raise ValueError('Input changed during filtering')
        receipt.update(status='passed',input_sha256=hashes,script_sha256=sha(__file__),duckdb_version=duckdb.__version__)
        for name in output_names[:-1]:
            receipt['outputs'][name]=dict(sha256=sha(out/name),bytes=(out/name).stat().st_size,rows=count(f'SELECT count(*) FROM read_parquet({lit(out/name)})'))
    except Exception as exc:
        receipt['error_type']=type(exc).__name__
        for name in output_names[:-1]: (out/name).unlink(missing_ok=True)
        raise
    finally:
        if con: con.close()
        receipt['wall_seconds']=time.perf_counter()-start
        (out/'receipt.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'status':receipt['status'],'unit_id':meta['unit_id'],'counts':receipt['counts']}))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['missense','lof-hc','sites','metadata','outdir']: p.add_argument('--'+name,required=True)
    p.add_argument('--threads',type=int,default=2);p.add_argument('--memory',default='3GB')
    try: run(p.parse_args())
    except Exception: raise SystemExit('Pre-carrier filtering failed; inspect task receipt. No protected records printed.')
