#!/usr/bin/env python3
"""Standalone paired extraction timing and exhaustive local product comparison."""
import argparse
from collections import Counter
import csv
import gzip
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import resource
import statistics
import subprocess
import sys
import time

ROOT=Path(__file__).resolve().parents[2]
PRODUCTS=['carriers.tsv.gz','candidates.tsv','unmatched.tsv','samples.tsv','sample_burden.tsv',
          'sample_gene_burden.tsv','gene_burden.tsv','sample_distinct_alleles.tsv']
SCIENCE=['samples_in_source','by_candidate_type','candidate_annotation_records','distinct_candidate_alleles',
         'overlapping_type_alleles','carrier_annotation_records','distinct_allele_audit_by_class','contig_alias']


def identity(p,hash_content=True):
    p=Path(p);s=p.stat();r=dict(path=str(p.resolve()),bytes=s.st_size,mtime_ns=s.st_mtime_ns)
    if hash_content:
        h=hashlib.sha256()
        with p.open('rb') as f:
            for b in iter(lambda:f.read(8388608),b''):h.update(b)
        r['sha256']=h.hexdigest()
    return r


def read_table(p):
    with (gzip.open(p,'rt',newline='') if str(p).endswith('.gz') else p.open(newline='')) as f:
        rows=list(csv.reader(f,delimiter='\t'))
    if not rows or len(rows[0])!=len(set(rows[0])) or any(len(r)!=len(rows[0]) for r in rows):
        raise ValueError('Malformed benchmark output table')
    return rows[0], [tuple(r) for r in rows[1:]]


def compare(left,right):
    report={};passed=True
    for name in PRODUCTS:
        lh,lr=read_table(left/name);rh,rr=read_table(right/name)
        lc,rc=Counter(lr),Counter(rr)
        # Full values and duplicate multiplicities, not only identity keys or totals.
        equal=lh==rh and lc==rc
        order=lr==rr
        report[name]=dict(left_rows=len(lr),right_rows=len(rr),schema_equal=lh==rh,
                          full_records_equal=equal,row_order_equal=order,
                          left_only_rows=sum((lc-rc).values()),right_only_rows=sum((rc-lc).values()))
        passed &= equal and order  # includes complete source sample order
    l=json.loads((left/'receipt.json').read_text());r=json.loads((right/'receipt.json').read_text())
    matched={k:l.get(k)==r.get(k) for k in SCIENCE}
    passed &= l['status']==r['status']=='passed' and all(matched.values())
    return dict(passed=bool(passed),products=report,receipt_fields_equal=matched)


def run(a):
    if a.pairs<2:raise ValueError('At least two reversed-order pairs required')
    out=Path(a.outdir).resolve();out.mkdir(parents=True,exist_ok=False)
    os.chmod(out,0o700)
    fields=['metadata','missense','lof_hc','vcf','index']
    paths={k:Path(getattr(a,k)).resolve() for k in fields}
    ids={k:identity(p,k!='vcf') for k,p in paths.items()}
    container=json.loads(Path(a.container_receipt).read_text()) if a.container_receipt else None
    report=dict(status='failed',versions={p:importlib.metadata.version(p) for p in ['pysam','cyvcf2','numpy','duckdb']},
                python=sys.version,inputs=ids,container=container,
                measurement='parent perf_counter and child user/system CPU; includes interpreter startup, input checks, extraction, compression and summary writes; excludes scheduler/container startup and comparison',
                cache_policy='No cache flush or cold-cache claim; consecutive pairs reverse order. First-run I/O and shared filesystem contention may affect results.',
                comparison='Full field values and duplicate multiplicity in every product; exact row/sample order; scientific receipt fields. Individual rows stay local.',
                runs=[],comparisons=[],code_identities={str(p.relative_to(ROOT)):identity(p) for p in [ROOT/'scripts/extract_filtered_carriers.py',ROOT/'scripts/extract_exact_carriers.py',Path(__file__),Path(__file__).with_name('run_backend.py'),Path(__file__).with_name('cyvcf2_backend.py')]})
    env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',NUMEXPR_NUM_THREADS='1')
    try:
        baseline=None
        for pair in range(a.pairs):
            order=['pysam','cyvcf2'] if pair%2==0 else ['cyvcf2','pysam']
            for backend in order:
                number=len(report['runs'])+1;directory=out/f'{number:02d}-{backend}'
                cmd=[sys.executable,str(Path(__file__).with_name('run_backend.py')),'--backend',backend,'--outdir',str(directory)]
                for k,p in paths.items():cmd+=['--'+k.replace('_','-'),str(p)]
                for k in ['expected_missense','expected_hc']:
                    if getattr(a,k) is not None:cmd+=['--'+k.replace('_','-'),str(getattr(a,k))]
                before=resource.getrusage(resource.RUSAGE_CHILDREN);start=time.perf_counter()
                with (out/f'{number:02d}-private-process.log').open('w') as log:
                    p=subprocess.run(cmd,stdout=log,stderr=log,env=env)
                wall=time.perf_counter()-start;after=resource.getrusage(resource.RUSAGE_CHILDREN)
                timing=dict(run=number,pair=pair+1,backend=backend,wall_seconds=wall,
                            cpu_user_seconds=after.ru_utime-before.ru_utime,cpu_system_seconds=after.ru_stime-before.ru_stime,exit_code=p.returncode)
                timing['cpu_seconds']=timing['cpu_user_seconds']+timing['cpu_system_seconds'];report['runs'].append(timing)
                if p.returncode:raise ValueError('Backend failed; inspect private local logs')
                if baseline is None:baseline=directory
                comparison=compare(baseline,directory);report['comparisons'].append(dict(run=number,**comparison))
                if not comparison['passed']:raise ValueError('Record/summary parity failed; outputs retained locally')
        if ids!={k:identity(p,k!='vcf') for k,p in paths.items()}:raise ValueError('Input identity changed during benchmark')
        report['timing_summary']={b:{m:dict(median=statistics.median(r[m] for r in report['runs'] if r['backend']==b),
                                                     minimum=min(r[m] for r in report['runs'] if r['backend']==b),
                                                     maximum=max(r[m] for r in report['runs'] if r['backend']==b)) for m in ['wall_seconds','cpu_seconds']} for b in ['pysam','cyvcf2']}
        receipt=json.loads((baseline/'receipt.json').read_text());report['aggregate_baseline']={k:receipt[k] for k in SCIENCE}
        report['status']='passed'
    except Exception as exc:
        report['error']=str(exc) if isinstance(exc,ValueError) else 'Benchmark failure; inspect local artifacts'
        raise
    finally:
        (out/'benchmark.json').write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
    # Only aggregate status/timing on stdout. No sample IDs, allele records or mismatches.
    print(json.dumps({k:report[k] for k in ['status','versions','runs','timing_summary','aggregate_baseline']},indent=2))


if __name__=='__main__':
    os.umask(0o077)
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['metadata','missense','lof-hc','vcf','index','outdir']:p.add_argument('--'+name,required=True)
    p.add_argument('--container-receipt')
    p.add_argument('--pairs',type=int,default=2,help='2 = pysam, cyvcf2, cyvcf2, pysam; use 4 for eight runs')
    p.add_argument('--expected-missense',type=int)
    p.add_argument('--expected-hc',type=int)
    try:run(p.parse_args())
    except Exception:sys.exit('Benchmark failed. Inspect benchmark.json/private artifacts locally; do not share individual records.')
