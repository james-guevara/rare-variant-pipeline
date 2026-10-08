#!/usr/bin/env python3
"""Execute one isolated backend; timed by the parent driver."""
import argparse
import contextlib
import importlib.metadata
import json
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
import extract_filtered_carriers as production


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--backend',choices=['pysam','cyvcf2'],required=True)
    for name in ['metadata','missense','lof-hc','vcf','index','outdir']:p.add_argument('--'+name,required=True)
    p.add_argument('--expected-missense',type=int)
    p.add_argument('--expected-hc',type=int)
    a=p.parse_args()
    if a.backend=='cyvcf2':
        import cyvcf2_backend
        cyvcf2_backend.install(production)
    out=Path(a.outdir);out.mkdir(parents=True,exist_ok=False)
    with open(out/'private-backend.log','w') as log,contextlib.redirect_stdout(log):
        production.extract(a)
    r=json.loads((out/'receipt.json').read_text())
    r['benchmark_backend']=a.backend
    r['benchmark_versions']={x:importlib.metadata.version(x) for x in ['pysam','cyvcf2','numpy','duckdb']}
    (out/'receipt.json').write_text(json.dumps(r,indent=2,sort_keys=True)+'\n')


if __name__=='__main__':
    try:main()
    except Exception:sys.exit('Backend failed; inspect local receipt/private logs. No individual records printed.')
