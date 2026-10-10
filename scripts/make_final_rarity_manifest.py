#!/usr/bin/env python3
"""Inventory saved passed extraction/frequency receipts; never launch upstream work."""
import argparse
import csv
import json
from pathlib import Path
import re
import sys


def build(roots):
    rows=[];seen=set()
    for root in roots:
        paths=sorted(Path(root).glob('*/receipt.json'))
        if not paths:raise ValueError('No extraction receipts in input root')
        for p in paths:
            r=json.loads(p.read_text());unit=r['unit_id'];c=r['chromosome']
            if not re.fullmatch('[A-Za-z0-9][A-Za-z0-9_.-]*',unit) or unit in seen:raise ValueError('Duplicate or invalid unit identifier')
            if r.get('status')!='passed' or r.get('frequencies',{}).get('status')!='passed':raise ValueError('Passed extraction/frequency receipt required')
            if c not in ['chr'+str(i) for i in range(1,23)]+['chrX','chrY']:raise ValueError('Unsupported chromosome')
            files=[p.parent/n for n in ['carriers.tsv.gz','samples.tsv','receipt.json','candidates.tsv','variant_frequencies.tsv']]
            if not all(f.is_file() for f in files):raise ValueError('Missing saved extraction product')
            rows.append([unit,c,*[str(f.resolve()) for f in files]]);seen.add(unit)
    return rows


if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--root',action='append',required=True);a.add_argument('--output',required=True);args=a.parse_args()
    try:
        rows=build(args.root)
        with open(args.output,'x') as f:
            w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(['unit_id','chromosome','carriers','samples','source_receipt','candidates','frequencies']);w.writerows(rows)
        print(json.dumps(dict(units=len(rows),autosomal=sum(r[1] not in ('chrX','chrY') for r in rows),sex_chromosome=sum(r[1] in ('chrX','chrY') for r in rows))))
    except Exception:sys.exit('Manifest generation failed; inspect passed receipts/input roots locally. Existing manifests are not overwritten.')
