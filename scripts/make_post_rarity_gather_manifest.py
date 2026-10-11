#!/usr/bin/env python3
"""Inventory saved passed post-rarity QC products, without invoking upstream work."""
import argparse
import csv
import json
from pathlib import Path
import re
import sys


def build(root):
    rows=[];seen=set();psam=None
    for p in sorted(Path(root).glob('*/receipt.json')):
        r=json.loads(p.read_text());unit=r['unit_id'];chromosome=r['chromosome']
        if r.get('status')!='passed' or r.get('stage')!='post_rarity_qc' or not r.get('reconciliation',{}).get('passed'):
            raise ValueError('Passed post-rarity QC receipt required')
        if not re.fullmatch('[A-Za-z0-9][A-Za-z0-9_.-]*',unit) or unit in seen:raise ValueError('Invalid/duplicate unit')
        if chromosome not in ['chr'+str(i) for i in range(1,23)]+['chrX','chrY']:raise ValueError('Invalid chromosome')
        identity=r.get('psam_identity')
        if not identity or r.get('input_identities',{}).get('psam')!=identity or (psam is not None and psam!=identity):raise ValueError('PSAM identities disagree')
        psam=identity
        files=[p.parent/'carriers.qc.tsv.gz',p.parent/'samples.tsv',p]
        if not all(f.is_file() for f in files):raise ValueError('Missing post-rarity QC products')
        rows.append([unit,chromosome,*[str(f.resolve()) for f in files]]);seen.add(unit)
    if not rows:raise ValueError('No post-rarity QC receipts')
    return rows


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',required=True);p.add_argument('--output',required=True);a=p.parse_args()
    try:
        rows=build(a.root)
        with open(a.output,'x') as f:
            w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(['unit_id','chromosome','carriers','samples','source_receipt']);w.writerows(rows)
        print(json.dumps(dict(units=len(rows),autosomal=sum(r[1] not in ('chrX','chrY') for r in rows),sex_chromosome=sum(r[1] in ('chrX','chrY') for r in rows))))
    except Exception:sys.exit('Post-rarity QC manifest failed; inspect saved receipts locally. Existing manifests are not overwritten.')
