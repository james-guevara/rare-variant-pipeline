#!/usr/bin/env python3
"""Hash the staged pre-carrier resources once, outside block tasks."""
import argparse
import json
from pathlib import Path
from annotation_resources import identity, LOF_SIF


def build(root, chromosomes, container):
    root=Path(root).resolve()
    def item(relative):
        result=identity(root/relative)
        if not Path(result['path']).is_relative_to(root): raise ValueError('Resource escapes read-only root')
        return result
    db={}
    for c in chromosomes:
        c=c if c.startswith('chr') else 'chr'+c
        if c not in ['chr'+str(i) for i in range(1,23)]+['chrX','chrY']:
            raise ValueError('Invalid chromosome')
        db[c]=item(f'dbNSFP/5.3.1a/parquet_scores_af/{c}.parquet')
    regions={t:item('problematic-regions/'+t+'.bed') for t in ['genomicSuperDups','simpleRepeat','rmsk']}
    image=identity(container)
    if image['sha256']!=LOF_SIF: raise ValueError('Use the existing validated LOFTEE SIF')
    return dict(schema=1,stage='pre_carrier',resource_root=str(root),popmax=db,regions=regions,container=image)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['resource-root','container','output']: p.add_argument('--'+name,required=True)
    p.add_argument('--chromosomes',default='chr22');a=p.parse_args()
    Path(a.output).write_text(json.dumps(build(a.resource_root,a.chromosomes.split(','),a.container),indent=2,sort_keys=True)+'\n')
    print('Staged resource hashes recorded; SIF matches validated image.')
