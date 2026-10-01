#!/usr/bin/env python3
"""Verify candidate-stage resources against the canonical v1 checksum inventory."""
import argparse
import json
from pathlib import Path
from annotation_resources import identity, LOF_SIF


def build(root, chromosomes, annotation_lock, canonical):
    root = Path(root).resolve()
    expected = {item['path']: item for item in json.loads(Path(canonical).read_text())['files']}
    def verified(relative, local):
        item = identity(root/local)
        if not Path(item['path']).is_relative_to(root):
            raise ValueError('Resource must be physically contained under read-only root')
        reference = expected[relative]
        if (item['bytes'], item['sha256']) != (reference['bytes'], reference['sha256']):
            raise ValueError('Resource differs from canonical v1: '+relative)
        return item
    db = {}
    for chromosome in chromosomes:
        chrom = chromosome if chromosome.startswith('chr') else 'chr'+chromosome
        relative = f'dbNSFP/5.3.1a/parquet_expanded_mane_select/{chrom}.parquet'
        db[chrom] = verified(relative, relative)
    gene = verified('targeted-annotation/GeneBayes.Supplementary_Table_1.tsv',
                    'GeneBayes/GeneBayes.Supplementary_Table_1.tsv')
    old = json.loads(Path(annotation_lock).read_text())['containers']['loftee']
    image = identity(old['path'])
    if image['sha256'] != LOF_SIF or image['sha256'] != old['sha256']:
        raise ValueError('Container differs from validated SIF')
    return dict(schema=1, resource_root=str(root), dbnsfp=db, genebayes=gene, container=image)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--resource-root', required=True)
    p.add_argument('--chromosomes', default='chr22')
    p.add_argument('--annotation-lock', required=True)
    p.add_argument('--output', required=True)
    p.add_argument('--canonical-inventory', default=str(Path(__file__).resolve().parents[1]/'docs/resources/repair-20260929/post-validation-hashes.json'))
    a = p.parse_args()
    result = build(a.resource_root, a.chromosomes.split(','), a.annotation_lock, a.canonical_inventory)
    Path(a.output).write_text(json.dumps(result, indent=2, sort_keys=True)+'\n')
    print('Candidate resources and SIF match canonical SHA-256 values.')
