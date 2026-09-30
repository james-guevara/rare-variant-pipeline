#!/usr/bin/env python3
"""Compare Nextflow products with user-reported NBDC pilot counts/hashes; print no records."""
import argparse
import json
from pathlib import Path
from annotation_resources import sha

EXPECTED = {
    'chr22_block0': (713983, 409),
    'chr22_block12': (48001, 332),
    'chr22_block19': (372375, 430),
}
HASHES = {
    'picked': 'cbf5c6bb2f469bbb5e6e7bee1f357f325fca6f24b6246f8ea86d262bd4ff6208',
    'loftee': 'b9fee4e3527371c8ff11b6caaf25de99088f215d20fc063463be737c99a7ca82',
}


def count(path):
    with open(path) as handle:
        if not handle.readline().strip():
            raise ValueError('Missing TSV header')
        return sum(1 for line in handle if line.strip())


def check(root, unit):
    picked_count, lof_count = EXPECTED[unit]
    results = {}
    for stage, filename, expected in [('fastvep-picker', 'picked', picked_count), ('loftee', 'loftee', lof_count)]:
        directory = Path(root) / stage / unit
        product = directory / (filename + '.tsv')
        receipt = json.loads((directory / 'receipt.json').read_text())
        digest = sha(product)
        if receipt['status'] != 'passed' or receipt['unit_id'] != unit:
            raise ValueError(f'{unit}/{stage}: receipt failed or wrong unit')
        if count(product) != expected or digest != receipt['output_sha256']:
            raise ValueError(f'{unit}/{stage}: count/hash differs from receipt or pilot')
        if stage == 'fastvep-picker' and receipt['input_records'] != picked_count:
            raise ValueError(f'{unit}: input count differs from pilot')
        if stage == 'loftee' and (receipt['input_picked_rows'] != picked_count or
                                  receipt['input_sha256'] != results['picked']['sha256']):
            raise ValueError(f'{unit}: LOFTEE input differs from picked product')
        if unit == 'chr22_block12' and digest != HASHES[filename]:
            raise ValueError(f'{unit}/{stage}: bytes differ from manual pilot')
        results[filename] = dict(rows=expected, sha256=digest)
    return dict(unit_id=unit, status='passed', outputs=results)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--outdir', required=True)
    p.add_argument('--units', default='chr22_block12')
    a = p.parse_args()
    for unit in a.units.split(','):
        print(json.dumps(check(a.outdir, unit), sort_keys=True))
