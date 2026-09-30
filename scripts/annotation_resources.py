#!/usr/bin/env python3
"""Create a checksum lock for immutable, shared annotation resources (no copying)."""
import argparse
import hashlib
import json
from pathlib import Path
import re

FAST_SIF = '6eb015d6cb41ae10c64d373b508963496f14a2f57f14d3294cfcd5e06d922d8f'
LOF_SIF = 'a1722375533ad25de58d01635d3cf36667b4ca909149c59404a8deb6491907cd'
FASTA = 'Homo_sapiens.GRCh38.dna.primary_assembly.fa'


def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def identity(path):
    path = Path(path).resolve()
    before = path.stat()
    if not path.is_file() or not before.st_size:
        raise ValueError(f'Missing/empty resource: {path}')
    digest = sha(path)
    after = path.stat()
    if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
        raise ValueError(f'Resource changed while hashing: {path}')
    return dict(path=str(path), bytes=after.st_size, mtime_ns=after.st_mtime_ns,
                mtime_ms=after.st_mtime_ns // 1000000, sha256=digest)


def build(annotation, loftee, chromosomes, fast_sif, lof_sif):
    annotation, loftee = Path(annotation).resolve(), Path(loftee).resolve()
    fast_shared = {k: identity(annotation / v) for k, v in {
        'fasta': FASTA, 'fai': FASTA + '.fai',
        'priority': 'vep115.transcript-priority.tsv', 'ranks': 'vep115.consequence-ranks.tsv'}.items()}
    lof = {k: identity(loftee / v) for k, v in {
        'transcripts': 'ensembl-115/transcripts.sqlite',
        'ancestor': 'loftee-grch38/human_ancestor.fa.gz',
        'ancestor_fai': 'loftee-grch38/human_ancestor.fa.gz.fai',
        'ancestor_gzi': 'loftee-grch38/human_ancestor.fa.gz.gzi',
        'gerp': 'loftee-grch38/gerp_conservation_scores.homo_sapiens.GRCh38.bw',
        'conservation': 'loftee-grch38/loftee.sql'}.items()}
    lof.update(reference=fast_shared['fasta'], reference_fai=fast_shared['fai'])
    chroms = {}
    for chrom in chromosomes:
        chrom = chrom if chrom.startswith('chr') else 'chr' + chrom
        if not re.fullmatch(r'chr(?:[1-9]|1[0-9]|2[0-2]|X|Y|M|MT)', chrom):
            raise ValueError('Invalid chromosome: ' + chrom)
        gff = annotation / f'Homo_sapiens.GRCh38.115.{chrom}.gff3'
        files = dict(fast_shared, gff3=identity(gff), cache=identity(str(gff) + '.fastvep.cache'))
        if files['cache']['mtime_ns'] <= files['gff3']['mtime_ns']:
            raise ValueError(f'{chrom}: cache must be newer than GFF3; resources were not modified')
        chroms[chrom] = files
    for files in [*chroms.values(), lof]:
        for item in files.values():
            path = Path(item['path'])
            if not (path.is_relative_to(annotation) or path.is_relative_to(loftee)):
                raise ValueError(f'Resource symlink escapes the read-only roots: {path}')
    containers = {'fastvep': identity(fast_sif), 'loftee': identity(lof_sif)}
    for name, expected in [('fastvep', FAST_SIF), ('loftee', LOF_SIF)]:
        if containers[name]['sha256'] != expected:
            raise ValueError(f'{name} SIF checksum differs from validated image')
    return dict(schema=1, annotation_root=str(annotation), loftee_root=str(loftee),
                chromosomes=chroms, loftee=lof, containers=containers)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for arg in ['annotation', 'loftee', 'fastvep-sif', 'loftee-sif', 'output']:
        p.add_argument('--' + arg, required=True)
    p.add_argument('--chromosomes', required=True, help='Comma-separated resource chromosomes, e.g. chr22,chr21')
    a = p.parse_args()
    data = build(a.annotation, a.loftee, a.chromosomes.split(','), a.fastvep_sif, a.loftee_sif)
    Path(a.output).write_text(json.dumps(data, indent=2, sort_keys=True) + '\n')
    print(f'Wrote resource lock for {len(data["chromosomes"])} chromosomes; SIF hashes verified.')


if __name__ == '__main__':
    main()
