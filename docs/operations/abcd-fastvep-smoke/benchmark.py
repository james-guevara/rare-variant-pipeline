#!/usr/bin/env python3
"""Operational NBDC command wrapper; invokes existing annotation binaries only."""
import argparse
import gzip
import hashlib
import json
import os
from pathlib import Path
import resource
import shutil
import subprocess
import tempfile
import time

p = argparse.ArgumentParser()
p.add_argument('--input', required=True)
p.add_argument('--resources', required=True)
p.add_argument('--chromosome', default='chr22', help='Resource filename component')
p.add_argument('--outdir', required=True)
p.add_argument('--sif-sha256', required=True)
a = p.parse_args()
root, source, out = Path(a.resources), Path(a.input), Path(a.outdir)
out.mkdir(parents=True, exist_ok=True)
picked, receipt = out / 'picked.tsv', out / 'receipt.json'
if picked.exists() or receipt.exists():
    raise SystemExit('Use a fresh output directory for each measurement.')
names = {
    'gff3': f'Homo_sapiens.GRCh38.115.{a.chromosome}.gff3',
    'cache': f'Homo_sapiens.GRCh38.115.{a.chromosome}.gff3.fastvep.cache',
    'fasta': 'Homo_sapiens.GRCh38.dna.primary_assembly.fa',
    'fai': 'Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai',
    'priority': 'vep115.transcript-priority.tsv',
    'ranks': 'vep115.consequence-ranks.tsv',
}
files = {k: root / v for k, v in names.items()}
for path in [source, *files.values()]:
    if not path.is_file() or path.stat().st_size == 0:
        raise SystemExit(f'Missing/empty required file: {path.name}')
bins = {name: shutil.which(name) for name in ['bcftools', 'fastvep', 'fastvep-picker']}
if not all(bins.values()):
    raise SystemExit('Container must provide bcftools, fastvep, and fastvep-picker.')

def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(8 * 1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()

def inspect(path):
    opener = gzip.open if str(path).endswith('.gz') else open
    count, contigs, csq, header = 0, set(), False, False
    with opener(path, 'rt') as f:
        for line in f:
            if line.startswith('##INFO=<ID=CSQ,'):
                csq = True
            if line.startswith('#CHROM\t'):
                header = True
                if len(line.rstrip('\n').split('\t')) != 8:
                    raise RuntimeError('Input must have zero samples and no FORMAT column.')
            if line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) != 8 or ',' in fields[4]:
                raise RuntimeError('This smoke test requires sites-only biallelic records.')
            count += 1
            contigs.add(fields[0])
    if not header or not count:
        raise RuntimeError('Missing VCF header or empty input.')
    return count, contigs, csq

log = out / 'private-failure.log'
try:
    prep_start = time.perf_counter()
    input_count, contigs, has_csq = inspect(source)
    with open(files['fai']) as f:
        reference_contigs = {line.split('\t')[0] for line in f}
    with open(files['gff3']) as f:
        annotation_contigs = {line.split('\t')[0] for line in f if not line.startswith('#') and line.strip()}
    valid = reference_contigs & annotation_contigs
    rename = {}
    for contig in contigs:
        if contig in valid:
            continue
        alias = contig[3:] if contig.startswith('chr') else 'chr' + contig
        if alias not in valid:
            raise RuntimeError('Input contig has no matching GFF3/reference contig.')
        rename[contig] = alias
    if len({rename.get(c, c) for c in contigs}) != len(contigs):
        raise RuntimeError('Ambiguous contig aliases in input.')
    with tempfile.TemporaryDirectory(prefix='.smoke-', dir=out) as temp, open(log, 'w') as errors:
        temp = Path(temp)
        clean = temp / 'input.vcf'
        cmd = [bins['bcftools'], 'annotate']
        if has_csq:
            cmd += ['-x', 'INFO/CSQ']
        if rename:
            mapping = temp / 'contigs.tsv'
            mapping.write_text(''.join(f'{k}\t{v}\n' for k, v in sorted(rename.items())))
            cmd += ['--rename-chrs', str(mapping)]
        if not has_csq and not rename:
            cmd = [bins['bcftools'], 'view']
        subprocess.run(cmd + ['-Ov', '-o', str(clean), str(source)], check=True, stdout=errors, stderr=errors)
        clean_count, _, remaining_csq = inspect(clean)
        if clean_count != input_count or remaining_csq:
            raise RuntimeError('Temporary input count/CSQ validation failed.')
        prep_seconds = time.perf_counter() - prep_start
        fast = [bins['fastvep'], 'annotate', '--input', str(clean), '--gff3', str(files['gff3']),
                '--fasta', str(files['fasta']), '--transcript-cache', str(files['cache']),
                '--hgvs', '--symbol', '--canonical', '--output-format', 'vcf', '--output', '-']
        pick = [bins['fastvep-picker'], '--fastvep', '-', '--transcript-priority', str(files['priority']),
                '--consequence-ranks', str(files['ranks']), '--output', str(picked)]
        # A fresh supervisor isolates child CPU/RSS from bcftools input preparation.
        supervisor = '''import json, resource, subprocess, sys, time
fast, pick, stats = json.loads(sys.argv[1])
t = time.perf_counter()
first = subprocess.Popen(fast, stdout=subprocess.PIPE)
try:
    second = subprocess.Popen(pick, stdin=first.stdout, stdout=subprocess.DEVNULL)
finally:
    first.stdout.close()
b = second.wait()
a = first.wait()
wall = time.perf_counter() - t
r = resource.getrusage(resource.RUSAGE_CHILDREN)
if a or b:
    sys.exit(1)
with open(stats, 'w') as f:
    json.dump(dict(wall_seconds=wall, cpu_seconds=r.ru_utime+r.ru_stime,
                   peak_child_rss_kib=r.ru_maxrss), f)
'''
        import sys
        stats = temp / 'timing.json'
        subprocess.run([sys.executable, '-c', supervisor, json.dumps([fast, pick, str(stats)])],
                       check=True, stdout=errors, stderr=errors)
        timing = json.loads(stats.read_text())
    with open(picked) as f:
        rows = sum(1 for _ in f) - 1
    if rows != input_count:
        raise RuntimeError('Picked row count differs from input record count.')
    # Hash resources after annotation so the timed run does not warm them by hashing.
    data = dict(status='passed', input_records=input_count, output_picked_rows=rows,
                output_bytes=picked.stat().st_size, input_preparation_seconds=prep_seconds,
                input_sha256=sha(source), output_sha256=sha(picked), source_vcf=str(source),
                removed_input_csq=has_csq, contig_renames=rename, timing=timing,
                rss_scope='largest individual child RSS; not summed concurrent pipeline RSS',
                documented_fastvep_commit='cb8113d7bab2db42cb06bb2b2a40c57b60ea2561',
                documented_picker_version='0.1.0',
                documented_picker_source_commit='1671c5a76c369da50c64320d7dc3c719ac2ab95a',
                documented_oci_digest='sha256:7d5b76a28e2427ca97af6ebec1d5e38aec63b93609419c0cc77c063e48e2917d',
                supplied_sif_sha256=a.sif_sha256,
                identity_note='Documented revisions require a SIF built from the pinned OCI image; hashes below identify actual files.',
                binaries={k: dict(path=v, sha256=sha(v)) for k, v in bins.items()},
                resources={k: dict(path=str(v), bytes=v.stat().st_size, sha256=sha(v)) for k, v in files.items()})
    receipt.write_text(json.dumps(data, indent=2) + '\n')
    log.unlink()
    print(json.dumps(dict(status='passed', input_records=input_count, output_picked_rows=rows,
                          output_bytes=data['output_bytes'], **timing), indent=2))
except Exception:
    picked.unlink(missing_ok=True)
    receipt.unlink(missing_ok=True)
    raise SystemExit('Smoke test failed. Inspect private-failure.log locally; do not share variant records.')
