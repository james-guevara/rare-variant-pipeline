#!/usr/bin/env python3
"""Orchestration/receipts only: reuse the validated benchmark and container LOFTEE CLI."""
import argparse
import hashlib
import json
from pathlib import Path
import resource
import subprocess
import sys
import time


def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def rows(path):
    with open(path) as handle:
        if not handle.readline().strip():
            raise ValueError('Missing TSV header')
        return sum(1 for line in handle if line.strip())


def check_resources(resources):
    for item in resources.values():
        stat = Path(item['path']).stat()
        if (stat.st_size, stat.st_mtime_ns) != (item['bytes'], item['mtime_ns']):
            raise ValueError('Resource changed since checksum lock creation')


def run(a):
    meta = json.loads(Path(a.metadata).read_text())
    identities = meta[a.stage]
    receipt = dict(status='failed', stage=a.stage, unit_id=meta['unit_id'], chromosome=meta['chromosome'],
                   source_vcf=meta['source_vcf'], resources=identities,
                   container=meta['containers'][a.stage],
                   resource_verification='SHA-256 at lock creation; size/mtime checked before and after execution; immutable read-only resources required')
    start = time.perf_counter()
    output = Path('picked.tsv' if a.stage == 'fastvep' else 'loftee.tsv')
    try:
        check_resources(identities)
        if a.stage == 'fastvep':
            # A failed/resumed task can reuse its work directory. The standalone
            # benchmark intentionally refuses preexisting success products.
            for name in ['picked.tsv', 'receipt.json']:
                Path('benchmark-output', name).unlink(missing_ok=True)
            Path('fastvep-resources.json').write_text(json.dumps(identities))
            command = [sys.executable, a.benchmark, '--input', a.input,
                       '--resources', meta['annotation_root'], '--chromosome', meta['chromosome'],
                       '--outdir', 'benchmark-output', '--sif-sha256', meta['containers']['fastvep']['sha256'],
                       '--resource-identities', 'fastvep-resources.json']
            with open('private-failure.log', 'w') as log:
                subprocess.run(command, check=True, stdout=log, stderr=log)
            data = json.loads(Path('benchmark-output/receipt.json').read_text())
            receipt.update(data)
            # Restore original published source path instead of the task's staged alias.
            receipt['source_vcf'] = meta['source_vcf']
            Path('benchmark-output/picked.tsv').replace(output)
        else:
            upstream = json.loads(Path(a.upstream).read_text())
            if upstream['status'] != 'passed' or sha(a.input) != upstream['output_sha256']:
                raise ValueError('Picked output does not match upstream receipt')
            script = '/opt/rvp/scripts/run_standalone_loftee.py'
            command = [sys.executable, script, '--input', a.input]
            for flag in ['transcripts', 'reference', 'ancestor', 'gerp', 'conservation']:
                command += ['--' + flag, identities[flag]['path']]
            command += ['--output', str(output)]
            before = resource.getrusage(resource.RUSAGE_CHILDREN)
            timing = time.perf_counter()
            with open('private-failure.log', 'w') as log:
                subprocess.run(command, check=True, stdout=log, stderr=log)
            wall = time.perf_counter() - timing
            usage = resource.getrusage(resource.RUSAGE_CHILDREN)
            receipt.update(input_picked_rows=rows(a.input), output_loftee_rows=rows(output),
                           input_sha256=sha(a.input), upstream_receipt_sha256=sha(a.upstream),
                           output_sha256=sha(output), output_bytes=output.stat().st_size,
                           implementation=dict(path=script, sha256=sha(script)),
                           timing=dict(wall_seconds=wall,
                                       cpu_seconds=usage.ru_utime+usage.ru_stime-before.ru_utime-before.ru_stime,
                                       peak_child_rss_kib=usage.ru_maxrss))
        check_resources(identities)
        receipt['status'] = 'passed'
    except Exception as exc:
        output.unlink(missing_ok=True)
        receipt['status'] = 'failed'
        receipt['error_type'] = type(exc).__name__
        # Never copy exception contents/variant-bearing stderr to public aggregate output.
        raise
    finally:
        receipt['stage_wall_seconds'] = time.perf_counter() - start
        Path('receipt.json').write_text(json.dumps(receipt, indent=2, sort_keys=True) + '\n')
    Path('private-failure.log').unlink(missing_ok=True)
    print(json.dumps({key: receipt[key] for key in ['status', 'stage', 'unit_id', 'stage_wall_seconds']}))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--stage', choices=['fastvep', 'loftee'], required=True)
    p.add_argument('--metadata', required=True)
    p.add_argument('--input', required=True)
    p.add_argument('--benchmark')
    p.add_argument('--upstream')
    a = p.parse_args()
    try:
        run(a)
    except Exception:
        raise SystemExit('Annotation failed; receipt.json and private-failure.log remain in the task work directory.')


if __name__ == '__main__':
    main()
