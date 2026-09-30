import csv
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

import pytest

ROOT = Path(__file__).resolve().parents[1]


def module(name):
    spec = importlib.util.spec_from_file_location(name, ROOT/'scripts'/f'{name}.py')
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


res = module('annotation_resources')
stage = module('run_annotation_stage')


def resources(tmp_path):
    ann, lof = tmp_path/'annotation', tmp_path/'loftee'
    ann.mkdir(); lof.mkdir()
    names = [res.FASTA, res.FASTA+'.fai', 'vep115.transcript-priority.tsv', 'vep115.consequence-ranks.tsv']
    for chrom in ['chr21', 'chr22']:
        names += [f'Homo_sapiens.GRCh38.115.{chrom}.gff3', f'Homo_sapiens.GRCh38.115.{chrom}.gff3.fastvep.cache']
    for name in names:
        (ann/name).write_text('synthetic resource\n')
    for name in names:
        if name.endswith('.cache'):
            gff = ann/name.removesuffix('.fastvep.cache')
            ns = gff.stat().st_mtime_ns + 1000000000
            os.utime(ann/name, ns=(ns, ns))
    for name in ['ensembl-115/transcripts.sqlite', 'loftee-grch38/human_ancestor.fa.gz',
                 'loftee-grch38/human_ancestor.fa.gz.fai', 'loftee-grch38/human_ancestor.fa.gz.gzi',
                 'loftee-grch38/gerp_conservation_scores.homo_sapiens.GRCh38.bw', 'loftee-grch38/loftee.sql']:
        path = lof/name; path.parent.mkdir(exist_ok=True); path.write_text('synthetic resource\n')
    fast, lof_sif = tmp_path/'fast.sif', tmp_path/'lof.sif'
    fast.write_text('synthetic fast image'); lof_sif.write_text('synthetic lof image')
    return ann, lof, fast, lof_sif


def lock(tmp_path, monkeypatch):
    ann, lof, fast, lof_sif = resources(tmp_path)
    monkeypatch.setattr(res, 'FAST_SIF', res.sha(fast))
    monkeypatch.setattr(res, 'LOF_SIF', res.sha(lof_sif))
    return res.build(ann, lof, ['chr22', '21'], fast, lof_sif)


def test_lock_requires_exact_images_and_newer_cache(tmp_path, monkeypatch):
    ann, lof, fast, lof_sif = resources(tmp_path)
    with pytest.raises(ValueError, match='SIF checksum'):
        res.build(ann, lof, ['22'], fast, lof_sif)
    monkeypatch.setattr(res, 'FAST_SIF', res.sha(fast)); monkeypatch.setattr(res, 'LOF_SIF', res.sha(lof_sif))
    cache = ann/'Homo_sapiens.GRCh38.115.chr22.gff3.fastvep.cache'
    gff = ann/'Homo_sapiens.GRCh38.115.chr22.gff3'
    ns = gff.stat().st_mtime_ns
    os.utime(cache, ns=(ns, ns))
    with pytest.raises(ValueError, match='cache must be newer'):
        res.build(ann, lof, ['22'], fast, lof_sif)
    assert cache.stat().st_mtime_ns == ns  # no repair/mutation


def test_changed_resource_fails_before_annotation(tmp_path, monkeypatch):
    data = lock(tmp_path, monkeypatch)
    stage.check_resources(data['loftee'])
    Path(data['loftee']['transcripts']['path']).write_text('changed')
    with pytest.raises(ValueError, match='Resource changed'):
        stage.check_resources(data['loftee'])


def test_loftee_uses_validated_cli_and_preserves_failure_receipt(tmp_path, monkeypatch):
    data = lock(tmp_path, monkeypatch)
    meta = dict(unit_id='explicit_id', chromosome='chr22', source_vcf='/published/input.vcf.gz',
                loftee=data['loftee'], containers=data['containers'])
    monkeypatch.chdir(tmp_path)
    Path('unit.json').write_text(json.dumps(meta)); Path('picked.tsv').write_text('header\nrow\nrow2\n')
    Path('upstream.json').write_text(json.dumps(dict(status='passed', output_sha256=stage.sha('picked.tsv'))))
    from argparse import Namespace
    args = Namespace(metadata='unit.json', stage='loftee', input='picked.tsv', upstream='upstream.json')
    calls = []
    def fake_run(command, **kwargs):
        calls.append(command)
        Path('loftee.tsv').write_text('header\nlof\n')
    real_sha = stage.sha
    monkeypatch.setattr(stage, 'sha', lambda p: 'a'*64 if str(p).startswith('/opt/') else real_sha(p))
    monkeypatch.setattr(stage.subprocess, 'run', fake_run)
    stage.run(args)
    command = calls[0]
    assert command[:2] == [sys.executable, '/opt/rvp/scripts/run_standalone_loftee.py']
    for flag in ['transcripts', 'reference', 'ancestor', 'gerp', 'conservation']:
        assert command[command.index('--'+flag)+1] == data['loftee'][flag]['path']
    receipt = json.loads(Path('receipt.json').read_text())
    assert receipt['status'] == 'passed' and receipt['input_picked_rows'] == 2 and receipt['output_loftee_rows'] == 1
    def fail(*args, **kwargs):
        raise RuntimeError('PROTECTED_RECORD_MUST_NOT_APPEAR_IN_RECEIPT')
    monkeypatch.setattr(stage.subprocess, 'run', fail)
    with pytest.raises(RuntimeError):
        stage.run(args)
    assert not Path('loftee.tsv').exists()
    assert json.loads(Path('receipt.json').read_text())['status'] == 'failed'
    assert 'PROTECTED_RECORD' not in Path('receipt.json').read_text()


@pytest.mark.skipif(not shutil.which('nextflow'), reason='Nextflow required')
def test_nextflow_selection_resume_and_resource_invalidation(tmp_path, monkeypatch):
    data = lock(tmp_path, monkeypatch)
    repo = tmp_path/'repo'
    repo.mkdir()
    for directory in ['conf', 'lib', 'modules', 'subworkflows', 'scripts', 'docs/operations/abcd-fastvep-smoke']:
        shutil.copytree(ROOT/directory, repo/directory, dirs_exist_ok=True)
    for name in ['annotation.nf', 'annotation.config']:
        shutil.copy(ROOT/name, repo/name)
    # Test the real Nextflow graph/cache/publication, replacing only the unavailable scientific executables.
    (repo/'scripts/run_annotation_stage.py').write_text('''import argparse,json,pathlib
p=argparse.ArgumentParser()
for n in ['stage','metadata','input','benchmark','upstream']: p.add_argument('--'+n)
a=p.parse_args(); m=json.loads(pathlib.Path(a.metadata).read_text())
pathlib.Path('picked.tsv' if a.stage=='fastvep' else 'loftee.tsv').write_text('header\\nsynthetic\\n')
pathlib.Path('receipt.json').write_text(json.dumps(dict(status='passed',unit_id=m['unit_id'])))
''')
    resource_lock = tmp_path/'resources.json'; resource_lock.write_text(json.dumps(data))
    manifest = tmp_path/'sites.tsv'
    manifest.write_text('unit_id\tchromosome\tvcf\n')
    for unit, chrom in [('block12', 'chr22'), ('block19', '22'), ('other_chrom', 'chr21')]:
        source = tmp_path/(unit+'.vcf.gz'); source.write_text('synthetic staged input for graph test')
        with manifest.open('a') as handle: handle.write(f'{unit}\t{chrom}\t{source}\n')
    base = ['nextflow', '-C', str(repo/'annotation.config'), 'run', str(repo/'annotation.nf'),
            '-ansi-log', 'false', '--sites_manifest', str(manifest), '--resource_lock', str(resource_lock),
            '--annotation_root', data['annotation_root'], '--loftee_root', data['loftee_root'],
            '--fastvep_container', data['containers']['fastvep']['path'],
            '--loftee_container', data['containers']['loftee']['path'],
            '--outdir', str(tmp_path/'out'), '--annotation_memory', '1 GB']
    # Local graph-only test; do not invoke SIFs or a Slurm scheduler.
    env = dict(os.environ)
    bin_dir = tmp_path/'bin'; bin_dir.mkdir()
    (bin_dir/'python').symlink_to(sys.executable)
    env['PATH'] = str(bin_dir)+os.pathsep+env['PATH']
    def run(label, units, resume=True, success=True):
        trace = tmp_path/(label+'.trace.tsv')
        command = base + ['--select_units', units, '-with-trace', str(trace)]
        if resume: command += ['-resume']
        result = subprocess.run(command, cwd=tmp_path, env=env, capture_output=True, text=True, timeout=120)
        if not success:
            assert result.returncode != 0
            return result
        assert result.returncode == 0, result.stdout+result.stderr
        return list(csv.DictReader(trace.open(), delimiter='\t'))
    first = run('one', 'block12', False)
    assert len(first) == 2 and all(r['status'] == 'COMPLETED' for r in first)
    second = run('three', 'all')
    assert len(second) == 6
    assert all(r['status'] == ('CACHED' if '(block12)' in r['name'] else 'COMPLETED') for r in second)
    assert all(r['status'] == 'CACHED' for r in run('subset', 'block19'))
    for invalid in ['missing', 'block12,block12', 'block12,']:
        assert 'Unknown, empty, or duplicate' in str(run('invalid', invalid, success=False))
    changed = Path(data['chromosomes']['chr22']['priority']['path']); changed.write_text('new priority table')
    assert 'Resource changed/missing' in str(run('stale', 'block12', success=False))
    data['chromosomes']['chr22']['priority'] = res.identity(changed)
    resource_lock.write_text(json.dumps(data))
    assert all(r['status'] == 'COMPLETED' for r in run('new_resource', 'block12'))
    products = {p: p.read_bytes() for p in (tmp_path/'out').rglob('*.tsv')}
    assert len(products) == 6
    shutil.rmtree(tmp_path/'work')
    assert all(not p.is_symlink() and p.read_bytes() == content for p, content in products.items())


@pytest.mark.skipif(not shutil.which('bcftools'), reason='Native bcftools required')
def test_existing_benchmark_wrapper_with_synthetic_tools(tmp_path, monkeypatch):
    """Exercise actual adapter/benchmark and real bcftools; annotation binaries are explicit stubs."""
    data = lock(tmp_path, monkeypatch)
    ann = Path(data['annotation_root'])
    (ann/res.FASTA).write_text('>22\nAAAA\n')
    (ann/(res.FASTA+'.fai')).write_text('22\t4\t4\t4\t5\n')
    gff = ann/'Homo_sapiens.GRCh38.115.chr22.gff3'
    gff.write_text('##gff-version 3\n22\tEnsembl\tgene\t1\t4\t.\t+\t.\tID=gene:x\n')
    cache = Path(str(gff)+'.fastvep.cache')
    ns = gff.stat().st_mtime_ns + 1000000000; os.utime(cache, ns=(ns,ns))
    for key, item in data['chromosomes']['chr22'].items():
        data['chromosomes']['chr22'][key] = res.identity(item['path'])
    bins = tmp_path/'bin'; bins.mkdir()
    fast = bins/'fastvep'
    fast.write_text('#!'+sys.executable+'''\nimport pathlib,sys
assert sys.argv[1]=='annotate'
for flag in ['--hgvs','--symbol','--canonical']: assert flag in sys.argv
text=pathlib.Path(sys.argv[sys.argv.index('--input')+1]).read_text()
assert 'CSQ=' not in text and 'ID=CSQ,' not in text
assert '\\n22\\t' in text
sys.stdout.write(text)
''')
    picker = bins/'fastvep-picker'
    picker.write_text('#!'+sys.executable+'''\nimport pathlib,sys
assert sys.argv[sys.argv.index('--fastvep')+1]=='-'
rows=[line for line in sys.stdin if not line.startswith('#')]
pathlib.Path(sys.argv[sys.argv.index('--output')+1]).write_text('header\\n'+''.join(rows))
''')
    fast.chmod(0o755); picker.chmod(0o755)
    monkeypatch.setenv('PATH', str(bins)+os.pathsep+os.environ['PATH'])
    source = tmp_path/'sites.vcf'
    source.write_text('##fileformat=VCFv4.2\n##contig=<ID=chr22,length=4>\n'
                      '##INFO=<ID=CSQ,Number=.,Type=String,Description="VEP112">\n'
                      '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
                      'chr22\t1\t.\tA\tG\t.\tPASS\tCSQ=old\n'
                      'chr22\t2\t.\tA\tC\t.\tPASS\tCSQ=old\n')
    before = source.read_bytes()
    meta = dict(unit_id='block12', chromosome='chr22', source_vcf=str(source),
                annotation_root=str(ann), fastvep=data['chromosomes']['chr22'], containers=data['containers'])
    (tmp_path/'unit.json').write_text(json.dumps(meta))
    command = [sys.executable, str(ROOT/'scripts/run_annotation_stage.py'), '--stage', 'fastvep',
               '--metadata', 'unit.json', '--input', str(source),
               '--benchmark', str(ROOT/'docs/operations/abcd-fastvep-smoke/benchmark.py')]
    result = subprocess.run(command, cwd=tmp_path, capture_output=True, text=True)
    assert result.returncode == 0, result.stdout+result.stderr
    receipt = json.loads((tmp_path/'receipt.json').read_text())
    assert receipt['input_records'] == receipt['output_picked_rows'] == 2
    assert receipt['removed_input_csq'] and receipt['contig_renames'] == {'chr22': '22'}
    assert receipt['resources'] == data['chromosomes']['chr22']
    assert receipt['source_vcf'] == str(source) and receipt['unit_id'] == 'block12'
    assert source.read_bytes() == before
    # Nonzero FastVEP must fail even when the downstream picker exits successfully.
    fast.write_text('#!'+sys.executable+'\nimport sys\nsys.exit(7)\n'); fast.chmod(0o755)
    # Retry in the same work directory, including prior benchmark products.
    result = subprocess.run(command, cwd=tmp_path, capture_output=True, text=True)
    assert result.returncode != 0
    assert json.loads((tmp_path/'receipt.json').read_text())['status'] == 'failed'
    assert not (tmp_path/'picked.tsv').exists()


def test_pilot_checker_rejects_changed_bytes_even_with_matching_receipt(tmp_path, monkeypatch):
    monkeypatch.syspath_prepend(str(ROOT/'scripts'))
    checker = module('check_annotation_pilot')
    unit = 'chr22_block12'
    monkeypatch.setitem(checker.EXPECTED, unit, (2, 1))
    for stage_name, filename, content in [('fastvep-picker', 'picked', 'header\na\nb\n'),
                                          ('loftee', 'loftee', 'header\nlof\n')]:
        directory = tmp_path/stage_name/unit; directory.mkdir(parents=True)
        product = directory/(filename+'.tsv'); product.write_text(content)
        digest = checker.sha(product); monkeypatch.setitem(checker.HASHES, filename, digest)
        receipt = dict(status='passed', unit_id=unit, output_sha256=digest, input_records=2,
                       input_picked_rows=2, input_sha256=checker.HASHES['picked'])
        (directory/'receipt.json').write_text(json.dumps(receipt))
    assert checker.check(tmp_path, unit)['status'] == 'passed'
    product = tmp_path/'loftee'/unit/'loftee.tsv'; product.write_text('header\nchanged\n')
    receipt_file = product.with_name('receipt.json'); receipt = json.loads(receipt_file.read_text())
    receipt['output_sha256'] = checker.sha(product); receipt_file.write_text(json.dumps(receipt))
    with pytest.raises(ValueError, match='bytes differ from manual pilot'):
        checker.check(tmp_path, unit)
