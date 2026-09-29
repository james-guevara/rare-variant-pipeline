import csv
import gzip
import json
from pathlib import Path
import shutil
import subprocess

import pytest

ROOT = Path(__file__).resolve().parents[1]
HAS_BCFTOOLS = shutil.which('bcftools') is not None
HAS_NEXTFLOW = shutil.which('nextflow') is not None


def fixture_vcf(path, csq=True, empty=False, chromosome='chr22', samples=3):
    header = ('##fileformat=VCFv4.2\n##contig=<ID=chr22,length=100000>\n'
              '##VEP="v112"\n##FILTER=<ID=LowQual,Description="Low quality">\n'
              '##INFO=<ID=AC,Number=A,Type=Integer,Description="Original AC">\n'
              '##INFO=<ID=AN,Number=1,Type=Integer,Description="Original AN">\n'
              '##INFO=<ID=NOTE,Number=1,Type=String,Description="Annotation">\n'
              '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    if csq:
        header += '##INFO=<ID=CSQ,Number=.,Type=String,Description="VEP112. Format: Allele|Consequence|Gene">\n'
    header += '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' + '\t'.join('s%d' % i for i in range(samples)) + '\n'
    records = []
    if not empty:
        for pos, ident, ref, alt, qual, filt in [(1, 'rs1', 'A', 'G', '99', 'PASS'),
                                               (2, '.', 'AT', 'A', '.', 'LowQual'),
                                               (3, 'custom', 'C', 'T', '12.5', '.')]:
            info = 'AC=9;AN=20;NOTE=keep_me'
            if csq:
                info += ';CSQ=%s|missense_variant|GENE1,%s|intron_variant|GENE2' % (alt, alt)
            records.append('\t'.join([chromosome, str(pos), ident, ref, alt, qual, filt, info, 'GT'] +
                                     ['0/1']*samples) + '\n')
    # Standard gzip input deliberately has no index; the adapter needs neither.
    with gzip.open(path, 'wt') as handle:
        handle.write(header + ''.join(records))
    return header, records


def direct_run(tmp_path, source, chromosome='chr22'):
    (tmp_path/'receipt-prefix.json').write_text(json.dumps(dict(
        unit_id='explicit_block', chromosome=chromosome, source_vcf=str(source)))[:-1]+',')
    return subprocess.run(['bash', str(ROOT/'scripts/make_sites_catalog.sh'), 'explicit_block',
        chromosome, str(source), str(ROOT/'scripts/check_sites_records.awk'), '1'],
        cwd=tmp_path, capture_output=True, text=True)


@pytest.mark.skipif(not HAS_BCFTOOLS, reason='Native bcftools required')
@pytest.mark.parametrize('csq,empty', [(True, False), (False, False), (True, True)])
def test_records_annotations_and_empty_inputs(tmp_path, csq, empty):
    source = tmp_path/'source.vcf.gz'
    header, records = fixture_vcf(source, csq=csq, empty=empty)
    result = direct_run(tmp_path, source)
    assert result.returncode == 0, result.stdout+result.stderr
    output = tmp_path/'explicit_block.sites.vcf.gz'
    with gzip.open(output, 'rt') as handle:
        lines = handle.readlines()
    assert [line for line in lines if not line.startswith('#')] == ['\t'.join(r.split('\t')[:8])+'\n' for r in records]
    assert len(next(line for line in lines if line.startswith('#CHROM')).rstrip().split('\t')) == 8
    assert '##VEP="v112"\n' in lines
    if csq:
        assert next(line for line in lines if line.startswith('##INFO=<ID=CSQ')) in header
    receipt = json.loads((tmp_path/'explicit_block.receipt.json').read_text())
    assert receipt['input_records'] == receipt['output_records'] == len(records)
    assert receipt['output_samples'] == 0 and receipt['csq_present'] == csq and receipt['status'] == 'PASS'
    assert (tmp_path/'explicit_block.sites.vcf.gz.csi').stat().st_size > 0


@pytest.mark.skipif(not HAS_BCFTOOLS, reason='Native bcftools required')
def test_bad_chromosome_never_emits_success_receipt(tmp_path):
    source = tmp_path/'source.vcf.gz'; fixture_vcf(source)
    result = direct_run(tmp_path, source, chromosome='chr21')
    assert result.returncode != 0 and 'chromosome disagrees' in result.stderr
    assert not (tmp_path/'explicit_block.receipt.json').exists()


@pytest.mark.skipif(not (HAS_NEXTFLOW and HAS_BCFTOOLS), reason='Nextflow and native bcftools required')
def test_block_subset_resume_and_durable_publication(tmp_path):
    source = tmp_path/"source with spaces'quote.vcf.gz"
    _, expected = fixture_vcf(source, samples=8877)
    source2 = tmp_path/'another-name.vcf.gz'; fixture_vcf(source2)
    manifest = tmp_path/'blocks.tsv'
    manifest.write_text('unit_id\tchromosome\tvcf\nchr22_block12\tchr22\t'+source.name+'\n')

    def run(name, resume=False):
        command = ['nextflow', '-C', str(ROOT/'sites_catalog.config'), 'run', str(ROOT/'sites_catalog.nf'),
            '-ansi-log', 'false', '--sites_manifest', str(manifest), '--outdir', str(tmp_path/'results'),
            '-with-trace', str(tmp_path/(name+'.trace.tsv'))]
        if resume:
            command.append('-resume')
        result = subprocess.run(command, cwd=tmp_path, capture_output=True, text=True, timeout=120)
        assert result.returncode == 0, result.stdout+result.stderr
        return list(csv.DictReader((tmp_path/(name+'.trace.tsv')).open(), delimiter='\t'))

    first = run('one')
    assert len(first) == 1 and first[0]['name'] == 'BLOCK_SITES_CATALOG (chr22_block12)'
    manifest.write_text(manifest.read_text()+'chr22_block13\t22\t'+source2.name+'\n')
    second = run('two', True)
    statuses = {r['name']: r['status'] for r in second}
    assert statuses == {'BLOCK_SITES_CATALOG (chr22_block12)': 'CACHED',
                        'BLOCK_SITES_CATALOG (chr22_block13)': 'COMPLETED'}
    sites = tmp_path/'results/sites'
    assert len(list(sites.iterdir())) == 6
    receipt = json.loads((sites/'chr22_block12.receipt.json').read_text())
    assert receipt['source_vcf'] == str(source)
    assert receipt['input_records'] == receipt['output_records'] == len(expected)
    assert receipt['output_samples'] == 0 and receipt['csq_present']
    manifest.write_text('unit_id\tchromosome\tvcf\nchr22_block13\t22\t'+source2.name+'\n')
    assert all(r['status'] == 'CACHED' for r in run('subset', True))
    saved = {p: p.read_bytes() for p in sites.iterdir()}
    shutil.rmtree(tmp_path/'work')
    assert all(not p.is_symlink() and p.read_bytes() == value for p, value in saved.items())
    result = subprocess.run(['bcftools', 'view', '-H', '-r', 'chr22:1-3', str(sites/'chr22_block12.sites.vcf.gz')],
                            capture_output=True, text=True, check=True)
    assert len(result.stdout.splitlines()) == 3


@pytest.mark.skipif(not HAS_NEXTFLOW, reason='Nextflow required')
def test_manifest_validation_without_vcf_inspection(tmp_path):
    shutil.copytree(ROOT/'lib', tmp_path/'lib')
    (tmp_path/'input.vcf.gz').write_text('intentionally not a VCF: manifest validation must not inspect it')
    (tmp_path/'main.nf').write_text(r'''
workflow {
    def m = file('blocks.tsv')
    def head = 'unit_id\tchromosome\tvcf\n'
    def valid = 'explicit_name\tchr22\tinput.vcf.gz\n'
    m.text = head + valid
    assert SitesCatalogManifest.load(m)[0].unit_id == 'explicit_name'
    [head, 'chromosome\tvcf\nchr22\tinput.vcf.gz\n', head+valid+valid,
     head+'\tchr22\tinput.vcf.gz\n', head+'../unsafe\tchr22\tinput.vcf.gz\n',
     head+'unit\tchr23\tinput.vcf.gz\n', head+'unit\tchr22\t\n',
     head+'unit\tchr22\tmissing.vcf.gz\n', head+'unit\tchr22\n',
     'unit_id\tchromosome\tvcf\tvcf\nu\t22\tx\ty\n'].each { invalid ->
        m.text = invalid
        try { SitesCatalogManifest.load(m); assert false: invalid }
        catch (IllegalArgumentException expected) { }
    }
    println 'MANIFEST_VALIDATION_PASS'
}
''')
    result = subprocess.run(['nextflow', 'run', 'main.nf', '-ansi-log', 'false'], cwd=tmp_path,
                            capture_output=True, text=True, timeout=120)
    assert result.returncode == 0, result.stdout+result.stderr
    assert 'MANIFEST_VALIDATION_PASS' in result.stdout
