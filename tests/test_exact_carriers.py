import csv
import gzip
import importlib.util
import json
from pathlib import Path
import shutil
import subprocess
import sys
from argparse import Namespace

import pysam
import pytest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('extract', ROOT/'scripts/extract_exact_carriers.py')
extract = importlib.util.module_from_spec(spec); spec.loader.exec_module(extract)


def fixture(tmp_path, extra='', index='tbi'):
    tmp_path.mkdir(exist_ok=True)
    source = tmp_path/'original.vcf'
    source.write_text('''##fileformat=VCFv4.2
##contig=<ID=chr22,length=1000>
##FILTER=<ID=q10,Description="Low quality">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Quality">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allele depths">
##FORMAT=<ID=FT,Number=1,Type=String,Description="Filter">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4
chr22\t10\t.\tA\tG\t.\tPASS\t.\tGT:GQ:DP:AD:FT\t0/1:50:20:10,10:PASS\t1|1:60:25:0,25:PASS\t0/0:40:18:18,0:PASS\t./.:.:.:.:.
chr22\t10\t.\tA\tT\t.\tPASS\t.\tGT\t0/0\t0/0\t0/1\t0/0
chr22\t12\t.\tAAA\tA\t.\tPASS\t.\tGT\t0/1\t0/0\t0/0\t0/0
chr22\t13\t.\tA\tC\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/0
chr22\t20\t.\tC\tT\t.\tq10\t.\tGT\t1/.\t1\t0\t./.
chr22\t30\t.\tG\tA\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0\t./.
chr22\t40\t.\tT\tG\t.\tPASS\t.\tGT\t0/1\t0/1\t0/1\t0/1
'''+extra)
    vcf = tmp_path/'source.vcf.gz'; pysam.tabix_compress(str(source), str(vcf), force=True)
    pysam.tabix_index(str(vcf), preset='vcf', force=True, csi=index=='csi')
    loftee = tmp_path/'loftee.tsv'
    loftee.write_text('CHROM\tPOS\tREF\tALT\tGene\tFeature\tSYMBOL\tLoF\n'
                      '22\t10\tA\tG\tGENE1\tTX1\tG1\tHC\n'
                      '22\t13\tA\tC\tGENE1\tTX1\tG1\tHC\n'
                      '22\t20\tC\tT\tGENE1\tTX1\tG1\tHC\n'
                      '22\t30\tG\tA\tGENE2\tTX2\tG2\tHC\n'
                      '22\t40\tT\tG\tGENE3\tTX3\tG3\tLC\n'
                      '22\t10\tC\tG\tGENE4\tTX4\tG4\tHC\n'
                      '22\t50\tG\tC\tGENE4\tTX4\tG4\tHC\n')
    idx = Path(str(vcf)+'.'+index)
    metadata = tmp_path/'unit.json'
    metadata.write_text(json.dumps(dict(unit_id='block12', chromosome='chr22', vcf=str(vcf), index=str(idx), loftee=str(loftee))))
    return Namespace(metadata=str(metadata), loftee=str(loftee), vcf=str(vcf), index=str(idx),
                     outdir=str(tmp_path/'out'), expected_hc=None)


@pytest.mark.parametrize('index', ['tbi', 'csi'])
def test_exact_alleles_genotypes_and_counts(tmp_path, index, capsys):
    args = fixture(tmp_path, index=index)
    before = {p: p.read_bytes() for p in [Path(args.vcf), Path(args.index), Path(args.loftee)]}
    extract.extract(args)
    out = Path(args.outdir)
    with gzip.open(out/'carriers.tsv.gz', 'rt') as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    assert len(rows) == 4  # excludes other ALT, overlapping deletion, LC, and hom-ref
    assert {(r['POS'], r['sample']) for r in rows} == {('10','S1'),('10','S2'),('20','S1'),('20','S2')}
    assert rows[1]['GT'] == '1|1' and rows[1]['alt_dosage'] == '2'
    assert rows[0]['GQ'] == '50' and rows[0]['DP'] == '20' and rows[0]['AD'] == '10,10' and rows[0]['FT'] == 'PASS'
    assert rows[2]['GT'] == '1/.' and rows[2]['GQ'] == '.' and rows[2]['site_FILTER'] == 'q10'
    assert rows[3]['GT'] == '1'
    receipt = json.loads((out/'receipt.json').read_text())
    assert receipt['candidate_hc_variants'] == 6 and receipt['non_hc_annotation_rows'] == 1
    assert receipt['matched_candidate_variants'] == 4 and receipt['unmatched_candidate_variants'] == 2
    assert receipt['matched_variants_without_carriers'] == 2
    assert receipt['carrier_records'] == 4 and receipt['samples_with_hc_plof'] == 2 and receipt['genes_with_hc_plof'] == 1
    assert receipt['partial_call_carrier_records'] == 1 and receipt['observed_alt_alleles'] == 5
    assert 'S1\tGENE1\t2\t0\n' in (out/'sample_gene_burden.tsv').read_text()
    assert 'S3\t0\t0\n' in (out/'sample_burden.tsv').read_text()
    assert len((out/'samples.tsv').read_text().splitlines()) == 5
    assert all(p.read_bytes() == value for p, value in before.items())
    assert 'S1' not in capsys.readouterr().out


@pytest.mark.parametrize('case', ['duplicate_candidate', 'duplicate_source', 'multiallelic', 'missing_index', 'wrong_hc'])
def test_invalid_inputs_fail_without_products(tmp_path, case):
    extra = ''
    if case == 'duplicate_source':
        extra = 'chr22\t50\t.\tG\tC\t.\tPASS\t.\tGT\t0/1\t0/0\t0/0\t0/0\n'*2
    if case == 'multiallelic':
        extra = 'chr22\t50\t.\tG\tC,T\t.\tPASS\t.\tGT\t0/2\t0/0\t0/0\t0/0\n'
    args = fixture(tmp_path, extra=extra)
    if case == 'duplicate_candidate':
        with open(args.loftee, 'a') as handle:
            handle.write('22\t10\tA\tG\tGENE1\tTX1\tG1\tHC\n')
    if case == 'missing_index': Path(args.index).unlink()
    if case == 'wrong_hc': args.expected_hc = 302
    with pytest.raises(Exception): extract.extract(args)
    receipt = json.loads((Path(args.outdir)/'receipt.json').read_text())
    assert receipt['status'] == 'failed'
    assert not (Path(args.outdir)/'carriers.tsv.gz').exists()


def test_empty_hc_set(tmp_path):
    args = fixture(tmp_path)
    path = Path(args.loftee); path.write_text(path.read_text().replace('\tHC\n','\tLC\n'))
    extract.extract(args)
    receipt = json.loads((Path(args.outdir)/'receipt.json').read_text())
    assert receipt['candidate_hc_variants'] == receipt['carrier_records'] == 0
    assert len((Path(args.outdir)/'sample_burden.tsv').read_text().splitlines()) == 5


@pytest.mark.skipif(not shutil.which('nextflow'), reason='Nextflow required')
def test_nextflow_real_extraction_subset_resume(tmp_path):
    a = star_fixture(tmp_path/'block12'); b = fixture(tmp_path/'block19', index='csi')
    manifest = tmp_path/'manifest.tsv'
    manifest.write_text('unit_id\tchromosome\tloftee\tvcf\tindex\n'+''.join(
        f'{unit}\tchr22\t{args.loftee}\t{args.vcf}\t{args.index}\n' for unit, args in [('block12',a),('block19',b)]))
    config = tmp_path/'python.config'
    # Execute the real extraction script using the test environment's installed pysam.
    config.write_text('process.beforeScript = "export PATH='+str(Path(sys.executable).parent)+':\\$PATH"\n')
    base = ['nextflow', '-C', str(ROOT/'carriers.config')+','+str(config), 'run', str(ROOT/'carriers.nf'),
            '-ansi-log', 'false', '--carrier_manifest', str(manifest), '--outdir', str(tmp_path/'published'),
            '--carrier_memory', '1 GB']
    def run(name, units, resume=False):
        command = base+['--select_units',units,'-with-trace',str(tmp_path/(name+'.trace'))]
        if resume: command += ['-resume']
        result = subprocess.run(command, cwd=tmp_path, capture_output=True, text=True, timeout=120)
        assert result.returncode == 0, result.stdout+result.stderr
        return list(csv.DictReader((tmp_path/(name+'.trace')).open(), delimiter='\t'))
    assert len(run('one','block12')) == 1
    second = run('two','all',True)
    assert {r['name']:r['status'] for r in second} == {'EXACT_HC_CARRIERS (block12)':'CACHED','EXACT_HC_CARRIERS (block19)':'COMPLETED'}
    assert all(r['status']=='CACHED' for r in run('subset','block19',True))
    out = tmp_path/'published/carriers/block12'
    assert json.loads((out/'receipt.json').read_text())['carrier_records'] == 7
    shutil.rmtree(tmp_path/'work')
    assert (out/'carriers.tsv.gz').is_file() and not (out/'carriers.tsv.gz').is_symlink()


def test_optional_format_headers_absent(tmp_path):
    args = fixture(tmp_path)
    source = tmp_path/'minimal.vcf'
    source.write_text('##fileformat=VCFv4.2\n##contig=<ID=chr22,length=1000>\n'
                      '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
                      '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n'
                      'chr22\t10\t.\tA\tG\t.\tPASS\t.\tGT\t0/1\n')
    pysam.tabix_compress(str(source), args.vcf, force=True)
    pysam.tabix_index(args.vcf, preset='vcf', force=True)
    extract.extract(args)
    with gzip.open(Path(args.outdir)/'carriers.tsv.gz','rt') as handle:
        row = next(csv.DictReader(handle, delimiter='\t'))
    assert all(row[k]=='.' for k in ['GQ','DP','AD','FT'])


def star_fixture(tmp_path, index='tbi'):
    args = fixture(tmp_path, index=index, extra=(
        'chr22\t59\t.\tAAA\tA\t.\tPASS\t.\tGT\t0/1\t0/1\t0/1\t0/1\n'
        'chr22\t60\t.\tAT\t*\t.\tPASS\t.\tGT:AD\t0/1:5,5\t1|1:0,10\t1/.:0,5\t./.:.\n'
        'chr22\t60\t.\tAT\tATT\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0\t1/1\n'
        'chr22\t60\t.\tA\t*\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/1\n'
        'chr22\t61\t.\tC\tT\t.\tPASS\t.\tGT\t0/1\t0/1\t0/1\t0/1\n'
        'chr22\t62\t.\tGG\t*\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0\t./.\n'))
    with open(args.loftee, 'a') as handle:
        handle.write('22\t60\tAT\t*\tGENE1\tTX1\tG1\tHC\n'
                     '22\t61\tC\t*\tGENE4\tTX4\tG4\tHC\n'
                     '22\t62\tGG\t*\tGENE2\tTX2\tG2\tHC\n')
    return args


@pytest.mark.parametrize('index', ['tbi','csi'])
def test_spanning_deletions_exact_matching_and_separate_burdens(tmp_path, index):
    args = star_fixture(tmp_path, index)
    args.expected_hc = 9  # guard includes both classes, without discarding stars
    extract.extract(args)
    out = Path(args.outdir)
    with gzip.open(out/'carriers.tsv.gz', 'rt') as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    stars = [r for r in rows if r['allele_class'] == 'spanning_deletion']
    assert len(rows) == 7 and len(stars) == 3
    assert {r['sample'] for r in stars} == {'S1','S2','S3'}
    assert all(r['POS']=='60' and r['REF']=='AT' and r['ALT']=='*' for r in stars)
    assert {r['GT'] for r in stars} == {'0/1','1|1','1/.'}
    assert stars[1]['alt_dosage'] == '2' and stars[1]['AD'] == '0,10'
    receipt = json.loads((out/'receipt.json').read_text())
    assert receipt['schema_version'] == 2
    assert receipt['candidate_hc_variants'] == 9 and receipt['carrier_records'] == 7
    assert receipt['samples_with_hc_plof'] == 3 and receipt['genes_with_hc_plof'] == 1
    seq, star = (receipt['by_allele_class'][k] for k in ['sequence','spanning_deletion'])
    assert seq['candidate_hc_variants'] == 6 and seq['carrier_records'] == 4
    assert star == dict(candidate_hc_variants=3, matched_candidate_variants=2,
                        unmatched_candidate_variants=1, matched_variants_without_carriers=1,
                        carrier_records=3, samples_with_hc_plof=3, genes_with_hc_plof=1,
                        partial_call_carrier_records=1, observed_alt_alleles=4)
    assert 'S1\tGENE1\t2\t1\n' in (out/'sample_gene_burden.tsv').read_text()
    assert 'S3\tGENE1\t0\t1\n' in (out/'sample_gene_burden.tsv').read_text()
    assert 'S4\t0\t0\n' in (out/'sample_burden.tsv').read_text()
    assert 'GENE1\t4\t2\t3\t3\n' in (out/'gene_burden.tsv').read_text()
    with (out/'unmatched.tsv').open() as handle:
        unmatched = list(csv.DictReader(handle, delimiter='\t'))
    assert [(r['POS'],r['ALT']) for r in unmatched if r['allele_class']=='spanning_deletion'] == [('61','*')]
    with (out/'candidates.tsv').open() as handle:
        audit = list(csv.DictReader(handle, delimiter='\t'))
    assert sum(r['allele_class']=='spanning_deletion' for r in audit) == 3


def test_multiallelic_source_with_star_still_fails(tmp_path):
    args = fixture(tmp_path, extra='chr22\t50\t.\tG\tC,*\t.\tPASS\t.\tGT\t0/2\t0/1\t0/0\t0/0\n')
    with open(args.loftee, 'a') as handle:
        handle.write('22\t50\tG\t*\tGENE4\tTX4\tG4\tHC\n')
    with pytest.raises(extract.ValidationError, match='Multiallelic source'):
        extract.extract(args)
    assert not (Path(args.outdir)/'carriers.tsv.gz').exists()


@pytest.mark.parametrize('alt', ['<*>','<DEL>','G,*','.',''])
def test_other_symbolic_or_multiple_candidate_alts_not_conflated_with_star(tmp_path, alt):
    args = fixture(tmp_path)
    with open(args.loftee,'a') as handle:
        handle.write(f'22\t60\tA\t{alt}\tGENE1\tTX1\tG1\tHC\n')
    with pytest.raises(extract.ValidationError):
        extract.extract(args)
