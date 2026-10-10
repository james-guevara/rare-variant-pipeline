import gzip
import json
from pathlib import Path
import shutil
import subprocess
import sys
from types import SimpleNamespace

import duckdb
import pysam
import pytest

from test_filtered_carriers import rows, extract, fixture as autosome_fixture
from test_carrier_frequencies import psam as autosome_psam
from carrier_frequencies import (PAR, SEX_POLICY, FrequencyCounter, frequency_region,
                                 expected_ploidy, normalize_sex, validate_chromosome)
from extract_exact_carriers import ValidationError

ROOT=Path(__file__).resolve().parents[1]


@pytest.mark.parametrize('chromosome,interval',[('X',p) for p in PAR['X']]+[('Y',p) for p in PAR['Y']])
def test_inclusive_boundaries(chromosome,interval):
    start,end=interval
    for alias in [chromosome,'chr'+chromosome]:
        assert frequency_region(alias,start-1)==chromosome+'_nonPAR'
        assert frequency_region(alias,end+1)==chromosome+'_nonPAR'
        for pos in [start,start+1,end-1,end]:
            if chromosome=='Y':
                with pytest.raises(ValidationError,match='Y-PAR'):frequency_region(alias,pos)
            else:
                assert frequency_region(alias,pos)=='X_PAR'


@pytest.mark.parametrize('region,sex,result',[
    ('X_PAR','male',(2,None)),('X_PAR','female',(2,None)),
    ('X_nonPAR','male',(1,None)),('X_nonPAR','female',(2,None)),
    ('Y_nonPAR','male',(1,None)),('Y_nonPAR','female',(None,'excluded_female_Y')),
    *[(r,'unknown',(None,'excluded_unknown_sex')) for r in ['X_PAR','X_nonPAR','Y_nonPAR']],
    ('autosome','unknown',(2,None)),
])
def test_expected_ploidy(region,sex,result):
    assert expected_ploidy(region,sex)==result


def test_policy_and_sex_encoding():
    for c in ['X','Y','chrX','chrY']:
        with pytest.raises(ValidationError):validate_chromosome(c)
        validate_chromosome(c,SEX_POLICY)
    with pytest.raises(ValidationError):validate_chromosome('X','guess')
    assert normalize_sex('1')=='male' and normalize_sex('2')=='female'
    assert all(normalize_sex(s)=='unknown' for s in ['0','-9','NA','.','M','unknown'])


def fixture(tmp,chromosome,index='tbi',par_candidate=False):
    tmp.mkdir(parents=True,exist_ok=True)
    positions=[10000,10001,10002,2781479,2781480,155701382,155701383,156030895,156030896] if chromosome=='X' else [10000,2781480,56887902,57217416]
    if par_candidate: positions=sorted(set(positions+[10001]))
    source=tmp/'source.vcf'
    header=f'''##fileformat=VCFv4.2
##contig=<ID=chr{chromosome},length=160000000>
##FILTER=<ID=q10,Description="Low quality">
##INFO=<ID=AC,Number=A,Type=Integer,Description="Source AC">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Source AN">
##INFO=<ID=AF,Number=A,Type=Float,Description="Source AF">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tM1\tF1\tU1\tD1\tM2\tM3
'''
    records=[]
    for pos in positions:
        # The Y PAR candidate is deliberately UNMATCHED to test candidate-level rejection.
        if chromosome=='Y' and pos==10001:continue
        gt='0|1' if chromosome=='X' and pos==10002 else '1|1'
        alt='*' if pos==positions[-1] else 'G'
        records.append(f'chr{chromosome}\t{pos}\t.\tA\t{alt}\t.\tq10\tAC=9;AN=12;AF=0.5\tGT\t{gt}\t0/1\t1/1\t1/1\t1/.\t0\n')
    source.write_text(header+''.join(records))
    vcf=tmp/'source.vcf.gz';pysam.tabix_compress(str(source),str(vcf),force=True);pysam.tabix_index(str(vcf),preset='vcf',force=True,csi=index=='csi')
    con=duckdb.connect()
    con.execute('CREATE TABLE c(CHROM VARCHAR,POS BIGINT,REF VARCHAR,ALT VARCHAR,Gene VARCHAR,Feature VARCHAR,SYMBOL VARCHAR,Consequence VARCHAR,LoF VARCHAR,tier VARCHAR,allele_class VARCHAR,pcf_retained BOOLEAN)')
    for kind in ['missense','lof_hc']:
        for pos in positions:
            alt='*' if pos==positions[-1] else 'G'
            if kind=='missense' and alt=='*':continue
            con.execute('INSERT INTO c VALUES (?,?,?,?,?,?,?,?,?,?,?,?)',[chromosome,pos,'A',alt,'GENE','TX','GENE','missense_variant' if kind=='missense' else 'frameshift_variant',None if kind=='missense' else 'HC','miss_t1' if kind=='missense' else None,'sequence' if alt=='G' else 'spanning_deletion',True])
        con.execute('COPY c TO ? (FORMAT PARQUET)',[str(tmp/(kind+'.parquet'))]);con.execute('DELETE FROM c')
    con.close()
    psam=tmp/'samples.psam';psam.write_text('#IID\tSEX\tparticipant_id\tfrequency_representative\tunrelated\n'
        'M1\t1\tP1\t1\t1\nF1\t2\tP2\t1\t1\nU1\t0\tP3\t1\t1\nD1\t1\tP1\t0\t1\nM2\t1\tP4\t1\t1\nM3\t1\tP5\t1\t0\n')
    meta=tmp/'meta.json';meta.write_text(json.dumps(dict(unit_id='chr'+chromosome+'_pilot',chromosome='chr'+chromosome)))
    return SimpleNamespace(metadata=str(meta),missense=str(tmp/'missense.parquet'),lof_hc=str(tmp/'lof_hc.parquet'),vcf=str(vcf),index=str(vcf)+'.'+index,psam=str(psam),outdir=str(tmp/'out'),expected_hc=None,expected_missense=None,compute_frequencies=True,sex_chromosome_policy=SEX_POLICY)


@pytest.mark.parametrize('chromosome',['X','Y'])
@pytest.mark.parametrize('index',['tbi','csi'])
def test_extraction_raw_parity_ploidy_counts_audits_and_stars(tmp_path,chromosome,index,capsys):
    a=fixture(tmp_path,chromosome,index)
    a.compute_frequencies=False;extract(a);raw=Path(a.outdir)
    a.outdir=str(tmp_path/'corrected');a.compute_frequencies=True;extract(a);out=Path(a.outdir)
    for p in raw.iterdir():
        if p.name=='receipt.json':continue
        other=out/p.name
        assert (gzip.decompress(p.read_bytes()) if p.suffix=='.gz' else p.read_bytes())==(gzip.decompress(other.read_bytes()) if other.suffix=='.gz' else other.read_bytes())
    data=rows(out/'variant_frequencies.tsv');bypos={int(r['POS']):r for r in data}
    if chromosome=='X':
        for pos in [10000,2781480,155701382,156030896]:
            r=bypos[pos]
            assert r['frequency_region']=='X_nonPAR'
            assert (r['cohort_ac'],r['cohort_an'],r['unrelated_ac'],r['unrelated_an'])==('2','4','2','3')
        for pos in [10001,2781479,155701383,156030895]:
            r=bypos[pos]
            assert r['frequency_region']=='X_PAR'
            assert (r['cohort_ac'],r['cohort_an'])==('4','5')
        assert (bypos[10002]['cohort_ac'],bypos[10002]['cohort_an'])==('3','5')  # male and female heterozygotes in PAR
    else:
        for r in data:
            assert r['frequency_region']=='Y_nonPAR'
            assert (r['cohort_ac'],r['cohort_an'],r['unrelated_ac'],r['unrelated_an'])==('1','2','1','1')
    assert all(r['source_info_ac']=='9' and r['source_info_an']=='12' and r['source_ac_an']=='0.75' and r['source_info_af']=='0.5' for r in data)
    assert data[-1]['allele_class']=='spanning_deletion'
    assert all(r['matched']=='True' for r in data)
    receipt=json.loads((out/'receipt.json').read_text());f=receipt['frequencies']
    assert receipt['samples_in_source']==6
    assert f['sample_selection']['cohort']==5 and f['sample_selection']['unrelated']==4
    assert f['sample_selection']['by_sex']['cohort']==dict(male=3,female=1,unknown=1)
    for cls in f['by_allele_class'].values():
        n=cls['matched_distinct_alleles']
        assert sum(cls['matched_by_region'].values())==n
        for group,c in cls['sample_sets'].items():
            assert sum(c['by_reason'].values())==n*f['sample_selection'][group]
            assert c['by_reason']['excluded_unknown_sex']==n
            assert c['by_reason']['excluded_female_Y']==(n if chromosome=='Y' else 0)
    assert f['by_allele_class']['spanning_deletion']['matched_distinct_alleles']==1
    assert any(r['sample']=='U1' for r in rows(out/'carriers.tsv.gz'))
    assert any(r['sample']=='F1' for r in rows(out/'carriers.tsv.gz'))
    assert len(bypos)==len(data)  # no duplicated frequencies for cross-type annotations
    assert 'M1' not in capsys.readouterr().out


def test_unmatched_y_par_candidate_fails_before_vcf_access(tmp_path,monkeypatch):
    a=fixture(tmp_path,'Y',par_candidate=True)
    def forbidden(*args,**kwargs):raise AssertionError('Must reject Y-PAR candidate before opening VCF')
    monkeypatch.setattr(pysam,'VariantFile',forbidden)
    with pytest.raises(ValidationError,match='Y-PAR'):extract(a)
    out=Path(a.outdir)
    assert json.loads((out/'receipt.json').read_text())['status']=='failed'
    assert not list(out.glob('*.tsv')) and not list(out.glob('*.gz'))


def test_autosomal_counts_and_policy_unchanged_by_flag(tmp_path):
    a=autosome_fixture(tmp_path);a.psam=str(autosome_psam(tmp_path/'samples.psam'));a.compute_frequencies=True
    extract(a);before=Path(a.outdir)
    a.sex_chromosome_policy=SEX_POLICY;a.outdir=str(tmp_path/'with-policy');extract(a)
    after=Path(a.outdir)
    for name in ['variant_frequencies.tsv','frequency_audit.tsv','source_frequencies.tsv']:
        assert (before/name).read_bytes()==(after/name).read_bytes()
    old=json.loads((before/'receipt.json').read_text());new=json.loads((after/'receipt.json').read_text())
    assert old['frequencies']==new['frequencies']
    assert old['frequency_status']==new['frequency_status']


@pytest.mark.skipif(not shutil.which('nextflow'),reason='Nextflow required')
def test_nextflow_xy_wiring_and_resume(tmp_path):
    x=fixture(tmp_path/'X','X');y=fixture(tmp_path/'Y','Y','csi')
    manifest=tmp_path/'manifest.tsv';manifest.write_text('unit_id\tchromosome\tmissense\tlof_hc\tvcf\tindex\n'+''.join(
        f'chr{c}_pilot\tchr{c}\t{a.missense}\t{a.lof_hc}\t{a.vcf}\t{a.index}\n' for c,a in [('X',x),('Y',y)]))
    config=tmp_path/'python.config';config.write_text('process.beforeScript = "export PATH='+str(Path(sys.executable).parent)+':\\$PATH"\n')
    cmd=['nextflow','-C',str(ROOT/'filtered_carriers.config')+','+str(config),'run',str(ROOT/'filtered_carriers.nf'),'-ansi-log','false',
         '--carrier_manifest',str(manifest),'--psam',x.psam,'--compute_frequencies','true','--select_units','all',
         '--outdir',str(tmp_path/'published'),'--carrier_memory','1 GB']
    fail=subprocess.run(cmd,cwd=tmp_path,capture_output=True,text=True,timeout=120)
    assert fail.returncode and '--sex_chromosome_policy grch38_x_only_par' in fail.stdout+fail.stderr
    cmd+=['--sex_chromosome_policy',SEX_POLICY]
    for n,status in [(1,'COMPLETED'),(2,'CACHED')]:
        result=subprocess.run(cmd+['-with-trace',str(tmp_path/f'{n}.trace')]+(['-resume'] if n==2 else []),cwd=tmp_path,capture_output=True,text=True,timeout=120)
        assert result.returncode==0,result.stdout+result.stderr
        assert [r['status'] for r in rows(tmp_path/f'{n}.trace')]==[status,status]
    shutil.rmtree(tmp_path/'work')
    for c in ['X','Y']:
        out=tmp_path/f'published/filtered-carriers/chr{c}_pilot'
        receipt=json.loads((out/'receipt.json').read_text())
        assert receipt['frequency_policy']['sex_chromosomes']==SEX_POLICY
        assert receipt['sources']['sex_chromosome_policy']==SEX_POLICY
        assert (out/'variant_frequencies.tsv').is_file() and not (out/'variant_frequencies.tsv').is_symlink()


@pytest.mark.parametrize('chromosome,position,sex,reason',[
    ('X',10001,'0','excluded_unknown_sex'),('X',2781480,'NA','excluded_unknown_sex'),
    ('Y',2781480,'0','excluded_unknown_sex'),('Y',2781480,'2','excluded_female_Y')])
def test_excluded_sex_never_contributes_or_reads_gt(chromosome,position,sex,reason):
    metadata=[{'#IID':'S','SEX':sex,'participant_id':'P','frequency_representative':'1','unrelated':'1'}]
    class NoGenotypes:
        def __getitem__(self,key):raise AssertionError('Excluded sex must not be counted')
    counter=FrequencyCounter(metadata,chromosome,SEX_POLICY)
    value,audit=counter.count(SimpleNamespace(pos=position,samples=NoGenotypes()),'spanning_deletion')
    for group in ['cohort','unrelated']:
        assert value[group+'_ac']==value[group+'_an']==0 and value[group+'_af'] is None
        assert value[group+'_excluded_genotypes']==1
        assert audit[group]=={reason:1}
        assert counter.receipt()['by_allele_class']['spanning_deletion']['sample_sets'][group]['zero_an_alleles']==1
