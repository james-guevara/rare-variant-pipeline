from argparse import Namespace
import json
from pathlib import Path
import shutil
import subprocess
import sys

import pytest
from test_post_carrier_qc import base_row
from qc_filtered_carriers import identity, write_table, run, evaluate, ValidationError
from test_filtered_carriers import rows
import qc_post_rarity as adapter
from carrier_frequencies import policy_for, SEX_POLICY
from filter_final_rarity import POLICY as RARITY_POLICY
from make_post_rarity_qc_manifest import build

ROOT=Path(__file__).resolve().parents[1]


def context(chromosome='chrX',sex='1'):
    c=object.__new__(adapter.Context);c.chromosome=chromosome;c.sex={'S1':{'1':'male','2':'female'}.get(sex,'unknown')}
    return c


@pytest.mark.parametrize('chromosome,pos,sex,gt,ad,passed,reason,effective',[
    ('chr22',10,'0','0/1','3,1',True,None,1),
    ('chr22',10,'0','1/1','1,9',True,None,2),
    ('chr22',10,'0','1','1,9',False,'unexpected_ploidy',None),
    ('chrX',10001,'1','0|1','3,1',True,None,1),
    ('chrX',2781479,'1','1|1','1,9',True,None,2),
    ('chrX',155701383,'2','0/1','1,3',True,None,1),
    ('chrX',156030895,'1','1','0,10',False,'unexpected_ploidy',None),
    ('chrX',10000,'1','1','1,9',True,None,1),
    ('chrX',2781480,'1','1|1','1,9',True,None,1),
    ('chrX',155701382,'1','0/1','5,5',False,'haploid_heterozygous',None),
    ('chrX',156030896,'1','1/.','0,10',False,'partial_call',None),
    ('chrX',10000,'2','0/1','5,5',True,None,1),
    ('chrX',10000,'2','1','0,10',False,'unexpected_ploidy',None),
    ('chrX',10001,'0','0/1','5,5',False,'unknown_sex',None),
    ('chrX',10000,'0','1/1','0,10',False,'unknown_sex',None),
    ('chrY',10000,'1','1','1,9',True,None,1),
    ('chrY',2781480,'1','1/1','1,9',True,None,1),
    ('chrY',56887902,'1','0|1','5,5',False,'haploid_heterozygous',None),
    ('chrY',57217416,'1','./1','0,10',False,'partial_call',None),
    ('chrY',10000,'2','1/1','0,10',False,'female_y',None),
    ('chrY',10000,'0','1','0,10',False,'unknown_sex',None),
    ('chrY',10000,'1','1/1/1','0,10',False,'unexpected_ploidy',None),
])
def test_ploidy_regions_and_dosage(chromosome,pos,sex,gt,ad,passed,reason,effective):
    r=base_row(CHROM=chromosome,POS=str(pos),GT=gt,AD=ad)
    ab,bad,extra=context(chromosome,sex).evaluate(r)
    assert (not bad)==passed
    if reason:assert reason in bad
    assert extra['qc_effective_alt_dosage']==effective
    assert r['GT']==gt and r['alt_dosage']=='1'  # evaluator never rewrites raw fields


@pytest.mark.parametrize('changes',[
    {},{'GQ':'19.9'},{'DP':'9.9'},{'AD':'3.001,1'},{'AD':'1,3.001'},
    {'AD':'3,1'},{'AD':'1,3'},{'FT':'FAIL'},{'site_FILTER':'.'},
    {'GQ':'nan'},{'DP':'-1'},{'AD':'0,0'},{'AD':'1,2,3'},
    {'GT':'1|1','AD':'1,9'},{'GT':'1/1','AD':'1,8.999'},
])
def test_autosomal_quality_rules_are_reused(changes):
    r=base_row(**changes)
    ab,bad,extra=context('chr22','0').evaluate(r)
    assert (ab,bad)==evaluate(r)


@pytest.mark.parametrize('pos',[10001,2781479,56887903,57217415])
def test_y_par_fails(pos):
    with pytest.raises(ValidationError,match='Y-PAR'):
        context('chrY').evaluate(base_row(CHROM='chrY',POS=str(pos)))


def fixture(tmp,chromosome='chrX',empty=False):
    tmp.mkdir(parents=True,exist_ok=True)
    raw=[base_row(CHROM=chromosome,POS='10000',GT='1|1',alt_dosage='2',AD='1,9'),
         base_row(CHROM=chromosome,POS='10000',GT='1|1',alt_dosage='2',AD='1,9',candidate_type='lof_hc',tier='untiered',Gene='G2',Feature='TX2'),
         base_row(CHROM=chromosome,POS='2781480',GT='0/1',AD='5,5'),
         base_row(CHROM=chromosome,POS='2781481',sample='S2',GT='1/1',alt_dosage='2',AD='1,9'),
         base_row(CHROM=chromosome,POS='2781482',sample='S3',GT='1',AD='1,9'),
         base_row(CHROM=chromosome,POS='2781483',GT='1',AD='1,9',ALT='*',allele_class='spanning_deletion',candidate_type='lof_hc',tier='lof_t1')]
    if empty:raw=[]
    carriers=tmp/'carriers.rare.tsv.gz';write_table(carriers,list(base_row()),raw)
    samples=tmp/'samples.tsv';write_table(samples,['sample'],[dict(sample=s) for s in ['S1','S2','S3','S4']])
    source=tmp/'receipt.json';source.write_text(json.dumps(dict(status='passed',schema_version=1,stage='final_rarity',unit_id='pilot',chromosome=chromosome,
        policy=RARITY_POLICY,reconciliation=dict(passed=True),samples_in_source=4,
        source_frequency_policy=policy_for(chromosome,SEX_POLICY),
        distinct_alleles=dict(input=5,pass_count=5,fail_count=0),candidate_annotations=dict(input=6,pass_count=6,fail_count=0),
        carrier_annotations=dict(input=len(raw),pass_count=len(raw),fail_count=0),
        outputs={'carriers.rare.tsv.gz':identity(carriers),'samples.tsv':identity(samples)})))
    psam=tmp/'samples.psam';psam.write_text('#IID\tSEX\tparticipant_id\tfrequency_representative\tunrelated\n'
        'S1\t1\tP1\t0\t0\nS2\t2\tP2\t1\t1\nS3\t0\tP3\t1\t1\nS4\t1\tP4\t1\t1\n')
    meta=tmp/'meta.json';meta.write_text(json.dumps(dict(unit_id='pilot',chromosome=chromosome)))
    return Namespace(metadata=str(meta),carriers=str(carriers),samples=str(samples),source_receipt=str(source),psam=str(psam),sex_chromosome_policy=SEX_POLICY,outdir=str(tmp/'qc'))


@pytest.mark.parametrize('chromosome',['chrX','chrY'])
def test_adapter_hashes_summaries_and_reconciliation(tmp_path,chromosome,capsys):
    a=fixture(tmp_path,chromosome);inputs=[Path(getattr(a,k)) for k in ['carriers','samples','source_receipt','psam']];before={p:identity(p) for p in inputs}
    run(a,adapter);out=Path(a.outdir);passed=rows(out/'carriers.qc.tsv.gz')
    assert len(passed)==(4 if chromosome=='chrX' else 3)
    hom=[r for r in passed if r['sample']=='S1' and r['POS']=='10000']
    assert len(hom)==2 and all(r['GT']=='1|1' and r['alt_dosage']=='2' and r['qc_effective_alt_dosage']=='1' for r in hom)
    summary=rows(out/'sample_burden.tsv')
    assert all(r['observed_alt_alleles']=='1' for r in summary if r['sample']=='S1' and r['carrier_variant_records']=='1')
    assert all(r['carrier_variant_records']=='0' for r in summary if r['sample']=='S4')
    distinct=rows(out/'sample_distinct_alleles.tsv');assert distinct[0]['sequence_distinct_variant_count']=='1'
    r=json.loads((out/'receipt.json').read_text());assert r['stage']=='post_rarity_qc' and r['reconciliation']['passed']
    assert r['input_rows']==r['pass_rows']+r['fail_rows']==sum(r['failure_combinations'].values())==6
    assert r['independent_failures']['unknown_sex']==1
    assert r['independent_failures']['haploid_heterozygous']==1
    assert r['independent_failures']['female_y']==(chromosome=='chrY')
    assert r['effective_dosage_by_class']['spanning_deletion']['effective_observed_alt_alleles']==1
    assert r['by_candidate_type']['lof_hc']['by_allele_class']['sequence']['by_tier']['untiered']['pass_rows']==1
    assert all(identity(p)==v for p,v in before.items())
    assert 'S1' not in capsys.readouterr().out


@pytest.mark.parametrize('case',['empty','all_fail','changed_carriers','changed_samples','wrong_stage','missing_psam_sample','dosage_mismatch'])
def test_empty_failures_and_provenance(tmp_path,case):
    a=fixture(tmp_path,empty=case=='empty')
    if case=='all_fail':
        data=rows(a.carriers)
        for r in data:r['GQ']='0'
        write_table(a.carriers,list(data[0]),data)
        p=Path(a.source_receipt);src=json.loads(p.read_text());src['outputs']['carriers.rare.tsv.gz']=identity(a.carriers);p.write_text(json.dumps(src))
    if case=='changed_carriers':Path(a.carriers).write_bytes(b'invalid')
    if case=='changed_samples':Path(a.samples).write_text('sample\nS1\n')
    if case=='wrong_stage':
        p=Path(a.source_receipt);s=json.loads(p.read_text());s['stage']='post_extraction_qc';p.write_text(json.dumps(s))
    if case=='missing_psam_sample':
        p=Path(a.psam);p.write_text(p.read_text().replace('S4\t1\tP4\t1\t1\n',''))
    if case=='dosage_mismatch':
        data=rows(a.carriers)
        for r in data:
            if r['POS']=='10000':r['alt_dosage']='1'
        write_table(a.carriers,list(data[0]),data)
        p=Path(a.source_receipt);s=json.loads(p.read_text());s['outputs']['carriers.rare.tsv.gz']=identity(a.carriers);p.write_text(json.dumps(s))
    if case in ['empty','all_fail']:
        run(a,adapter);assert rows(Path(a.outdir)/'carriers.qc.tsv.gz')==[]
        assert len(rows(Path(a.outdir)/'samples.tsv'))==4
    else:
        with pytest.raises(ValidationError):run(a,adapter)
        assert json.loads((Path(a.outdir)/'receipt.json').read_text())['status']=='failed'
        assert not (Path(a.outdir)/'carriers.qc.tsv.gz').exists()


@pytest.mark.skipif(not shutil.which('nextflow'),reason='Nextflow required')
def test_nextflow_post_rarity_resume(tmp_path):
    a=fixture(tmp_path/'pilot');manifest=tmp_path/'manifest.tsv'
    data=build(tmp_path);assert len(data)==1
    import csv
    with manifest.open('w') as f:
        w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(['unit_id','chromosome','carriers','samples','source_receipt']);w.writerows(data)
    cfg=tmp_path/'python.config';cfg.write_text('process.beforeScript = "export PATH='+str(Path(sys.executable).parent)+':\\$PATH"\n')
    cmd=['nextflow','-C',str(ROOT/'post_rarity_qc.config')+','+str(cfg),'run',str(ROOT/'post_rarity_qc.nf'),'-ansi-log','false',
        '--post_rarity_qc_manifest',str(manifest),'--psam',a.psam,'--sex_chromosome_policy',SEX_POLICY,'--select_units','all',
        '--outdir',str(tmp_path/'published'),'--post_rarity_qc_memory','1 GB','-with-trace',str(tmp_path/'trace.tsv')]
    for n,status in [(1,'COMPLETED'),(2,'CACHED')]:
        r=subprocess.run(cmd+(['-resume'] if n==2 else []),cwd=tmp_path,capture_output=True,text=True,timeout=120)
        assert r.returncode==0,r.stdout+r.stderr
        assert [x['status'] for x in rows(tmp_path/'trace.tsv')]==[status]
    shutil.rmtree(tmp_path/'work')
    out=tmp_path/'published/post-rarity-qc/pilot';assert len(rows(out/'carriers.qc.tsv.gz'))==4
    assert not (out/'carriers.qc.tsv.gz').is_symlink()


def test_legacy_gather_rejects_new_policy(tmp_path):
    from gather_post_qc import run as gather
    a=fixture(tmp_path/'pilot');run(a,adapter);out=Path(a.outdir)
    meta=tmp_path/'gather.json'
    meta.write_text(json.dumps(dict(units=[dict(unit_id='pilot',chromosome='chrX')],
        chromosome='chrX',gather_id='pilot',manifest_units=['pilot'])))
    g=Namespace(carriers=[str(out/'carriers.qc.tsv.gz')],samples=[str(out/'samples.tsv')],
        receipts=[str(out/'receipt.json')],metadata=str(meta),outdir=str(tmp_path/'gather'))
    with pytest.raises(ValidationError):gather(g)
    assert json.loads((Path(g.outdir)/'receipt.json').read_text())['status']=='failed'
