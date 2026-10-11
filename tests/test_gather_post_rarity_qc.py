from argparse import Namespace
import csv
import json
from pathlib import Path
import shutil
import subprocess
import sys
import pytest
from test_post_rarity_qc import fixture
from qc_filtered_carriers import run as qc, identity, write_table, ValidationError
from test_filtered_carriers import rows
from test_gather_post_qc import rewrite
from gather_post_qc import run as gather
import qc_post_rarity
import gather_post_rarity_qc as adapter
from make_post_rarity_gather_manifest import build

ROOT=Path(__file__).resolve().parents[1]


def make(tmp,empty=False,chrom='chrX'):
    units=[];carriers=[];samples=[];receipts=[]
    for i in range(2):
        a=fixture(tmp/str(i),chromosome=chrom,empty=empty)
        raw=rows(a.carriers)
        for r in raw:r['POS']=str(int(r['POS'])+i*3000000)
        if raw:write_table(a.carriers,list(raw[0]),raw)
        p=Path(a.source_receipt);s=json.loads(p.read_text());s['unit_id']='block'+str(i);s['outputs']['carriers.rare.tsv.gz']=identity(a.carriers);p.write_text(json.dumps(s))
        p=Path(a.metadata);s=json.loads(p.read_text());s['unit_id']='block'+str(i);p.write_text(json.dumps(s));qc(a,qc_post_rarity)
        d=Path(a.outdir);units.append(dict(unit_id='block'+str(i),chromosome=chrom))
        carriers.append(str(d/'carriers.qc.tsv.gz'));samples.append(str(d/'samples.tsv'));receipts.append(str(d/'receipt.json'))
    meta=tmp/'gather.json';meta.write_text(json.dumps(dict(units=units,chromosome=chrom,gather_id='all',manifest_units=['block0','block1'])))
    return Namespace(carriers=carriers,samples=samples,receipts=receipts,metadata=str(meta),outdir=str(tmp/'gather'))


def test_effective_dosage_unions_and_raw_preservation(tmp_path,capsys):
    a=make(tmp_path);before=[identity(p) for p in a.carriers+a.samples+a.receipts];capsys.readouterr();gather(a,adapter)
    out=Path(a.outdir);r=json.loads((out/'receipt.json').read_text())
    assert r['stage']=='post_rarity_qc_gather' and r['carrier_annotation_records']==8
    assert r['psam_identity']==json.loads(Path(a.receipts[0]).read_text())['psam_identity']
    assert r['by_allele_class']['sequence']['carrier_annotation_records']==6
    assert r['by_allele_class']['sequence']['distinct_allele_sample_records']==4
    assert r['by_allele_class']['sequence']['samples_with_carriers']==2
    assert r['by_allele_class']['sequence']['annotation_observed_alt_alleles']==8
    assert r['by_allele_class']['sequence']['distinct_allele_sample_observed_alt_alleles']==6
    assert sum(int(x['observed_alt_alleles']) for x in rows(out/'gene_burden.tsv'))==10
    assert r['by_allele_class']['spanning_deletion']['distinct_allele_sample_records']==2
    assert r['by_allele_class']['spanning_deletion']['samples_with_carriers']==1
    rr=rows(out/'carriers.qc.tsv.gz')
    assert any(x['GT']=='1|1' and x['alt_dosage']=='2' and x['qc_effective_alt_dosage']=='1' for x in rr)
    burden=rows(out/'sample_burden.tsv')
    assert {x['sample'] for x in burden}=={'S1','S2','S3','S4'}
    assert all(x['carrier_variant_records']=='0' for x in burden if x['sample'] in ('S3','S4'))
    assert next(x for x in burden if x['sample']=='S1' and x['candidate_type']=='missense' and x['tier']=='miss_t1')['observed_alt_alleles']=='2'
    assert next(x for x in burden if x['sample']=='S2' and x['candidate_type']=='missense' and x['tier']=='miss_t1')['observed_alt_alleles']=='4'
    assert any(x['tier']=='untiered' and x['carrier_variant_records']=='2' for x in rows(out/'sample_gene_burden.tsv'))
    assert max(int(x['carrier_samples']) for x in rows(out/'gene_burden.tsv'))==2
    assert before==[identity(p) for p in a.carriers+a.samples+a.receipts]
    assert 'S1' not in capsys.readouterr().out


@pytest.mark.parametrize('case',['psam','psam_lineage','policy','stage','hash','count','duplicate','overlap_other_sample','cross_type_sex','effective_dosage','raw_dosage','region','ploidy','samples'])
def test_fail_closed(tmp_path,case):
    a=make(tmp_path)
    if case in ('duplicate','overlap_other_sample'):
        def mutate(rr):
            for r in rr:
                r['POS']=str(int(r['POS'])-3000000)
                if case=='overlap_other_sample':r['sample']='S4'
        rewrite(a,1,'carriers',mutate)
    elif case in ('cross_type_sex','effective_dosage','raw_dosage','region','ploidy'):
        field,value={'cross_type_sex':('qc_sex','female'),'effective_dosage':('qc_effective_alt_dosage','2'),
            'raw_dosage':('alt_dosage','1'),'region':('qc_frequency_region','autosome'),'ploidy':('qc_expected_ploidy','2')}[case]
        def mutate(rr):
            rr[0][field]=value
            if case=='cross_type_sex':rr[0].update(qc_expected_ploidy='2',qc_effective_alt_dosage='2')
        rewrite(a,0,'carriers',mutate)
    elif case=='samples':rewrite(a,1,'samples',lambda rr:rr.reverse())
    elif case=='hash':Path(a.carriers[0]).write_bytes(b'changed')
    else:
        p=Path(a.receipts[1]);r=json.loads(p.read_text())
        if case in ('psam','psam_lineage'):
            r['psam_identity']['sha256']='f'*64
            if case=='psam':r['input_identities']['psam']=r['psam_identity']
        if case=='policy':r['policy']['min_gq']=21
        if case=='stage':r['stage']='post_extraction_qc'
        if case=='count':r['pass_rows']+=1
        p.write_text(json.dumps(r))
    with pytest.raises(ValidationError):gather(a,adapter)
    out=Path(a.outdir)
    assert json.loads((out/'receipt.json').read_text())['status']=='failed'
    assert not (out/'sample_burden.tsv').exists()


def test_8877_zero_samples_and_empty_blocks(tmp_path):
    a=make(tmp_path,empty=True)
    for i in range(2):
        sp=Path(a.samples[i]);write_table(sp,['sample'],(dict(sample=f'S{n}') for n in range(8877)))
        rp=Path(a.receipts[i]);r=json.loads(rp.read_text());r['samples_in_source']=8877;r['outputs']['samples.tsv']=identity(sp);rp.write_text(json.dumps(r))
    gather(a,adapter);out=Path(a.outdir)
    assert len(rows(out/'samples.tsv'))==len(rows(out/'sample_distinct_alleles.tsv'))==8877
    assert len({r['sample'] for r in rows(out/'sample_burden.tsv')})==8877
    assert rows(out/'sample_gene_burden.tsv')==[]


def test_same_gene_cross_type_union(tmp_path):
    a=make(tmp_path)
    for i in range(2):rewrite(a,i,'carriers',lambda rr:rr[0].update(Gene=rr[1]['Gene']))
    gather(a,adapter)
    rr=rows(Path(a.outdir)/'sample_gene_distinct_alleles.tsv')
    assert next(r for r in rr if r['sample']=='S1' and r['allele_class']=='sequence')['distinct_allele_sample_records']=='2'


@pytest.mark.skipif(not shutil.which('nextflow'),reason='Nextflow required')
def test_nextflow_pilot_resume_mixed_full_and_global_psam(tmp_path):
    units=[]
    for chrom in ['chr22','chrX','chrY']:
        a=make(tmp_path/chrom,chrom=chrom)
        for i in range(2):
            rp=Path(a.receipts[i]);r=json.loads(rp.read_text());unit=chrom+'_block'+str(i);r['unit_id']=unit;rp.write_text(json.dumps(r))
            units.append([unit,chrom,a.carriers[i],a.samples[i],a.receipts[i]])
    manifest=tmp_path/'manifest.tsv'
    with manifest.open('w') as f:
        w=csv.writer(f,delimiter='\t');w.writerow(['unit_id','chromosome','carriers','samples','source_receipt']);w.writerows(units)
    # Manifest builder inventories only direct child receipt products.
    inventory=tmp_path/'inventory';inventory.mkdir()
    for unit,chrom,cp,sp,rp in units:
        d=inventory/unit;d.mkdir()
        for p in [cp,sp,rp]:shutil.copy(p,d/Path(p).name)
    assert len(build(inventory))==6
    cfg=tmp_path/'python.config';cfg.write_text('process.beforeScript = "export PATH='+str(Path(sys.executable).parent)+':\\$PATH"\n')
    trace=tmp_path/'trace.tsv'
    cmd=['nextflow','-C',str(ROOT/'gather_post_rarity_qc.config')+','+str(cfg),'run',str(ROOT/'gather_post_rarity_qc.nf'),'-ansi-log','false',
        '--gather_manifest',str(manifest),'--outdir',str(tmp_path/'published'),'--gather_memory','1 GB','-with-trace',str(trace)]
    for selection,gid,expected in [('chr22_block0','pilot',['COMPLETED']),('chr22_block0','pilot',['CACHED']),
            ('chr22_block0,chrX_block0,chrY_block0','mixed',['COMPLETED']*3),('all','all',['COMPLETED']*3),('all','all',['CACHED']*3)]:
        p=subprocess.run(cmd+['--select_units',selection,'--gather_id',gid,'-resume'],cwd=tmp_path,capture_output=True,text=True,timeout=120)
        assert p.returncode==0,p.stdout+p.stderr
        assert sorted(x['status'] for x in rows(trace))==expected
    rp=Path(units[-1][-1]);r=json.loads(rp.read_text());r['psam_identity']['sha256']='f'*64;r['input_identities']['psam']=r['psam_identity'];rp.write_text(json.dumps(r))
    p=subprocess.run(cmd+['--select_units','all','--gather_id','bad'],cwd=tmp_path,capture_output=True,text=True,timeout=120)
    assert p.returncode!=0 and 'PSAM identities differ' in p.stdout+p.stderr
    shutil.rmtree(tmp_path/'work')
    assert json.loads((tmp_path/'published/post-rarity-qc-gather/chrX/all/receipt.json').read_text())['full_manifest_selected']
