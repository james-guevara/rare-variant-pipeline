import csv
import json
from pathlib import Path
import shutil
import subprocess
import sys
from argparse import Namespace
import pytest
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from test_post_carrier_qc import fixture,rows
from qc_filtered_carriers import run as qc,identity,write_table,ValidationError
from gather_post_qc import run


def make(tmp,empty=False):
    units=[];carriers=[];samples=[];receipts=[]
    for i in range(2):
        a=fixture(tmp/str(i),empty=empty);raw=rows(Path(a.carriers))
        for r in raw:r['POS']=str(int(r['POS'])+i*100)
        if raw:write_table(a.carriers,list(raw[0]),raw)
        p=Path(a.source_receipt);s=json.loads(p.read_text());s['unit_id']='block'+str(i);s['outputs']['carriers.tsv.gz']=identity(a.carriers);p.write_text(json.dumps(s))
        p=Path(a.metadata);s=json.loads(p.read_text());s['unit_id']='block'+str(i);p.write_text(json.dumps(s));qc(a)
        d=Path(a.outdir);units.append(dict(unit_id='block'+str(i),chromosome='chr22'))
        carriers.append(str(d/'carriers.qc.tsv.gz'));samples.append(str(d/'samples.tsv'));receipts.append(str(d/'receipt.json'))
    meta=tmp/'gather.json';meta.write_text(json.dumps(dict(units=units,chromosome='chr22',gather_id='all',manifest_units=['block0','block1'])))
    return Namespace(carriers=carriers,samples=samples,receipts=receipts,metadata=str(meta),outdir=str(tmp/'gather'))


def rewrite(a,index,kind,mutate):
    p=Path(getattr(a,kind)[index]);rr=rows(p);fields=list(rr[0]);mutate(rr);write_table(p,fields,rr)
    rp=Path(a.receipts[index]);r=json.loads(rp.read_text());r['outputs'][p.name]=identity(p);rp.write_text(json.dumps(r))


def test_gather_associations_distinct_genes_stars_zeros(tmp_path,capsys):
    a=make(tmp_path);before=[identity(p) for p in a.carriers+a.samples+a.receipts];capsys.readouterr();run(a)
    out=Path(a.outdir);r=json.loads((out/'receipt.json').read_text())
    assert r['carrier_annotation_records']==6 and r['samples_in_source']==3 and r['reconciliation']['passed']
    assert r['by_allele_class']['sequence']['carrier_annotation_records']==4
    assert r['by_allele_class']['sequence']['distinct_allele_sample_records']==2
    assert r['by_allele_class']['sequence']['samples_with_carriers']==1 # not sum of block counts
    assert r['by_allele_class']['spanning_deletion']['distinct_allele_sample_records']==2
    assert {x['sample'] for x in rows(out/'sample_burden.tsv')}=={'S1','S2','S3'}
    assert all(x['carrier_variant_records']=='0' for x in rows(out/'sample_burden.tsv') if x['sample']!='S1')
    assert {x['Feature'] for x in rows(out/'carriers.qc.tsv.gz')}=={'TX1','TX_OTHER'}
    assert any(x['tier']=='untiered' and x['allele_class']=='spanning_deletion' for x in rows(out/'sample_gene_burden.tsv'))
    assert all(x['carrier_samples']=='1' for x in rows(out/'gene_burden.tsv'))
    assert [identity(p) for p in a.carriers+a.samples+a.receipts]==before
    assert 'S1' not in capsys.readouterr().out


def test_empty_all_samples_and_zero_strata(tmp_path):
    a=make(tmp_path,True);run(a);out=Path(a.outdir)
    assert len(rows(out/'samples.tsv'))==3 and len(rows(out/'sample_distinct_alleles.tsv'))==3
    assert rows(out/'sample_gene_burden.tsv')==[]
    assert all(x['carrier_variant_records']=='0' for x in rows(out/'sample_burden.tsv'))


@pytest.mark.parametrize('case',['duplicate','overlap_other_sample','samples','hash','policy','count','chromosome','within','failed'])
def test_fail_closed(tmp_path,case):
    a=make(tmp_path)
    if case in ('duplicate','overlap_other_sample'):
        def mutate(rr):
            for r in rr:
                r['POS']=str(int(r['POS'])-100);r['CHROM']='22'
                if case=='overlap_other_sample':r['sample']='S2'
        rewrite(a,1,'carriers',mutate)
    elif case=='within':rewrite(a,1,'carriers',lambda rr:rr.__setitem__(2,dict(rr[0])))
    elif case=='samples':rewrite(a,1,'samples',lambda rr:rr.reverse())
    elif case=='hash':Path(a.carriers[0]).write_bytes(b'bad')
    else:
        p=Path(a.receipts[0]);r=json.loads(p.read_text())
        if case=='policy':r['policy']['min_gq']=21
        if case=='count':r['pass_rows']+=1
        if case=='chromosome':r['sources']['chromosome']='chr21'
        if case=='failed':r['status']='failed'
        p.write_text(json.dumps(r))
    with pytest.raises(ValidationError):run(a)
    assert not (Path(a.outdir)/'sample_burden.tsv').exists()
    r=json.loads((Path(a.outdir)/'receipt.json').read_text());assert r['status']=='failed'
    if case.startswith('duplicate') or case=='overlap_other_sample':assert r['duplicate_checks']['cross_block_exact_alleles']==2


def test_same_gene_cross_type_distinct_union(tmp_path):
    a=make(tmp_path)
    for i in range(2):rewrite(a,i,'carriers',lambda rr:rr[1].update(Gene=rr[0]['Gene']))
    run(a);rr=rows(Path(a.outdir)/'sample_gene_distinct_alleles.tsv')
    assert next(r for r in rr if r['allele_class']=='sequence')['distinct_allele_sample_records']=='2'


@pytest.mark.skipif(not shutil.which('nextflow'),reason='Nextflow required')
def test_nextflow_subset_resume_and_full_gather(tmp_path):
    a=make(tmp_path);manifest=tmp_path/'manifest.tsv'
    manifest.write_text('unit_id\tchromosome\tcarriers\tsamples\tsource_receipt\n'+''.join(f'block{i}\tchr22\t{a.carriers[i]}\t{a.samples[i]}\t{a.receipts[i]}\n' for i in range(2)))
    config=tmp_path/'python.config';config.write_text('process.beforeScript = "export PATH='+str(Path(sys.executable).parent)+':\\$PATH"\n')
    cmd=['nextflow','-C',str(ROOT/'gather_post_qc.config')+','+str(config),'run',str(ROOT/'gather_post_qc.nf'),'-ansi-log','false','--gather_manifest',str(manifest),'--outdir',str(tmp_path/'published'),'--gather_memory','1 GB']
    for label,units,gid,status in [('one','block0','pilot','COMPLETED'),('resume','block0','pilot','CACHED'),('all','all','all','COMPLETED')]:
        p=subprocess.run(cmd+['--select_units',units,'--gather_id',gid,'-with-trace',str(tmp_path/(label+'.trace'))]+(['-resume'] if label!='one' else []),cwd=tmp_path,capture_output=True,text=True,timeout=120)
        assert p.returncode==0,p.stdout+p.stderr
        assert [x['status'] for x in rows(tmp_path/(label+'.trace'))]==[status]
    shutil.rmtree(tmp_path/'work');out=tmp_path/'published/post-qc-gather/chr22/all'
    assert len(rows(out/'carriers.qc.tsv.gz'))==6 and not (out/'carriers.qc.tsv.gz').is_symlink()
    assert json.loads((out/'receipt.json').read_text())['full_manifest_selected']


def test_all_8877_samples_retained(tmp_path):
    a=make(tmp_path,True)
    for i in range(2):
        sp=Path(a.samples[i]);write_table(sp,['sample'],(dict(sample=f'S{n}') for n in range(8877)))
        rp=Path(a.receipts[i]);r=json.loads(rp.read_text());r['samples_in_source']=8877;r['outputs']['samples.tsv']=identity(sp);rp.write_text(json.dumps(r))
    run(a);out=Path(a.outdir)
    assert len(rows(out/'samples.tsv'))==len(rows(out/'sample_distinct_alleles.tsv'))==8877
    assert len({r['sample'] for r in rows(out/'sample_burden.tsv')})==8877


def test_output_cannot_overwrite_input(tmp_path):
    a=make(tmp_path);a.outdir=str(Path(a.carriers[0]).parent)
    before=identity(a.carriers[0])
    with pytest.raises(ValidationError):run(a)
    assert identity(a.carriers[0])==before
