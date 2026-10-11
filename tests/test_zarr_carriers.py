"""Backend parity: the same records must give the same policy products."""
import copy
import gzip
import json
from pathlib import Path
import numpy as np
import pysam
import pytest
import zarr
from test_filtered_carriers import fixture,extract,rows
from test_carrier_frequencies import psam
from zarr_frequencies import ZarrFrequencyCounter
from carrier_frequencies import FrequencyCounter
from types import SimpleNamespace


def to_zarr(vcf,path):
    with pysam.VariantFile(vcf) as f:
        records=list(f);samples=list(f.header.samples);contigs=list(f.header.contigs);filters=list(f.header.filters)
    g=zarr.open_group(str(path),mode='w')
    def arr(name,data,dtype=None):
        a=np.asarray(data,dtype=dtype);g.create_array(name,data=a,chunks=(min(2,len(a)),)+a.shape[1:])
    arr('sample_id',samples);arr('contig_id',contigs);arr('filter_id',filters)
    arr('variant_position',[r.pos for r in records],'i4');arr('variant_contig',[contigs.index(r.contig) for r in records],'i4')
    arr('variant_allele',[[r.ref,*r.alts] for r in records])
    arr('variant_filter',[[f in r.filter for f in filters] for r in records],bool)
    gt=np.full((len(records),len(samples),2),-2,dtype='i1')
    for j,r in enumerate(records):
        for i,s in enumerate(samples):
            v=r.samples[s].get('GT') or ();gt[j,i,:len(v)]=[-1 if x is None else x for x in v]
    arr('call_genotype',gt);arr('call_genotype_mask',gt<0)
    arr('call_genotype_phased',[[r.samples[s].phased for s in samples] for r in records],bool)
    for name in ['GQ','DP']:
        arr('call_'+name,[[r.samples[s].get(name) if r.samples[s].get(name) is not None else -1 for s in samples] for r in records],'i4')
    arr('call_AD',[[[-1 if x is None else x for x in (r.samples[s].get('AD') or (-1,))] + [-2]*(2-len(r.samples[s].get('AD') or (-1,))) for s in samples] for r in records],'i4')
    arr('call_FT',[[r.samples[s].get('FT') or '.' for s in samples] for r in records])
    return g


def test_products_and_frequencies_match_vcf(tmp_path):
    a=fixture(tmp_path);a.psam=str(psam(tmp_path/'samples.psam'));a.compute_frequencies=True
    extract(a);old=Path(a.outdir);store=tmp_path/'source.zarr';to_zarr(a.vcf,store)
    b=copy.copy(a);b.zarr=str(store);b.outdir=str(tmp_path/'zarr-out');extract(b);new=Path(b.outdir)
    for f in old.iterdir():
        if f.name=='receipt.json':continue
        expected=rows(f);actual=rows(new/f.name)
        if f.name=='carriers.tsv.gz':actual=[{k:v for k,v in r.items() if not k.startswith('source_')} for r in actual]
        assert expected==actual,f.name
    ra=json.loads((old/'receipt.json').read_text());rb=json.loads((new/'receipt.json').read_text())
    assert ra['frequencies']==rb['frequencies']
    assert rb['source_backend']=='zarr'
    assert rb['genotype_chunk_reads']<=4


def test_vector_counts_match_scalar_randomized():
    rng=np.random.default_rng(25);samples=['s'+str(i) for i in range(80)]
    meta=[{'#IID':s,'SEX':'0','participant_id':s,'frequency_representative':str(i%3!=0 and 1 or 0),'unrelated':str(i%2)} for i,s in enumerate(samples)]
    a=FrequencyCounter(meta,'22');b=ZarrFrequencyCounter(meta,'22',samples)
    for _ in range(40):
        gt=rng.choice([-1,0,1],size=(80,2));gt[::11,1]=-2
        calls={s:{'GT':tuple(None if x==-1 else int(x) for x in row if x!=-2)} for s,row in zip(samples,gt)}
        r=SimpleNamespace(gt=gt,samples=calls,pos=10)
        assert a.count(r,'sequence')==b.count(r,'sequence')
    assert a.receipt()==b.receipt()


def test_empty_candidates_keep_all_samples(tmp_path):
    import duckdb
    a=fixture(tmp_path);store=tmp_path/'source.zarr';to_zarr(a.vcf,store)
    con=duckdb.connect()
    for name in ['missense','lof_hc']:
        f=Path(getattr(a,name));con.execute('CREATE OR REPLACE TABLE empty AS SELECT * FROM read_parquet(?) WHERE FALSE',[str(f)]);f.unlink();con.execute('COPY empty TO ? (FORMAT PARQUET)',[str(f)])
    a.expected_missense=0;a.expected_hc=0;a.zarr=str(store);a.psam=str(psam(tmp_path/'samples.psam'));a.compute_frequencies=True
    extract(a)
    assert len(rows(Path(a.outdir)/'samples.tsv'))==4
    assert rows(Path(a.outdir)/'carriers.tsv.gz')==[]
    assert rows(Path(a.outdir)/'variant_frequencies.tsv')==[]


def test_sites_preserve_source_frequencies_without_genotypes(tmp_path):
    from zarr_sites import export
    a=fixture(tmp_path);store=tmp_path/'sites.zarr';g=to_zarr(a.vcf,store)
    n=g['variant_position'].shape[0]
    g.create_array('variant_AC',data=np.full((n,1),5,dtype='i4'))
    g.create_array('variant_AN',data=np.full(n,1000,dtype='i4'))
    g.create_array('variant_AF',data=np.full((n,1),0.4,dtype='f4'))
    for name in list(g.array_keys()):
        if name.startswith('call_'):del g[name]
    out=tmp_path/'sites.vcf.gz';export(store,'chr22',out,tmp_path/'sites.json')
    with pysam.VariantFile(out) as f:
        assert not list(f.header.samples)
        rs=list(f);assert len(rs)==n
        assert all(r.info['AC']==(5,) and r.info['AN']==1000 for r in rs)
        assert all(abs(r.info['AF'][0]-0.4)<1e-6 for r in rs)
    assert json.loads((tmp_path/'sites.json').read_text())['genotype_reads']==0


def test_nextflow_zarr_backend_and_resume(tmp_path):
    import subprocess,sys,shutil
    if not shutil.which('nextflow'):pytest.skip('Nextflow required')
    root=Path(__file__).resolve().parents[1]
    a=fixture(tmp_path);store=tmp_path/'source.zarr';to_zarr(a.vcf,store)
    sample=psam(tmp_path/'samples.psam')
    manifest=tmp_path/'zarr.tsv';manifest.write_text('unit_id\tchromosome\tmissense\tlof_hc\tzarr\n'+f'block12\tchr22\t{a.missense}\t{a.lof_hc}\t{store}\n')
    config=tmp_path/'python.config';config.write_text('process.beforeScript = "export PATH='+str(Path(sys.executable).parent)+':\\$PATH"\n')
    cmd=['nextflow','-C',str(root/'zarr_carriers.config')+','+str(config),'run',str(root/'zarr_carriers.nf'),'-ansi-log','false','--carrier_manifest',str(manifest),'--psam',str(sample),'--carrier_memory','1 GB','--outdir',str(tmp_path/'published')]
    for name,extra in [('first',[]),('resume',['-resume'])]:
        trace=tmp_path/(name+'.trace');r=subprocess.run(cmd+extra+['-with-trace',str(trace)],cwd=tmp_path,capture_output=True,text=True,timeout=120)
        assert r.returncode==0,r.stdout+r.stderr
        assert rows(trace)[0]['status']==('CACHED' if extra else 'COMPLETED')
    receipt=json.loads((tmp_path/'published/filtered-carriers/block12/receipt.json').read_text())
    assert receipt['source_backend']=='zarr' and receipt['frequencies']['status']=='passed'


def test_pedigree_selection_is_order_independent_and_retains_source_order(tmp_path):
    from make_family_frequency_psam import build
    meta=tmp_path/'meta.tsv';meta.write_text('sample\tperson\tfamily\tsex\ns3\tp2\tf1\tMale\ns2\tp1\tf1\tFemale\ns1\tp1\tf1\tFemale\ns4\tp4\tf2\tMale\n')
    source=tmp_path/'source.psam';source.write_text('#IID\n'+'\n'.join(['s4','s2','s3','s1'])+'\n')
    out=tmp_path/'out.psam';r=build(meta,source,out,'sample','person','family','sex')
    data=rows(out);assert [x['#IID'] for x in data]==['s4','s2','s3','s1']
    assert {x['#IID'] for x in data if x['unrelated']=='1'}=={'s1','s4'}
    assert r['participants']==3 and r['families']==2


@pytest.mark.parametrize('problem',['duplicate','mask'])
def test_invalid_source_fails_closed(tmp_path,problem):
    from extract_exact_carriers import ValidationError
    a=fixture(tmp_path);store=tmp_path/'source.zarr';g=to_zarr(a.vcf,store);a.zarr=str(store)
    if problem=='mask':g['call_genotype_mask'][0,0,0]=True
    elif problem=='duplicate':
        g['variant_position'][1]=g['variant_position'][0]
        g['variant_allele'][1]=g['variant_allele'][0]
    else:
        alleles=g['variant_allele'][:];del g['variant_allele']
        extended=np.column_stack([alleles,np.full(len(alleles),'')]);extended[0,2]='T'
        g.create_array('variant_allele',data=extended)
    with pytest.raises(ValidationError):extract(a)
    assert json.loads((Path(a.outdir)/'receipt.json').read_text())['status']=='failed'
    assert not (Path(a.outdir)/'carriers.tsv.gz').exists()


def test_zarr_outputs_feed_final_rarity(tmp_path):
    from argparse import Namespace
    from filter_final_rarity import run
    a=fixture(tmp_path);a.psam=str(psam(tmp_path/'samples.psam'));a.compute_frequencies=True
    store=tmp_path/'source.zarr';to_zarr(a.vcf,store);a.zarr=str(store);extract(a)
    out=Path(a.outdir)
    b=Namespace(metadata=a.metadata,carriers=str(out/'carriers.tsv.gz'),samples=str(out/'samples.tsv'),candidates=str(out/'candidates.tsv'),frequencies=str(out/'variant_frequencies.tsv'),source_receipt=str(out/'receipt.json'),outdir=str(tmp_path/'rarity'))
    run(b)
    assert json.loads((Path(b.outdir)/'receipt.json').read_text())['status']=='passed'
    assert len(rows(Path(b.outdir)/'samples.tsv'))==4


def test_nextflow_sites_entrypoint(tmp_path):
    import subprocess,sys,shutil
    if not shutil.which('nextflow'):pytest.skip('Nextflow required')
    root=Path(__file__).resolve().parents[1];a=fixture(tmp_path);store=tmp_path/'source.zarr';g=to_zarr(a.vcf,store)
    n=g['variant_position'].shape[0];g.create_array('variant_AC',data=np.ones((n,1),dtype='i4'));g.create_array('variant_AN',data=np.full(n,1000,dtype='i4'))
    manifest=tmp_path/'sites.tsv';manifest.write_text('unit_id\tchromosome\tzarr\n'+f'b1\tchr22\t{store}\n')
    config=tmp_path/'python.config';config.write_text('process.beforeScript = "export PATH='+str(Path(sys.executable).parent)+':\\$PATH"\nprocess.memory="1 GB"\n')
    r=subprocess.run(['nextflow','-C',str(root/'zarr_sites.config')+','+str(config),'run',str(root/'zarr_sites.nf'),'-ansi-log','false','--zarr_manifest',str(manifest),'--outdir',str(tmp_path/'published')],cwd=tmp_path,capture_output=True,text=True,timeout=120)
    assert r.returncode==0,r.stdout+r.stderr
    assert json.loads((tmp_path/'published/sites/b1/receipt.json').read_text())['output_records']==n


def test_multiallelic_projection_frequencies_and_source_provenance(tmp_path):
    from zarr_carrier_source import Source
    from zarr_sites import export
    a=fixture(tmp_path);store=tmp_path/'multi.zarr';g=to_zarr(a.vcf,store)
    alleles=g['variant_allele'][:];del g['variant_allele']
    als=np.column_stack([alleles,np.full(len(alleles),'')]);als[0,2]='T';als[1,1]='C';g.create_array('variant_allele',data=als)
    g['call_genotype'][0]=np.array([[1,2],[2,2],[0,2],[-1,2]])
    g['call_genotype_mask'][0]=g['call_genotype'][0]<0
    ad=g['call_AD'][:];del g['call_AD'];new=np.full((*ad.shape[:2],3),-2,dtype='i4');new[:,:,:2]=ad;new[0]=[[2,8,10],[0,0,20],[10,0,10],[-1,-1,-1]];g.create_array('call_AD',data=new)
    metadata=[{'#IID':'S'+str(i),'participant_id':'P'+str(i),'frequency_representative':'1','unrelated':'1','SEX':'0'} for i in range(1,5)]
    selected={('22',10,'A','G'):[],('22',10,'A','T'):[]}
    with Source(store,'chr22') as source:
        records=list(source.records(selected));assert len(records)==2 and source.chunk_reads==1
        counter=ZarrFrequencyCounter(metadata,'22',source.samples)
        first=counter.count(records[0][1],'sequence')[0];second=counter.count(records[1][1],'sequence')[0]
        assert (first['cohort_ac'],first['cohort_an'])==(1,7)
        assert (second['cohort_ac'],second['cohort_an'])==(5,7)
        calls=list(source.calls(records[1][1]));r=calls[0][0]
        assert r['GT']=='0/1' and r['source_GT']=='1/2'
        assert r['AD']=='2,10' and r['source_AD']=='2,8,10'
        assert r['source_alt_index']==2 and r['source_variant_index']==0
    n=len(alleles);ac=np.ones((n,2),dtype='i4');ac[:,1]=-2;ac[0]=[1,5]
    g.create_array('variant_AC',data=ac);g.create_array('variant_AN',data=np.full(n,7,dtype='i4'))
    path=tmp_path/'sites.vcf.gz';export(store,'chr22',path,tmp_path/'sites.json')
    with pysam.VariantFile(path) as f:
        at10=[r for r in f if r.pos==10 and r.alts in [('G',),('T',)]]
        assert [r.alts for r in at10]==[('G',),('T',)]
        assert [r.info['AC'] for r in at10]==[(1,),(5,)]
        assert [r.info['ZARR_ALT_INDEX'] for r in at10]==[1,2]


def test_sex_count_does_not_hide_other_alt_heterozygosity():
    meta=[{'#IID':'s','participant_id':'p','frequency_representative':'1','unrelated':'1','SEX':'1'}]
    counter=ZarrFrequencyCounter(meta,'X',['s'],'grch38_x_only_par')
    record=SimpleNamespace(pos=5000000,source_gt=np.array([[2,3]]),samples={'s':{'GT':(0,0)}})
    values,audit=counter.count(record,'sequence')
    assert values['unrelated_an']==0
    assert audit['unrelated']['excluded_haploid_heterozygous']==1
