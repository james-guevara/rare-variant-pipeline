import csv
import gzip
import json
from pathlib import Path
import shutil
import subprocess
import sys
from argparse import Namespace
import pytest

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from qc_filtered_carriers import evaluate,run,identity,write_table,ValidationError


def base_row(**changes):
    return dict(CHROM='chr22',POS='10',REF='A',ALT='G',sample='S1',candidate_type='missense',tier='miss_t1',
                Gene='G1',Feature='TX1',SYMBOL='SYM',allele_class='sequence',GT='0/1',GQ='20',DP='10',
                AD='3,1',FT='.',alt_dosage='1',site_FILTER='PASS',**{}) | changes


@pytest.mark.parametrize('changes,expected',[
    ({},True),({'AD':'1,3'},True),({'AD':'3.0001,1'},False),({'AD':'1,3.0001'},False),
    ({'GT':'0|1'},True),({'GT':'1|0'},True),({'GT':'1/0'},True),
    ({'GT':'1/1','AD':'1,9'},True),({'GT':'1|1','AD':'1,8.999'},False),
    ({'GT':'1','AD':'1,9'},True),({'GT':'1','AD':'1,8.999'},False),
    ({'GT':'1/.','AD':'0,10'},False),({'GT':'./1'},False),({'GT':'0/0'},False),({'GT':'0/2'},False),
    ({'GT':'1/1/1'},False),({'GT':'.'},False),({'GQ':'19.999'},False),({'DP':'9.999'},False),
    ({'site_FILTER':'.'},False),({'site_FILTER':'PASS;q10'},False),({'FT':'FAIL'},True)
])
def test_gt_and_threshold_boundaries(changes,expected):
    ab,fail=evaluate(base_row(**changes));assert (not fail)==expected


@pytest.mark.parametrize('field',['GQ','DP'])
@pytest.mark.parametrize('value',['.','','garbage','NaN','inf','-inf','-1','1_000','1e999'])
def test_bad_quality_fails_closed(field,value):
    ab,fail=evaluate(base_row(**{field:value}));assert field.lower()+'_invalid' in fail


@pytest.mark.parametrize('ad',['.','','1','1,2,3','-1,10','10,-1','NaN,10','1,inf','0,0','x,1','1_0,10','1e308,1e308'])
def test_bad_ad_fails_closed(ad):
    ab,fail=evaluate(base_row(AD=ad));assert ab is None and 'ad_invalid' in fail and 'ab_out_of_range' not in fail


def fixture(tmp,empty=False):
    tmp.mkdir(exist_ok=True)
    raw=[base_row(),base_row(candidate_type='lof_hc',tier='lof_t1',Gene='G_OTHER',Feature='TX_OTHER'),
         base_row(POS='20',ALT='*',allele_class='spanning_deletion',candidate_type='lof_hc',tier='untiered',GT='1',AD='1,9'),
         base_row(POS='30',sample='S2',candidate_type='lof_hc',tier='lof_t2',GQ='99',DP='40',AD='39,1'),
         base_row(POS='40',ALT='*',allele_class='spanning_deletion',candidate_type='lof_hc',tier='lof_t1',GQ='NaN',DP='-1',AD='-1,10',site_FILTER='q10'),
         base_row(POS='50',GT='1/.',AD='0,10')]
    if empty:raw=[]
    carriers=tmp/'carriers.tsv.gz';write_table(carriers,list(base_row()),raw)
    samples=tmp/'samples.tsv';write_table(samples,['sample'],[{'sample':s} for s in ['S1','S2','S3']])
    source=tmp/'extraction-receipt.json';source.write_text(json.dumps(dict(status='passed',schema_version=1,unit_id='block12',samples_in_source=3,carrier_annotation_records=len(raw),outputs={'carriers.tsv.gz':identity(carriers),'samples.tsv':identity(samples)})))
    metadata=tmp/'metadata.json';metadata.write_text(json.dumps(dict(unit_id='block12',chromosome='chr22',carriers=str(carriers),samples=str(samples),source_receipt=str(source))))
    return Namespace(carriers=str(carriers),samples=str(samples),source_receipt=str(source),metadata=str(metadata),outdir=str(tmp/'out'))


def rows(path):
    with (gzip.open(path,'rt') if str(path).endswith('.gz') else open(path)) as f:return list(csv.DictReader(f,delimiter='\t'))


def test_summaries_overlap_stars_untiered_zero_samples_and_reconciliation(tmp_path,capsys):
    a=fixture(tmp_path);before={p:identity(p) for p in [a.carriers,a.samples,a.source_receipt]};run(a)
    out=Path(a.outdir);r=json.loads((out/'receipt.json').read_text())
    assert (r['input_rows'],r['pass_rows'],r['fail_rows'])==(6,3,3)
    assert r['independent_failures']==dict(site_filter=1,gq_invalid=1,gq_below_min=0,dp_invalid=1,dp_below_min=0,gt_unsupported=1,ad_invalid=1,ab_out_of_range=1)
    assert r['failure_combinations']=={'PASS':3,'ab_out_of_range':1,'gt_unsupported':1,'site_filter|gq_invalid|dp_invalid|ad_invalid':1}
    assert sum(r['failure_combinations'].values())==6 and r['reconciliation']['passed']
    qc=rows(out/'carriers.qc.tsv.gz')
    assert {(x['candidate_type'],x['tier'],x['Gene'],x['Feature']) for x in qc if x['POS']=='10'}=={('missense','miss_t1','G1','TX1'),('lof_hc','lof_t1','G_OTHER','TX_OTHER')}
    assert next(x for x in qc if x['ALT']=='*')['tier']=='untiered'
    distinct=rows(out/'sample_distinct_alleles.tsv')
    assert distinct==[dict(sample=s,sequence_distinct_variant_count='1' if s=='S1' else '0',spanning_deletion_distinct_record_count='1' if s=='S1' else '0') for s in ['S1','S2','S3']]
    summaries=rows(out/'sample_burden.tsv')
    assert {x['sample'] for x in summaries}=={'S1','S2','S3'}
    assert all(x['carrier_variant_records']=='0' for x in summaries if x['sample'] in ('S2','S3'))
    assert r['by_candidate_type']['lof_hc']['by_allele_class']['spanning_deletion']['by_tier']['untiered']['pass_rows']==1
    assert all(identity(p)==v for p,v in before.items())
    assert 'S1' not in capsys.readouterr().out


def test_empty_carriers_retain_sample_universe(tmp_path):
    a=fixture(tmp_path,empty=True);run(a);out=Path(a.outdir)
    assert rows(out/'carriers.qc.tsv.gz')==[]
    assert len(rows(out/'sample_burden.tsv'))==12
    assert len(rows(out/'samples.tsv'))==len(rows(out/'sample_distinct_alleles.tsv'))==3
    assert json.loads((out/'receipt.json').read_text())['reconciliation']['passed']


@pytest.mark.parametrize('case',['hash','duplicate','conflicting_overlap','unknown_sample','overwrite'])
def test_invalid_provenance_or_associations_fail(tmp_path,case):
    a=fixture(tmp_path)
    if case=='overwrite':a.outdir=str(tmp_path)
    else:
        data=rows(Path(a.carriers));fields=list(data[0])
        if case=='duplicate':data.append(dict(data[0]))
        elif case=='conflicting_overlap':data[1]['GQ']='99'
        elif case=='unknown_sample':data[0]['sample']='UNKNOWN'
        else:data[0]['GQ']='21'
        write_table(a.carriers,fields,data)
        if case!='hash':
            p=Path(a.source_receipt);r=json.loads(p.read_text());r['outputs']['carriers.tsv.gz']=identity(a.carriers);r['carrier_annotation_records']=len(data);p.write_text(json.dumps(r))
    with pytest.raises(ValidationError):run(a)
    assert not (Path(a.outdir)/'carriers.qc.tsv.gz').exists()


def test_existing_qc_genotype_agrees_for_supported_well_formed_values(tmp_path):
    import duckdb
    data=[base_row(**x) for x in [{},{'AD':'1,3'},{'AD':'9,1'},{'GT':'1','AD':'1,9'},
                                  {'GT':'1/1','AD':'0,10'},{'GT':'1|0'},{'GQ':'19'},{'DP':'9'},{'GT':'1/.'}]]
    con=duckdb.connect();con.execute('CREATE TABLE x(rowid INTEGER,GT VARCHAR,GQ VARCHAR,DP VARCHAR,AD VARCHAR)')
    con.executemany('INSERT INTO x VALUES (?,?,?,?,?)',[(i,r['GT'],r['GQ'],r['DP'],r['AD']) for i,r in enumerate(data)])
    inp=tmp_path/'in.parquet';out=tmp_path/'out.parquet';con.execute('COPY x TO ? (FORMAT PARQUET)',[str(inp)])
    config=tmp_path/'config.json';config.write_text(json.dumps(dict(cohorts={'synthetic':{}},output_base=str(tmp_path),qc=dict(min_gq=20,min_dp=10,het_ab_min=0.25,het_ab_max=0.75,hom_ab_min=0.9))))
    p=subprocess.run([sys.executable,str(ROOT/'scripts/postprocess/qc_genotype.py'),'--cohort','synthetic','--chrom','chr22','--resources',str(config),'--input',str(inp),'--output',str(out)],capture_output=True,text=True)
    assert p.returncode==0,p.stderr
    assert [x[0] for x in con.execute('SELECT rowid FROM read_parquet(?) ORDER BY rowid',[str(out)]).fetchall()]==[i for i,r in enumerate(data) if not evaluate(r)[1]]


@pytest.mark.skipif(not shutil.which('nextflow'),reason='Nextflow required')
def test_nextflow_subset_resume_durable_outputs(tmp_path):
    a=fixture(tmp_path/'a');b=fixture(tmp_path/'b')
    p=Path(b.source_receipt);data=json.loads(p.read_text());data['unit_id']='block19';p.write_text(json.dumps(data))
    manifest=tmp_path/'manifest.tsv';manifest.write_text('unit_id\tchromosome\tcarriers\tsamples\tsource_receipt\n'+''.join(
        f'{unit}\tchr22\t{x.carriers}\t{x.samples}\t{x.source_receipt}\n' for unit,x in [('block12',a),('block19',b)]))
    config=tmp_path/'python.config';config.write_text('process.beforeScript = "export PATH='+str(Path(sys.executable).parent)+':\\$PATH"\n')
    cmd=['nextflow','-C',str(ROOT/'post_carrier_qc.config')+','+str(config),'run',str(ROOT/'post_carrier_qc.nf'),'-ansi-log','false','--qc_manifest',str(manifest),'--outdir',str(tmp_path/'published'),'--qc_memory','1 GB']
    for label,units,expected in [('one','block12',{'POST_CARRIER_QC (block12)':'COMPLETED'}),('two','all',{'POST_CARRIER_QC (block12)':'CACHED','POST_CARRIER_QC (block19)':'COMPLETED'})]:
        r=subprocess.run(cmd+['--select_units',units,'-with-trace',str(tmp_path/(label+'.trace'))]+(['-resume'] if label=='two' else []),cwd=tmp_path,capture_output=True,text=True,timeout=120)
        assert r.returncode==0,r.stdout+r.stderr
        assert {x['name']:x['status'] for x in rows(tmp_path/(label+'.trace'))}==expected
    shutil.rmtree(tmp_path/'work');p=tmp_path/'published/post-carrier-qc/block12/carriers.qc.tsv.gz'
    assert p.is_file() and not p.is_symlink() and len(rows(p))==3


def test_reads_actual_filtered_extractor_products(tmp_path):
    from test_filtered_carriers import fixture as carrier_fixture
    from extract_filtered_carriers import extract
    raw=carrier_fixture(tmp_path/'raw');extract(raw)
    directory=Path(raw.outdir)
    metadata=tmp_path/'metadata.json';metadata.write_text(json.dumps(dict(unit_id='block12',chromosome='chr22')))
    a=Namespace(carriers=str(directory/'carriers.tsv.gz'),samples=str(directory/'samples.tsv'),source_receipt=str(directory/'receipt.json'),metadata=str(metadata),outdir=str(tmp_path/'qc'))
    run(a);receipt=json.loads((Path(a.outdir)/'receipt.json').read_text())
    assert receipt['input_rows']==9 and receipt['pass_rows']==4 and receipt['fail_rows']==5
    assert receipt['samples_in_source']==4
    assert receipt['post_qc_distinct_alleles_by_class']['sequence']['distinct_variant_sample_records']==2
    assert receipt['post_qc_distinct_alleles_by_class']['spanning_deletion']['distinct_variant_sample_records']==0
