import csv
import gzip
import json
from pathlib import Path
import shutil
import subprocess
import sys
from argparse import Namespace
import duckdb
import pytest
import pysam

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from extract_filtered_carriers import extract,ValidationError
from extract_exact_carriers import sha
from test_exact_carriers import star_fixture


def fixture(tmp,index='tbi'):
    old=star_fixture(tmp,index)
    con=duckdb.connect()
    con.execute('CREATE TABLE c(CHROM VARCHAR,POS BIGINT,REF VARCHAR,ALT VARCHAR,Gene VARCHAR,Feature VARCHAR,SYMBOL VARCHAR,Consequence VARCHAR,LoF VARCHAR,tier VARCHAR,allele_class VARCHAR,pcf_retained BOOLEAN)')
    for kind,rows in {
        'missense':[(10,'A','G','GENE_M','TX_M','miss_t1'),(20,'C','T','GENE_M','TX_M','miss_t2')],
        'lof_hc':[(10,'A','G','GENE1','TX1','lof_t1'),(13,'A','C','GENE1','TX1',None),(30,'G','A','GENE2','TX2',None),(60,'AT','*','GENE1','TX1','lof_t2'),(61,'C','*','GENE4','TX4',None)]
    }.items():
        con.execute('DELETE FROM c')
        for pos,ref,alt,gene,tx,tier in rows:
            con.execute('INSERT INTO c VALUES (?,?,?,?,?,?,?,?,?,?,?,?)',['22',pos,ref,alt,gene,tx,'SYMBOL','missense_variant' if kind=='missense' else 'frameshift_variant',None if kind=='missense' else 'HC',tier,'spanning_deletion' if alt=='*' else 'sequence',True])
        con.execute('COPY c TO ? (FORMAT PARQUET)',[str(tmp/(kind+'.parquet'))])
    metadata=tmp/'filtered-meta.json'
    metadata.write_text(json.dumps(dict(unit_id='block12',chromosome='chr22',vcf=old.vcf,index=old.index,missense=str(tmp/'missense.parquet'),lof_hc=str(tmp/'lof_hc.parquet'))))
    return Namespace(metadata=str(metadata),missense=str(tmp/'missense.parquet'),lof_hc=str(tmp/'lof_hc.parquet'),vcf=old.vcf,index=old.index,outdir=str(tmp/'filtered-out'),expected_missense=2,expected_hc=5)


def rows(path):
    opener=gzip.open if str(path).endswith('.gz') else open
    with opener(path,'rt') as f:return list(csv.DictReader(f,delimiter='\t'))


@pytest.mark.parametrize('index',['tbi','csi'])
def test_filtered_selection_overlap_tiers_stars_and_zero_samples(tmp_path,index,capsys):
    a=fixture(tmp_path,index);before={p:sha(p) for p in [a.missense,a.lof_hc,a.vcf,a.index]}
    extract(a);out=Path(a.outdir);data=rows(out/'carriers.tsv.gz')
    assert len(data)==9
    assert not any(r['POS'] in ('40','59') for r in data)
    assert {r['sample'] for r in data}=={'S1','S2','S3'}
    shared=[r for r in data if r['POS']=='10' and r['sample']=='S2']
    assert {(r['candidate_type'],r['tier'],r['Gene'],r['Feature']) for r in shared}=={('missense','miss_t1','GENE_M','TX_M'),('lof_hc','lof_t1','GENE1','TX1')}
    assert all(r['GT']=='1|1' and r['alt_dosage']=='2' for r in shared)
    assert shared[0]['GQ']=='60' and shared[0]['DP']=='25' and shared[0]['AD']=='0,25' and shared[0]['FT']=='PASS'
    assert any(r['GT']=='1/.' and r['site_FILTER']=='q10' for r in data)
    stars=[r for r in data if r['allele_class']=='spanning_deletion']
    assert len(stars)==3 and all(r['tier']=='lof_t2' and r['ALT']=='*' for r in stars)
    r=json.loads((out/'receipt.json').read_text())
    assert r['candidate_annotation_records']==7 and r['distinct_candidate_alleles']==6 and r['overlapping_type_alleles']==1
    assert r['carrier_annotation_records']==9
    assert r['distinct_allele_audit_by_class']['sequence']['carrier_variant_sample_records']==4
    assert r['distinct_allele_audit_by_class']['spanning_deletion']['carrier_variant_sample_records']==3
    hc=r['by_candidate_type']['lof_hc'];seq=hc['by_allele_class']['sequence'];star=hc['by_allele_class']['spanning_deletion']
    assert hc['candidate_records']==5 and hc['matched_candidates']==4 and hc['unmatched_candidates']==1
    assert seq['by_tier']['untiered']['candidate_records']==2
    assert seq['by_tier']['untiered']['matched_candidates_without_carriers']==2
    assert star['by_tier']['untiered']['unmatched_candidates']==1
    assert star['by_tier']['lof_t2']['carrier_records']==3
    sample=rows(out/'sample_burden.tsv')
    assert {r['sample'] for r in sample}=={'S1','S2','S3','S4'}
    assert all(r['carrier_variant_records']=='0' for r in sample if r['sample']=='S4')
    assert next(r for r in sample if r['sample']=='S2' and r['candidate_type']=='missense' and r['tier']=='miss_t1')['carrier_variant_records']=='1'
    assert rows(out/'sample_distinct_alleles.tsv')[-1]==dict(sample='S4',sequence_distinct_variant_count='0',spanning_deletion_distinct_record_count='0')
    assert all(sha(p)==h for p,h in before.items())
    assert 'S1' not in capsys.readouterr().out


def mutate(path,sql):
    c=duckdb.connect();c.execute('CREATE TABLE c AS SELECT * FROM read_parquet(?)',[str(path)]);c.execute(sql);Path(path).unlink();c.execute('COPY c TO ? (FORMAT PARQUET)',[str(path)])


@pytest.mark.parametrize('case',['empty','not_filtered','duplicate','wrong_tier','wrong_class','missing_gene','wrong_expected','bad_index','multiallelic'])
def test_empty_and_invalid_inputs(tmp_path,case):
    a=fixture(tmp_path)
    if case=='empty':
        for p in [a.missense,a.lof_hc]:mutate(p,'DELETE FROM c')
        a.expected_missense=a.expected_hc=0;extract(a)
        out=Path(a.outdir);r=json.loads((out/'receipt.json').read_text())
        assert r['candidate_annotation_records']==r['carrier_annotation_records']==0
        assert len(rows(out/'samples.tsv'))==4 and len(rows(out/'sample_distinct_alleles.tsv'))==4
        assert len(rows(out/'sample_burden.tsv'))==16
        assert all(x['carrier_variant_records']=='0' for x in rows(out/'sample_burden.tsv'))
        return
    mutations={'not_filtered':'UPDATE c SET pcf_retained=false','duplicate':'INSERT INTO c SELECT * FROM c LIMIT 1',
               'wrong_tier':"UPDATE c SET tier='lof_t1'",'wrong_class':"UPDATE c SET allele_class='spanning_deletion'",'missing_gene':"UPDATE c SET Gene=NULL"}
    if case in mutations:mutate(a.missense,mutations[case])
    if case=='wrong_expected':a.expected_missense=99
    if case=='bad_index':Path(a.index).unlink()
    if case=='multiallelic':
        source=tmp_path/'original.vcf';source.write_text(source.read_text().replace('chr22\t60\t.\tAT\t*\t','chr22\t60\t.\tAT\t*,ATT\t'))
        pysam.tabix_compress(str(source),a.vcf,force=True);pysam.tabix_index(a.vcf,preset='vcf',force=True)
    with pytest.raises(Exception):extract(a)
    assert json.loads((Path(a.outdir)/'receipt.json').read_text())['status']=='failed'
    assert not (Path(a.outdir)/'carriers.tsv.gz').exists()


@pytest.mark.skipif(not shutil.which('nextflow'),reason='Nextflow required')
@pytest.mark.parametrize('with_psam',[False,True])
def test_nextflow_wiring_subset_resume(tmp_path,with_psam):
    a=fixture(tmp_path/'a');b=fixture(tmp_path/'b',index='csi')
    manifest=tmp_path/'manifest.tsv';manifest.write_text('unit_id\tchromosome\tmissense\tlof_hc\tvcf\tindex\n'+''.join(
        f'{unit}\tchr22\t{x.missense}\t{x.lof_hc}\t{x.vcf}\t{x.index}\n' for unit,x in [('block12',a),('block19',b)]))
    config=tmp_path/'python.config';config.write_text('process.beforeScript = "export PATH='+str(Path(sys.executable).parent)+':\\$PATH"\n')
    cmd=['nextflow','-C',str(ROOT/'filtered_carriers.config')+','+str(config),'run',str(ROOT/'filtered_carriers.nf'),'-ansi-log','false',
         '--carrier_manifest',str(manifest),'--outdir',str(tmp_path/'published'),'--carrier_memory','1 GB','--expected_missense','2','--expected_hc','5']
    if with_psam:
        metadata=tmp_path/'samples.psam'
        metadata.write_text('#IID\tSEX\tparticipant_id\tfrequency_representative\tunrelated\n'+''.join(f'S{i}\t0\tP{i}\t1\t1\n' for i in range(1,5)))
        cmd+=['--psam',str(metadata)]
    for name,units,expected in [('one','block12',{'FILTERED_CARRIERS (block12)':'COMPLETED'}),
                               ('two','all',{'FILTERED_CARRIERS (block12)':'CACHED','FILTERED_CARRIERS (block19)':'COMPLETED'})]:
        r=subprocess.run(cmd+['--select_units',units,'-with-trace',str(tmp_path/(name+'.trace'))]+(['-resume'] if name=='two' else []),cwd=tmp_path,capture_output=True,text=True,timeout=120)
        assert r.returncode==0,r.stdout+r.stderr
        assert {x['name']:x['status'] for x in rows(tmp_path/(name+'.trace'))}==expected
    if with_psam:
        assert [r['#IID'] for r in rows(tmp_path/'published/filtered-carriers/block12/sample_metadata.tsv')]==['S1','S2','S3','S4']
    shutil.rmtree(tmp_path/'work')
    p=tmp_path/'published/filtered-carriers/block12/carriers.tsv.gz'
    assert p.is_file() and not p.is_symlink() and len(rows(p))==9
