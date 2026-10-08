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

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from annotation_resources import identity
from filter_pre_carrier import run,sha,lit,AF


def fixture(tmp):
    tmp.mkdir(exist_ok=True)
    root=tmp/'resources';root.mkdir()
    con=duckdb.connect()
    con.execute('CREATE TABLE candidates(CHROM VARCHAR,POS BIGINT,REF VARCHAR,ALT VARCHAR,allele_class VARCHAR,Gene VARCHAR,tier VARCHAR,preserved_score DOUBLE)')
    con.executemany('INSERT INTO candidates VALUES (?,?,?,?,?,?,?,?)',[
        ('chr22',i,'A','*' if i==20 else 'G','spanning_deletion' if i==20 else 'sequence','GENE','tier',0.9) for i in range(1,21)])
    for name in ['missense','lof_hc']:
        con.execute(f'COPY candidates TO {lit(tmp/(name+".parquet"))} (FORMAT PARQUET)')
    sites=tmp/'sites.vcf.gz'
    with gzip.open(sites,'wt') as f:
        f.write('##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
        for i in range(1,21):
            info='AC=1;AN=2000;AF=0.99;MLEAF=0.99'
            if i==2:info='AC=1;AN=200;AF=0;MLEAF=0'
            if i==3:info='AC=2;AN=399'
            if i==8:info='AN=2000;AF=0'
            if i==9:info='AC=0;AN=0'
            alt='T' if i==7 else '*' if i==20 else 'G'
            f.write(f'22\t{i}\t.\tA\t{alt}\t.\tPASS\t{info}\n')
    con.execute(f'CREATE TABLE pop("#chr" VARCHAR,"pos(1-based)" VARCHAR,ref VARCHAR,alt VARCHAR,"{AF}" VARCHAR,"gnomAD4.1_joint_AF" VARCHAR)')
    con.executemany('INSERT INTO pop VALUES (?,?,?,?,?,?)',[
        ('22',str(i),'A','T' if i==7 else 'G',None if i==1 else '0.001' if i==4 else '0.002' if i==5 else '0.000999' if i==6 else '0','0.99') for i in range(1,20)])
    db=root/'chr22.parquet';con.execute(f'COPY pop TO {lit(db)} (FORMAT PARQUET)')
    regions={}
    for t,body in {
        'genomicSuperDups':'chr22\t10\t12\nchr22\t17\t18\nchr21\t18\t19\n',
        'simpleRepeat':'22\t13\t14\n22\t17\t18\n',
        'rmsk':'22\t14\t15\tx\t0\t+\tSimple_repeat\n22\t15\t16\tx\t0\t+\tLow_complexity\n22\t16\t17\tx\t0\t+\tLINE\n22\t17\t18\tx\t0\t+\tSimple_repeat\n'
    }.items():
        p=root/(t+'.bed');p.write_text(body);regions[t]=identity(p)
    image=tmp/'image.sif';image.write_text('synthetic')
    meta=dict(unit_id='block12',chromosome='chr22',popmax=identity(db),regions=regions,container=identity(image))
    metadata=tmp/'meta.json';metadata.write_text(json.dumps(meta))
    return Namespace(missense=str(tmp/'missense.parquet'),lof_hc=str(tmp/'lof_hc.parquet'),sites=str(sites),metadata=str(metadata),outdir=str(tmp/'out'),threads=1,memory='512MB')


def test_policies_both_types_exact_joins_boundaries_and_unchanged_inputs(tmp_path,capsys):
    a=fixture(tmp_path)
    before={p:sha(p) for p in [a.missense,a.lof_hc,a.sites]}
    run(a)
    assert before=={p:sha(p) for p in before}
    r=json.loads((Path(a.outdir)/'receipt.json').read_text())
    for kind in ['missense','lof_hc']:
        c=r['counts'][kind]
        assert c['input_candidates']==20
        assert c['independent']=={'cohort':{'pass_count':15,'fail_count':5},'gnomad':{'pass_count':18,'fail_count':2},'region':{'pass_count':14,'fail_count':6}}
        assert c['sequential']=={'after_cohort':15,'after_gnomad':13,'final_retained':7}
        assert c['missing_site_matches']==1 and c['missing_cohort_af']==3
        assert c['matched_missing_or_invalid_cohort_af']==2
        assert c['missing_gnomad_popmax']==3 and c['missing_gnomad_matches']==2
        assert c['overlaps']=={'genomicSuperDups':3,'simpleRepeat':2,'rmsk':3}
        assert c['overlap_union']==6
        assert c['allele_classes']['spanning_deletion']=={'input_candidates':1,'retained':1}
        rows=duckdb.sql(f'SELECT POS,preserved_score FROM read_parquet({lit(Path(a.outdir)/(kind+".filtered.parquet"))}) ORDER BY POS').fetchall()
        assert rows==[(i,0.9) for i in [1,6,10,13,17,19,20]]
    assert r['outputs']['filter_audit.parquet']['rows']==40
    assert 'GENE' not in capsys.readouterr().out


def test_genotyped_vcf_rejected_before_payload_or_hash(tmp_path,monkeypatch):
    import filter_pre_carrier as script
    a=fixture(tmp_path)
    with gzip.open(a.sites,'wt') as f:
        f.write('##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tPROTECTED_SAMPLE\nINVALID_PAYLOAD\n')
    realsha=script.sha
    def guarded(path):
        assert Path(path)!=Path(a.sites), 'Genotyped payload must not be hashed'
        return realsha(path)
    monkeypatch.setattr(script,'sha',guarded)
    with pytest.raises(ValueError,match='genotype-free'):run(a)
    assert not list(Path(a.outdir).glob('*.parquet'))


@pytest.mark.parametrize('kind',['conflicting_popmax','invalid_popmax','changed_resource','empty'])
def test_resource_failures_and_empty_candidates(tmp_path,kind):
    a=fixture(tmp_path);meta=json.loads(Path(a.metadata).read_text())
    if kind=='empty':
        for name in ['missense','lof_hc']:
            p=Path(getattr(a,name));con=duckdb.connect();con.execute(f'CREATE TABLE x AS SELECT * FROM read_parquet({lit(p)}) WHERE false');p.unlink();con.execute(f'COPY x TO {lit(p)} (FORMAT PARQUET)')
        run(a)
        r=json.loads((Path(a.outdir)/'receipt.json').read_text())
        assert all(c['input_candidates']==c['final_retained']==0 for c in r['counts'].values())
        return
    if kind=='changed_resource':Path(meta['regions']['rmsk']['path']).write_text('changed')
    else:
        p=Path(meta['popmax']['path']);con=duckdb.connect();con.execute(f'CREATE TABLE x AS SELECT * FROM read_parquet({lit(p)})')
        if kind=='conflicting_popmax':con.execute(f'INSERT INTO x SELECT * REPLACE (\'0.3\' AS "{AF}") FROM x WHERE "pos(1-based)"=\'6\'')
        else:con.execute(f'UPDATE x SET "{AF}"=\'invalid\' WHERE "pos(1-based)"=\'6\'')
        p.unlink();con.execute(f'COPY x TO {lit(p)} (FORMAT PARQUET)');meta['popmax']=identity(p);Path(a.metadata).write_text(json.dumps(meta))
    with pytest.raises(ValueError):run(a)
    assert json.loads((Path(a.outdir)/'receipt.json').read_text())['status']=='failed'
    assert not list(Path(a.outdir).glob('*.parquet'))


def test_lock_records_all_resources_and_pins_container(tmp_path,monkeypatch):
    import lock_pre_carrier_resources as lock
    a=fixture(tmp_path);meta=json.loads(Path(a.metadata).read_text());root=tmp_path/'bundle'
    db=root/'dbNSFP/5.3.1a/parquet_scores_af';db.mkdir(parents=True)
    beds=root/'problematic-regions';beds.mkdir()
    shutil.copy(meta['popmax']['path'],db/'chr22.parquet')
    for t,v in meta['regions'].items():shutil.copy(v['path'],beds/(t+'.bed'))
    monkeypatch.setattr(lock,'LOF_SIF',meta['container']['sha256'])
    result=lock.build(root,['22'],meta['container']['path'])
    assert result['popmax']['chr22']['sha256']==meta['popmax']['sha256']
    assert set(result['regions'])=={'genomicSuperDups','simpleRepeat','rmsk'}
    monkeypatch.setattr(lock,'LOF_SIF','wrong')
    with pytest.raises(ValueError):lock.build(root,['22'],meta['container']['path'])


@pytest.mark.skipif(not shutil.which('nextflow'),reason='Nextflow required')
def test_nextflow_resume_and_publication(tmp_path):
    a=fixture(tmp_path);meta=json.loads(Path(a.metadata).read_text())
    lock=dict(schema=1,stage='pre_carrier',resource_root=str(tmp_path/'resources'),popmax={'chr22':meta['popmax']},regions=meta['regions'],container=meta['container'])
    (tmp_path/'lock.json').write_text(json.dumps(lock))
    (tmp_path/'manifest.tsv').write_text(f'unit_id\tchromosome\tmissense\tlof_hc\tsites\nblock12\tchr22\t{a.missense}\t{a.lof_hc}\t{a.sites}\n')
    config=tmp_path/'python.config';config.write_text('process.beforeScript = "export PATH='+str(Path(sys.executable).parent)+':\\$PATH"\n')
    cmd=['nextflow','-C',str(ROOT/'pre_carrier.config')+','+str(config),'run',str(ROOT/'pre_carrier.nf'),'-ansi-log','false',
         '--filter_manifest',str(tmp_path/'manifest.tsv'),'--filter_resource_lock',str(tmp_path/'lock.json'),'--filter_resource_root',lock['resource_root'],
         '--select_units','block12','--outdir',str(tmp_path/'published'),'--filter_memory','1 GB']
    for n,expected in [(1,'COMPLETED'),(2,'CACHED')]:
        result=subprocess.run(cmd+['-with-trace',str(tmp_path/f'trace{n}.tsv')]+(['-resume'] if n==2 else []),cwd=tmp_path,capture_output=True,text=True,timeout=120)
        assert result.returncode==0,result.stdout+result.stderr
        rows=list(csv.DictReader((tmp_path/f'trace{n}.tsv').open(),delimiter='\t'))
        assert len(rows)==1 and rows[0]['status']==expected
    shutil.rmtree(tmp_path/'work')
    assert (tmp_path/'published/pre-carrier/block12/missense.filtered.parquet').is_file()
    assert not (tmp_path/'published/pre-carrier/block12/missense.filtered.parquet').is_symlink()

@pytest.mark.parametrize('field,value',[('REF','C'),('CHROM','21'),('ALT','T')])
def test_site_and_population_require_exact_allele(tmp_path,field,value):
    a=fixture(tmp_path)
    with gzip.open(a.sites,'rt') as f:lines=f.readlines()
    for n,line in enumerate(lines):
        if line.startswith('22\t6\t'):
            row=line.rstrip('\n').split('\t');row[{'CHROM':0,'REF':3,'ALT':4}[field]]=value
            lines[n]='\t'.join(row)+'\n'
    with gzip.open(a.sites,'wt') as f:f.writelines(lines)
    meta=json.loads(Path(a.metadata).read_text());p=Path(meta['popmax']['path']);con=duckdb.connect()
    con.execute(f'CREATE TABLE x AS SELECT * FROM read_parquet({lit(p)})')
    column={'CHROM':'#chr','REF':'ref','ALT':'alt'}[field]
    con.execute(f'UPDATE x SET "{column}"=? WHERE "pos(1-based)"=\'6\'',[value])
    p.unlink();con.execute(f'COPY x TO {lit(p)} (FORMAT PARQUET)')
    meta['popmax']=identity(p);Path(a.metadata).write_text(json.dumps(meta));run(a)
    audit=con.execute('SELECT pcf_site_matched,pcf_gnomad_matched,pcf_cohort_pass,pcf_gnomad_pass,pcf_retained FROM read_parquet(?) WHERE POS=6',
                      [str(Path(a.outdir)/'filter_audit.parquet')]).fetchall()
    assert audit==[(False,False,False,True,False)]*2


def test_duplicate_identical_population_rows_do_not_multiply_candidates(tmp_path):
    a=fixture(tmp_path);meta=json.loads(Path(a.metadata).read_text());p=Path(meta['popmax']['path'])
    con=duckdb.connect();con.execute(f'CREATE TABLE x AS SELECT * FROM read_parquet({lit(p)})')
    con.execute('INSERT INTO x SELECT * FROM x');p.unlink();con.execute(f'COPY x TO {lit(p)} (FORMAT PARQUET)')
    meta['popmax']=identity(p);Path(a.metadata).write_text(json.dumps(meta));run(a)
    receipt=json.loads((Path(a.outdir)/'receipt.json').read_text())
    assert all(v['input_candidates']==20 and v['final_retained']==7 for v in receipt['counts'].values())


def test_relaxed_screen_preserves_original_info_and_is_not_final_rarity(tmp_path):
    a=fixture(tmp_path)
    with gzip.open(a.sites,'rt') as f:s=f.read()
    s=s.replace('22\t1\t.\tA\tG\t.\tPASS\tAC=1;AN=2000;AF=0.99;MLEAF=0.99',
                '22\t1\t.\tA\tG\t.\tPASS\tAC=4;AN=1000;AF=0.123;MLEAF=0.9')
    with gzip.open(a.sites,'wt') as f:f.write(s)
    run(a)
    con=duckdb.connect()
    for name in ['missense','lof_hc']:
        found=con.execute('SELECT pcf_cohort_af,pcf_source_info_ac,pcf_source_info_an,pcf_source_info_af FROM read_parquet(?) WHERE POS=1',[str(Path(a.outdir)/(name+'.filtered.parquet'))]).fetchone()
        assert found==(0.004,'4','1000','0.123')
        assert not con.execute('SELECT * FROM read_parquet(?) WHERE POS=2',[str(Path(a.outdir)/(name+'.filtered.parquet'))]).fetchall()
