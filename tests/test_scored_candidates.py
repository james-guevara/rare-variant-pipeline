import csv
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
from argparse import Namespace

import duckdb
import pytest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
sys.path.insert(0,str(ROOT/'scripts/postprocess'))
from annotation_resources import identity
from join_scores import SCORES
from tier_variants import T_STARS
from select_scored_candidates import run, lit


def fixture(tmp_path):
    tmp_path.mkdir(exist_ok=True)
    root=tmp_path/'resources';root.mkdir()
    fields=['CHROM','POS','REF','ALT','Gene','Feature','SYMBOL','Consequence']
    rows=[['chr22',str(i),'A','G','ENSG1','TX1','SAME','missense_variant'] for i in range(1,10)]
    rows += [['chr22','10','A','G','ENSG1','TX1','SAME','synonymous_variant'],
             ['chr22','12','A','*','ENSG1','TX1','SAME','missense_variant']]
    lof=[]
    for pos,gene,alt,status in [(20,'ENSG1','G','HC'),(21,'ENSG2','*','HC'),(22,'ENSG3','G','HC'),
                               (23,'ENSG_MISSING','G','HC'),(24,'ENSG1','G','LC'),(25,'ENSG4','G','HC')]:
        row=['22',str(pos),'A',alt,gene,'TX'+str(pos),'SAME','frameshift_variant']
        rows.append(row);lof.append(row+[status])
    for name,head,data in [('picked.tsv',fields,rows),('loftee.tsv',fields+['LoF'],lof)]:
        with (tmp_path/name).open('w',newline='') as f:
            w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(head);w.writerows(data)
    score_fields=[entry[0] for entry in SCORES]
    records=[]
    for i in [1,2,3,4,5,7,8,9]:
        row={'#chr':'22','pos(1-based)':str(i),'ref':'A','alt':'G',**{name:'.' for name in score_fields}}
        for j,(name,threshold) in enumerate(T_STARS.items()):
            row[name]=str(threshold if j < max(0,5-i) else threshold-0.00001)
        row['MPC_score']='.;1.2;3.4';row['popEVE_score']='-2.0;.;-8.0'
        if i==7:row['alt']='T'
        if i==8:row['ref']='C'
        if i==9:row['#chr']='21'
        records.append(row)
    records.append(dict(records[0])) # duplicate exact key must not multiply output
    db_tsv=root/'scores.tsv'
    with db_tsv.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(records[0]),delimiter='\t');w.writeheader();w.writerows(records)
    db=root/'chr22.parquet'
    con=duckdb.connect();con.execute(f"COPY (SELECT * FROM read_csv({lit(db_tsv)},delim='\t',header=true,all_varchar=true)) TO {lit(db)} (FORMAT PARQUET)");con.close()
    gene=root/'GeneBayes.tsv'
    gene.write_text('ensg\thgnc\tobs_lof\texp_lof\tprior_mean\tpost_mean\tpost_lower_95\tpost_upper_95\n'+''.join(
        f'{g}\tSAME\t1\t2\t0.1\t{value}\t0.01\t0.5\n' for g,value in [('ENSG1',0.18),('ENSG2',0.03),('ENSG3',0.02999),('ENSG4',0.17999)]))
    sif=tmp_path/'image.sif';sif.write_text('synthetic image')
    meta=dict(unit_id='block12',chromosome='chr22',picked=str(tmp_path/'picked.tsv'),loftee=str(tmp_path/'loftee.tsv'),
              dbnsfp=identity(db),genebayes=identity(gene),container=identity(sif))
    metadata=tmp_path/'unit.json';metadata.write_text(json.dumps(meta))
    return Namespace(picked=meta['picked'],loftee=meta['loftee'],metadata=str(metadata),
                     postprocess_dir=str(ROOT/'scripts/postprocess'),outdir=str(tmp_path/'out'),threads=1,memory='512MB')


def test_exact_scores_missing_annotations_thresholds_and_gene_join(tmp_path,capsys):
    a=fixture(tmp_path);run(a)
    out=Path(a.outdir);r=json.loads((out/'receipt.json').read_text())
    assert r['missense_input_rows']==10 and r['hc_rows']==5
    seq=r['by_allele_class']['sequence'];star=r['by_allele_class']['spanning_deletion']
    assert seq['dbnsfp_matches']==5 and seq['dbnsfp_nonmatches']==4
    assert seq['n_flag_counts']=={'0':5,'1':1,'2':1,'3':1,'4':1}
    assert seq['n_scored_counts']=={'0':4,'1':0,'2':0,'3':0,'4':5}
    assert seq['fully_scored_below_thresholds']==1
    assert seq['matched_without_rankscores']==0
    assert star['dbnsfp_nonmatches']==1 and star['missense_selected_rows']==0
    assert star['hc_rows']==1 and star['lof_tier_counts']['lof_t2']==1
    assert seq['genebayes_matches']==3 and seq['genebayes_nonmatches']==1
    assert seq['lof_tier_counts']=={'lof_t1':1,'lof_t2':1,'untiered':2}
    assert r['duplicate_dbnsfp_keys']==1 and r['outputs']['missense.parquet']['rows']==4
    con=duckdb.connect()
    miss=con.execute('SELECT POS,tier,n_flag,MPC_score,popEVE_score FROM read_parquet(?) ORDER BY POS',[str(out/'missense.parquet')]).fetchall()
    assert miss==[(i,f'miss_t{i}',5-i,3.4,-8.0) for i in range(1,5)]
    hc=con.execute('SELECT POS,tier,allele_class,genebayes_matched FROM read_parquet(?) ORDER BY POS',[str(out/'lof_hc.parquet')]).fetchall()
    assert hc==[(20,'lof_t1','sequence',True),(21,'lof_t2','spanning_deletion',True),(22,None,'sequence',True),(23,None,'sequence',False),(25,'lof_t2','sequence',True)]
    assert not list(out.glob('*.tsv'))
    assert 'ENSG1' not in capsys.readouterr().out


def test_mpc_non_mane_value_and_missing_rankscore_policy(tmp_path):
    a=fixture(tmp_path)
    meta=json.loads(Path(a.metadata).read_text());db=Path(meta['dbnsfp']['path'])
    con=duckdb.connect()
    # One scored non-MANE transcript and a missing MANE value: historical
    # list-max behavior must not be replaced by positional MANE extraction.
    con.execute('CREATE TABLE source AS SELECT *, ? AS MANE, ? AS Ensembl_transcriptid FROM read_parquet(?)',
                ['.;Select','OLD_TX;MANE_TX',str(db)])
    con.execute("UPDATE source SET MPC_score='3.4;.'")
    # POS 2 stays selected on three other predictors, with unavailable MPC.
    con.execute("UPDATE source SET MPC_rankscore='.' WHERE \"pos(1-based)\"='2'")
    for col in T_STARS:
        con.execute(f'UPDATE source SET "{col}"=\'.\' WHERE "pos(1-based)"=\'5\'')
    db.unlink();con.execute(f'COPY source TO {lit(db)} (FORMAT PARQUET)')
    meta['dbnsfp']=identity(db);Path(a.metadata).write_text(json.dumps(meta))
    run(a)
    rows=con.execute('SELECT POS,MPC_score,MPC_rankscore,n_scored,n_flag,tier FROM read_parquet(?) ORDER BY POS',
                     [str(Path(a.outdir)/'missense.parquet')]).fetchall()
    assert rows[0]==(1,3.4,T_STARS['MPC_rankscore'],4,4,'miss_t1')
    assert rows[1]==(2,3.4,None,3,3,'miss_t2')
    r=json.loads((Path(a.outdir)/'receipt.json').read_text())['by_allele_class']['sequence']
    assert r['matched_without_rankscores']==1
    assert r['dbnsfp_nonmatches']==4
    assert r['rankscore_missing_counts']['MPC_rankscore']==6
    assert r['fully_scored_below_thresholds']==0
    assert r['n_scored_counts']=={'0':5,'1':0,'2':0,'3':1,'4':3}


@pytest.mark.parametrize('kind',['empty','duplicate_gene','wrong_pair','changed_resource'])
def test_edge_inputs(tmp_path,kind):
    a=fixture(tmp_path)
    if kind=='empty':
        for name in ['picked','loftee']:
            p=Path(getattr(a,name));p.write_text(p.read_text().splitlines()[0]+'\n')
        run(a)
        r=json.loads((Path(a.outdir)/'receipt.json').read_text())
        assert r['missense_input_rows']==r['hc_rows']==0
        assert all(v['rows']==0 for v in r['outputs'].values())
        return
    meta=json.loads(Path(a.metadata).read_text())
    if kind=='duplicate_gene':
        p=Path(meta['genebayes']['path']);p.write_text(p.read_text()+p.read_text().splitlines()[1]+'\n')
        meta['genebayes']=identity(p);Path(a.metadata).write_text(json.dumps(meta))
    if kind=='changed_resource':
        Path(meta['genebayes']['path']).write_text('changed')
    if kind=='wrong_pair':
        p=Path(a.loftee);p.write_text(p.read_text().replace('TX20','DIFFERENT'))
    with pytest.raises(ValueError):run(a)
    assert json.loads((Path(a.outdir)/'receipt.json').read_text())['status']=='failed'
    assert not list(Path(a.outdir).glob('*.parquet'))


@pytest.mark.skipif(not shutil.which('nextflow'),reason='Nextflow required')
def test_nextflow_pilot_resume_publication(tmp_path):
    a=fixture(tmp_path/'a');b=fixture(tmp_path/'b')
    # Share immutable resources but keep separate block inputs.
    meta=json.loads(Path(a.metadata).read_text())
    lock=dict(schema=1,resource_root=str(tmp_path/'a/resources'),dbnsfp={'chr22':meta['dbnsfp']},genebayes=meta['genebayes'],container=meta['container'])
    resource_lock=tmp_path/'lock.json';resource_lock.write_text(json.dumps(lock))
    manifest=tmp_path/'blocks.tsv';manifest.write_text('unit_id\tchromosome\tpicked\tloftee\n'+''.join(
        f'{unit}\tchr22\t{x.picked}\t{x.loftee}\n' for unit,x in [('block12',a),('block19',b)]))
    config=tmp_path/'python.config';config.write_text('process.beforeScript = "export PATH='+str(Path(sys.executable).parent)+':\\$PATH"\n')
    base=['nextflow','-C',str(ROOT/'candidates.config')+','+str(config),'run',str(ROOT/'candidates.nf'),'-ansi-log','false',
          '--candidate_manifest',str(manifest),'--candidate_resource_lock',str(resource_lock),
          '--candidate_resource_root',lock['resource_root'],'--outdir',str(tmp_path/'published'),'--candidate_memory','1 GB']
    def execute(label,units,resume=False):
        cmd=base+['--select_units',units,'-with-trace',str(tmp_path/(label+'.trace'))]
        if resume:cmd+=['-resume']
        r=subprocess.run(cmd,cwd=tmp_path,capture_output=True,text=True,timeout=120)
        assert r.returncode==0,r.stdout+r.stderr
        return list(csv.DictReader((tmp_path/(label+'.trace')).open(),delimiter='\t'))
    assert len(execute('one','block12'))==1
    second=execute('two','all',True)
    assert {r['name']:r['status'] for r in second}=={'SCORE_CANDIDATES (block12)':'CACHED','SCORE_CANDIDATES (block19)':'COMPLETED'}
    assert all(r['status']=='CACHED' for r in execute('subset','block12',True))
    shutil.rmtree(tmp_path/'work')
    out=tmp_path/'published/candidates/block12'
    assert not (out/'missense.parquet').is_symlink()
    assert json.loads((out/'receipt.json').read_text())['outputs']['missense.parquet']['rows']==4


def test_tiers_match_existing_tier_script(tmp_path):
    a=fixture(tmp_path);run(a)
    con=duckdb.connect()
    old_input=tmp_path/'existing-tier-input.parquet'
    con.execute(f'''COPY (
        SELECT * EXCLUDE(tier,miss_n_flag,n_flag), CAST(NULL AS VARCHAR) AS LoF,
               CAST(NULL AS DOUBLE) AS genebayes_post_mean
        FROM read_parquet({lit(Path(a.outdir)/'missense.parquet')})
        UNION ALL BY NAME
        SELECT * EXCLUDE(tier) FROM read_parquet({lit(Path(a.outdir)/'lof_hc.parquet')})
        ) TO {lit(old_input)} (FORMAT PARQUET)''')
    config=tmp_path/'resources.json';config.write_text(json.dumps({'output_base':str(tmp_path)}))
    output=tmp_path/'existing-tier-output.parquet'
    result=subprocess.run([sys.executable,str(ROOT/'scripts/postprocess/tier_variants.py'),
        '--cohort','synthetic','--chrom','chr22','--resources',str(config),
        '--input',str(old_input),'--output',str(output)],capture_output=True,text=True)
    assert result.returncode==0,result.stderr
    old=con.execute('SELECT POS,tier FROM read_parquet(?) ORDER BY POS',[str(output)]).fetchall()
    new=con.execute('SELECT POS,tier FROM read_parquet(?) UNION ALL SELECT POS,tier FROM read_parquet(?) ORDER BY POS',
                    [str(Path(a.outdir)/'missense.parquet'),str(Path(a.outdir)/'lof_hc.parquet')]).fetchall()
    assert old==new


def test_resource_lock_canonical_hash_enforcement(tmp_path,monkeypatch):
    import lock_candidate_resources as lock
    a=fixture(tmp_path)
    meta=json.loads(Path(a.metadata).read_text())
    root=tmp_path/'bundle';dbdir=root/'dbNSFP/5.3.1a/parquet_expanded_mane_select';dbdir.mkdir(parents=True)
    genedir=root/'GeneBayes';genedir.mkdir()
    shutil.copy(meta['dbnsfp']['path'],dbdir/'chr22.parquet')
    shutil.copy(meta['genebayes']['path'],genedir/'GeneBayes.Supplementary_Table_1.tsv')
    old_lock=tmp_path/'annotation-lock.json';old_lock.write_text(json.dumps({'containers':{'loftee':meta['container']}}))
    canonical=tmp_path/'canonical.json'
    canonical.write_text(json.dumps({'files':[
        {**meta['dbnsfp'],'path':'dbNSFP/5.3.1a/parquet_expanded_mane_select/chr22.parquet'},
        {**meta['genebayes'],'path':'targeted-annotation/GeneBayes.Supplementary_Table_1.tsv'}]}))
    monkeypatch.setattr(lock,'LOF_SIF',meta['container']['sha256'])
    data=lock.build(root,['22'],old_lock,canonical)
    assert data['dbnsfp']['chr22']['sha256']==meta['dbnsfp']['sha256']
    (genedir/'GeneBayes.Supplementary_Table_1.tsv').write_text('changed')
    with pytest.raises(ValueError,match='canonical v1'):
        lock.build(root,['22'],old_lock,canonical)
