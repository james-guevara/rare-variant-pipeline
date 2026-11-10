import csv
import gzip
import importlib.util
import json
from pathlib import Path
import sqlite3
import pytest

spec=importlib.util.spec_from_file_location('genome',Path(__file__).parents[1]/'scripts/gather_genome_burdens.py')
g=importlib.util.module_from_spec(spec);spec.loader.exec_module(g)
STRATA=[('lof_hc','sequence','untiered'),('missense','sequence','miss_t1'),('lof_hc','spanning_deletion','untiered')]


def fixture(tmp_path):
    entries=[]
    for chrom in ['1','2']:
        p=tmp_path/('chr'+chrom);p.mkdir()
        (p/'samples.tsv').write_text('sample\ns1\ns2\nzero\n')
        fields=['CHROM','POS','REF','ALT','sample','candidate_type','allele_class','tier','Gene','qc_effective_alt_dosage','source_unit_id']
        data=[['chr'+chrom,'10','A','C','s1','lof_hc','sequence','untiered','G','1','u'+chrom],
              ['chr'+chrom,'10','A','C','s1','missense','sequence','miss_t1','G','1','u'+chrom],
              ['chr'+chrom,'20','A','*','s2','lof_hc','spanning_deletion','untiered','G','2','u'+chrom]]
        with gzip.open(p/'carriers.qc.tsv.gz','wt') as f:
            w=csv.writer(f,delimiter='\t');w.writerow(fields);w.writerows(data)
        receipt=dict(stage='post_rarity_qc_gather',status='passed',reconciliation={'passed':True},full_manifest_selected=True,
            chromosome='chr'+chrom,summary_dosage_field='qc_effective_alt_dosage',policy={'version':'test'},psam_identity={'sha256':'same'},
            duplicate_checks=dict(cross_block_annotation_associations=0,cross_block_exact_alleles=0),selected_units=['u'+chrom],manifest_units=['u'+chrom],
            by_stratum=[dict(candidate_type=t,allele_class=c,tier=tier,carrier_annotation_records=1,observed_alt_alleles=2 if c=='spanning_deletion' else 1) for t,c,tier in STRATA],
            samples_in_source=3,carrier_annotation_records=3,source_counts={'pass_rows':3},
            outputs={n:g.identity(p/n) for n in ['samples.tsv','carriers.qc.tsv.gz']})
        (p/'receipt.json').write_text(json.dumps(receipt));entries.append('chr'+chrom+'\t'+str(p))
    manifest=tmp_path/'manifest.tsv';manifest.write_text('chromosome\tdirectory\n'+'\n'.join(entries)+'\n')
    return manifest


def test_union_dosage_zeros_and_gene_across_chromosomes(tmp_path):
    m=fixture(tmp_path);out=tmp_path/'out';r=g.run(m,['1','2'],out)
    assert r['carrier_annotation_records']==6
    assert r['by_allele_class']['sequence']['samples_with_carriers']==1
    assert r['by_allele_class']['sequence']['distinct_allele_sample_records']==2
    assert r['by_allele_class']['sequence']['carrier_annotation_records']==4
    assert r['by_allele_class']['spanning_deletion']['distinct_allele_sample_observed_alt_alleles']==4
    assert all(x['samples_with_carriers']==1 for x in r['by_stratum'])
    genes=list(g.table(out/'gene_burden.tsv'))
    assert all(x['carrier_samples']=='1' for x in genes)
    zero=[x for x in g.table(out/'sample_burden.tsv') if x['sample']=='zero']
    assert len(zero)==3 and all(x['carrier_variant_records']=='0' for x in zero)
    assert list(g.table(out/'sample_distinct_alleles.tsv'))[0]['sequence_distinct_variant_count']=='2'


@pytest.mark.parametrize('change',['missing','duplicate','policy','psam','incomplete','checksum','stratum'])
def test_fail_closed(tmp_path,change):
    m=fixture(tmp_path)
    if change=='missing':m.write_text('\n'.join(m.read_text().splitlines()[:2])+'\n')
    elif change=='duplicate':m.write_text(m.read_text().replace('chr2\t','chr1\t'))
    elif change=='checksum':(tmp_path/'chr2/samples.tsv').write_text('sample\ns2\ns1\nzero\n')
    else:
        p=tmp_path/'chr2/receipt.json';r=json.loads(p.read_text())
        if change=='policy':r['policy']={'version':'different'}
        if change=='psam':r['psam_identity']={'sha256':'different'}
        if change=='incomplete':r['full_manifest_selected']=False
        if change=='stratum':r['by_stratum'][0]['carrier_annotation_records']=9
        p.write_text(json.dumps(r))
    with pytest.raises((ValueError,sqlite3.IntegrityError)):
        g.run(m,['1','2'],tmp_path/'out')
    assert not (tmp_path/'out/receipt.json').exists()


def test_nextflow_wrapper_and_resume(tmp_path):
    import os
    import shutil
    import subprocess
    if not shutil.which('nextflow'):
        pytest.skip('Nextflow unavailable')
    m=fixture(tmp_path)
    root=Path(__file__).parents[1]
    cfg=tmp_path/'test.config';cfg.write_text("process.executor='local'\nprocess.memory='1 GB'\n")
    cmd=['nextflow','-C',str(cfg),'run',str(root/'gather_genome_burdens.nf'),'-ansi-log','false',
         '--chromosome_manifest',str(m),'--expected_chromosomes','1,2','--outdir',str(tmp_path/'published'),
         '-work-dir',str(tmp_path/'work'),'-with-trace',str(tmp_path/'first.trace')]
    result=subprocess.run(cmd,cwd=tmp_path,capture_output=True,text=True,env=os.environ)
    assert result.returncode==0,result.stdout+result.stderr
    r=json.loads((tmp_path/'published/genome/receipt.json').read_text())
    assert r['carrier_annotation_records']==6 and r['samples_in_source']==3
    cmd[-1]=str(tmp_path/'resume.trace')
    result=subprocess.run(cmd+['-resume'],cwd=tmp_path,capture_output=True,text=True,env=os.environ)
    assert result.returncode==0,result.stdout+result.stderr
    assert 'CACHED' in (tmp_path/'resume.trace').read_text()
