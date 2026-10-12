import json,csv,gzip,sys,subprocess,shutil
from pathlib import Path
from argparse import Namespace
import numpy as np
import duckdb
from test_pre_carrier_filter import fixture
from zarr_preliminary_frequencies import run as frequencies
from filter_pre_carrier import run as filtering
import zarr

def setup(tmp):
 a=fixture(tmp)
 # Original INFO remains missing throughout.
 with gzip.open(a.sites,'wt') as f:
  f.write('##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
  for i in range(1,21):f.write('22\t%d\t.\tA\t%s\t.\tPASS\tAC=.;AN=.\n'%(i,'*' if i==20 else 'G'))
 g=zarr.open_group(str(tmp/'source.zarr'),mode='w')
 def arr(k,x):g.create_array(k,data=np.asarray(x))
 arr('sample_id',['s'+str(i) for i in range(200)]);arr('contig_id',['chr22']);arr('filter_id',['PASS'])
 arr('variant_position',np.arange(1,21));arr('variant_contig',np.zeros(20,dtype='i4'))
 arr('variant_allele',[['A','*' if i==20 else 'G'] for i in range(1,21)])
 arr('variant_filter',np.ones((20,1),bool));gt=np.zeros((20,200,2),dtype='i1');gt[:,0,1]=1;gt[1,1,1]=1
 arr('call_genotype',gt);arr('call_genotype_mask',gt<0)
 ps=tmp/'samples.psam';ps.write_text('#IID\tSEX\tparticipant_id\tfrequency_representative\tunrelated\n'+''.join('s%d\t1\tp%d\t0\t0\n'%(i,i) for i in range(200)))
 f=Namespace(zarr=str(tmp/'source.zarr'),chromosome='chr22',missense=a.missense,lof_hc=a.lof_hc,psam=str(ps),output=str(tmp/'freq.tsv'),receipt=str(tmp/'freq.json'),sex_chromosome_policy=None)
 return a,f

def test_all_samples_preliminary_and_source_preserved(tmp_path):
 a,f=setup(tmp_path);frequencies(f)
 rr=list(csv.DictReader(open(f.output),delimiter='\t'));assert rr[0]['AC']=='1' and rr[0]['AN']=='400'
 a.cohort_frequencies=f.output;a.frequency_receipt=f.receipt;filtering(a)
 data=duckdb.sql("select POS,pcf_cohort_ac,pcf_cohort_an,pcf_source_info_ac,pcf_cohort_pass from read_parquet(?) where candidate_type='missense' order by POS",params=[str(Path(a.outdir)/'filter_audit.parquet')]).fetchall()
 assert data[0]==(1,1,400,'.',True)
 assert data[1]==(2,2,400,'.',False) # equality at .005 fails
 r=json.load(open(Path(a.outdir)/'receipt.json'));assert r['preliminary_frequency_provenance']['samples']==200

def test_nextflow_preliminary(tmp_path):
 if not shutil.which('nextflow'):
  import pytest;pytest.skip('Nextflow required')
 a,f=setup(tmp_path);root=Path(__file__).resolve().parents[1];meta=json.load(open(a.metadata));res=tmp_path/'resources'
 lock={'schema':1,'stage':'pre_carrier','resource_root':str(res),'popmax':{'chr22':meta['popmax']},'regions':meta['regions'],'container':meta['container']}
 lp=tmp_path/'lock.json';lp.write_text(json.dumps(lock))
 manifest=tmp_path/'manifest.tsv';manifest.write_text('unit_id\tchromosome\tmissense\tlof_hc\tsites\tzarr\nb1\tchr22\t'+ '\t'.join([a.missense,a.lof_hc,a.sites,f.zarr])+'\n')
 config=tmp_path/'python.config';config.write_text('process.beforeScript="export PATH='+str(Path(sys.executable).parent)+':\\$PATH"\nprocess.memory="1 GB"\n')
 r=subprocess.run(['nextflow','-C',str(root/'zarr_pre_carrier.config')+','+str(config),'run',str(root/'zarr_pre_carrier.nf'),'--filter_manifest',str(manifest),'--filter_resource_lock',str(lp),'--filter_resource_root',str(res),'--psam',f.psam,'--select_units','all','--outdir',str(tmp_path/'published')],cwd=tmp_path,capture_output=True,text=True,timeout=120)
 assert r.returncode==0,r.stdout+r.stderr
 assert json.load(open(tmp_path/'published/pre-carrier/b1/receipt.json'))['status']=='passed'
