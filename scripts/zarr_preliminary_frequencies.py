"""All-source-sample candidate AC/AN for an explicit preliminary rarity screen."""
import argparse,csv,json
from pathlib import Path
import duckdb
from zarr_carrier_source import Source,identity
from zarr_frequencies import ZarrFrequencyCounter
from carrier_psam import load_psam
from carrier_frequencies import validate_chromosome,validate_candidate_regions
from extract_exact_carriers import chrom,file_identity

def run(a):
    validate_chromosome(a.chromosome,a.sex_chromosome_policy)
    con=duckdb.connect();keys=set()
    for path in [a.missense,a.lof_hc]:
        for c,p,r,t in con.execute('SELECT CHROM,POS,REF,ALT FROM read_parquet(?)',[str(path)]).fetchall():
            if chrom(c)!=chrom(a.chromosome):raise ValueError('Candidate chromosome mismatch')
            keys.add((chrom(c),int(p),r,t))
    con.close();validate_candidate_regions(keys,a.chromosome,a.sex_chromosome_policy)
    before=identity(a.zarr);inputs={k:file_identity(getattr(a,k),True) for k in ['missense','lof_hc','psam']}
    with Source(a.zarr,a.chromosome,genotype_only=True) as src:
        _,metadata,_=load_psam(a.psam,src.samples)
        # Every source sample counts here, independent of later participant/family selection.
        meta=[dict(r,participant_id=r['#IID'],frequency_representative='1',unrelated='1') for r in metadata]
        counter=ZarrFrequencyCounter(meta,a.chromosome,src.samples,a.sex_chromosome_policy);found={}
        for key,record in src.records(keys):
            v,_=counter.count(record,'spanning_deletion' if key[3]=='*' else 'sequence')
            found[key]=(v['cohort_ac'],v['cohort_an'])
        with open(a.output,'w') as f:
            w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(['CHROM','POS','REF','ALT','matched','AC','AN'])
            for k in sorted(keys):w.writerow([*k,int(k in found),*found.get(k,('.','.'))])
        if identity(a.zarr)!=before:raise ValueError('Source metadata changed')
        if inputs!={k:file_identity(getattr(a,k),True) for k in inputs}:raise ValueError('Input changed')
        receipt=dict(status='passed',frequency_source='zarr_genotypes_all_source_samples',source=before,input_identities=inputs,samples=len(src.samples),candidate_alleles=len(keys),matched_alleles=len(found),genotype_chunk_reads=src.chunk_reads,selection='all source samples; no participant or unrelated subsetting; no genotype QC',sex_chromosome_policy=a.sex_chromosome_policy,counting_audit=counter.receipt(),output=file_identity(a.output,True))
        Path(a.receipt).write_text(json.dumps(receipt,indent=2)+'\n')
if __name__=='__main__':
    p=argparse.ArgumentParser()
    for k in ['zarr','chromosome','missense','lof-hc','psam','output','receipt']:p.add_argument('--'+k,required=True)
    p.add_argument('--sex-chromosome-policy',default=None);run(p.parse_args())
