"""Explicit, deterministic pedigree-based frequency selection; not KING unrelateds."""
import argparse,csv,json,hashlib
from collections import defaultdict
from pathlib import Path


def build(metadata,source_psam,output,sample_col,participant_col,family_col,sex_col):
    with open(metadata) as f:rows=list(csv.DictReader(f,delimiter='\t'))
    with open(source_psam) as f:
        reader=csv.DictReader(f,delimiter='\t');samples=[r.get('#IID',r.get('IID')) for r in reader]
    if not samples or None in samples or len(samples)!=len(set(samples)):raise ValueError('Invalid source sample IDs')
    lookup={r[sample_col]:r for r in rows}
    if len(lookup)!=len(rows):raise ValueError('Duplicate metadata sample IDs')
    selected=[lookup[s] for s in samples]
    bad={'','0','.','NA','N/A'}
    if any(r[participant_col] in bad or r[family_col] in bad for r in selected):raise ValueError('Missing participant/family identity; explicit resolution required')
    by_participant=defaultdict(list)
    for s,r in zip(samples,selected):by_participant[r[participant_col]].append((s,r))
    representatives={};families=defaultdict(list)
    for participant,rs in by_participant.items():
        if len({r[family_col] for _,r in rs})!=1:raise ValueError('Participant has conflicting family IDs')
        s,r=min(rs,key=lambda pair:pair[0]);representatives[participant]=s;families[r[family_col]].append(participant)
    unrelated={representatives[min(ps)] for ps in families.values()}
    fields=['#IID','SEX','participant_id','frequency_representative','unrelated']
    with open(output,'w') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader()
        for s,r in zip(samples,selected):
            sex={'Male':'1','Female':'2','1':'1','2':'2'}.get(r[sex_col],'0')
            w.writerow({'#IID':s,'SEX':sex,'participant_id':r[participant_col],
              'frequency_representative':int(representatives[r[participant_col]]==s),'unrelated':int(s in unrelated)})
    sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
    receipt=dict(policy='one-lexicographically-first-participant-per-family-v1',
        representative_rule='lexicographically first source sample per participant',
        unrelated_rule='representative of lexicographically first participant per family',
        limitation='pedigree-based proxy; no genetic relatedness verification or cross-family relationship exclusion',
        samples=len(samples),participants=len(representatives),families=len(families),unrelated=len(unrelated),
        metadata_sha256=sha(metadata),source_psam_sha256=sha(source_psam),output_sha256=sha(output))
    Path(str(output)+'.selection.json').write_text(json.dumps(receipt,indent=2)+'\n')
    return receipt

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for n in ['metadata','source-psam','output','sample-col','participant-col','family-col','sex-col']:p.add_argument('--'+n,required=True)
    a=p.parse_args();print(json.dumps(build(a.metadata,a.source_psam,a.output,a.sample_col,a.participant_col,a.family_col,a.sex_col)))
