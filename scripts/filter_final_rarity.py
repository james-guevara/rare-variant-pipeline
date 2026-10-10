#!/usr/bin/env python3
"""Final exact-allele rarity from saved frequencies; no VCF reads or genotype QC."""
import argparse
from collections import Counter, defaultdict
import csv
from decimal import Decimal, InvalidOperation
import gzip
import json
from pathlib import Path
import re
import sys
import time
from qc_filtered_carriers import (ValidationError, identity, write_table, summaries,
                                 KEYS, TYPES, CLASSES, group, row_key)

INPUT_NAMES = dict(carriers='carriers.tsv.gz', samples='samples.tsv',
                   candidates='candidates.tsv', frequencies='variant_frequencies.tsv')
PRODUCTS = ['carriers.rare.tsv.gz','candidates.rare.tsv','variant_frequencies.rare.tsv',
            'rarity_audit.tsv','samples.tsv','sample_burden.tsv','sample_gene_burden.tsv',
            'gene_burden.tsv','sample_distinct_alleles.tsv']
POLICY = dict(unrelated_af='saved unrelated_af < 0.001', denominator='saved unrelated_an > 0',
              missing='missing/invalid/nonfinite AF or AN fails closed; unmatched fails',
              annotation_scope='one exact-allele decision applies to all candidate types/transcripts/carriers',
              genotype_qc='none; raw genotypes/dosages unchanged',
              stars='separate spanning-deletion audit records; never merged with sequence biological burdens')
MISSING = ('', '.', 'NA', 'null', 'None')


def table(path, required):
    with (gzip.open(path,'rt',newline='') if str(path).endswith('.gz') else open(path,newline='')) as f:
        reader=csv.DictReader(f,delimiter='\t');fields=reader.fieldnames or []
        if not set(required)<=set(fields) or len(fields)!=len(set(fields)):
            raise ValidationError('Invalid saved input schema')
        rows=list(reader)
    if any(None in r or any(v is None for v in r.values()) for r in rows):
        raise ValidationError('Malformed saved input row')
    return fields,rows


def key(row, chromosome):
    if row['CHROM'].removeprefix('chr')!=chromosome.removeprefix('chr') or not re.fullmatch('[0-9]+',row['POS']) or int(row['POS'])<1:
        raise ValidationError('Invalid saved allele locus')
    if not re.fullmatch('[ACGTN]+',row['REF']) or not (row['ALT']=='*' or re.fullmatch('[ACGTN]+',row['ALT'])) or row['REF']==row['ALT']:
        raise ValidationError('Invalid saved allele')
    if row['allele_class']!=('spanning_deletion' if row['ALT']=='*' else 'sequence'):
        raise ValidationError('Allele class mismatch')
    return row['CHROM'].removeprefix('chr'),int(row['POS']),row['REF'],row['ALT']


def decision(row):
    if row['matched'] not in ('True','False'):raise ValidationError('Invalid matched flag')
    if row['matched']=='False':return 'unmatched'
    an=row['unrelated_an'];af=row['unrelated_af']
    if an in MISSING:return 'missing_an'
    if not re.fullmatch('[0-9]+',an):return 'invalid_an'
    if int(an)==0:return 'zero_an'
    if af in MISSING:return 'missing_af'
    try:value=Decimal(af)
    except InvalidOperation:return 'invalid_af'
    if not value.is_finite() or value<0 or value>1:return 'invalid_af'
    return 'PASS' if value<Decimal('0.001') else 'af_at_or_above_threshold'


def annotation(row):
    if row['candidate_type'] not in TYPES or row['allele_class'] not in CLASSES:
        raise ValidationError('Invalid candidate type/class')
    allowed=['miss_t'+str(i) for i in range(1,5)] if row['candidate_type']=='missense' else ['lof_t1','lof_t2']
    if row['tier'] not in allowed+['untiered'] or any(row[k] in ('','.','-') for k in ['Gene','Feature']):
        raise ValidationError('Invalid candidate tier/gene/transcript')
    return tuple(row[k] for k in ['candidate_type','tier','Gene','Feature','SYMBOL','allele_class'])


def counts(total,passed):return dict(input=total,pass_count=passed,fail_count=total-passed)


def run(a):
    out=Path(a.outdir).resolve();out.mkdir(parents=True,exist_ok=True)
    meta=json.loads(Path(a.metadata).read_text());start=time.perf_counter()
    inputs={k:Path(getattr(a,k)).resolve() for k in [*INPUT_NAMES,'source_receipt']}
    if any(out/n in inputs.values() for n in PRODUCTS+['receipt.json']):
        raise ValidationError('Separate final-rarity output directory required')
    receipt=dict(schema_version=1,status='failed',stage='final_rarity',unit_id=meta['unit_id'],chromosome=meta['chromosome'],sources=meta,policy=POLICY)
    try:
        ids={k:identity(p) for k,p in inputs.items()}
        source=json.loads(inputs['source_receipt'].read_text())
        if source.get('status')!='passed' or source.get('unit_id')!=meta['unit_id'] or source.get('chromosome','').removeprefix('chr')!=meta['chromosome'].removeprefix('chr'):
            raise ValidationError('Source extraction receipt does not match passed unit')
        if source.get('frequencies',{}).get('status')!='passed':raise ValidationError('Passed corrected-frequency receipt required')
        for name,product in INPUT_NAMES.items():
            if source['outputs'][product]['sha256']!=ids[name]['sha256']:
                raise ValidationError('Saved input differs from extraction receipt')
        sf,sample_rows=table(inputs['samples'],['sample']);samples=[r['sample'] for r in sample_rows]
        if sf!=['sample'] or not samples or len(set(samples))!=len(samples) or any(s in ('','.') for s in samples) or len(samples)!=source['samples_in_source']:
            raise ValidationError('Invalid source sample universe')
        sample_set=set(samples)
        ff,freq=table(inputs['frequencies'],KEYS+['allele_class','matched','unrelated_af','unrelated_an'])
        frequency={};reasons={}
        for r in freq:
            k=key(r,meta['chromosome'])
            if k in frequency:raise ValidationError('Duplicate exact frequency allele')
            frequency[k]=r;reasons[k]=decision(r)
        cf,candidates=table(inputs['candidates'],KEYS+['candidate_type','tier','Gene','Feature','SYMBOL','allele_class','matched','carrier_records'])
        associations={};candidate_keys=set()
        for r in candidates:
            k=key(r,meta['chromosome']);sig=annotation(r);assoc=(k,r['candidate_type'])
            if assoc in associations:raise ValidationError('Duplicate candidate association')
            if k not in frequency:raise ValidationError('Candidate lacks exact frequency allele')
            if r['matched']!=frequency[k]['matched']:raise ValidationError('Candidate/frequency matched status differs')
            associations[assoc]=sig;candidate_keys.add(k)
        if candidate_keys!=set(frequency):raise ValidationError('Frequency/candidate allele sets differ')
        if len(candidates)!=source['candidate_annotation_records'] or len(frequency)!=source['distinct_candidate_alleles']:
            raise ValidationError('Candidate count differs from extraction receipt')
        fields,carriers=table(inputs['carriers'],KEYS+['sample','candidate_type','tier','Gene','Feature','SYMBOL','allele_class','GT','GQ','DP','AD','FT','site_FILTER','alt_dosage'])
        seen=set();genotypes={};observed=Counter();passed=[];carrier_counts=defaultdict(Counter);candidate_counts=defaultdict(Counter)
        for r in carriers:
            k=key(r,meta['chromosome']);sig=annotation(r);assoc=(k,r['candidate_type'])
            if associations.get(assoc)!=sig:raise ValidationError('Carrier lacks matching candidate annotation')
            if r['sample'] not in sample_set:raise ValidationError('Carrier outside source sample universe')
            if frequency[k]['matched']!='True':raise ValidationError('Carrier for unmatched allele')
            record=(assoc,r['sample'])
            if record in seen:raise ValidationError('Duplicate carrier association')
            seen.add(record)
            gt=tuple(r[n] for n in ['GT','GQ','DP','AD','FT','site_FILTER','alt_dosage'])
            gkey=(k,r['sample'])
            if gkey in genotypes and genotypes[gkey]!=gt:raise ValidationError('Conflicting cross-type carrier genotype')
            genotypes[gkey]=gt
            if not re.fullmatch('[1-9][0-9]*',r['alt_dosage']):raise ValidationError('Invalid saved ALT dosage')
            observed[assoc]+=1;g=group(r);carrier_counts[g]['input']+=1
            keep=reasons[k]=='PASS';carrier_counts[g]['pass']+=keep
            if keep:passed.append(r)
        if len(carriers)!=source['carrier_annotation_records']:raise ValidationError('Carrier count differs from extraction receipt')
        for r in candidates:
            k=key(r,meta['chromosome']);g=group(r)
            if str(observed[k,r['candidate_type']])!=r['carrier_records']:raise ValidationError('Candidate/carrier counts do not reconcile')
            candidate_counts[g]['input']+=1;candidate_counts[g]['pass']+=reasons[k]=='PASS'
        retained={k for k,r in reasons.items() if r=='PASS'}
        write_table(out/'carriers.rare.tsv.gz',fields,passed)
        write_table(out/'candidates.rare.tsv',cf,(r for r in candidates if key(r,meta['chromosome']) in retained))
        write_table(out/'variant_frequencies.rare.tsv',ff,(r for r in freq if key(r,meta['chromosome']) in retained))
        write_table(out/'rarity_audit.tsv',ff+['rarity_pass','rarity_reason'],({**r,'rarity_pass':key(r,meta['chromosome']) in retained,'rarity_reason':reasons[key(r,meta['chromosome'])]} for r in freq))
        write_table(out/'samples.tsv',['sample'],sample_rows)
        distinct=summaries(out,candidates,passed,samples)
        strata=sorted(set(candidate_counts)|set(carrier_counts))
        by_type={}
        for t in TYPES:
            by_type[t]={'by_allele_class':{}}
            for c in CLASSES:
                tiers=sorted({'untiered'}|{tr for typ,cl,tr in strata if (typ,cl)==(t,c)})
                sub={}
                for tr in tiers:
                    sub[tr]={name:counts(d[t,c,tr]['input'],d[t,c,tr]['pass']) for name,d in [('candidate_annotations',candidate_counts),('carrier_annotations',carrier_counts)]}
                by_type[t]['by_allele_class'][c]={**{name:counts(sum(sub[tr][name]['input'] for tr in tiers),sum(sub[tr][name]['pass_count'] for tr in tiers)) for name in ['candidate_annotations','carrier_annotations']},'by_tier':sub}
        for t in by_type.values():
            for name in ['candidate_annotations','carrier_annotations']:
                t[name]=counts(sum(c[name]['input'] for c in t['by_allele_class'].values()),sum(c[name]['pass_count'] for c in t['by_allele_class'].values()))
        by_class={}
        for c in CLASSES:
            keys={k for k,r in frequency.items() if r['allele_class']==c}
            raw={(row_key(r)) for r in carriers if r['allele_class']==c}
            kept={k for k in raw if k[:-1] in retained}
            by_class[c]=dict(distinct_alleles=counts(len(keys),len(keys&retained)),
                distinct_allele_sample_records=counts(len(raw),len(kept)),
                reasons=dict(Counter(reasons[k] for k in keys)))
        if ids!={k:identity(p) for k,p in inputs.items()}:raise ValidationError('Input changed during final rarity filtering')
        if sum(v['carrier_annotations']['pass_count'] for t in by_type.values() for v in t['by_allele_class'].values())!=len(passed):raise ValidationError('Summary count reconciliation failed')
        receipt.update(status='passed',input_identities=ids,source_frequency_policy=source['frequencies']['policy'],
            samples_in_source=len(samples),distinct_alleles=counts(len(freq),len(retained)),
            candidate_annotations=counts(len(candidates),sum(reasons[key(r,meta['chromosome'])]=='PASS' for r in candidates)),
            carrier_annotations=counts(len(carriers),len(passed)),by_candidate_type=by_type,by_allele_class=by_class,
            retained_distinct_alleles_by_class=distinct,reconciliation=dict(passed=True,source_hashes_match=True,source_counts_match=True,exact_allele_annotation_consistency=True,input_equals_pass_plus_fail=True),
            code_identities={Path(p).name:identity(p) for p in [__file__,Path(__file__).with_name('qc_filtered_carriers.py')]},
            outputs={n:identity(out/n) for n in PRODUCTS})
    except Exception as exc:
        receipt.update(error_type=type(exc).__name__,error=str(exc) if isinstance(exc,ValidationError) else 'Final rarity input/processing failure; inspect locally')
        for n in PRODUCTS:(out/n).unlink(missing_ok=True)
        raise
    finally:
        receipt['wall_seconds']=time.perf_counter()-start
        (out/'receipt.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n')
    print(json.dumps({k:receipt[k] for k in ['status','unit_id','distinct_alleles','candidate_annotations','carrier_annotations','by_allele_class']}))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in [*INPUT_NAMES,'source_receipt','metadata','outdir']:p.add_argument('--'+name.replace('_','-'),required=True)
    try:run(p.parse_args())
    except Exception:sys.exit('Final rarity failed; inspect local receipt. No protected records printed.')
