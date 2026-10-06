#!/usr/bin/env python3
"""Fixed post-extraction per-row QC; no upstream execution or family propagation."""
import argparse
from collections import Counter, defaultdict
import csv
import gzip
import hashlib
import io
import json
import math
from pathlib import Path
import re
import sys
import time

TYPES=('missense','lof_hc')
CLASSES=('sequence','spanning_deletion')
KEYS=['CHROM','POS','REF','ALT']
REASONS=('site_filter','gq_invalid','gq_below_min','dp_invalid','dp_below_min','gt_unsupported','ad_invalid','ab_out_of_range')
HET={'0/1','1/0','0|1','1|0'}
HOM={'1/1','1|1','1'}
POLICY=dict(site_FILTER='PASS',min_gq=20,min_dp=10,het_ab_min=0.25,het_ab_max=0.75,
            hom_alt_ab_min=0.90,haploid_alt='GT=1 uses AB >= 0.90; no chromosome-specific policy',
            ab='AD_alt / (AD_ref + AD_alt); exactly two finite nonnegative values with finite positive sum',
            malformed='missing/malformed/nonfinite/negative GQ/DP/AD fail closed; unsupported GT fails',
            FT='preserved but not filtered',family_propagation=False,
            ab_failure='ab_out_of_range only when GT and AD are evaluable; otherwise gt_unsupported/ad_invalid')


class ValidationError(ValueError):pass


def sha(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
    return h.hexdigest()


def identity(path):
    p=Path(path);s=p.stat()
    return dict(bytes=s.st_size,sha256=sha(p))


def number(value):
    try:
        if value is None or not re.fullmatch(r'[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?', str(value).strip()):return None
        n=float(value)
        return n if math.isfinite(n) and n>=0 else None
    except (ValueError,TypeError):return None


def evaluate(row):
    bad=[]
    if row['site_FILTER']!='PASS':bad.append('site_filter')
    for field,minimum in [('GQ',20),('DP',10)]:
        n=number(row[field]);name=field.lower()
        if n is None:bad.append(name+'_invalid')
        elif n<minimum:bad.append(name+'_below_min')
    gt=row['GT'];supported=gt in HET|HOM
    if not supported:bad.append('gt_unsupported')
    parts=row['AD'].split(',') if isinstance(row['AD'],str) else []
    values=[number(x) for x in parts]
    valid=len(values)==2 and all(x is not None for x in values)
    total=sum(values) if valid else None
    valid=valid and math.isfinite(total) and total>0
    ab=values[1]/total if valid else None
    if not valid:bad.append('ad_invalid')
    if valid and supported:
        if not (0.25<=ab<=0.75 if gt in HET else ab>=0.90):bad.append('ab_out_of_range')
    return ab,tuple(bad)


def write_table(path,fields,rows):
    # Reproducible gzip with no original filename/timestamp in the header.
    raw=open(path,'wb')
    zipped=gzip.GzipFile(filename='',mode='wb',fileobj=raw,mtime=0) if str(path).endswith('.gz') else raw
    with io.TextIOWrapper(zipped,newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n',extrasaction='ignore');w.writeheader()
        for r in rows:w.writerow({k:'.' if r.get(k) is None else r[k] for k in fields})
    if not raw.closed:raw.close()


def row_key(r):return (r['CHROM'].removeprefix('chr'),int(r['POS']),r['REF'],r['ALT'],r['sample'])
def group(r):return (r['candidate_type'],r['allele_class'],r['tier'])


def summaries(out,all_rows,passed,samples):
    strata=sorted({group(r) for r in all_rows}|{(t,c,'untiered') for t in TYPES for c in CLASSES})
    cells=Counter();dosages=Counter();distinct=set()
    for r in passed:
        cell=(r['sample'],*group(r),r['Gene']);cells[cell]+=1;dosages[cell]+=int(r['alt_dosage'])
        distinct.add((row_key(r),r['allele_class']))
    totals=Counter();total_dosage=Counter();gene_counts=Counter();gene_samples=Counter()
    for (s,t,c,tier,gene),n in cells.items():
        totals[s,t,c,tier]+=n;total_dosage[s,t,c,tier]+=dosages[s,t,c,tier,gene]
        gene_counts[t,c,tier,gene]+=n;gene_samples[t,c,tier,gene]+=1
    fields=['sample','candidate_type','allele_class','tier','carrier_variant_records','observed_alt_alleles']
    write_table(out/'sample_burden.tsv',fields,(dict(zip(fields,[s,*g,totals[(s,*g)],total_dosage[(s,*g)]])) for s in samples for g in strata))
    fields=['sample','candidate_type','allele_class','tier','Gene','carrier_variant_records','observed_alt_alleles']
    write_table(out/'sample_gene_burden.tsv',fields,(dict(zip(fields,[*cell,n,dosages[cell]])) for cell,n in sorted(cells.items())))
    fields=['candidate_type','allele_class','tier','Gene','carrier_variant_records','carrier_samples']
    genes=sorted({(*group(r),r['Gene']) for r in all_rows})
    write_table(out/'gene_burden.tsv',fields,(dict(zip(fields,[*g,gene_counts[g],gene_samples[g]])) for g in genes))
    counts=Counter((key[-1],c) for key,c in distinct)
    write_table(out/'sample_distinct_alleles.tsv',['sample','sequence_distinct_variant_count','spanning_deletion_distinct_record_count'],(
        dict(sample=s,sequence_distinct_variant_count=counts[s,'sequence'],spanning_deletion_distinct_record_count=counts[s,'spanning_deletion']) for s in samples))
    return {c:dict(distinct_variant_sample_records=sum(cl==c for key,cl in distinct),samples_with_carriers=sum(n>0 for (s,cl),n in counts.items() if cl==c)) for c in CLASSES}


def run(a):
    out=Path(a.outdir).resolve();out.mkdir(parents=True,exist_ok=True)
    meta=json.loads(Path(a.metadata).read_text());start=time.perf_counter()
    products=['carriers.qc.tsv.gz','qc_audit.tsv.gz','samples.tsv','sample_burden.tsv','sample_gene_burden.tsv','gene_burden.tsv','sample_distinct_alleles.tsv']
    inputs={k:Path(getattr(a,k)).resolve() for k in ['carriers','samples','source_receipt']}
    if any(out/n in inputs.values() for n in products+['receipt.json']):raise ValidationError('Separate QC output directory required')
    receipt=dict(schema_version=1,status='failed',stage='post_extraction_qc',unit_id=meta['unit_id'],sources=meta,policy=POLICY,
                 summary_policy='candidate-type associations kept separately; distinct allele/sample union within allele class; never combine star records with sequence burdens',
                 failure_count_policy='independent reason flags may overlap; mutually exclusive combinations partition failed rows; PASS combination contains passing rows')
    try:
        ids={k:identity(p) for k,p in inputs.items()}
        source=json.loads(inputs['source_receipt'].read_text())
        if source.get('status')!='passed' or source.get('unit_id')!=meta['unit_id']:raise ValidationError('Source extraction receipt does not match passed unit')
        for arg,name in [('carriers','carriers.tsv.gz'),('samples','samples.tsv')]:
            if source['outputs'][name]['sha256']!=ids[arg]['sha256']:raise ValidationError('Input differs from extraction receipt')
        with open(inputs['samples']) as f:
            reader=csv.DictReader(f,delimiter='\t')
            if reader.fieldnames!=['sample']:raise ValidationError('Invalid source samples schema')
            samples=[r['sample'] for r in reader]
        sample_set=set(samples)
        if not samples or any(not s or s=='.' for s in samples) or len(samples)!=len(sample_set):raise ValidationError('Expected unique nonempty source samples')
        if source['samples_in_source']!=len(samples):raise ValidationError('Source sample count disagrees with receipt')
        all_rows=[];passed=[];audits=[];counts=defaultdict(Counter);independent=Counter();combinations=Counter()
        seen=set();genotypes={}
        with gzip.open(inputs['carriers'],'rt',newline='') as f:
            reader=csv.DictReader(f,delimiter='\t');fields=reader.fieldnames or []
            required=set(KEYS+['Gene','Feature','SYMBOL','sample','GT','GQ','DP','AD','FT','alt_dosage','site_FILTER','allele_class','candidate_type','tier'])
            if not required<=set(fields) or len(fields)!=len(set(fields)) or any(x.startswith('qc_') for x in fields):raise ValidationError('Invalid raw carrier schema')
            for r in reader:
                if None in r or any(v is None for v in r.values()):raise ValidationError('Malformed carrier TSV row')
                if r['sample'] not in sample_set:raise ValidationError('Carrier sample absent from samples.tsv')
                if r['candidate_type'] not in TYPES or r['allele_class'] not in CLASSES:raise ValidationError('Invalid carrier type/class')
                if not re.fullmatch('[0-9]+',r['POS']) or int(r['POS'])<1 or r['CHROM'].removeprefix('chr')!=meta['chromosome'].removeprefix('chr'):raise ValidationError('Invalid carrier locus')
                if not re.fullmatch('[ACGTN]+',r['REF']) or not (r['ALT']=='*' or re.fullmatch('[ACGTN]+',r['ALT'])):raise ValidationError('Invalid carrier allele')
                if r['allele_class']!=('spanning_deletion' if r['ALT']=='*' else 'sequence'):raise ValidationError('Carrier class/allele mismatch')
                if any(r[k] in ('','.','-') for k in ['Gene','Feature']):raise ValidationError('Missing carrier gene/transcript')
                allowed=['miss_t'+str(i) for i in range(1,5)] if r['candidate_type']=='missense' else ['lof_t1','lof_t2']
                if r['tier'] not in allowed+['untiered']:raise ValidationError('Invalid carrier tier')
                key=row_key(r);association=(key,r['candidate_type'])
                if association in seen:raise ValidationError('Duplicate within-type allele/sample association')
                seen.add(association)
                signature=tuple(r[k] for k in ['GT','GQ','DP','AD','FT','site_FILTER','alt_dosage'])
                if key in genotypes and genotypes[key]!=signature:raise ValidationError('Conflicting genotypes across candidate types')
                genotypes[key]=signature
                ab,bad=evaluate(r);g=group(r);counts[g]['input']+=1;counts[g]['fail' if bad else 'pass']+=1
                combo='|'.join(bad) if bad else 'PASS';combinations[combo]+=1
                counts[g]['combination:'+combo]+=1
                for reason in bad:independent[reason]+=1;counts[g]['reason:'+reason]+=1
                if not bad:
                    # Validate bookkeeping only for supported passing GT; genotype
                    # failures (including partial GT) remain row-level QC failures.
                    if r['alt_dosage']!=str(r['GT'].replace('|','/').split('/').count('1')):raise ValidationError('ALT dosage disagrees with supported genotype')
                    passed.append({**r,'qc_AB':ab})
                all_rows.append(r)
                audits.append({**r,'qc_AB':ab,'qc_pass':not bad,'qc_failure_combination':combo})
        if source['carrier_annotation_records']!=len(all_rows):raise ValidationError('Carrier row count disagrees with extraction receipt')
        def report(counter):
            return dict(input_rows=counter['input'],pass_rows=counter['pass'],fail_rows=counter['fail'],
                        independent_failures={r:counter['reason:'+r] for r in REASONS},
                        failure_combinations={k.removeprefix('combination:'):v for k,v in sorted(counter.items()) if k.startswith('combination:')})
        by_type={}
        for t in TYPES:
            total=Counter()
            classes={}
            for c in CLASSES:
                tiers=sorted({'untiered'}|{tier for typ,cl,tier in counts if typ==t and cl==c})
                subtotal=Counter()
                for tier in tiers:subtotal.update(counts[t,c,tier])
                classes[c]={**report(subtotal),'by_tier':{tier:report(counts[t,c,tier]) for tier in tiers}}
                total.update(subtotal)
            by_type[t]={**report(total),'by_allele_class':classes}
        write_table(out/'carriers.qc.tsv.gz',fields+['qc_AB'],passed)
        audit_fields=KEYS+['sample','candidate_type','allele_class','tier','Gene','Feature','qc_AB','qc_pass','qc_failure_combination']
        write_table(out/'qc_audit.tsv.gz',audit_fields,audits)
        write_table(out/'samples.tsv',['sample'],(dict(sample=s) for s in samples))
        distinct=summaries(out,all_rows,passed,samples)
        if ids!={k:identity(p) for k,p in inputs.items()}:raise ValidationError('Input changed during QC')
        reconciled=(len(all_rows)==len(passed)+sum(v for k,v in combinations.items() if k!='PASS')==sum(combinations.values())
                    and all(c['input']==c['pass']+c['fail'] for c in counts.values()))
        if not reconciled:raise ValidationError('QC count reconciliation failed')
        receipt.update(status='passed',input_identities=ids,source_stage=source.get('schema_version'),samples_in_source=len(samples),
                       input_rows=len(all_rows),pass_rows=len(passed),fail_rows=len(all_rows)-len(passed),
                       independent_failures={r:independent[r] for r in REASONS},failure_combinations=dict(sorted(combinations.items())),
                       by_candidate_type=by_type,post_qc_distinct_alleles_by_class=distinct,
                       reconciliation=dict(passed=True,input_equals_pass_plus_fail=True,combinations_partition_input=True),
                       script_identity=identity(__file__),outputs={n:identity(out/n) for n in products})
    except Exception as exc:
        receipt.update(error_type=type(exc).__name__,error=str(exc) if isinstance(exc,ValidationError) else 'QC input or processing failure; inspect locally')
        for n in products:(out/n).unlink(missing_ok=True)
        raise
    finally:
        receipt['wall_seconds']=time.perf_counter()-start
        (out/'receipt.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n')
    print(json.dumps({k:receipt[k] for k in ['status','unit_id','input_rows','pass_rows','fail_rows','independent_failures','failure_combinations','by_candidate_type']}))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for n in ['carriers','samples','source_receipt','metadata','outdir']:p.add_argument('--'+n.replace('_','-'),required=True)
    try:run(p.parse_args())
    except Exception:sys.exit('Post-extraction QC failed; inspect task receipt. No protected records printed.')
