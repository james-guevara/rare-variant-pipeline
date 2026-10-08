#!/usr/bin/env python3
"""Filtered-candidate adapter using the existing indexed exact carrier core."""
import argparse
from collections import Counter, defaultdict
import csv
import json
from pathlib import Path
import re
import sys
import time
import duckdb
import pysam
from carrier_psam import load_psam
from extract_exact_carriers import (ValidationError, KEY_FIELDS, FIELDS, ALLELE_CLASSES,
                                   chrom, text, file_identity, source_header, exact_records, carrier_calls)

TYPES=('missense','lof_hc')
ANNOTATION_FIELDS=KEY_FIELDS+['candidate_type','tier','Gene','Feature','SYMBOL','allele_class']
CARRIER_FIELDS=FIELDS+['candidate_type','tier']


def load_candidates(paths, chromosome):
    selected=defaultdict(list)
    con=duckdb.connect()
    try:
        for kind,path in paths.items():
            cursor=con.execute('SELECT * FROM read_parquet(?)',[str(path)])
            names=[d[0] for d in cursor.description]
            required=set(KEY_FIELDS+['Gene','Feature','allele_class','tier','pcf_retained','Consequence'])
            if kind=='lof_hc':required.add('LoF')
            if not required<=set(names):raise ValidationError('Filtered candidate Parquet lacks required columns')
            seen=set()
            for values in cursor.fetchall():
                row=dict(zip(names,values))
                if row['pcf_retained'] is not True:raise ValidationError('Input includes candidates not retained by pre-carrier filtering')
                if not all(isinstance(row[k],str) and row[k] not in ('','.','-') for k in ['CHROM','REF','ALT','Gene','Feature']):
                    raise ValidationError('Missing candidate allele/gene/transcript identity')
                if chrom(row['CHROM'])!=chrom(chromosome):raise ValidationError('Candidate chromosome disagrees with manifest')
                if not re.fullmatch('[ACGTN]+',row['REF']) or not (row['ALT']=='*' or re.fullmatch('[ACGTN]+',row['ALT'])) or row['REF']==row['ALT']:
                    raise ValidationError('Candidate must have one sequence ALT or spanning-deletion *')
                if not re.fullmatch('[0-9]+',str(row['POS'])) or int(row['POS'])<1:raise ValidationError('Invalid candidate position')
                row['POS']=int(row['POS'])
                expected='spanning_deletion' if row['ALT']=='*' else 'sequence'
                if row['allele_class']!=expected:raise ValidationError('Candidate allele_class disagrees with ALT')
                if kind=='lof_hc' and row['LoF']!='HC':raise ValidationError('Non-HC row in filtered HC input')
                if kind=='missense' and 'missense_variant' not in str(row['Consequence']).split('&'):
                    raise ValidationError('Non-missense row in filtered missense input')
                tier=row['tier']
                if tier in (None,'','.','untiered'):tier='untiered'
                allowed=['miss_t'+str(i) for i in range(1,5)] if kind=='missense' else ['lof_t1','lof_t2']
                if tier not in allowed+['untiered']:raise ValidationError('Tier incompatible with candidate type')
                if 'candidate_type' in row and row['candidate_type']!=kind:raise ValidationError('Conflicting candidate type')
                row.update(candidate_type=kind,tier=tier)
                key=(chrom(row['CHROM']),row['POS'],row['REF'],row['ALT'])
                if key in seen:raise ValidationError('Duplicate exact allele within candidate type')
                seen.add(key);selected[key].append(row)
    finally:con.close()
    return selected


def write_rows(path,fields,rows):
    with open(path,'w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n',extrasaction='ignore');w.writeheader()
        for row in rows:w.writerow({k:text(row.get(k)) for k in fields})


def extract(a):
    start=time.perf_counter();meta=json.loads(Path(a.metadata).read_text());out=Path(a.outdir);out.mkdir(parents=True,exist_ok=True)
    products=['carriers.tsv.gz','candidates.tsv','unmatched.tsv','samples.tsv','sample_burden.tsv','sample_gene_burden.tsv','gene_burden.tsv','sample_distinct_alleles.tsv','source_frequencies.tsv']
    receipt=dict(schema_version=1,status='failed',unit_id=meta['unit_id'],chromosome=meta['chromosome'],sources=meta,
        pysam_version=pysam.__version__,duckdb_version=duckdb.__version__,
        definition='exact allele lookup; one carrier association per candidate type and variant/sample; ALT index 1; no genotype QC',
        overlap_policy='retain both annotations for cross-type exact-key overlaps; type counts are not additive; distinct-allele summary counts each allele/sample once per allele class',
        spanning_deletion_policy='star record audits are separate from sequence variant counts in all summaries including tiers; no deletion-event mapping or combined biological burden',
        genotype_policy='preserve GT/GQ/DP/AD/FT/site FILTER and ALT dosage; include partial and haploid ALT calls; homozygous ALT is one carrier variant',
        vcf_identity_method='path/size/mtime and index SHA-256; no full genotype VCF hash')
    psam=getattr(a,'psam',None)
    if psam: products.append('sample_metadata.tsv')
    receipt['frequency_status']='deferred; source frequencies are uncorrected; no final rarity decision'
    receipt['query_batch_bp']=getattr(a,'batch_bp',10000)
    try:
        paths={k:getattr(a,k) for k in TYPES}
        ids={k:file_identity(getattr(a,k),k!='vcf') for k in ['vcf','index',*TYPES]}
        if psam: ids['psam']=file_identity(psam,True)
        selected=load_candidates(paths,meta['chromosome'])
        for kind,expected in [('missense',a.expected_missense),('lof_hc',a.expected_hc)]:
            if expected is not None and sum(r['candidate_type']==kind for rs in selected.values() for r in rs)!=expected:
                raise ValidationError('Candidate count differs from expected count for '+kind)
        source_frequencies={}
        found=set();carrier_keys=set();records=Counter();dosages=Counter();partials=Counter()
        burdens=Counter();burden_dosages=Counter();distinct=Counter();samples_by_group=defaultdict(set)
        carrier_counts=Counter()
        annotations=[(k,r) for k,rs in selected.items() for r in rs]
        def group(row):return (row['candidate_type'],row['allele_class'],row['tier'])
        strata=sorted({group(r) for k,r in annotations} | {(t,c,'untiered') for t in TYPES for c in ALLELE_CLASSES})
        with pysam.VariantFile(a.vcf,index_filename=a.index) as vcf:
            contig,samples=source_header(vcf,meta['chromosome'])
            if psam:
                header,metadata,extra=load_psam(psam,samples)
                write_rows(out/'sample_metadata.tsv',header,metadata)
                receipt['psam']=dict(mode='metadata only; no sample exclusion or frequency calculation',
                    source_samples_annotated=len(metadata),extra_psam_samples=extra,
                    participant_count=len({r['participant_id'] for r in metadata}),
                    representative_flags=sum(r['frequency_representative']=='1' for r in metadata),
                    unrelated_flags=sum(r['unrelated']=='1' for r in metadata),
                    representative_selection_status='not applied; policy deferred')
            with pysam.BGZFile(str(out/'carriers.tsv.gz'),'w') as target:
                target.write(('\t'.join(CARRIER_FIELDS)+'\n').encode())
                for key,record in exact_records(vcf,contig,selected,batch_bp=getattr(a,'batch_bp',10000)):
                    found.add(key)
                    info={name:record.info.get(name) for name in ['AC','AN','AF'] if name in record.header.info}
                    ac,an=info.get('AC'),info.get('AN')
                    if isinstance(ac,tuple): ac=ac[0] if len(ac)==1 else None
                    if isinstance(an,tuple): an=an[0] if len(an)==1 else None
                    valid=isinstance(ac,int) and isinstance(an,int) and an>0 and 0<=ac<=an
                    source_frequencies[key]=dict(source_info_ac=text(info.get('AC')),
                        source_info_an=text(info.get('AN')),source_info_af=text(info.get('AF')),
                        source_ac_an=ac/an if valid else None)
                    for call,is_partial in carrier_calls(record):
                        carrier_keys.add(key);carrier_counts[key]+=1
                        kind=selected[key][0]['allele_class']
                        distinct[call['sample'],kind]+=1
                        for annotation in selected[key]:
                            g=group(annotation);row={**call,**{k:annotation.get(k,'.') for k in ['Gene','Feature','SYMBOL','allele_class','candidate_type','tier']}}
                            target.write(('\t'.join(text(row[k]) for k in CARRIER_FIELDS)+'\n').encode())
                            records[g]+=1;dosages[g]+=call['alt_dosage'];partials[g]+=is_partial
                            samples_by_group[g].add(call['sample'])
                            b=(call['sample'],*g,annotation['Gene'])
                            burdens[b]+=1;burden_dosages[b]+=call['alt_dosage']
        frequency_fields=KEY_FIELDS+['allele_class','matched','source_info_ac','source_info_an','source_info_af','source_ac_an']
        write_rows(out/'source_frequencies.tsv',frequency_fields,(
            {**dict(zip(KEY_FIELDS,key)),'allele_class':selected[key][0]['allele_class'],
             'matched':key in found,**source_frequencies.get(key,{})} for key in sorted(selected)))
        audit=[{**r,'matched':k in found,'carrier_records':carrier_counts[k]} for k,r in sorted(annotations,key=lambda x:(x[0],x[1]['candidate_type']))]
        write_rows(out/'candidates.tsv',ANNOTATION_FIELDS+['matched','carrier_records'],audit)
        write_rows(out/'unmatched.tsv',ANNOTATION_FIELDS,[r for r in audit if not r['matched']])
        write_rows(out/'samples.tsv',['sample'],[dict(sample=s) for s in samples])
        sample_totals=Counter();sample_dosages=Counter();gene_records=Counter();gene_samples=Counter()
        for (sample,t,c,tier,gene),n in burdens.items():
            sample_totals[sample,t,c,tier]+=n;sample_dosages[sample,t,c,tier]+=burden_dosages[sample,t,c,tier,gene]
            gene_records[t,c,tier,gene]+=n;gene_samples[t,c,tier,gene]+=1
        fields=['sample','candidate_type','allele_class','tier','carrier_variant_records','observed_alt_alleles']
        write_rows(out/'sample_burden.tsv',fields,(
            dict(zip(fields,[s,*g,sample_totals[(s,*g)],sample_dosages[(s,*g)]])) for s in samples for g in strata))
        fields=['sample','candidate_type','allele_class','tier','Gene','carrier_variant_records','observed_alt_alleles']
        write_rows(out/'sample_gene_burden.tsv',fields,(dict(zip(fields,[*key,n,burden_dosages[key]])) for key,n in sorted(burdens.items())))
        fields=['candidate_type','allele_class','tier','Gene','carrier_variant_records','carrier_samples']
        gene_groups=sorted({(*group(r),r['Gene']) for k,r in annotations})
        write_rows(out/'gene_burden.tsv',fields,(dict(zip(fields,[*g,gene_records[g],gene_samples[g]])) for g in gene_groups))
        write_rows(out/'sample_distinct_alleles.tsv',['sample','sequence_distinct_variant_count','spanning_deletion_distinct_record_count'],(
            dict(sample=s,sequence_distinct_variant_count=distinct[s,'sequence'],spanning_deletion_distinct_record_count=distinct[s,'spanning_deletion']) for s in samples))
        def summarize(subset):
            keys={k for k,r in subset};groups={group(r) for k,r in subset}
            return dict(candidate_records=len(subset),matched_candidates=len(keys&found),unmatched_candidates=len(keys-found),
                        matched_candidates_without_carriers=len((keys&found)-carrier_keys),carrier_records=sum(records[g] for g in groups),
                        samples_with_carriers=len(set().union(*(samples_by_group[g] for g in groups))),
                        observed_alt_alleles=sum(dosages[g] for g in groups),partial_call_carrier_records=sum(partials[g] for g in groups))
        by_type={}
        for t in TYPES:
            sub=[(k,r) for k,r in annotations if r['candidate_type']==t]
            by_type[t]={**summarize(sub),'by_allele_class':{}}
            for c in ALLELE_CLASSES:
                cls=[(k,r) for k,r in sub if r['allele_class']==c]
                by_type[t]['by_allele_class'][c]={**summarize(cls),'by_tier':{tier:summarize([(k,r) for k,r in cls if r['tier']==tier]) for tier in sorted({'untiered'}|{r['tier'] for k,r in cls})}}
        for name,old in ids.items():
            current=file_identity(getattr(a,name),name!='vcf')
            if current!=old:raise ValidationError('Input changed during extraction')
        receipt.update(status='passed',input_identities=ids,samples_in_source=len(samples),by_candidate_type=by_type,
            candidate_annotation_records=len(annotations),distinct_candidate_alleles=len(selected),overlapping_type_alleles=sum(len(rs)>1 for rs in selected.values()),
            carrier_annotation_records=sum(records.values()),
            distinct_allele_audit_by_class={c:dict(candidate_alleles=sum(rs[0]['allele_class']==c for rs in selected.values()),
                carrier_variant_sample_records=sum(n for (s,cl),n in distinct.items() if cl==c)) for c in ALLELE_CLASSES},
            contig_alias=dict(annotation=meta['chromosome'],source=contig),
            code_identities={Path(p).name:file_identity(p,True) for p in [__file__,Path(__file__).with_name('extract_exact_carriers.py'),Path(__file__).with_name('carrier_psam.py')]},
            outputs={name:file_identity(out/name,True) for name in products})
    except Exception as exc:
        receipt.update(error_type=type(exc).__name__,error=str(exc) if isinstance(exc,ValidationError) else 'Input or extraction failure; inspect locally')
        for name in products:(out/name).unlink(missing_ok=True)
        raise
    finally:
        receipt['wall_seconds']=time.perf_counter()-start
        (out/'receipt.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n')
    print(json.dumps({k:receipt[k] for k in ['status','unit_id','samples_in_source','overlapping_type_alleles','by_candidate_type','distinct_allele_audit_by_class']}))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['metadata','missense','lof_hc','vcf','index','outdir']:p.add_argument('--'+name.replace('_','-'),required=True)
    p.add_argument('--psam')
    p.add_argument('--batch-bp',type=int,default=10000)
    p.add_argument('--expected-missense',type=int);p.add_argument('--expected-hc',type=int)
    try:extract(p.parse_args())
    except Exception:sys.exit('Filtered carrier extraction failed; inspect task receipt. No individual records printed.')
