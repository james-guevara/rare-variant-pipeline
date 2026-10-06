#!/usr/bin/env python3
"""Gather passed post-extraction QC products; never query genotypes or reapply QC."""
import argparse
from collections import Counter, defaultdict
import csv
import gzip
import json
from pathlib import Path
import sys
import time
from qc_filtered_carriers import CLASSES, TYPES, POLICY, KEYS, ValidationError, identity, row_key, group, write_table

PRODUCTS=['carriers.qc.tsv.gz','samples.tsv','sample_burden.tsv','sample_gene_burden.tsv',
          'gene_burden.tsv','sample_distinct_alleles.tsv','sample_gene_distinct_alleles.tsv']
REQUIRED=KEYS+['sample','candidate_type','allele_class','tier','Gene','Feature','SYMBOL',
               'GT','GQ','DP','AD','FT','site_FILTER','alt_dosage','qc_AB']


def table(path,required=None):
    with (gzip.open(path,'rt') if str(path).endswith('.gz') else open(path)) as f:
        reader=csv.DictReader(f,delimiter='\t');fields=reader.fieldnames or []
        if len(fields)!=len(set(fields)) or (required and not set(required)<=set(fields)):
            raise ValidationError('Invalid input table schema')
        rows=list(reader)
        if any(None in r or any(v is None for v in r.values()) for r in rows):
            raise ValidationError('Malformed input table')
        return fields,rows


def aggregate(out,rows,samples,strata):
    cells=Counter();dosage=Counter();distinct=defaultdict(set);gene_distinct=defaultdict(set)
    group_samples=defaultdict(set);group_alleles=defaultdict(set)
    for r in rows:
        g=group(r);key=row_key(r);cell=(r['sample'],*g,r['Gene'])
        cells[cell]+=1;dosage[cell]+=int(r['alt_dosage'])
        distinct[r['sample'],r['allele_class']].add(key[:-1])
        gene_distinct[r['sample'],r['Gene'],r['allele_class']].add(key[:-1])
        group_samples[g].add(r['sample']);group_alleles[g].add(key[:-1])
    totals=Counter();total_dosage=Counter();genes=Counter();gene_samples=defaultdict(set)
    for (s,t,c,tier,gene),n in cells.items():
        totals[s,t,c,tier]+=n;total_dosage[s,t,c,tier]+=dosage[s,t,c,tier,gene]
        genes[t,c,tier,gene]+=n;gene_samples[t,c,tier,gene].add(s)
    fields=['sample','candidate_type','allele_class','tier','carrier_variant_records','observed_alt_alleles']
    write_table(out/'sample_burden.tsv',fields,(dict(zip(fields,[s,*g,totals[s,*g],total_dosage[s,*g]])) for s in samples for g in sorted(strata)))
    fields=['sample','candidate_type','allele_class','tier','Gene','carrier_variant_records','observed_alt_alleles']
    write_table(out/'sample_gene_burden.tsv',fields,(dict(zip(fields,[*k,n,dosage[k]])) for k,n in sorted(cells.items())))
    fields=['candidate_type','allele_class','tier','Gene','carrier_variant_records','carrier_samples']
    write_table(out/'gene_burden.tsv',fields,(dict(zip(fields,[*k,n,len(gene_samples[k])])) for k,n in sorted(genes.items())))
    write_table(out/'sample_distinct_alleles.tsv',['sample','sequence_distinct_variant_count','spanning_deletion_distinct_record_count'],(
        dict(sample=s,sequence_distinct_variant_count=len(distinct[s,'sequence']),spanning_deletion_distinct_record_count=len(distinct[s,'spanning_deletion'])) for s in samples))
    fields=['sample','Gene','allele_class','distinct_allele_sample_records']
    write_table(out/'sample_gene_distinct_alleles.tsv',fields,(dict(zip(fields,[*k,len(v)])) for k,v in sorted(gene_distinct.items())))
    by_stratum=[dict(candidate_type=t,allele_class=c,tier=tier,carrier_annotation_records=sum(n for (s,typ,cl,tr),n in totals.items() if (typ,cl,tr)==(t,c,tier)),
                     samples_with_carriers=len(group_samples[t,c,tier]),distinct_alleles=len(group_alleles[t,c,tier])) for t,c,tier in sorted(strata)]
    by_class={c:dict(carrier_annotation_records=sum(r['allele_class']==c for r in rows),
                    distinct_allele_sample_records=sum(len(distinct[s,c]) for s in samples),
                    samples_with_carriers=sum(bool(distinct[s,c]) for s in samples),
                    distinct_alleles=len(set().union(*(distinct[s,c] for s in samples))),
                    genes_with_carriers=len({r['Gene'] for r in rows if r['allele_class']==c})) for c in CLASSES}
    if sum(totals.values())!=len(rows) or sum(cells.values())!=len(rows) or sum(genes.values())!=len(rows):
        raise ValidationError('Summary reconciliation failed')
    return by_stratum,by_class


def run(a):
    start=time.perf_counter();out=Path(a.outdir).resolve();meta=json.loads(Path(a.metadata).read_text())
    units=meta['units'];chrom=meta['chromosome'].removeprefix('chr')
    if not units or len({u['unit_id'] for u in units})!=len(units) or not len(units)==len(a.carriers)==len(a.samples)==len(a.receipts):
        raise ValidationError('Unique units and matching input lists required')
    paths=[Path(x).resolve() for x in a.carriers+a.samples+a.receipts]
    if any((out/n).resolve() in paths for n in PRODUCTS+['receipt.json']):raise ValidationError('Separate gather output directory required')
    out.mkdir(parents=True,exist_ok=True)
    receipt=dict(schema_version=1,stage='post_qc_gather',status='failed',chromosome='chr'+chrom,gather_id=meta['gather_id'],
                 selected_units=[u['unit_id'] for u in units],manifest_units=meta['manifest_units'],
                 container=meta.get('container'),container_identity=meta.get('container_identity'),
                 full_manifest_selected=set(meta['manifest_units'])=={u['unit_id'] for u in units},
                 policy=POLICY,duplicate_policy='Reject exact allele present in multiple selected blocks; reject within-type duplicate allele/sample associations; preserve cross-type associations within a block',
                 burden_policy='Annotation associations and distinct allele/sample counts are separate; sequence and spanning-deletion counts never combined as biological burden',
                 sample_gene_policy='Sparse tables: absent cells are zero; samples.tsv and dense sample tables retain every source sample')
    try:
        ids=[identity(p) for p in paths];lineage=[];rows=[];samples=None;fields=None;strata=set();expected=Counter()
        allele_units=defaultdict(set);association_units=defaultdict(set);seen=set();genotypes={};source_counts=Counter()
        for u,cp,sp,rp in zip(units,a.carriers,a.samples,a.receipts):
            source=json.loads(Path(rp).read_text())
            if source.get('status')!='passed' or source.get('stage')!='post_extraction_qc' or source.get('unit_id')!=u['unit_id'] or source.get('policy')!=POLICY or not source.get('reconciliation',{}).get('passed'):
                raise ValidationError('Passed matching post-QC receipt and fixed policy required')
            if source['sources']['chromosome'].removeprefix('chr')!=chrom:raise ValidationError('Receipt chromosome mismatch')
            for p,name in [(cp,'carriers.qc.tsv.gz'),(sp,'samples.tsv')]:
                if identity(p)!= {k:source['outputs'][name][k] for k in ['bytes','sha256']}:raise ValidationError('QC product differs from its receipt')
            sf,sr=table(sp,['sample']);ss=[r['sample'] for r in sr]
            if sf!=['sample'] or not ss or len(ss)!=len(set(ss)) or any(s in ('','.') for s in ss):raise ValidationError('Invalid sample universe')
            if samples is not None and samples!=ss:raise ValidationError('Sample universe/order differs across blocks')
            samples=ss;sample_set=set(ss)
            if source['samples_in_source']!=len(ss):raise ValidationError('Sample count mismatch')
            ff,rr=table(cp,REQUIRED)
            if 'source_unit_id' in ff or (fields is not None and fields!=ff):raise ValidationError('Carrier schemas differ across blocks')
            fields=ff
            if len(rr)!=source['pass_rows'] or source['input_rows']!=source['pass_rows']+source['fail_rows']:raise ValidationError('Source row reconciliation failed')
            observed=Counter()
            for r in rr:
                key=row_key(r);g=group(r);typ,cl,tier=g
                if key[0]!=chrom or key[1]<1 or r['sample'] not in sample_set:raise ValidationError('Invalid carrier locus/sample')
                if typ not in TYPES or cl not in CLASSES or cl!=('spanning_deletion' if r['ALT']=='*' else 'sequence'):raise ValidationError('Invalid type or allele class')
                if tier not in (['miss_t1','miss_t2','miss_t3','miss_t4','untiered'] if typ=='missense' else ['lof_t1','lof_t2','untiered']):raise ValidationError('Invalid tier')
                if r['alt_dosage'] not in ('1','2'):raise ValidationError('Invalid ALT dosage')
                association=(key,typ)
                if (u['unit_id'],association) in seen:raise ValidationError('Duplicate within-unit annotation association')
                seen.add((u['unit_id'],association));allele_units[key[:-1]].add(u['unit_id']);association_units[association].add(u['unit_id'])
                signature=tuple(r[k] for k in ['GT','GQ','DP','AD','FT','site_FILTER','alt_dosage','qc_AB'])
                if key in genotypes and signature!=genotypes[key]:raise ValidationError('Conflicting cross-type genotype associations')
                genotypes[key]=signature;observed[g]+=1
                rows.append(dict(r,CHROM='chr'+chrom,source_unit_id=u['unit_id']))
            declared=Counter()
            for t,td in source['by_candidate_type'].items():
                for c,cd in td['by_allele_class'].items():
                    for tier,counts in cd['by_tier'].items():declared[t,c,tier]=counts['pass_rows'];strata.add((t,c,tier))
            if +observed!=+declared:raise ValidationError('Stratum counts differ from source receipt')
            expected.update(declared)
            for k in ['input_rows','pass_rows','fail_rows']:source_counts[k]+=source[k]
            lineage.append(dict(unit=u,receipt_identity=identity(rp),input_rows=source['input_rows'],pass_rows=source['pass_rows'],fail_rows=source['fail_rows']))
        receipt['duplicate_checks']=dict(cross_block_exact_alleles=sum(len(v)>1 for v in allele_units.values()),cross_block_annotation_associations=sum(len(v)>1 for v in association_units.values()))
        if any(receipt['duplicate_checks'].values()):raise ValidationError('Cross-block duplicates detected; no burdens published')
        rows.sort(key=lambda r:(row_key(r),group(r),r['Gene'],r['Feature']))
        write_table(out/'carriers.qc.tsv.gz',fields+['source_unit_id'],rows)
        write_table(out/'samples.tsv',['sample'],(dict(sample=s) for s in samples))
        strata_report,classes=aggregate(out,rows,samples,strata)
        if len(rows)!=source_counts['pass_rows'] or sum(expected.values())!=len(rows):raise ValidationError('Gather count reconciliation failed')
        if ids!=[identity(p) for p in paths]:raise ValidationError('Inputs changed during gathering')
        receipt.update(status='passed',samples_in_source=len(samples),carrier_annotation_records=len(rows),source_counts=dict(source_counts),
                       by_stratum=strata_report,by_allele_class=classes,lineage=lineage,input_identities=ids,
                       reconciliation=dict(passed=True,source_pass_equals_output=True,stratum_counts_match=True,summary_associations_match=True),
                       script_identity=identity(__file__),helper_identity=identity(Path(__file__).with_name('qc_filtered_carriers.py')),
                       outputs={n:identity(out/n) for n in PRODUCTS})
    except Exception as e:
        receipt.update(error_type=type(e).__name__,error=str(e) if isinstance(e,ValidationError) else 'Gather failed; inspect inputs locally')
        for n in PRODUCTS:(out/n).unlink(missing_ok=True)
        raise
    finally:
        receipt['wall_seconds']=time.perf_counter()-start
        (out/'receipt.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n')
    print(json.dumps({k:receipt[k] for k in ['status','chromosome','samples_in_source','carrier_annotation_records','by_allele_class','duplicate_checks']}))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['carriers','samples','receipts']:p.add_argument('--'+name,nargs='+',required=True)
    for name in ['metadata','outdir']:p.add_argument('--'+name,required=True)
    try:run(p.parse_args())
    except Exception:sys.exit('Post-QC gathering failed; inspect receipt. No protected records printed.')
