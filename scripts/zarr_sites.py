"""Export genotype-free site keys and source AC/AN/AF for the existing policy stages.

No GT reads, recalculated frequencies, variant filtering or normalization.
Multiallelic source rows emit one exact REF/ALT row per ALT, with that ALT AC/AF.
"""
import argparse,json
from pathlib import Path
import numpy as np
import pysam
import zarr
from zarr_carrier_source import value,identity
from extract_exact_carriers import file_identity,ValidationError,chrom


def export(store,chromosome,output,receipt):
    before=identity(store);g=zarr.open_group(str(store),mode='r')
    for name in ['variant_AC','variant_AN','variant_position','variant_allele','variant_contig','variant_filter','filter_id']:
        if name not in g:raise ValidationError('Missing source array '+name)
    contigs=[str(x) for x in g['contig_id'][:]];filters=[str(x) for x in g['filter_id'][:]]
    h=pysam.VariantHeader()
    for c in contigs:h.contigs.add(c)
    for f in filters:
        if f not in ('.','PASS'):h.add_meta('FILTER',items=[('ID',f),('Description','Source FILTER')])
    for key,number,kind in [('AC','A','Integer'),('AN',1,'Integer'),('AF','A','Float')]:
        if 'variant_'+key in g:h.add_meta('INFO',items=[('ID',key),('Number',number),('Type',kind),('Description','Source INFO; not corrected frequency')])
    h.add_meta('INFO',items=[('ID','ZARR_ROW'),('Number',1),('Type','Integer'),('Description','Zero-based original Zarr variant row')])
    h.add_meta('INFO',items=[('ID','ZARR_ALT_INDEX'),('Number',1),('Type','Integer'),('Description','One-based original ALT index')])
    n=g['variant_position'].shape[0];size=g['variant_position'].chunks[0];count=0
    with pysam.VariantFile(str(output),'wz',header=h) as out:
        for start in range(0,n,size):
            stop=min(start+size,n)
            names=['variant_position','variant_allele','variant_contig','variant_filter','variant_AC','variant_AN','variant_AF','variant_id','variant_quality']
            block={k:np.asarray(g[k][start:stop]) for k in names if k in g}
            for j in range(stop-start):
                contig=contigs[int(block['variant_contig'][j])]
                if chrom(contig)!=chrom(chromosome):raise ValidationError('Store contains unexpected chromosome')
                als=tuple(str(x) for x in block['variant_allele'][j] if str(x))
                if len(als)<2:raise ValidationError('Source row lacks an ALT')
                for alt_index,alt in enumerate(als[1:],1):
                    pos=int(block['variant_position'][j]);r=out.new_record(contig=contig,start=pos-1,stop=pos-1+len(als[0]),alleles=(als[0],alt))
                    if 'variant_id' in block:r.id=str(block['variant_id'][j])
                    if 'variant_quality' in block:r.qual=value(block['variant_quality'][j])
                    for f,on in zip(filters,block['variant_filter'][j]):
                        if on and f!='.':r.filter.add(f)
                    for key in ['AC','AN','AF']:
                        if 'variant_'+key not in block:continue
                        v=value(block['variant_'+key][j])
                        if key in ['AC','AF']:
                            if not isinstance(v,tuple):v=(v,)
                            if len(v)!=len(als)-1:raise ValidationError('INFO Number=A mismatch')
                            v=(v[alt_index-1],)
                        elif isinstance(v,tuple):
                            if len(v)!=1:raise ValidationError('INFO AN must be scalar')
                            v=v[0]
                        r.info[key]=v
                    r.info['ZARR_ROW']=start+j;r.info['ZARR_ALT_INDEX']=alt_index
                    out.write(r);count+=1
    if identity(store)!=before:raise ValidationError('Source metadata changed during export')
    Path(receipt).write_text(json.dumps(dict(status='passed',source=before,source_records=n,output_records=count,record_expansion="one site row per source ALT with original row/ALT pointers",
        genotype_reads=0,output_samples=0,source_AC_AN_preserved=True,sites=file_identity(output,True)),indent=2)+'\n')

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for n in ['zarr','chromosome','output','receipt']:p.add_argument('--'+n,required=True)
    a=p.parse_args();export(a.zarr,a.chromosome,a.output,a.receipt)
