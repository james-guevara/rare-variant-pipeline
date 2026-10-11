"""Read-only VCZ backend for the shared filtered-carrier policy.

Each candidate ALT is projected to a biallelic view; source GT/AD are retained. -1 is a
missing allele and -2 is a padded (absent) allele, not an additional missing call.
"""
from collections.abc import Mapping
from types import SimpleNamespace
from pathlib import Path
import numpy as np
import zarr
from extract_exact_carriers import ValidationError, chrom, carrier_calls, file_identity


def identity(path):
    root=Path(path)
    if not (root/'zarr.json').is_file():
        raise ValidationError('Expected a Zarr v3 group')
    # Metadata identity, not an exhaustive value-integrity claim.
    files=sorted(root.glob('*/zarr.json'))+[root/'zarr.json']
    return {'path':str(root.resolve()),'identity_method':'array metadata hashes; source must remain immutable',
            'metadata':{str(p.relative_to(root)):file_identity(p,True) for p in files}}


def value(raw):
    a=np.asarray(raw)
    if a.ndim:
        return tuple(value(x) for x in a if not (np.issubdtype(a.dtype,np.integer) and x == -2))
    v=a.item()
    if isinstance(v,bytes):v=v.decode()
    if isinstance(v,float) and not np.isfinite(v):return None
    if isinstance(v,int) and v < 0:return None
    return v


def allele_info(raw, n_alt):
    """Decode Number=A by allele count; only trailing VCZ fill is padding."""
    a=np.asarray(raw).reshape(-1)
    if len(a)<n_alt:
        raise ValidationError('Source INFO Number=A shorter than ALT count')
    tail=a[n_alt:]
    if np.issubdtype(a.dtype,np.integer):
        valid=np.all(tail == -2)
        if np.any(a[:n_alt] == -2):
            raise ValidationError('Source INFO Number=A padding inside active ALTs')
    elif np.issubdtype(a.dtype,np.floating):
        valid=np.all(np.isnan(tail))
    else:
        valid=not len(tail)
    if not valid:
        raise ValidationError('Source INFO Number=A has non-padding values beyond ALT count')
    return tuple(value(x) for x in a[:n_alt])


class Calls(Mapping):
    def __init__(self,record):self.record=record
    def __iter__(self):return iter(self.record.source.samples)
    def __len__(self):return len(self.record.source.samples)
    def __getitem__(self,sample):
        r=self.record;i=r.source.sample_index[sample];raw=r.gt[i]
        gt=tuple(None if x == -1 else int(x) for x in raw if x != -2)
        fields={'GT':gt}
        for name in ['GQ','DP','AD','FT']:
            if 'call_'+name in r.block:
                v=value(r.block['call_'+name][r.offset,i])
                if name=='AD':
                    v=(v[0],v[r.alt_index]) if isinstance(v,tuple) and len(v)>r.alt_index else None
                fields[name]=v
        result=Call(fields)
        phased=r.block.get('call_genotype_phased')
        result.phased=bool(phased[r.offset,i]) if phased is not None else False
        return result


class Call(dict):pass


class Source:
    def __init__(self,path,chromosome):
        self.root=zarr.open_group(str(path),mode='r');self.chromosome=chrom(chromosome)
        required=['sample_id','variant_position','variant_allele','variant_contig','contig_id','call_genotype','call_genotype_mask','variant_filter','filter_id']
        if any(k not in self.root for k in required):raise ValidationError('Zarr lacks required genotype/site arrays')
        self.samples=[str(x) for x in self.root['sample_id'][:]]
        if not self.samples or len(set(self.samples))!=len(self.samples):raise ValidationError('Expected unique source samples')
        self.sample_index={s:i for i,s in enumerate(self.samples)}
        contigs=[str(x) for x in self.root['contig_id'][:]]
        aliases=[i for i,c in enumerate(contigs) if chrom(c)==self.chromosome]
        if len(aliases)!=1:raise ValidationError('Ambiguous source chromosome')
        self.contig=contigs[aliases[0]];self.contig_index=aliases[0]
        self.filters=[str(x) for x in self.root['filter_id'][:]]
        self.chunk_reads=0
    def __enter__(self):return self
    def __exit__(self,*args):pass
    def records(self,selected):
        positions={k[1] for k in selected};seen=set();g=self.root
        width=g['call_genotype'].chunks[0]
        n=g['variant_position'].shape[0]
        for start in range(0,n,width):
            stop=min(start+width,n);pos=np.asarray(g['variant_position'][start:stop]);contigs=np.asarray(g['variant_contig'][start:stop])
            offsets=np.flatnonzero(np.isin(pos,list(positions)) & (contigs==self.contig_index))
            if not len(offsets):continue
            alleles=g['variant_allele'][start:stop];matches=[]
            for j in offsets:
                als=[str(x) for x in alleles[j] if str(x)]
                if len(als)<2:raise ValidationError('Source record lacks ALT')
                for alt_index,alt in enumerate(als[1:],1):
                    key=(self.chromosome,int(pos[j]),als[0],alt)
                    if key not in selected:continue
                    if key in seen:raise ValidationError('Duplicate exact source allele record')
                    seen.add(key);matches.append((j,key,alt_index,len(als)))
            if not matches:continue
            fields=['call_genotype','call_genotype_mask','call_genotype_phased','call_GQ','call_DP','call_AD','call_FT']
            block={f:np.asarray(g[f][start:stop]) for f in fields if f in g};self.chunk_reads+=1
            for j,key,alt_index,n_alleles in matches:
                gt=block['call_genotype'][j];mask=block['call_genotype_mask'][j]
                if gt.shape!=mask.shape or not np.array_equal(mask,gt<0):raise ValidationError('GT and missingness mask disagree')
                if np.any(gt < -2) or np.any(gt >= n_alleles):raise ValidationError('Non-biallelic genotype allele index')
                # Padding must be trailing; preserve haploid versus partial calls.
                if gt.shape[1]>1 and np.any((gt[:,:-1]==-2)&(gt[:,1:]!=-2)):raise ValidationError('Non-trailing genotype padding')
                info={name:(allele_info(g['variant_'+name][start+j],n_alleles-1) if name in ['AC','AF'] else value(g['variant_'+name][start+j])) for name in ['AC','AN','AF'] if 'variant_'+name in g}
                for name in ['AC','AF']:
                    if name in info and isinstance(info[name],tuple):
                        if len(info[name])!=n_alleles-1:raise ValidationError('Source INFO Number=A mismatch')
                        info[name]=(info[name][alt_index-1],)
                fl=np.asarray(g['variant_filter'][start+j])
                projected=np.where(gt<0,gt,(gt==alt_index).astype(gt.dtype))
                r=SimpleNamespace(source=self,offset=j,block=block,gt=projected,source_gt=gt,alt_index=alt_index,variant_index=start+j,contig=self.contig,pos=key[1],ref=key[2],alts=(key[3],),
                    info=info,header=SimpleNamespace(info=info),filter=[f for f,on in zip(self.filters,fl) if on])
                r.samples=Calls(r)
                yield key,r
    def calls(self,record):
        # Materialize per-sample dictionaries only for ALT carriers.
        indexes=np.flatnonzero(np.any(record.gt==1,axis=1))
        samples={self.samples[i]:record.samples[self.samples[i]] for i in indexes}
        sparse=SimpleNamespace(**{**vars(record),'samples':samples})
        for row,partial in carrier_calls(sparse):
            i=self.sample_index[row['sample']]
            sep='|' if samples[row['sample']].phased else '/'
            row['source_GT']=sep.join('.' if x==-1 else str(int(x)) for x in record.source_gt[i] if x!=-2)
            raw_ad=record.block.get('call_AD')
            from extract_exact_carriers import text
            row['source_AD']=text(value(raw_ad[record.offset,i])) if raw_ad is not None else '.'
            row['source_variant_index']=record.variant_index
            row['source_alt_index']=record.alt_index
            yield row,partial
