"""Vectorized autosomal counts with the shared VCF policy and receipts."""
from collections import Counter
import numpy as np
from carrier_frequencies import FrequencyCounter, SETS, frequency_region, expected_ploidy
from types import SimpleNamespace


class ZarrFrequencyCounter(FrequencyCounter):
    def __init__(self,metadata,chromosome,samples,sex_policy=None):
        super().__init__(metadata,chromosome,sex_policy)
        lookup={s:i for i,s in enumerate(samples)}
        self.sample_lookup=lookup
        self.indexes={group:np.asarray([lookup[s] for s,u,_ in self.samples if group=='cohort' or u],dtype=np.int64) for group in SETS}

    def count(self,record,allele_class):
        if self.sex_chromosome:
            # Same explicit PAR/ploidy implementation. No pysam genotype reads.
            region=frequency_region(self.chromosome,record.pos)
            calls={}
            for sample,_,sex in self.samples:
                call=dict(record.samples[sample])
                raw=record.source_gt[self.sample_lookup[sample]]
                expected,_=expected_ploidy(region,sex)
                # A multiallelic heterozygote remains invalid in a haploid region,
                # even when both alleles are non-target ALTs projecting to 0/0.
                if expected==1 and np.all(raw>=0) and len(set(raw.tolist()))>1:
                    call['GT']=(0,1)
                calls[sample]=call
            return super().count(SimpleNamespace(pos=record.pos,samples=calls),allele_class)
        result={};audit={}
        for group in SETS:
            gt=record.gt[self.indexes[group]]
            length=np.sum(gt!=-2,axis=1)
            invalid=np.any((gt < -2)|(gt > 1),axis=1)
            called=np.sum(gt>=0,axis=1)
            complete=(~invalid)&(length==2)&(called==2)
            partial=(~invalid)&(length==2)&(called==1)
            eligible=complete|partial
            acs=np.sum(gt==1,axis=1)*eligible
            ans=called*eligible
            reasons=Counter()
            reasons['excluded_missing']=int(np.sum((length==0)|((~invalid)&(length==2)&(called==0))))
            reasons['excluded_invalid_allele']=int(np.sum(invalid & (length>0)))
            reasons['excluded_unexpected_ploidy']=int(np.sum((~invalid)&(length!=2)&(length>0)))
            reasons['counted_diploid_complete']=int(np.sum(complete))
            reasons['counted_diploid_partial']=int(np.sum(partial))
            reasons=Counter({k:v for k,v in reasons.items() if v})
            ac=int(acs.sum());an=int(ans.sum());counted=int(eligible.sum())
            value=dict(ac=ac,an=an,af=ac/an if an else None,counted_genotypes=counted,
                       reference_genotypes=int(np.sum(eligible & (acs==0))),excluded_genotypes=self.sizes[group]-counted)
            assert sum(reasons.values())==self.sizes[group]
            result.update({group+'_'+k:v for k,v in value.items()});audit[group]=reasons
            for reason,n in reasons.items():self.audit[allele_class,group,reason]+=n
            self.zero_an[allele_class,group]+=an==0
        self.matched[allele_class]+=1
        return result,audit
