#!/usr/bin/env python3
"""Post-rarity receipt adapter and explicit sex/PAR-aware genotype QC."""
import argparse
from collections import Counter
from pathlib import Path
import re
import sys
import qc_filtered_carriers as core
from carrier_psam import load_psam
from carrier_frequencies import (SEX_POLICY, PAR, normalize_sex, expected_ploidy,
                                 frequency_region, validate_chromosome)
from filter_final_rarity import POLICY as RARITY_POLICY

REASONS=core.REASONS+('unknown_sex','female_y','unexpected_ploidy','partial_call','haploid_heterozygous')
EXTRA_FIELDS=['qc_expected_ploidy','qc_frequency_region','qc_sex','qc_effective_alt_dosage']
POLICY=dict(core.POLICY,version='post-rarity-grch38-x-only-par-qc-v1',
    upstream='passed final_rarity; saved unrelated_af <0.001 and unrelated_an >0',
    sex_chromosome_policy=SEX_POLICY,par_coordinates_1based_inclusive={k:[list(v) for v in intervals] for k,intervals in PAR.items()},
    sex_encoding='PSAM SEX 1=male, 2=female, other nonempty values=unknown',
    ploidy='autosomes diploid; X PAR diploid for known sexes; X non-PAR male haploid/female diploid; Y non-PAR male haploid/female excluded; unknown sex excluded on X/Y',
    haploid_alt='literal GT=1 or diploid 1/1 (also phased), AB >=0.90, effective ALT dosage 1',
    exclusions='partial calls, unexpected ploidy and haploid-region heterozygotes fail; Y-PAR candidates fail the stage',
    dosage='raw GT and alt_dosage preserved; observed_alt_alleles summaries use qc_effective_alt_dosage',
    sample_selection='QC all source samples regardless of representative/unrelated flags; retain zeros')


def policy_for_run(a):
    mode=getattr(a,'site_filter_policy','pass')
    if mode not in ('pass','pass_or_missing'):
        raise core.ValidationError('Unsupported site FILTER policy')
    if mode=='pass':return POLICY
    return dict(POLICY,version='post-rarity-grch38-x-only-par-qc-v2',
                site_FILTER='PASS or .',site_filter_policy=mode)


def validate_source(source,meta):
    if source.get('stage')!='final_rarity' or source.get('policy')!=RARITY_POLICY or not source.get('reconciliation',{}).get('passed'):
        raise core.ValidationError('Passed final-rarity receipt with fixed rarity policy required')
    if source.get('chromosome','').removeprefix('chr')!=meta['chromosome'].removeprefix('chr'):
        raise core.ValidationError('Final-rarity chromosome differs from manifest')
    for name in ['distinct_alleles','candidate_annotations','carrier_annotations']:
        c=source[name]
        if c['input']!=c['pass_count']+c['fail_count']:
            raise core.ValidationError('Final-rarity source counts do not reconcile')


class Context:
    def __init__(self,a,source,samples,meta):
        self.site_filter_policy=getattr(a,'site_filter_policy','pass')
        policy_for_run(a)
        self.chromosome=meta['chromosome']
        policy=getattr(a,'sex_chromosome_policy',None)
        try:
            validate_chromosome(self.chromosome,policy)
            _,rows,_=load_psam(a.psam,samples)
        except ValueError as exc:raise core.ValidationError(str(exc)) from exc
        self.sex={r['#IID']:normalize_sex(r['SEX']) for r in rows}
        if self.chromosome.removeprefix('chr') in ('X','Y') and source['source_frequency_policy'].get('sex_chromosomes')!=SEX_POLICY:
            raise core.ValidationError('Final-rarity source does not declare GRCh38 X-only PAR frequencies')
        self.psam=a.psam

    def evaluate(self,row):
        try:region=frequency_region(self.chromosome,int(row['POS']))
        except ValueError as exc:raise core.ValidationError(str(exc)) from exc
        sex=self.sex[row['sample']];ploidy,excluded=expected_ploidy(region,sex)
        ab,bad=core.evaluate(row);bad=list(bad);ploidy_bad=[]
        if getattr(self,'site_filter_policy','pass')=='pass_or_missing' and row['site_FILTER']=='.':
            bad=[x for x in bad if x!='site_filter']
        if excluded:ploidy_bad.append('unknown_sex' if excluded=='excluded_unknown_sex' else 'female_y')
        gt=row['GT'];alleles=re.split(r'[/|]',gt)
        if '.' in alleles:ploidy_bad.append('partial_call')
        if ploidy==2 and len(alleles)!=2:ploidy_bad.append('unexpected_ploidy')
        if ploidy==1:
            if len(alleles) not in (1,2):ploidy_bad.append('unexpected_ploidy')
            if len(alleles)==2 and set(alleles)=={'0','1'}:ploidy_bad.append('haploid_heterozygous')
        # AB interpretation is inapplicable to sex/ploidy-ineligible calls.
        if ploidy_bad:bad=[x for x in bad if x!='ab_out_of_range']
        bad+=ploidy_bad
        effective=None if bad else 1 if ploidy==1 else alleles.count('1')
        return ab,tuple(bad),dict(qc_expected_ploidy=ploidy,qc_frequency_region=region,
            qc_sex=sex,qc_effective_alt_dosage=effective)

    def receipt(self,passed):
        dosage={c:dict(raw_observed_alt_alleles=sum(int(r['alt_dosage']) for r in passed if r['allele_class']==c),
                      effective_observed_alt_alleles=sum(r['qc_effective_alt_dosage'] for r in passed if r['allele_class']==c)) for c in core.CLASSES}
        return dict(upstream_stage='final_rarity',chromosome=self.chromosome,psam_identity=core.identity(self.psam),
            source_sample_sex_counts=dict(Counter(self.sex.values())),effective_dosage_by_class=dosage,
            dosage_scope='carrier annotation associations; cross-type totals are not distinct allele/sample dosage',
            summary_dosage_field='qc_effective_alt_dosage',
            adapter_identity=core.identity(__file__),
            helper_identities={n:core.identity(Path(__file__).with_name(n)) for n in
                ['carrier_frequencies.py','carrier_psam.py','extract_exact_carriers.py','filter_final_rarity.py']},
            downstream_gather='Legacy gather_post_qc rejects this stage/policy. A post-rarity adapter must validate this policy and sum qc_effective_alt_dosage.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for n in ['carriers','samples','source-receipt','psam','metadata','outdir']:p.add_argument('--'+n,required=True)
    p.add_argument('--site-filter-policy',choices=['pass','pass_or_missing'],default='pass')
    p.add_argument('--sex-chromosome-policy',choices=[SEX_POLICY])
    try:core.run(p.parse_args(),adapter=sys.modules[__name__])
    except Exception:sys.exit('Post-rarity QC failed; inspect local receipt. No protected records printed.')
