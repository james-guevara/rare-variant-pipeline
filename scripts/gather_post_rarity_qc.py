#!/usr/bin/env python3
"""Gather validated post-rarity QC without repeating QC or reading genotype VCFs."""
import argparse
from pathlib import Path
import re
import sys
import gather_post_qc as core
from qc_post_rarity import POLICY, EXTRA_FIELDS
from carrier_frequencies import frequency_region, expected_ploidy


class Context:
    def __init__(self):self.psam=None

    def validate_source(self,source,meta):
        p=source.get('psam_identity',{})
        if (set(p)!={'bytes','sha256'} or type(p['bytes']) is not int or p['bytes']<=0
                or not re.fullmatch('[0-9a-f]{64}',str(p['sha256']))):
            raise core.ValidationError('Valid PSAM content identity required')
        if source.get('input_identities',{}).get('psam')!=p:
            raise core.ValidationError('QC PSAM provenance disagrees')
        if self.psam is not None and self.psam!=p:
            raise core.ValidationError('PSAM identities differ across blocks')
        if meta.get('psam_identity',p)!=p:
            raise core.ValidationError('PSAM identity differs across chromosome selections')
        self.psam=p
        if source.get('summary_dosage_field')!='qc_effective_alt_dosage':
            raise core.ValidationError('Effective QC dosage summary policy required')
        if sum(source['failure_combinations'].values())!=source['input_rows']:
            raise core.ValidationError('QC failure combinations do not reconcile')
        def reconcile(counts,children=None):
            if counts['input_rows']!=counts['pass_rows']+counts['fail_rows']:
                raise core.ValidationError('QC count reconciliation failed')
            if children is not None and any(sum(c[k] for c in children)!=counts[k] for k in ('input_rows','pass_rows','fail_rows')):
                raise core.ValidationError('QC nested counts disagree')
        types=list(source['by_candidate_type'].values());reconcile(source,types)
        for t in types:
            classes=list(t['by_allele_class'].values());reconcile(t,classes)
            for c in classes:reconcile(c,list(c['by_tier'].values()))
        if source.get('chromosome','').removeprefix('chr')!=meta['chromosome'].removeprefix('chr'):
            raise core.ValidationError('QC chromosome mismatch')

    def validate_row(self,r):
        # Validate the saved ploidy/dosage contract; no quality or frequency recount.
        try:region=frequency_region(r['CHROM'],int(r['POS']))
        except ValueError as e:raise core.ValidationError('Invalid saved QC region') from e
        if r['qc_sex'] not in ('male','female','unknown'):
            raise core.ValidationError('Invalid saved QC sex')
        ploidy,excluded=expected_ploidy(region,r['qc_sex'])
        if excluded or r['qc_frequency_region']!=region or r['qc_expected_ploidy']!=str(ploidy):
            raise core.ValidationError('Inconsistent saved QC ploidy metadata')
        gt=r['GT']
        allowed={'1','1/1','1|1'} if ploidy==1 else {'0/1','1/0','0|1','1|0','1/1','1|1'}
        if gt not in allowed:raise core.ValidationError('GT incompatible with saved passing QC ploidy')
        raw=re.split(r'[/|]',gt).count('1');effective=1 if ploidy==1 else raw
        if r['alt_dosage']!=str(raw) or r['qc_effective_alt_dosage']!=str(effective):
            raise core.ValidationError('Saved raw/effective dosage disagrees with GT/ploidy')

    def receipt(self):
        return dict(psam_identity=self.psam,summary_dosage_field='qc_effective_alt_dosage',
            dosage_policy='Annotation-stratified allele totals use effective dosage; raw GT/dosage preserved; distinct allele/sample unions never sum cross-type associations',
            adapter_identity=core.identity(__file__),helper_identities={n:core.identity(Path(__file__).with_name(n)) for n in
            ['qc_post_rarity.py','carrier_frequencies.py','carrier_psam.py','extract_exact_carriers.py','filter_final_rarity.py']})


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['carriers','samples','receipts']:p.add_argument('--'+name,nargs='+',required=True)
    for name in ['metadata','outdir']:p.add_argument('--'+name,required=True)
    try:core.run(p.parse_args(),adapter=sys.modules[__name__])
    except Exception:sys.exit('Post-rarity QC gathering failed; inspect local receipt. No protected records printed.')
