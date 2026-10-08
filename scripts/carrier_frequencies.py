"""Representative genotype allele counts, independent of carrier genotype QC.

Only autosomal counting is exposed until the source X/Y PAR representation is
agreed. The pure counting primitive implements the agreed haploid rules so they
can be tested without guessing which source records have that expected ploidy.
"""
from collections import Counter
from extract_exact_carriers import ValidationError, chrom

SETS = ('cohort', 'unrelated')
REASONS = (
    'counted_diploid_complete', 'counted_diploid_partial', 'counted_haploid',
    'counted_diploid_encoded_haploid', 'excluded_missing',
    'excluded_invalid_allele', 'excluded_unexpected_ploidy',
    'excluded_haploid_partial', 'excluded_haploid_heterozygous',
    'excluded_unknown_ploidy',
)
POLICY = dict(
    version='representative-autosomal-no-genotype-qc-v1',
    cohort='frequency_representative == 1',
    unrelated='frequency_representative == 1 AND unrelated == 1',
    participant_validation='at most one selected representative per participant in source samples; no automatic selection',
    genotype_qc='none: no GQ/DP/AD/AB, FT, or site FILTER exclusions',
    diploid_partial='count each called allele; missing alleles do not contribute to AN',
    autosomal_ploidy='expect two GT alleles; unexpected ploidy excluded and audited',
    unknown_sex='included on autosomes; X/Y counting blocked pending PAR policy',
    sex_chromosomes='blocked pending source PAR representation; no sex or PAR defaults',
    haploid_rule='when expected ploidy is one: count literal 0 or 1; collapse 0/0 and 1/1; exclude partial and heterozygous calls; not yet applied to X/Y',
    zero_an='AF missing; cannot establish final rarity',
    unmatched='AC/AN/AF missing, distinct from a matched site with AN=0',
    allele_class='sequence and spanning_deletion frequencies/audits kept separate; no deletion-event inference',
    final_rarity='corrected unrelated AF < 0.001, separate downstream PR; not applied here',
)


def validate_chromosome(chromosome):
    if chrom(chromosome) not in {str(i) for i in range(1, 23)}:
        raise ValidationError('Corrected frequency counting supports autosomes only until X/Y PAR policy is agreed; raw extraction remains available')


def count_call(gt, expected_ploidy):
    """Return AC, AN, and one mutually exclusive audit reason. No quality inputs."""
    if expected_ploidy not in (1, 2):
        return 0, 0, 'excluded_unknown_ploidy'
    if gt is None or len(gt) == 0:
        return 0, 0, 'excluded_missing'
    if any(a not in (None, 0, 1) for a in gt):
        return 0, 0, 'excluded_invalid_allele'
    if expected_ploidy == 2:
        if len(gt) != 2:
            return 0, 0, 'excluded_unexpected_ploidy'
        called = [a for a in gt if a is not None]
        if not called:
            return 0, 0, 'excluded_missing'
        return sum(called), len(called), 'counted_diploid_partial' if len(called) == 1 else 'counted_diploid_complete'
    if len(gt) not in (1, 2):
        return 0, 0, 'excluded_unexpected_ploidy'
    if all(a is None for a in gt):
        return 0, 0, 'excluded_missing'
    if None in gt:
        return 0, 0, 'excluded_haploid_partial'
    if len(set(gt)) != 1:
        return 0, 0, 'excluded_haploid_heterozygous'
    return gt[0], 1, 'counted_haploid' if len(gt) == 1 else 'counted_diploid_encoded_haploid'


class FrequencyCounter:
    def __init__(self, metadata, chromosome):
        validate_chromosome(chromosome)
        self.samples = []
        selected_participants = set()
        for row in metadata:
            if row['frequency_representative'] != '1':
                continue
            participant = row['participant_id']
            if participant in selected_participants:
                raise ValidationError('Multiple source frequency representatives for one participant')
            selected_participants.add(participant)
            self.samples.append((row['#IID'], row['unrelated'] == '1'))
        self.sizes = dict(cohort=len(self.samples), unrelated=sum(u for _, u in self.samples))
        self.selection = dict(
            **self.sizes,
            source_participants=len({r['participant_id'] for r in metadata}),
            participants_without_representative=len({r['participant_id'] for r in metadata} - selected_participants),
            unrelated_flags_without_representative=sum(r['unrelated'] == '1' and r['frequency_representative'] != '1' for r in metadata),
        )
        self.audit = Counter()
        self.matched = Counter()
        self.zero_an = Counter()

    def count(self, record, allele_class):
        """Called once per exact source allele, including alleles with no carriers."""
        values = {s: dict(ac=0, an=0, counted_genotypes=0, reference_genotypes=0) for s in SETS}
        audit = {s: Counter() for s in SETS}
        for sample, unrelated in self.samples:
            ac, an, reason = count_call(record.samples[sample].get('GT'), 2)
            for group in SETS if unrelated else ('cohort',):
                audit[group][reason] += 1
                value = values[group]
                value['ac'] += ac
                value['an'] += an
                value['counted_genotypes'] += an > 0
                value['reference_genotypes'] += an > 0 and ac == 0
        result = {}
        for group in SETS:
            value = values[group]
            value['af'] = value['ac'] / value['an'] if value['an'] else None
            value['excluded_genotypes'] = self.sizes[group] - value['counted_genotypes']
            result.update({group + '_' + k: v for k, v in value.items()})
            if sum(audit[group].values()) != self.sizes[group]:
                raise ValidationError('Frequency genotype audit failed reconciliation')
            for reason, n in audit[group].items():
                self.audit[allele_class, group, reason] += n
            self.zero_an[allele_class, group] += value['an'] == 0
        self.matched[allele_class] += 1
        return result, audit

    def receipt(self):
        by_class = {}
        for kind in ('sequence', 'spanning_deletion'):
            by_class[kind] = dict(matched_distinct_alleles=self.matched[kind], sample_sets={})
            for group in SETS:
                reasons = {r: self.audit[kind, group, r] for r in REASONS}
                expected = self.matched[kind] * self.sizes[group]
                if sum(reasons.values()) != expected:
                    raise ValidationError('Frequency aggregate audit failed reconciliation')
                by_class[kind]['sample_sets'][group] = dict(
                    eligible_samples=self.sizes[group], evaluated_genotypes=sum(reasons.values()),
                    expected_evaluated_genotypes=expected, by_reason=reasons,
                    zero_an_alleles=self.zero_an[kind, group],
                )
        return dict(status='passed', policy=POLICY, sample_selection=self.selection, by_allele_class=by_class)
