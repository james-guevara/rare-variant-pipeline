"""Representative genotype allele counts, independent of carrier genotype QC.

Sex chromosome counting requires the explicit GRCh38 X-only PAR policy.
"""
from collections import Counter
from extract_exact_carriers import ValidationError, chrom

SETS = ('cohort', 'unrelated')
SEX_POLICY = 'grch38_x_only_par'
# Source-policy coordinates: 1-based inclusive, classified by variant POS.
PAR = {'X': ((10001, 2781479), (155701383, 156030895)),
       'Y': ((10001, 2781479), (56887903, 57217415))}
SEX_REASONS = ('excluded_unknown_sex', 'excluded_female_Y')
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


def validate_chromosome(chromosome, sex_policy=None):
    if sex_policy not in (None, SEX_POLICY):
        raise ValidationError('Unsupported sex chromosome frequency policy')
    if chrom(chromosome) in {str(i) for i in range(1, 23)}:
        return
    if chrom(chromosome) in ('X', 'Y') and sex_policy == SEX_POLICY:
        return
    raise ValidationError('Corrected X/Y frequencies require explicit --sex-chromosome-policy grch38_x_only_par; supported chromosomes are 1-22, X, Y')


def policy_for(chromosome, sex_policy=None):
    if chrom(chromosome) not in ('X', 'Y') or sex_policy != SEX_POLICY:
        return POLICY
    return dict(POLICY,
        version='representative-grch38-x-only-par-no-genotype-qc-v1',
        sex_chromosomes=SEX_POLICY, genome_build='GRCh38',
        sex_encoding='PSAM SEX: 1=male; 2=female; all other nonempty values=unknown',
        unknown_sex='exclude from X/Y frequency counts and audit, including X PAR; include on autosomes',
        par_coordinates_1based_inclusive=PAR,
        par_assignment='variant POS; X PAR diploid in both known sexes; Y PAR candidates forbidden even if unmatched',
        nonpar_ploidy='X male haploid/female diploid; Y male haploid/female excluded',
        haploid_rule='literal GT 0 or 1 counts one allele; collapse 0/0 or 1/1 to one; exclude heterozygous and partial calls',
        sample_eligibility='receipt eligible_samples means selected representatives before per-site sex/ploidy exclusions',
    )


def frequency_region(chromosome, position):
    chromosome = chrom(chromosome)
    if chromosome not in ('X', 'Y'):
        return 'autosome'
    in_par = any(start <= position <= end for start, end in PAR[chromosome])
    if chromosome == 'Y' and in_par:
        raise ValidationError('Unexpected Y-PAR candidate under GRCh38 X-only PAR policy; no corrected output produced')
    return chromosome + ('_PAR' if in_par else '_nonPAR')


def validate_candidate_regions(selected, chromosome, sex_policy=None):
    validate_chromosome(chromosome, sex_policy)
    if chrom(chromosome) == 'Y':
        for key in selected:
            frequency_region(chromosome, key[1])


def normalize_sex(value):
    return {'1': 'male', '2': 'female'}.get(str(value).strip(), 'unknown')


def expected_ploidy(region, sex):
    if region == 'autosome':
        return 2, None
    if sex == 'unknown':
        return None, 'excluded_unknown_sex'
    if region == 'X_PAR':
        return 2, None
    if region == 'X_nonPAR':
        return (1 if sex == 'male' else 2), None
    if region == 'Y_nonPAR':
        return (1, None) if sex == 'male' else (None, 'excluded_female_Y')
    raise ValidationError('Unsupported frequency region')


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
    def __init__(self, metadata, chromosome, sex_policy=None):
        validate_chromosome(chromosome, sex_policy)
        self.chromosome = chrom(chromosome)
        self.sex_chromosome = self.chromosome in ('X', 'Y')
        self.policy = policy_for(chromosome, sex_policy)
        self.reasons = REASONS + SEX_REASONS if self.sex_chromosome else REASONS
        self.regions = Counter()
        self.samples = []
        selected_participants = set()
        for row in metadata:
            if row['frequency_representative'] != '1':
                continue
            participant = row['participant_id']
            if participant in selected_participants:
                raise ValidationError('Multiple source frequency representatives for one participant')
            selected_participants.add(participant)
            self.samples.append((row['#IID'], row['unrelated'] == '1', normalize_sex(row['SEX'])))
        self.sizes = dict(cohort=len(self.samples), unrelated=sum(u for _, u, _ in self.samples))
        self.selection = dict(
            **self.sizes,
            source_participants=len({r['participant_id'] for r in metadata}),
            participants_without_representative=len({r['participant_id'] for r in metadata} - selected_participants),
            unrelated_flags_without_representative=sum(r['unrelated'] == '1' and r['frequency_representative'] != '1' for r in metadata),
        )
        if self.sex_chromosome:
            self.selection['by_sex'] = {
                group: {sex: sum(s == sex and (group == 'cohort' or u) for _, u, s in self.samples)
                        for sex in ('male', 'female', 'unknown')} for group in SETS}
        self.audit = Counter()
        self.matched = Counter()
        self.zero_an = Counter()

    def count(self, record, allele_class):
        """Called once per exact source allele, including alleles with no carriers."""
        region = frequency_region(self.chromosome, record.pos) if self.sex_chromosome else 'autosome'
        values = {s: dict(ac=0, an=0, counted_genotypes=0, reference_genotypes=0) for s in SETS}
        audit = {s: Counter() for s in SETS}
        for sample, unrelated, sex in self.samples:
            ploidy, excluded = expected_ploidy(region, sex)
            ac, an, reason = (0, 0, excluded) if excluded else count_call(record.samples[sample].get('GT'), ploidy)
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
        if self.sex_chromosome:
            result['frequency_region'] = region
            self.regions[allele_class, region] += 1
        return result, audit

    def receipt(self):
        by_class = {}
        for kind in ('sequence', 'spanning_deletion'):
            by_class[kind] = dict(matched_distinct_alleles=self.matched[kind], sample_sets={})
            for group in SETS:
                reasons = {r: self.audit[kind, group, r] for r in self.reasons}
                expected = self.matched[kind] * self.sizes[group]
                if sum(reasons.values()) != expected:
                    raise ValidationError('Frequency aggregate audit failed reconciliation')
                by_class[kind]['sample_sets'][group] = dict(
                    eligible_samples=self.sizes[group], evaluated_genotypes=sum(reasons.values()),
                    expected_evaluated_genotypes=expected, by_reason=reasons,
                    zero_an_alleles=self.zero_an[kind, group],
                )
        if self.sex_chromosome:
            for kind in by_class:
                by_class[kind]['matched_by_region'] = {
                    region: self.regions[kind, region] for region in
                    (('X_PAR', 'X_nonPAR') if self.chromosome == 'X' else ('Y_nonPAR',))}
        return dict(status='passed', policy=self.policy, sample_selection=self.selection, by_allele_class=by_class)
