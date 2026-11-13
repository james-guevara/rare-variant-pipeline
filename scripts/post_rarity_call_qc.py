#!/usr/bin/env python3
"""Call-level QC of one gathered chromosome (after GATHER_POST_RARITY_QC); same products, fewer calls, extra annotation.

Rules (sequence alleles; spanning-deletion '*' rows pass through unchanged):
  strand     tiered HC LoF (lof_t1, lof_t2) removed when its source FORMAT/SB is informative and shows a single called ALT with
             ALT reads on exactly one strand and REF reads on both (SB = REF-fwd, REF-rev, ALT-fwd, ALT-rev). Missing,
             uninformative or ambiguous multiallelic SB is retained. Requires --zarr-dir; without it the rule is not applied.
  site AB    every call of an allele whose site mean allele balance (mean qc_AB over the allele's carriers in this chromosome,
             one value per sample) is below --snv-site-ab-min (SNVs) or --indel-site-ab-min (indels) is removed.
  exemption  calls in the DNM rescue list (sample + raw or normalized allele) are never removed by these rules.
Annotation columns added to every row: callqc_site_mean_ab, callqc_site_carriers, callqc_n_alt_alleles (ALT alleles at the
source record; needs --zarr-dir), callqc_strand_status, callqc_dnm_exempt. Removed rows go to callqc_audit.tsv.gz with reasons.
Thresholds were chosen on SPARK iWGS/iWES and SSC (2026-11-12): clogit, trio labels and parent-to-child transmission.
"""
import argparse
from collections import Counter
import csv
import gzip
import json
from pathlib import Path
import sys
import time
from gather_post_qc import aggregate, write_parquet, PRODUCTS as GATHER_PRODUCTS
from qc_filtered_carriers import ValidationError, identity, write_table, row_key, group
from qc_post_rarity import load_rescue, RESCUE_FIELDS

PRODUCTS = GATHER_PRODUCTS + ['callqc_audit.tsv.gz']
ANNOTATION = ['callqc_site_mean_ab', 'callqc_site_carriers', 'callqc_n_alt_alleles', 'callqc_strand_status', 'callqc_dnm_exempt']
LOF_TIERS = ('lof_t1', 'lof_t2')


def policy(a):
    return dict(version='call-qc-v1', applies_to='sequence alleles; spanning-deletion rows unchanged',
                strand=('tiered HC LoF removed if informative SB, single called ALT, ALT reads on exactly one strand and REF reads on both; '
                        'missing/uninformative/ambiguous SB retained') if a.zarr_dir else 'not applied (no Zarr source)',
                site_ab=f'remove all calls of an allele with site mean qc_AB < {a.snv_site_ab_min} (SNV) / < {a.indel_site_ab_min} (indel); '
                        'site mean over one value per carrier sample in this chromosome',
                snv_site_ab_min=a.snv_site_ab_min, indel_site_ab_min=a.indel_site_ab_min,
                dnm_exemption='calls in the DNM rescue list are never removed by call QC' if a.dnm_rescue else 'none (no rescue list)')


def strand_status(sb, gt, ai):
    if sb is None or len(sb) != 4 or min(sb) < 0: return 'missing'
    if set(int(x) for x in gt) - {0} != {ai}: return 'ambiguous_multiallelic'
    if sum(sb[:2]) == 0 or sum(sb[2:]) == 0: return 'uninformative'
    return 'one_strand_alt_ref_both' if min(sb[2:]) == 0 and min(sb[:2]) > 0 else 'informative_retained'


def read_zarr(zarr_dir, rows):
    """SB (tiered LoF sequence calls) and ALT-allele counts per source record, read from <zarr_dir>/<unit>.no-pl.zarr."""
    import numpy as np
    import zarr
    sb, nalt, reads = {}, {}, 0
    by_unit = {}
    for r in rows:
        by_unit.setdefault(r['source_unit_id'], []).append(r)
    for unit, rr in by_unit.items():
        g = zarr.open_group(str(Path(zarr_dir) / (unit + '.no-pl.zarr')), mode='r')
        index = {str(s): i for i, s in enumerate(g['sample_id'][:])}
        vis = sorted({int(r['source_variant_index']) for r in rr})
        alleles = g['variant_allele'].oindex[vis, :]
        for v, al in zip(vis, alleles):
            nalt[unit, v] = int(sum(1 for x in list(al)[1:] if str(x) not in ('', 'None')))
        need = sorted({(int(r['source_variant_index']), r['sample'], int(r['source_alt_index'])) for r in rr
                       if r['tier'] in LOF_TIERS and r['allele_class'] == 'sequence'})
        if not need: continue
        if any(s not in index for _, s, _ in need): raise ValidationError('Carrier sample absent from Zarr store')
        vi = np.array([v for v, _, _ in need]); si = np.array([index[s] for _, s, _ in need])
        gt = np.asarray(g['call_genotype'].vindex[vi[:, None], si[:, None], np.arange(2)[None, :]])
        has_sb = 'call_SB' in g
        s = np.asarray(g['call_SB'].vindex[vi[:, None], si[:, None], np.arange(4)[None, :]]) if has_sb else None
        ids = np.asarray(g['sample_id'].vindex[si])
        for j, (v, smp, a) in enumerate(need):
            if str(ids[j]) != smp: raise ValidationError('Zarr sample mismatch')
            sb[unit, v, smp, a] = (strand_status(None if s is None else [int(x) for x in s[j]], gt[j].tolist(), a), gt[j].tolist())
        reads += len(need)
    return sb, nalt, reads


def run(a):
    start = time.perf_counter(); src = Path(a.gather_dir).resolve(); out = Path(a.outdir).resolve()
    if out == src or src in out.parents or out in src.parents: raise ValidationError('Separate call-QC output directory required')
    out.mkdir(parents=True, exist_ok=True)
    receipt = dict(schema_version=1, stage='post_rarity_call_qc', status='failed', call_qc_policy=policy(a))
    try:
        g = json.loads((src / 'receipt.json').read_text())
        if g.get('stage') != 'post_rarity_qc_gather' or g.get('status') != 'passed' or not g.get('reconciliation', {}).get('passed'):
            raise ValidationError('Passed post-rarity gather receipt required')
        ids = {n: identity(src / n) for n in ['carriers.qc.tsv.gz', 'samples.tsv', 'receipt.json']}
        if {k: ids['carriers.qc.tsv.gz'][k] for k in ['bytes', 'sha256']} != {k: g['outputs']['carriers.qc.tsv.gz'][k] for k in ['bytes', 'sha256']}:
            raise ValidationError('Gathered carriers differ from their receipt')
        with gzip.open(src / 'carriers.qc.tsv.gz', 'rt', newline='') as f:
            reader = csv.DictReader(f, delimiter='\t'); fields = list(reader.fieldnames or []); rows = list(reader)
        if len(rows) != g['carrier_annotation_records'] or any(x in fields for x in ANNOTATION): raise ValidationError('Gathered carrier table does not match receipt')
        with open(src / 'samples.tsv') as f: samples = [r['sample'] for r in csv.DictReader(f, delimiter='\t')]
        chrom = g['chromosome']
        rescue = load_rescue(a.dnm_rescue, chrom) if a.dnm_rescue else set()
        sb, nalt, reads = read_zarr(a.zarr_dir, rows) if a.zarr_dir else ({}, {}, 0)
        per_sample = {}
        for r in rows:
            if r['allele_class'] == 'sequence': per_sample[row_key(r)[:4], r['sample']] = float(r['qc_AB'])
        site = {}
        for (allele, _), ab in per_sample.items(): site.setdefault(allele, []).append(ab)
        site = {k: (sum(v) / len(v), len(v)) for k, v in site.items()}
        kept, removed, reasons = [], [], Counter()
        for r in rows:
            k = row_key(r); allele = k[:4]
            exempt = (r['sample'], k[0], k[1], k[2], k[3]) in rescue
            st = sb.get((r['source_unit_id'], int(r['source_variant_index']), r['sample'], int(r['source_alt_index'])))
            if st is not None and [int(x) for x in st[1]] != [int(x) for x in r['source_GT'].replace('|', '/').split('/')]:
                raise ValidationError('Zarr genotype differs from saved source genotype')
            mean, n = site.get(allele, (None, None))
            r.update(callqc_site_mean_ab='.' if mean is None else f'{mean:.6g}', callqc_site_carriers='.' if n is None else str(n),
                     callqc_n_alt_alleles=str(nalt.get((r['source_unit_id'], int(r['source_variant_index'])), '.')),
                     callqc_strand_status=st[0] if st else ('not_applicable' if r['tier'] not in LOF_TIERS or r['allele_class'] != 'sequence' else 'not_read'),
                     callqc_dnm_exempt=str(exempt))
            why = []
            if r['allele_class'] == 'sequence' and not exempt:
                if st and st[0] == 'one_strand_alt_ref_both': why.append('strand')
                snv = len(r['REF']) == 1 and len(r['ALT']) == 1
                if mean is not None and mean < (a.snv_site_ab_min if snv else a.indel_site_ab_min): why.append('site_ab')
            if why:
                reasons['|'.join(why)] += 1; removed.append(dict(r, callqc_reason='|'.join(why)))
            else:
                kept.append(r)
        if len(kept) + len(removed) != len(rows): raise ValidationError('Call-QC reconciliation failed')
        fields_out = fields + ANNOTATION
        write_table(out / 'carriers.qc.tsv.gz', fields_out, kept)
        write_parquet(out / 'carriers.qc.tsv.gz', out / 'carriers.qc.parquet', fields_out, len(kept))
        write_table(out / 'callqc_audit.tsv.gz', fields_out + ['callqc_reason'], removed)
        write_table(out / 'samples.tsv', ['sample'], (dict(sample=s) for s in samples))
        strata = {(x['candidate_type'], x['allele_class'], x['tier']) for x in g['by_stratum']}
        strata_report, classes = aggregate(out, kept, samples, strata, 'qc_effective_alt_dosage')
        if ids != {n: identity(src / n) for n in ids}: raise ValidationError('Gather inputs changed during call QC')
        for k in ['chromosome', 'gather_id', 'selected_units', 'manifest_units', 'full_manifest_selected', 'policy', 'psam_identity',
                  'duplicate_checks', 'summary_dosage_field', 'samples_in_source']:
            if k in g: receipt[k] = g[k]
        receipt.update(status='passed', source_gather=dict(directory=str(src), identities=ids),
                       input_rows=len(rows), carrier_annotation_records=len(kept), removed_rows=len(removed),
                       removal_reasons=dict(reasons), dnm_exempt_rows=sum(r['callqc_dnm_exempt'] == 'True' for r in rows),
                       dnm_rescue_identity=identity(a.dnm_rescue) if a.dnm_rescue else None, zarr_sb_reads=reads,
                       strand_status=dict(Counter(r['callqc_strand_status'] for r in rows)),
                       by_stratum=strata_report, by_allele_class=classes,
                       reconciliation=dict(passed=True, input_equals_kept_plus_removed=True),
                       script_identity=identity(__file__), outputs={n: identity(out / n) for n in PRODUCTS})
    except Exception as e:
        receipt.update(error_type=type(e).__name__, error=str(e) if isinstance(e, ValidationError) else 'Call QC failed; inspect inputs locally')
        for n in PRODUCTS: (out / n).unlink(missing_ok=True)
        raise
    finally:
        receipt['wall_seconds'] = time.perf_counter() - start
        (out / 'receipt.json').write_text(json.dumps(receipt, indent=2, sort_keys=True) + '\n')
    print(json.dumps({k: receipt.get(k) for k in ['status', 'chromosome', 'input_rows', 'carrier_annotation_records', 'removed_rows', 'removal_reasons']}))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--gather-dir', required=True); p.add_argument('--outdir', required=True)
    p.add_argument('--zarr-dir', help='directory holding <unit>.no-pl.zarr stores (strand rule and ALT-allele counts)')
    p.add_argument('--dnm-rescue', help='TSV(.gz) ' + ','.join(RESCUE_FIELDS) + ': calls exempt from call QC')
    p.add_argument('--snv-site-ab-min', type=float, default=0.35); p.add_argument('--indel-site-ab-min', type=float, default=0.375)
    try: run(p.parse_args())
    except Exception: sys.exit('Call QC failed; inspect receipt. No protected records printed.')
