#!/usr/bin/env python3
"""Indexed, exact biallelic HC-LOFTEE carrier extraction. No genotype QC filtering."""
import argparse
from collections import Counter
import csv
import hashlib
import json
from pathlib import Path
import re
import sys
import time

import pysam

FIELDS = ['CHROM', 'POS', 'REF', 'ALT', 'Gene', 'Feature', 'SYMBOL', 'sample',
          'GT', 'GQ', 'DP', 'AD', 'FT', 'alt_dosage', 'site_FILTER', 'allele_class']
KEY_FIELDS = ['CHROM', 'POS', 'REF', 'ALT']
ALLELE_CLASSES = ('sequence', 'spanning_deletion')


class ValidationError(ValueError):
    pass


def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def chrom(value):
    return value[3:] if value.startswith('chr') else value


def candidates(path, chromosome):
    counts = Counter()
    result = {}
    with open(path, newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        required = set(KEY_FIELDS + ['Gene', 'Feature', 'LoF'])
        if not required <= set(reader.fieldnames or ()):
            raise ValidationError('LOFTEE TSV lacks required allele/gene/LoF columns')
        for row in reader:
            counts['annotation_rows'] += 1
            counts['HC' if row['LoF'] == 'HC' else 'non_HC'] += 1
            if row['LoF'] != 'HC':
                continue
            if chrom(row['CHROM']) != chrom(chromosome):
                raise ValidationError('HC candidate chromosome disagrees with manifest')
            if len(row['ALT'].split(',')) != 1:
                raise ValidationError('HC candidate must have exactly one ALT')
            if not re.fullmatch('[ACGTN]+', row['REF']) or not (
                    row['ALT'] == '*' or re.fullmatch('[ACGTN]+', row['ALT'])):
                raise ValidationError('HC candidate requires sequence REF and sequence ALT or spanning-deletion *')
            row['allele_class'] = 'spanning_deletion' if row['ALT'] == '*' else 'sequence'
            counts[row['allele_class']] += 1
            if row['REF'] == row['ALT'] or int(row['POS']) < 1:
                raise ValidationError('Invalid HC candidate allele or position')
            if row['Gene'] in ('', '.', '-'):
                raise ValidationError('HC candidate lacks gene identity')
            key = (chrom(row['CHROM']), int(row['POS']), row['REF'], row['ALT'])
            if key in result:
                raise ValidationError('Duplicate exact HC candidate; expected one picked transcript per allele')
            result[key] = row
    return result, counts


def text(value):
    if value is None:
        return '.'
    if isinstance(value, (tuple, list)):
        return ','.join(text(x) for x in value)
    return str(value)


def file_identity(path, hashed=False):
    stat = Path(path).stat()
    data = dict(bytes=stat.st_size, mtime_ns=stat.st_mtime_ns)
    if hashed:
        data['sha256'] = sha(path)
    return data


def source_header(vcf, chromosome):
    aliases = [c for c in vcf.header.contigs if chrom(c) == chrom(chromosome)]
    if len(aliases) != 1:
        raise ValidationError('Source VCF must have one unambiguous declared chromosome alias')
    contig = aliases[0]
    samples = list(vcf.header.samples)
    if not samples or len(samples) != len(set(samples)):
        raise ValidationError('Expected genotype-bearing VCF with unique samples')
    next(vcf.fetch(contig, 0, 1), None)  # Require index even with no candidates.
    return contig, samples


def exact_records(vcf, contig, selected):
    """Shared indexed exact-allele lookup; region overlap alone never matches."""
    found = set()
    for position in sorted({key[1] for key in selected}):
        for record in vcf.fetch(contig, position-1, position):
            if record.pos != position:
                continue
            if len(record.alts or ()) != 1:
                raise ValidationError('Multiallelic source record at candidate position; normalize upstream')
            key = (chrom(record.contig), record.pos, record.ref, record.alts[0])
            if key not in selected:
                continue
            if key in found:
                raise ValidationError('Duplicate exact source allele record')
            if 'GT' not in record.format:
                raise ValidationError('Matching source record lacks FORMAT/GT')
            found.add(key)
            yield key, record


def carrier_calls(record):
    """Shared no-QC genotype semantics, including haploid and partial ALT calls."""
    for sample, call in record.samples.items():
        gt = call.get('GT')
        if gt is None:
            continue
        if any(allele not in (None, 0, 1) for allele in gt):
            raise ValidationError('Non-biallelic genotype allele index')
        if 1 not in gt:
            continue
        row = dict(CHROM=record.contig, POS=record.pos, REF=record.ref, ALT=record.alts[0],
                   sample=sample, GT=('|' if call.phased else '/').join(text(x) for x in gt),
                   alt_dosage=gt.count(1), site_FILTER=';'.join(record.filter) or '.')
        row.update({name: text(call.get(name)) for name in ['GQ', 'DP', 'AD', 'FT']})
        yield row, None in gt


def extract(a):
    started = time.perf_counter()
    meta = json.loads(Path(a.metadata).read_text())
    out = Path(a.outdir); out.mkdir(parents=True, exist_ok=True)
    receipt = dict(schema_version=2, status='failed', unit_id=meta['unit_id'], chromosome=meta['chromosome'],
                   sources=meta, pysam_version=pysam.__version__,
                   definition='LoF == HC; exact CHROM/POS/REF/ALT; one row per variant/sample carrying ALT index 1',
                   burden_definition='separate sequence HC variant counts and spanning-deletion HC record counts per sample/gene; homozygous ALT counts once',
                   aggregate_scope='top-level totals include both allele classes as a raw audit; use by_allele_class for separated counts',
                   spanning_deletion_policy='retain exact star records separately; no upstream deletion mapping, event deduplication, or combined biological burden',
                   genotype_filter='none; partial calls with a called ALT are included; quality fields preserved',
                   vcf_identity_method='path/size/mtime plus index SHA-256; no full genotype VCF hashing or scanning')
    products = ['carriers.tsv.gz', 'candidates.tsv', 'unmatched.tsv', 'samples.tsv',
                'sample_gene_burden.tsv', 'sample_burden.tsv', 'gene_burden.tsv']
    try:
        identities = {name: file_identity(getattr(a, name), name != 'vcf') for name in ['vcf', 'index', 'loftee']}
        selected, counts = candidates(a.loftee, meta['chromosome'])
        receipt.update(annotation_rows=counts['annotation_rows'], hc_annotation_rows=counts['HC'],
                       non_hc_annotation_rows=counts['non_HC'], candidate_hc_variants=len(selected),
                       by_allele_class={kind: dict(candidate_hc_variants=counts[kind]) for kind in ALLELE_CLASSES})
        if a.expected_hc is not None and len(selected) != a.expected_hc:
            raise ValidationError('HC candidate count does not match expected count')
        found = set()
        burdens = Counter()
        gene_records = Counter()
        carrier_variants = set()
        partial = dosage_total = carrier_records = 0
        class_partials, class_dosages = Counter(), Counter()
        with pysam.VariantFile(a.vcf, index_filename=a.index) as vcf:
            contig, samples = source_header(vcf, meta['chromosome'])
            with open(out/'samples.tsv', 'w') as handle:
                handle.write('sample\n' + ''.join(s+'\n' for s in samples))
            with pysam.BGZFile(str(out/'carriers.tsv.gz'), 'w') as target:
                target.write(('\t'.join(FIELDS)+'\n').encode())
                for key, record in exact_records(vcf, contig, selected):
                    found.add(key)
                    annotation = selected[key]
                    kind = annotation['allele_class']
                    for row, is_partial in carrier_calls(record):
                        row.update(Gene=annotation['Gene'], Feature=annotation['Feature'],
                                   SYMBOL=annotation.get('SYMBOL', '.'), allele_class=kind)
                        target.write(('\t'.join(text(row[name]) for name in FIELDS)+'\n').encode())
                        burdens[row['sample'], annotation['Gene'], kind] += 1
                        gene_records[annotation['Gene'], kind] += 1
                        class_partials[kind] += is_partial
                        class_dosages[kind] += row['alt_dosage']
                        carrier_variants.add(key)
                        carrier_records += 1
                        partial += is_partial
                        dosage_total += row['alt_dosage']
        for name, keys in [('candidates.tsv', selected), ('unmatched.tsv', selected.keys()-found)]:
            with open(out/name, 'w', newline='') as handle:
                fields = KEY_FIELDS + ['Gene', 'Feature', 'SYMBOL', 'allele_class']
                writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', lineterminator='\n', extrasaction='ignore')
                writer.writeheader()
                for key in sorted(keys):
                    writer.writerow({k: selected[key].get(k, '.') for k in fields})
        burden_columns = 'sequence_hc_variant_count\tspanning_deletion_hc_record_count'
        with open(out/'sample_gene_burden.tsv', 'w') as handle:
            handle.write('sample\tGene\t'+burden_columns+'\n')
            for sample, gene in sorted({(s, g) for s, g, kind in burdens}):
                handle.write(f'{sample}\t{gene}\t{burdens[sample, gene, "sequence"]}\t{burdens[sample, gene, "spanning_deletion"]}\n')
        totals = Counter()
        for (sample, gene, kind), count in burdens.items():
            totals[sample, kind] += count
        with open(out/'sample_burden.tsv', 'w') as handle:
            handle.write('sample\t'+burden_columns+'\n')
            for sample in samples:
                handle.write(f'{sample}\t{totals[sample, "sequence"]}\t{totals[sample, "spanning_deletion"]}\n')
        gene_sample_counts = Counter((gene, kind) for sample, gene, kind in burdens)
        with open(out/'gene_burden.tsv', 'w') as handle:
            handle.write('Gene\tsequence_carrier_records\tsequence_carrier_samples\tspanning_deletion_carrier_records\tspanning_deletion_carrier_samples\n')
            for gene in sorted({r['Gene'] for r in selected.values()}):
                values = [value for kind in ALLELE_CLASSES
                          for value in (gene_records[gene, kind], gene_sample_counts[gene, kind])]
                handle.write(gene+'\t'+'\t'.join(map(str, values))+'\n')
        by_class = {}
        for kind in ALLELE_CLASSES:
            keys = {key for key, row in selected.items() if row['allele_class'] == kind}
            by_class[kind] = dict(candidate_hc_variants=len(keys),
                                 matched_candidate_variants=len(keys & found),
                                 unmatched_candidate_variants=len(keys - found),
                                 matched_variants_without_carriers=len((keys & found) - carrier_variants),
                                 carrier_records=sum(n for (gene, c), n in gene_records.items() if c == kind),
                                 samples_with_hc_plof=len({s for s, g, c in burdens if c == kind}),
                                 genes_with_hc_plof=len({g for s, g, c in burdens if c == kind}),
                                 partial_call_carrier_records=class_partials[kind],
                                 observed_alt_alleles=class_dosages[kind])
        for name in identities:
            now = file_identity(getattr(a, name), False)
            if any(now[k] != identities[name][k] for k in now):
                raise ValidationError('Input changed during extraction')
        receipt.update(status='passed', input_identities=identities, annotation_rows=counts['annotation_rows'],
                       hc_annotation_rows=counts['HC'], non_hc_annotation_rows=counts['non_HC'],
                       candidate_hc_variants=len(selected), matched_candidate_variants=len(found),
                       unmatched_candidate_variants=len(selected)-len(found),
                       matched_variants_without_carriers=len(found-carrier_variants),
                       carrier_records=carrier_records, samples_in_source=len(samples),
                       samples_with_hc_plof=len({s for s, kind in totals}),
                       genes_with_hc_plof=len({g for g, kind in gene_records}), by_allele_class=by_class,
                       partial_call_carrier_records=partial, observed_alt_alleles=dosage_total,
                       contig_alias=dict(annotation=meta['chromosome'], source=contig),
                       outputs={name: file_identity(out/name, True) for name in products})
    except Exception as exc:
        receipt['error_type'] = type(exc).__name__
        # Controlled validation messages contain no source records/sample IDs.
        receipt['error'] = str(exc) if isinstance(exc, ValidationError) else 'Input or extraction failure; inspect inputs locally'
        for name in products:
            (out/name).unlink(missing_ok=True)
        raise
    finally:
        receipt['wall_seconds'] = time.perf_counter()-started
        (out/'receipt.json').write_text(json.dumps(receipt, indent=2, sort_keys=True)+'\n')
    print(json.dumps({k: receipt[k] for k in ['status', 'unit_id', 'candidate_hc_variants',
          'unmatched_candidate_variants', 'carrier_records', 'samples_with_hc_plof', 'genes_with_hc_plof', 'by_allele_class']}))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    for name in ['metadata', 'loftee', 'vcf', 'index', 'outdir']:
        p.add_argument('--'+name, required=True)
    p.add_argument('--expected-hc', type=int)
    try:
        extract(p.parse_args())
    except Exception:
        sys.exit('Carrier extraction failed; inspect the task receipt. No individual records are printed.')
