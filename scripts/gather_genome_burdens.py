#!/usr/bin/env python3
"""Recompute genome summaries from receipt-verified chromosome carrier products."""
import argparse
from collections import Counter
import csv
import gzip
import hashlib
import json
from pathlib import Path
import sqlite3
import tempfile
import time


def identity(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b''):
            h.update(block)
    return dict(bytes=Path(path).stat().st_size, sha256=h.hexdigest())


def require(ok, message):
    if not ok:
        raise ValueError(message)


def chrom(value):
    return value[3:] if value.startswith('chr') else value


def table(path):
    with open(path) as f:
        yield from csv.DictReader(f, delimiter='\t')


def export(db, path, fields, query):
    with open(path, 'w', newline='') as f:
        w = csv.writer(f, delimiter='\t', lineterminator='\n')
        w.writerow(fields.split(','))
        w.writerows(db.execute(query))


def run(manifest, expected, outdir):
    start = time.monotonic()
    sources = list(table(manifest))
    expected = [chrom(c) for c in expected]
    require(expected and len(set(expected)) == len(expected), 'Unique expected chromosomes required')
    require(len(sources) == len(expected) and {chrom(s['chromosome']) for s in sources} == set(expected),
            'Manifest must contain exactly the expected chromosomes, once each')
    out = Path(outdir).resolve()
    require(not out.exists(), 'Use a new output directory')
    parents = [Path(s['directory']).resolve() for s in sources]
    require(all(out != p and p not in out.parents and out not in p.parents for p in parents), 'Separate outputs required')
    lineage = []
    samples = None
    policy = None
    psam = None
    declared_strata = None
    expected_counts = Counter()
    total = 0
    with tempfile.TemporaryDirectory(prefix='genome-burdens-') as tmp:
        db = sqlite3.connect(str(Path(tmp) / 'burdens.sqlite'))
        db.execute('PRAGMA cache_size=-65536')
        db.execute('PRAGMA temp_store=FILE')
        db.execute('PRAGMA foreign_keys=ON')
        db.executescript('''
          CREATE TABLE samples(sample TEXT PRIMARY KEY, ordinal INTEGER UNIQUE);
          CREATE TABLE strata(typ TEXT, cl TEXT, tier TEXT, PRIMARY KEY(typ,cl,tier));
          CREATE TABLE carriers(chrom TEXT, pos INTEGER, ref TEXT, alt TEXT,
            sample TEXT REFERENCES samples(sample), typ TEXT, cl TEXT, tier TEXT, gene TEXT,
            dosage INTEGER CHECK(dosage IN (1,2)),
            PRIMARY KEY(chrom,pos,ref,alt,sample,typ),
            FOREIGN KEY(typ,cl,tier) REFERENCES strata(typ,cl,tier));
        ''')
        for s, directory in zip(sources, parents):
            rp = directory / 'receipt.json'
            rid = identity(rp)
            r = json.loads(rp.read_text())
            c = chrom(s['chromosome'])
            require(r.get('stage') == 'post_rarity_qc_gather' and r.get('status') == 'passed' and
                    r.get('reconciliation', {}).get('passed') and r.get('full_manifest_selected'),
                    'Passed complete chromosome gather required')
            require(chrom(r['chromosome']) == c, 'Receipt chromosome mismatch')
            require(r.get('summary_dosage_field') == 'qc_effective_alt_dosage', 'Effective dosage required')
            require(r.get('psam_identity') and r.get('policy'), 'Policy and PSAM provenance required')
            require(r.get('duplicate_checks') == dict(cross_block_annotation_associations=0, cross_block_exact_alleles=0),
                    'Chromosome duplicate checks must pass')
            require(r['selected_units'] and len(set(r['selected_units'])) == len(r['selected_units']) and
                    set(r['selected_units']) == set(r['manifest_units']), 'Incomplete unit selection')
            strata = {(x['candidate_type'], x['allele_class'], x['tier']) for x in r['by_stratum']}
            require(len(strata) == len(r['by_stratum']), 'Duplicate receipt strata')
            if policy is None:
                policy, psam, declared_strata = r['policy'], r['psam_identity'], strata
                db.executemany('INSERT INTO strata VALUES (?,?,?)', sorted(strata))
            require((policy, psam, declared_strata) == (r['policy'], r['psam_identity'], strata),
                    'Chromosome policies, PSAM identities or strata differ')
            ids = {}
            for name in ['samples.tsv', 'carriers.qc.tsv.gz']:
                ids[name] = identity(directory / name)
                require(ids[name] == {k: r['outputs'][name][k] for k in ['bytes', 'sha256']}, 'Product checksum mismatch: ' + name)
            ss = [x['sample'] for x in table(directory / 'samples.tsv')]
            require(ss and len(ss) == len(set(ss)) and all(x not in ('', '.') for x in ss), 'Invalid sample universe')
            require(len(ss) == r['samples_in_source'], 'Sample count mismatch')
            if samples is None:
                samples = ss
                db.executemany('INSERT INTO samples VALUES (?,?)', ((v, i) for i, v in enumerate(ss)))
            require(samples == ss, 'Sample identity/order mismatch')
            observed, doses = Counter(), Counter()
            with gzip.open(directory / 'carriers.qc.tsv.gz', 'rt') as f:
                reader = csv.DictReader(f, delimiter='\t')
                required = ['CHROM','POS','REF','ALT','sample','candidate_type','allele_class','tier','Gene','qc_effective_alt_dosage','source_unit_id']
                require(set(required) <= set(reader.fieldnames or []) and len(reader.fieldnames)==len(set(reader.fieldnames)), 'Invalid carrier schema')
                batch = []
                for row in reader:
                    require(None not in row and all(v is not None for v in row.values()), 'Malformed carrier row')
                    require(chrom(row['CHROM']) == c and int(row['POS']) > 0, 'Carrier chromosome/position mismatch')
                    require(row['source_unit_id'] in r['selected_units'], 'Unknown source unit')
                    require(row['allele_class'] == ('spanning_deletion' if row['ALT'] == '*' else 'sequence'), 'Allele class mismatch')
                    key = (row['candidate_type'], row['allele_class'], row['tier'])
                    d = int(row['qc_effective_alt_dosage'])
                    observed[key] += 1
                    doses[key] += d
                    batch.append((c, int(row['POS']), row['REF'], row['ALT'], row['sample'], *key, row['Gene'], d))
                    if len(batch) == 10000:
                        db.executemany('INSERT INTO carriers VALUES (?,?,?,?,?,?,?,?,?,?)', batch)
                        batch.clear()
                db.executemany('INSERT INTO carriers VALUES (?,?,?,?,?,?,?,?,?,?)', batch)
            declared = Counter({(x['candidate_type'],x['allele_class'],x['tier']):x['carrier_annotation_records'] for x in r['by_stratum']})
            declared_doses = Counter({(x['candidate_type'],x['allele_class'],x['tier']):x['observed_alt_alleles'] for x in r['by_stratum']})
            require(+observed == +declared and +doses == +declared_doses, 'Chromosome stratum reconciliation failed')
            require(sum(observed.values()) == r['carrier_annotation_records'] == r['source_counts']['pass_rows'], 'Chromosome total mismatch')
            total += sum(observed.values())
            expected_counts.update(observed)
            for name, value in ids.items():
                require(identity(directory / name) == value, 'Input changed during read')
            require(identity(rp) == rid, 'Receipt changed during read')
            lineage.append(dict(chromosome='chr'+c, directory=str(directory), receipt_identity=rid, outputs=ids,
                                units=len(r['selected_units']), carrier_annotation_records=sum(observed.values())))
            db.commit()
            print('Imported chr%s: %d carrier associations' % (c, sum(observed.values())), flush=True)
        require(db.execute('SELECT count(*) FROM carriers').fetchone()[0] == total, 'Total reconciliation failed')
        require(not db.execute('SELECT 1 FROM carriers GROUP BY chrom,pos,ref,alt,sample HAVING min(dosage)!=max(dosage) LIMIT 1').fetchone(),
                'Conflicting cross-type dosage')
        db.executescript('''
          CREATE TABLE cells AS SELECT sample,typ,cl,tier,gene,count(*) AS n,sum(dosage) AS dose
            FROM carriers GROUP BY sample,typ,cl,tier,gene;
          CREATE TABLE totals AS SELECT sample,typ,cl,tier,sum(n) AS n,sum(dose) AS dose
            FROM cells GROUP BY sample,typ,cl,tier;
          CREATE UNIQUE INDEX total_key ON totals(sample,typ,cl,tier);
          CREATE TABLE distinct_calls AS SELECT chrom,pos,ref,alt,sample,cl,max(dosage) AS dose
            FROM carriers GROUP BY chrom,pos,ref,alt,sample,cl;
          CREATE INDEX distinct_sample ON distinct_calls(sample,cl);
        ''')
        for t in ['cells', 'totals']:
            require(db.execute('SELECT coalesce(sum(n),0) FROM '+t).fetchone()[0] == total, 'Summary total mismatch')
        out.mkdir(parents=True)
        export(db,out/'samples.tsv','sample','SELECT sample FROM samples ORDER BY ordinal')
        export(db,out/'sample_burden.tsv','sample,candidate_type,allele_class,tier,carrier_variant_records,observed_alt_alleles',
               'SELECT s.sample,t.typ,t.cl,t.tier,coalesce(b.n,0),coalesce(b.dose,0) FROM samples s CROSS JOIN strata t LEFT JOIN totals b ON s.sample=b.sample AND t.typ=b.typ AND t.cl=b.cl AND t.tier=b.tier ORDER BY s.ordinal,t.typ,t.cl,t.tier')
        export(db,out/'sample_gene_burden.tsv','sample,candidate_type,allele_class,tier,Gene,carrier_variant_records,observed_alt_alleles',
               'SELECT sample,typ,cl,tier,gene,n,dose FROM cells ORDER BY sample,typ,cl,tier,gene')
        export(db,out/'gene_burden.tsv','candidate_type,allele_class,tier,Gene,carrier_variant_records,carrier_samples,observed_alt_alleles',
               'SELECT typ,cl,tier,gene,sum(n),count(DISTINCT sample),sum(dose) FROM cells GROUP BY typ,cl,tier,gene ORDER BY typ,cl,tier,gene')
        export(db,out/'sample_distinct_alleles.tsv','sample,sequence_distinct_variant_count,spanning_deletion_distinct_record_count',
               "SELECT s.sample,coalesce(sum(d.cl='sequence'),0),coalesce(sum(d.cl='spanning_deletion'),0) FROM samples s LEFT JOIN distinct_calls d ON s.sample=d.sample GROUP BY s.sample ORDER BY s.ordinal")
        export(db,out/'sample_gene_distinct_alleles.tsv','sample,Gene,allele_class,distinct_allele_sample_records',
               'SELECT sample,gene,cl,count(*) FROM (SELECT DISTINCT sample,gene,cl,chrom,pos,ref,alt FROM carriers) GROUP BY sample,gene,cl ORDER BY sample,gene,cl')
        stratum_report = []
        for typ, cl, tier in sorted(declared_strata):
            args = (typ,cl,tier)
            n, ns, dose = db.execute('SELECT coalesce(sum(n),0),count(DISTINCT sample),coalesce(sum(dose),0) FROM totals WHERE typ=? AND cl=? AND tier=?',args).fetchone()
            na = db.execute('SELECT count(*) FROM (SELECT DISTINCT chrom,pos,ref,alt FROM carriers WHERE typ=? AND cl=? AND tier=?)',args).fetchone()[0]
            require(n == expected_counts[args], 'Genome stratum mismatch')
            stratum_report.append(dict(candidate_type=typ,allele_class=cl,tier=tier,carrier_annotation_records=n,samples_with_carriers=ns,distinct_alleles=na,observed_alt_alleles=dose))
        classes = {}
        for cl in ['sequence','spanning_deletion']:
            n,dose,ng=db.execute('SELECT count(*),coalesce(sum(dosage),0),count(DISTINCT gene) FROM carriers WHERE cl=?',(cl,)).fetchone()
            nd,ns,dd=db.execute('SELECT count(*),count(DISTINCT sample),coalesce(sum(dose),0) FROM distinct_calls WHERE cl=?',(cl,)).fetchone()
            na=db.execute('SELECT count(*) FROM (SELECT DISTINCT chrom,pos,ref,alt FROM distinct_calls WHERE cl=?)',(cl,)).fetchone()[0]
            classes[cl]=dict(carrier_annotation_records=n,annotation_observed_alt_alleles=dose,genes_with_carriers=ng,
                             distinct_allele_sample_records=nd,samples_with_carriers=ns,distinct_alleles=na,distinct_allele_sample_observed_alt_alleles=dd)
        with (out/'stratum_summary.tsv').open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=list(stratum_report[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(stratum_report)
        # Retain one manifest referencing chromosome carriers, rather than duplicating them.
        with (out/'chromosome_sources.tsv').open('w',newline='') as f:
            w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(['chromosome','directory'])
            w.writerows((x['chromosome'],x['directory']) for x in lineage)
        receipt=dict(schema_version=1,stage='post_rarity_genome_burdens',status='passed',chromosomes=['chr'+c for c in expected],
                     policy=policy,psam_identity=psam,samples_in_source=len(samples),carrier_annotation_records=total,
                     by_stratum=stratum_report,by_allele_class=classes,lineage=lineage,summary_dosage_field='qc_effective_alt_dosage',
                     reconciliation=dict(passed=True,source_counts_match=True,unique_samples_recomputed=True,zero_carrier_samples_retained=True),
                     sample_gene_policy='Sparse tables: absent cells are zero; samples.tsv and dense sample tables retain every source sample',
                     burden_policy='Annotation associations distinct from allele/sample counts; sequence and spanning-deletion counts separate',
                     script_identity=identity(__file__),manifest_identity=identity(manifest),wall_seconds=time.monotonic()-start,
                     outputs={p.name:identity(p) for p in sorted(out.glob('*.tsv'))})
        (out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
        db.close()
        print('PASS: %d samples, %d chromosomes, %d carrier associations' % (len(samples),len(sources),total),flush=True)
        return receipt


if __name__ == '__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--manifest',required=True)
    p.add_argument('--chromosomes',required=True,help='Explicit comma-separated expected chromosome set')
    p.add_argument('--outdir',required=True)
    a=p.parse_args()
    run(a.manifest,a.chromosomes.split(','),a.outdir)
