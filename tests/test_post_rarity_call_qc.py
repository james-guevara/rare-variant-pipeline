from argparse import Namespace
import gzip
import json
from pathlib import Path
import sys
import pytest
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from qc_filtered_carriers import identity, write_table, ValidationError
from test_filtered_carriers import rows
import post_rarity_call_qc as cq
import gather_genome_burdens as gg

FIELDS = ['CHROM', 'POS', 'REF', 'ALT', 'sample', 'candidate_type', 'allele_class', 'tier', 'Gene', 'Feature', 'SYMBOL', 'GT',
          'alt_dosage', 'qc_AB', 'qc_effective_alt_dosage', 'source_GT', 'source_variant_index', 'source_alt_index', 'source_unit_id']


def row(pos, sample, ab, ref='A', alt='G', tier='lof_t1', typ='lof_hc', vi=None, gene='G1'):
    return dict(CHROM='chr22', POS=str(pos), REF=ref, ALT=alt, sample=sample, candidate_type=typ, allele_class='sequence' if alt != '*' else 'spanning_deletion',
                tier=tier, Gene=gene, Feature='TX1', SYMBOL='S', GT='0/1', alt_dosage='1', qc_AB=str(ab), qc_effective_alt_dosage='1',
                source_GT='0/1', source_variant_index=str(vi if vi is not None else pos), source_alt_index='1', source_unit_id='chr22')


def gather_dir(tmp, data):
    d = tmp / 'gather'; d.mkdir(parents=True)
    write_table(d / 'carriers.qc.tsv.gz', FIELDS, data)
    write_table(d / 'samples.tsv', ['sample'], [dict(sample=s) for s in ['S1', 'S2', 'S3', 'S4']])
    strata = sorted({(r['candidate_type'], r['allele_class'], r['tier']) for r in data} | {('missense', 'sequence', 'miss_t1')})
    receipt = dict(stage='post_rarity_qc_gather', status='passed', reconciliation=dict(passed=True), chromosome='chr22', gather_id='autosomes',
                   selected_units=['chr22'], manifest_units=['chr22'], full_manifest_selected=True, policy=dict(version='test'),
                   psam_identity=dict(bytes=1, sha256='0' * 64), duplicate_checks=dict(cross_block_annotation_associations=0, cross_block_exact_alleles=0),
                   summary_dosage_field='qc_effective_alt_dosage', samples_in_source=4, carrier_annotation_records=len(data),
                   source_counts=dict(pass_rows=len(data)), by_stratum=[dict(candidate_type=t, allele_class=c, tier=x) for t, c, x in strata],
                   outputs={'carriers.qc.tsv.gz': identity(d / 'carriers.qc.tsv.gz'), 'samples.tsv': identity(d / 'samples.tsv')})
    (d / 'receipt.json').write_text(json.dumps(receipt))
    return d


def args(tmp, d, **kw):
    return Namespace(gather_dir=str(d), outdir=str(tmp / 'callqc'), zarr_dir=None, dnm_rescue=None, snv_site_ab_min=0.35, indel_site_ab_min=0.375) if not kw else \
        Namespace(**dict(dict(gather_dir=str(d), outdir=str(tmp / 'callqc'), zarr_dir=None, dnm_rescue=None, snv_site_ab_min=0.35, indel_site_ab_min=0.375), **kw))


DATA = [row(100, 'S1', 0.5), row(100, 'S2', 0.45),                     # SNV site mean 0.475 -> kept
        row(200, 'S1', 0.30), row(200, 'S2', 0.36),                    # SNV site mean 0.33 -> both removed
        row(300, 'S3', 0.36, ref='A', alt='AT'),                       # indel 0.36 < 0.375 -> removed
        row(400, 'S3', 0.38, ref='AT', alt='A'),                       # indel 0.38 -> kept
        row(500, 'S4', 0.26),                                          # SNV 0.26 but listed de novo -> exempt, kept
        row(600, 'S4', 0.10, alt='*', tier='untiered'),                # spanning deletion -> untouched
        row(100, 'S1', 0.5, typ='missense', tier='miss_t1', gene='G2')]  # second annotation of a kept allele


def test_site_ab_rule_exemption_and_products(tmp_path):
    d = gather_dir(tmp_path, DATA)
    rescue = tmp_path / 'rescue.tsv.gz'
    with gzip.open(rescue, 'wt') as f: f.write('sample\tCHROM\tPOS\tREF\tALT\nS4\tchr22\t500\tA\tG\n')
    a = args(tmp_path, d, dnm_rescue=str(rescue)); cq.run(a)
    out = Path(a.outdir); r = json.loads((out / 'receipt.json').read_text())
    kept = rows(out / 'carriers.qc.tsv.gz'); removed = rows(out / 'callqc_audit.tsv.gz')
    assert r['status'] == 'passed' and r['stage'] == 'post_rarity_call_qc'
    assert r['input_rows'] == 9 and r['carrier_annotation_records'] == len(kept) == 6 and r['removed_rows'] == len(removed) == 3
    assert {(x['POS'], x['sample']) for x in removed} == {('200', 'S1'), ('200', 'S2'), ('300', 'S3')} and r['removal_reasons'] == {'site_ab': 3}
    assert next(x for x in kept if x['POS'] == '500')['callqc_dnm_exempt'] == 'True' and r['dnm_exempt_rows'] == 1
    assert next(x for x in kept if x['POS'] == '100')['callqc_site_mean_ab'] == '0.475' and next(x for x in kept if x['POS'] == '100')['callqc_site_carriers'] == '2'
    assert next(x for x in kept if x['ALT'] == '*')['callqc_site_mean_ab'] == '.'
    assert all(x['callqc_strand_status'] == 'not_read' for x in kept if x['tier'] == 'lof_t1' and x['allele_class'] == 'sequence')
    assert r['call_qc_policy']['strand'].startswith('not applied')
    assert (out / 'carriers.qc.parquet').exists() and r['outputs']['carriers.qc.parquet'] == identity(out / 'carriers.qc.parquet')
    burden = {(x['sample'], x['tier']): x['carrier_variant_records'] for x in rows(out / 'sample_burden.tsv') if x['allele_class'] == 'sequence'}
    assert burden[('S1', 'lof_t1')] == '1' and burden[('S2', 'lof_t1')] == '1' and burden[('S3', 'lof_t1')] == '1' and burden[('S4', 'lof_t1')] == '1'
    assert {x['sample'] for x in rows(out / 'samples.tsv')} == {'S1', 'S2', 'S3', 'S4'}


def test_genome_gather_accepts_call_qc(tmp_path):
    d = gather_dir(tmp_path, DATA); a = args(tmp_path, d); cq.run(a)
    m = tmp_path / 'm.tsv'; m.write_text('chromosome\tdirectory\nchr22\t' + a.outdir + '\n')
    rec = gg.run(m, ['22'], tmp_path / 'genome')
    assert rec['source_stage'] == 'post_rarity_call_qc' and rec['carrier_annotation_records'] == json.loads((Path(a.outdir) / 'receipt.json').read_text())['carrier_annotation_records']
    assert rec['call_qc_policy']['snv_site_ab_min'] == 0.35


@pytest.mark.parametrize('case', ['stage', 'checksum', 'count'])
def test_fail_closed(tmp_path, case):
    d = gather_dir(tmp_path, DATA); r = json.loads((d / 'receipt.json').read_text())
    if case == 'stage': r['stage'] = 'post_rarity_qc'
    if case == 'checksum': r['outputs']['carriers.qc.tsv.gz']['sha256'] = '1' * 64
    if case == 'count': r['carrier_annotation_records'] = 1
    (d / 'receipt.json').write_text(json.dumps(r))
    a = args(tmp_path, d)
    with pytest.raises(ValidationError): cq.run(a)
    assert json.loads((Path(a.outdir) / 'receipt.json').read_text())['status'] == 'failed'
    assert not (Path(a.outdir) / 'carriers.qc.tsv.gz').exists()


@pytest.mark.parametrize('sb,gt,ai,status', [
    ([10, 10, 0, 8], [0, 1], 1, 'one_strand_alt_ref_both'), ([10, 10, 4, 4], [0, 1], 1, 'informative_retained'),
    ([-1, -1, -1, -1], [0, 1], 1, 'missing'), ([10, 10, 0, 8], [1, 2], 1, 'ambiguous_multiallelic'),
    ([0, 0, 3, 4], [0, 1], 1, 'uninformative'), ([10, 0, 0, 8], [0, 1], 1, 'informative_retained')])
def test_strand_status(sb, gt, ai, status):
    assert cq.strand_status(sb, gt, ai) == status


def test_strand_rule_from_zarr(tmp_path):
    zarr = pytest.importorskip('zarr'); np = pytest.importorskip('numpy')
    store = tmp_path / 'zarr'; store.mkdir()
    g = zarr.open_group(str(store / 'chr22.no-pl.zarr'), mode='w')
    n_var = 700
    g.create_array('sample_id', data=np.asarray(['S1', 'S2', 'S3', 'S4']))
    al = np.full((n_var, 3), '', dtype=object); al[:, 0] = 'A'; al[:, 1] = 'G'; al[300, 1] = 'AT'; al[400, 1] = 'A'; al[100, 2] = 'C'
    g.create_array('variant_allele', data=np.asarray(al.tolist()))
    gt = np.zeros((n_var, 4, 2), dtype='i1'); sb = np.full((n_var, 4, 4), 5, dtype='i2')
    for r in DATA: gt[int(r['POS']), int(r['sample'][1]) - 1] = [0, 1]
    sb[500, 3] = [10, 10, 0, 8]   # one-strand ALT on the exempt de novo -> still kept (exempt)
    sb[100, 0] = [10, 10, 0, 8]   # one-strand ALT, S1 at 100 -> removed by strand rule
    g.create_array('call_genotype', data=gt); g.create_array('call_SB', data=sb)
    d = gather_dir(tmp_path, DATA); a = args(tmp_path, d, zarr_dir=str(store)); cq.run(a)
    out = Path(a.outdir); r = json.loads((out / 'receipt.json').read_text()); removed = rows(out / 'callqc_audit.tsv.gz')
    s1 = [x for x in removed if (x['POS'], x['sample']) == ('100', 'S1')]
    assert s1 and all('strand' in x['callqc_reason'] for x in s1) and len(s1) == 2   # both annotations of the allele/sample
    assert r['strand_status']['one_strand_alt_ref_both'] >= 1 and r['zarr_sb_reads'] > 0
    assert next(x for x in rows(out / 'carriers.qc.tsv.gz') if x['POS'] == '100' and x['sample'] == 'S2')['callqc_n_alt_alleles'] == '2'


def test_nextflow_call_qc_workflow(tmp_path):
    import shutil, subprocess
    if not shutil.which('nextflow'): pytest.skip('nextflow not available')
    root = Path(__file__).resolve().parents[1]
    d = gather_dir(tmp_path, DATA)
    m = tmp_path / 'manifest.tsv'; m.write_text('chromosome\tdirectory\nchr22\t' + str(d) + '\n')
    cfg = tmp_path / 'python.config'; cfg.write_text('process.beforeScript = "export PATH=' + str(Path(sys.executable).parent) + ':\\$PATH"\n')
    r = subprocess.run(['nextflow', '-C', str(root / 'post_rarity_call_qc.config') + ',' + str(cfg), 'run', str(root / 'post_rarity_call_qc.nf'), '-ansi-log', 'false',
                        '--call_qc_manifest', str(m), '--outdir', str(tmp_path / 'published'), '--call_qc_memory', '1 GB'],
                       cwd=tmp_path, capture_output=True, text=True, timeout=180)
    assert r.returncode == 0, r.stdout + r.stderr
    rec = json.loads((tmp_path / 'published/post-rarity-call-qc/chr22/out/receipt.json').read_text())
    assert rec['status'] == 'passed' and rec['removed_rows'] == 4 and rec['dnm_rescue_identity'] is None
