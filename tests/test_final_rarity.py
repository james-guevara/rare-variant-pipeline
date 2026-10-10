from argparse import Namespace
import builtins
import json
from pathlib import Path
import shutil
import subprocess
import sys

import pytest
from test_filtered_carriers import fixture as extraction_fixture, extract, rows
from test_carrier_frequencies import psam
from qc_filtered_carriers import identity, write_table
from filter_final_rarity import run, decision, ValidationError
from make_final_rarity_manifest import build

ROOT=Path(__file__).resolve().parents[1]


def fixture(tmp):
    a=extraction_fixture(tmp);a.psam=str(psam(tmp/'samples.psam'));a.compute_frequencies=True
    extract(a);out=Path(a.outdir)
    freq=rows(out/'variant_frequencies.tsv')
    for r in freq:
        if r['matched']=='True':
            r['unrelated_an']='1000000'
            r['unrelated_af']={'10':'0.000999','13':'0','20':'0.001','30':'.','60':'0'}[r['POS']]
    write_table(out/'variant_frequencies.tsv',list(freq[0]),freq)
    refresh(out,'variant_frequencies.tsv')
    return Namespace(metadata=a.metadata,carriers=str(out/'carriers.tsv.gz'),samples=str(out/'samples.tsv'),
        candidates=str(out/'candidates.tsv'),frequencies=str(out/'variant_frequencies.tsv'),source_receipt=str(out/'receipt.json'),outdir=str(tmp/'rare'))


def refresh(out,name):
    p=out/'receipt.json';r=json.loads(p.read_text());r['outputs'][name]=identity(out/name);p.write_text(json.dumps(r))


@pytest.mark.parametrize('af,an,matched,reason',[
    ('0.000999999999999999','10','True','PASS'),('0.001','10','True','af_at_or_above_threshold'),
    ('0.001000000000000001','10','True','af_at_or_above_threshold'),('0','10','True','PASS'),
    ('.','10','True','missing_af'),('nan','10','True','invalid_af'),('inf','10','True','invalid_af'),
    ('-0.01','10','True','invalid_af'),('1.1','10','True','invalid_af'),('broken','10','True','invalid_af'),
    ('0','0','True','zero_an'),('.','0','True','zero_an'),('0','.','True','missing_an'),
    ('0','-1','True','invalid_an'),('0','1.5','True','invalid_an'),('0','10','False','unmatched'),
])
def test_strict_threshold_and_missing_fail_closed(af,an,matched,reason):
    assert decision(dict(unrelated_af=af,unrelated_an=an,matched=matched))==reason


def test_complete_annotation_consistency_stars_tiers_zeros_and_no_genotype_reads(tmp_path,monkeypatch,capsys):
    a=fixture(tmp_path);inputs=[Path(getattr(a,n)) for n in ['carriers','samples','candidates','frequencies','source_receipt']]
    before={p:identity(p) for p in inputs}
    real_open=builtins.open
    def guarded(file,*args,**kwargs):
        assert str(file) not in [str(tmp_path/'source.vcf.gz'),str(tmp_path/'original.vcf')], 'No genotype VCF reads'
        return real_open(file,*args,**kwargs)
    monkeypatch.setattr(builtins,'open',guarded)
    import qc_filtered_carriers
    monkeypatch.setattr(qc_filtered_carriers,'evaluate',lambda *args:pytest.fail('No genotype QC'))
    capsys.readouterr();run(a);out=Path(a.outdir)
    carriers=rows(out/'carriers.rare.tsv.gz');assert len(carriers)==7
    assert len([r for r in carriers if r['POS']=='10'])==4
    assert {r['candidate_type'] for r in carriers if r['POS']=='10'}=={'missense','lof_hc'}
    assert not any(r['POS']=='20' for r in carriers)
    assert any(r['GT']=='1/.' and r['allele_class']=='spanning_deletion' for r in carriers)
    assert len(rows(out/'candidates.rare.tsv'))==4
    assert any(r['tier']=='untiered' for r in rows(out/'candidates.rare.tsv'))
    assert [r['sample'] for r in rows(out/'samples.tsv')]==['S1','S2','S3','S4']
    summary=rows(out/'sample_burden.tsv')
    assert all(r['carrier_variant_records']=='0' for r in summary if r['sample']=='S4')
    assert rows(out/'sample_distinct_alleles.tsv')[1]['sequence_distinct_variant_count']=='1'
    r=json.loads((out/'receipt.json').read_text())
    assert r['distinct_alleles']==dict(input=6,pass_count=3,fail_count=3)
    assert r['candidate_annotations']==dict(input=7,pass_count=4,fail_count=3)
    assert r['carrier_annotations']==dict(input=9,pass_count=7,fail_count=2)
    assert r['by_allele_class']['sequence']['distinct_allele_sample_records']['pass_count']==2
    assert r['by_allele_class']['spanning_deletion']['distinct_allele_sample_records']['pass_count']==3
    assert r['reconciliation']['passed']
    assert all(identity(p)==v for p,v in before.items())
    assert 'S1' not in capsys.readouterr().out


@pytest.mark.parametrize('case',['all_fail','duplicate_frequency','different_alt','changed_input','missing_frequency','conflicting_annotation'])
def test_failures_and_all_fail(tmp_path,case):
    a=fixture(tmp_path);out=Path(a.frequencies).parent
    if case in ['all_fail','duplicate_frequency','different_alt','missing_frequency']:
        data=rows(a.frequencies);fields=list(data[0])
        if case=='all_fail':
            for r in data:r['unrelated_af']='0.01'
        if case=='duplicate_frequency':data.append(data[0])
        if case=='different_alt':data[0]['ALT']='T'
        if case=='missing_frequency':data.pop(0)
        write_table(a.frequencies,fields,data);refresh(out,'variant_frequencies.tsv')
    if case=='changed_input':Path(a.samples).write_text('sample\nS1\n')
    if case=='conflicting_annotation':
        data=rows(a.carriers);data[0]['Gene']='DIFFERENT';write_table(a.carriers,list(data[0]),data);refresh(out,'carriers.tsv.gz')
    if case=='all_fail':
        run(a);assert rows(Path(a.outdir)/'carriers.rare.tsv.gz')==[]
        assert len(rows(Path(a.outdir)/'samples.tsv'))==4
        assert all(r['carrier_variant_records']=='0' for r in rows(Path(a.outdir)/'sample_burden.tsv'))
    else:
        with pytest.raises(ValidationError):run(a)
        assert json.loads((Path(a.outdir)/'receipt.json').read_text())['status']=='failed'
        assert not (Path(a.outdir)/'carriers.rare.tsv.gz').exists()


def test_manifest_and_duplicate_roots(tmp_path):
    a=fixture(tmp_path/'unit')
    data=build([tmp_path/'unit']);assert len(data)==1 and data[0][0]=='block12'
    with pytest.raises(ValueError):build([tmp_path/'unit',tmp_path/'unit'])


@pytest.mark.skipif(not shutil.which('nextflow'),reason='Nextflow required')
def test_nextflow_saved_products_only_resume_and_trace_overwrite(tmp_path):
    a=fixture(tmp_path/'unit');p=Path(a.source_receipt).parent
    manifest=tmp_path/'manifest.tsv';manifest.write_text('unit_id\tchromosome\tcarriers\tsamples\tsource_receipt\tcandidates\tfrequencies\n'+
        '\t'.join(['block12','chr22',a.carriers,a.samples,a.source_receipt,a.candidates,a.frequencies])+'\n')
    # Remove genotype inputs entirely before Nextflow; saved outputs suffice.
    for f in (tmp_path/'unit').glob('source.vcf.gz*'):f.unlink()
    cfg=tmp_path/'python.config';cfg.write_text('process.beforeScript = "export PATH='+str(Path(sys.executable).parent)+':\\$PATH"\n')
    cmd=['nextflow','-C',str(ROOT/'final_rarity.config')+','+str(cfg),'run',str(ROOT/'final_rarity.nf'),'-ansi-log','false',
         '--rarity_manifest',str(manifest),'--select_units','all','--outdir',str(tmp_path/'published'),
         '--rarity_memory','1 GB','-with-trace',str(tmp_path/'trace.tsv')]
    for resume,status in [(False,'COMPLETED'),(True,'CACHED')]:
        result=subprocess.run(cmd+(['-resume'] if resume else []),cwd=tmp_path,capture_output=True,text=True,timeout=120)
        assert result.returncode==0,result.stdout+result.stderr
        trace=rows(tmp_path/'trace.tsv');assert len(trace)==1 and trace[0]['status']==status
    shutil.rmtree(tmp_path/'work')
    out=tmp_path/'published/final-rarity/block12'
    assert len(rows(out/'carriers.rare.tsv.gz'))==7 and not (out/'carriers.rare.tsv.gz').is_symlink()


def test_empty_candidates_keep_source_samples(tmp_path):
    from test_filtered_carriers import mutate
    a=extraction_fixture(tmp_path);a.psam=str(psam(tmp_path/'samples.psam'));a.compute_frequencies=True
    for p in [a.missense,a.lof_hc]:mutate(p,'DELETE FROM c')
    a.expected_hc=a.expected_missense=0;extract(a);out=Path(a.outdir)
    final=Namespace(metadata=a.metadata,carriers=str(out/'carriers.tsv.gz'),samples=str(out/'samples.tsv'),
        candidates=str(out/'candidates.tsv'),frequencies=str(out/'variant_frequencies.tsv'),source_receipt=str(out/'receipt.json'),outdir=str(tmp_path/'rare'))
    run(final)
    r=json.loads((Path(final.outdir)/'receipt.json').read_text())
    assert r['distinct_alleles']==dict(input=0,pass_count=0,fail_count=0)
    assert len(rows(Path(final.outdir)/'sample_distinct_alleles.tsv'))==4
    assert all(x['carrier_variant_records']=='0' for x in rows(Path(final.outdir)/'sample_burden.tsv'))


@pytest.mark.parametrize('chromosome',['X','Y'])
def test_xy_rarity_does_not_apply_frequency_genotype_exclusions(tmp_path,chromosome):
    from test_carrier_sex_frequencies import fixture as xy_fixture
    a=xy_fixture(tmp_path,chromosome);extract(a);out=Path(a.outdir)
    data=rows(out/'variant_frequencies.tsv')
    for r in data:r['unrelated_af']='0';r['unrelated_an']='10000'
    write_table(out/'variant_frequencies.tsv',list(data[0]),data);refresh(out,'variant_frequencies.tsv')
    final=Namespace(metadata=a.metadata,carriers=str(out/'carriers.tsv.gz'),samples=str(out/'samples.tsv'),
        candidates=str(out/'candidates.tsv'),frequencies=str(out/'variant_frequencies.tsv'),source_receipt=str(out/'receipt.json'),outdir=str(tmp_path/'rare'))
    run(final)
    before=rows(final.carriers);after=rows(Path(final.outdir)/'carriers.rare.tsv.gz')
    assert before==after
    assert any(r['sample']=='U1' for r in after)
    assert any(r['sample']=='F1' for r in after)
    assert any(r['GT']=='1/.' for r in after)
    assert any(r['alt_dosage']=='2' and r['sample']=='M1' for r in after)
