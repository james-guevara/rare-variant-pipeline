import gzip
import json
from pathlib import Path
import shutil
import subprocess
import sys

import pytest

from test_filtered_carriers import fixture, extract, rows, mutate
from extract_exact_carriers import ValidationError
from carrier_frequencies import count_call, FrequencyCounter

ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.parametrize('gt,ploidy,expected', [
    ((0,0),2,(0,2,'counted_diploid_complete')),
    ((0,1),2,(1,2,'counted_diploid_complete')),
    ((1,1),2,(2,2,'counted_diploid_complete')),
    ((0,None),2,(0,1,'counted_diploid_partial')),
    ((None,1),2,(1,1,'counted_diploid_partial')),
    ((None,None),2,(0,0,'excluded_missing')),
    ((1,),2,(0,0,'excluded_unexpected_ploidy')),
    ((0,),2,(0,0,'excluded_unexpected_ploidy')),
    ((0,1,1),2,(0,0,'excluded_unexpected_ploidy')),
    ((0,2),2,(0,0,'excluded_invalid_allele')),
    ((0,),1,(0,1,'counted_haploid')),
    ((1,),1,(1,1,'counted_haploid')),
    ((0,0),1,(0,1,'counted_diploid_encoded_haploid')),
    ((1,1),1,(1,1,'counted_diploid_encoded_haploid')),
    ((0,1),1,(0,0,'excluded_haploid_heterozygous')),
    ((1,0),1,(0,0,'excluded_haploid_heterozygous')),
    ((1,None),1,(0,0,'excluded_haploid_partial')),
    ((None,0),1,(0,0,'excluded_haploid_partial')),
    ((None,),1,(0,0,'excluded_missing')),
    ((None,None),1,(0,0,'excluded_missing')),
    ((1,1,1),1,(0,0,'excluded_unexpected_ploidy')),
    ((0,0),None,(0,0,'excluded_unknown_ploidy')),
    (None,2,(0,0,'excluded_missing')),
])
def test_agreed_genotype_counting_rules(gt,ploidy,expected):
    assert count_call(gt,ploidy)==expected


def psam(path):
    path.write_text('#IID\tSEX\tparticipant_id\tfrequency_representative\tunrelated\n'
        'S4\t0\tP4\t1\t1\nS2\t2\tP1\t0\t1\nS1\t1\tP1\t1\t0\nS3\tNA\tP3\t1\t1\n')
    return path


@pytest.mark.parametrize('index',['tbi','csi'])
def test_all_eligible_calls_reference_partial_no_qc_and_intersection(tmp_path,index,capsys):
    a=fixture(tmp_path,index);a.psam=str(psam(tmp_path/'samples.psam'))
    extract(a)
    raw=Path(a.outdir)
    a.compute_frequencies=True;a.outdir=str(tmp_path/'with-frequencies')
    extract(a);out=Path(a.outdir)
    # All existing output records stay identical, including nonrepresentative carriers.
    for p in raw.iterdir():
        if p.name=='receipt.json':continue
        q=out/p.name
        assert (gzip.decompress(p.read_bytes()) if p.suffix=='.gz' else p.read_bytes()) == (gzip.decompress(q.read_bytes()) if q.suffix=='.gz' else q.read_bytes())
    freqs={int(r['POS']):r for r in rows(out/'variant_frequencies.tsv')}
    # The allele at 10 has two annotations, but one frequency row.
    assert len(freqs)==6
    r=freqs[10]
    assert (r['cohort_ac'],r['cohort_an'],r['cohort_af'])==('1','4','0.25')
    assert (r['unrelated_ac'],r['unrelated_an'],r['unrelated_af'])==('0','2','0.0')
    assert r['cohort_reference_genotypes']=='1'
    assert freqs[13]['cohort_ac']=='0' and freqs[13]['cohort_an']=='6'  # reference-only, no GQ/DP/AD fields
    assert freqs[20]['cohort_ac']==freqs[20]['cohort_an']=='1'  # q10 site, partial ALT, no GQ/DP/AD
    assert freqs[20]['unrelated_an']=='0' and freqs[20]['unrelated_af']=='.'
    assert freqs[60]['allele_class']=='spanning_deletion'
    assert freqs[60]['cohort_ac']=='2' and freqs[60]['cohort_an']=='3'
    assert freqs[61]['matched']=='False' and freqs[61]['cohort_ac']==freqs[61]['cohort_an']=='.'
    receipt=json.loads((out/'receipt.json').read_text());f=receipt['frequencies']
    assert f['sample_selection']['cohort']==3 and f['sample_selection']['unrelated']==2
    assert f['sample_selection']['unrelated_flags_without_representative']==1
    assert f['by_allele_class']['sequence']['matched_distinct_alleles']==4
    assert f['by_allele_class']['spanning_deletion']['matched_distinct_alleles']==1
    seq=f['by_allele_class']['sequence']['sample_sets']
    assert seq['cohort']['evaluated_genotypes']==12 and seq['unrelated']['evaluated_genotypes']==8
    assert seq['cohort']['by_reason']['excluded_unexpected_ploidy']==1
    for cls in f['by_allele_class'].values():
        for v in cls['sample_sets'].values():
            assert sum(v['by_reason'].values())==v['expected_evaluated_genotypes']
    audit=rows(out/'frequency_audit.tsv')
    assert sum(int(x['genotypes']) for x in audit if x['sample_set']=='cohort')==15
    assert len(rows(out/'samples.tsv'))==4
    assert 'S1' not in capsys.readouterr().out


@pytest.mark.parametrize('case',['empty','no_representatives','duplicate_representative','missing_psam','X','Y'])
def test_empty_counts_and_explicit_errors(tmp_path,case):
    a=fixture(tmp_path);a.compute_frequencies=True;a.psam=str(psam(tmp_path/'samples.psam'))
    if case=='empty':
        for p in [a.missense,a.lof_hc]:mutate(p,'DELETE FROM c')
        a.expected_missense=a.expected_hc=0
        extract(a)
        assert rows(Path(a.outdir)/'variant_frequencies.tsv')==[]
        assert len(rows(Path(a.outdir)/'samples.tsv'))==4
        return
    if case=='no_representatives':
        p=Path(a.psam);p.write_text(p.read_text().replace('\t1\t1','\t0\t1').replace('P1\t1\t0','P1\t0\t0'))
        extract(a)
        for r in rows(Path(a.outdir)/'variant_frequencies.tsv'):
            if r['matched']=='True':assert r['cohort_an']=='0' and r['cohort_af']=='.'
        return
    if case=='duplicate_representative':
        p=Path(a.psam);p.write_text(p.read_text().replace('S2\t2\tP1\t0','S2\t2\tP1\t1'))
    if case=='missing_psam':a.psam=None
    if case in ('X','Y'):
        m=json.loads(Path(a.metadata).read_text());m['chromosome']='chr'+case;Path(a.metadata).write_text(json.dumps(m))
    with pytest.raises(ValidationError):extract(a)
    out=Path(a.outdir)
    assert json.loads((out/'receipt.json').read_text())['status']=='failed'
    assert not (out/'carriers.tsv.gz').exists() and not (out/'variant_frequencies.tsv').exists()


@pytest.mark.skipif(not shutil.which('nextflow'),reason='Nextflow required')
def test_nextflow_frequency_wiring_resume_and_psam_change(tmp_path):
    a=fixture(tmp_path/'data');p=psam(tmp_path/'samples.psam')
    manifest=tmp_path/'manifest.tsv'
    manifest.write_text('unit_id\tchromosome\tmissense\tlof_hc\tvcf\tindex\n'+f'block12\tchr22\t{a.missense}\t{a.lof_hc}\t{a.vcf}\t{a.index}\n')
    cfg=tmp_path/'python.config';cfg.write_text('process.beforeScript = "export PATH='+str(Path(sys.executable).parent)+':\\$PATH"\n')
    cmd=['nextflow','-C',str(ROOT/'filtered_carriers.config')+','+str(cfg),'run',str(ROOT/'filtered_carriers.nf'),'-ansi-log','false',
        '--carrier_manifest',str(manifest),'--outdir',str(tmp_path/'published'),'--select_units','all',
        '--psam',str(p),'--compute_frequencies','true','--carrier_memory','1 GB']
    for n,expected in [(1,'COMPLETED'),(2,'CACHED'),(3,'COMPLETED')]:
        if n==3:
            p.write_text(p.read_text().replace('S1\t1\tP1\t1\t0','S1\t1\tP1\t1\t1'))
        run=subprocess.run(cmd+['-with-trace',str(tmp_path/f'{n}.trace')]+(['-resume'] if n>1 else []),cwd=tmp_path,capture_output=True,text=True,timeout=120)
        assert run.returncode==0,run.stdout+run.stderr
        assert [r['status'] for r in rows(tmp_path/f'{n}.trace')]==[expected]
    out=tmp_path/'published/filtered-carriers/block12'
    r=next(x for x in rows(out/'variant_frequencies.tsv') if x['POS']=='10')
    assert r['unrelated_ac']=='1' and r['unrelated_an']=='4'
    receipt=json.loads((out/'receipt.json').read_text());assert receipt['query_batch_bp']==10000
    assert receipt['frequencies']['sample_selection']['unrelated']==3
    shutil.rmtree(tmp_path/'work')
    assert (out/'variant_frequencies.tsv').is_file() and not (out/'variant_frequencies.tsv').is_symlink()


def test_frequency_counter_uses_gt_only_and_counts_phased_hom_alt(tmp_path):
    a=fixture(tmp_path)
    import pysam
    # Independent participants, all representative; raw source includes phased 1|1.
    metadata=[{'#IID':f'S{i}','participant_id':f'P{i}','frequency_representative':'1','unrelated':'1','SEX':'0'} for i in range(1,5)]
    counter=FrequencyCounter(metadata,'22')
    with pysam.VariantFile(a.vcf,index_filename=a.index) as source:
        record=next(source.fetch('chr22',9,10))
        value,_=counter.count(record,'sequence')
    assert value['cohort_ac']==3 and value['cohort_an']==6
    class GTOnly(dict):
        def get(self,key):
            assert key=='GT', 'Frequency counting must not inspect quality fields'
            return (0,0)
    class Record:
        samples={f'S{i}':GTOnly() for i in range(1,5)}
    value,_=counter.count(Record(),'sequence')
    assert value['cohort_ac']==0 and value['cohort_an']==8
