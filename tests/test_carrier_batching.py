import gzip
import json
from pathlib import Path
from types import SimpleNamespace
import pytest
from test_filtered_carriers import fixture, extract, rows
from extract_exact_carriers import candidate_windows, exact_records, ValidationError
from carrier_psam import load_psam


def test_batched_full_output_parity(tmp_path):
    a=fixture(tmp_path)
    a.batch_bp=1
    extract(a)
    original=Path(a.outdir)
    a.outdir=str(tmp_path/'batched');a.batch_bp=10000
    extract(a)
    for p in original.iterdir():
        if p.name=='receipt.json':continue
        other=Path(a.outdir)/p.name
        assert (gzip.decompress(p.read_bytes()) if p.suffix=='.gz' else p.read_bytes()) == (gzip.decompress(other.read_bytes()) if p.suffix=='.gz' else other.read_bytes())
    assert json.loads((Path(a.outdir)/'receipt.json').read_text())['query_batch_bp']==10000


def test_windows_exact_matching_duplicate_and_boundary_guards():
    keys={('22',p,'A','G'):{} for p in [1,9999,10000,10001,20000,30001]}
    assert list(candidate_windows(keys,10000))==[(0,10000),(10000,20000),(30000,30001)]
    def rec(pos,alt='G'):
        return SimpleNamespace(contig='chr22',pos=pos,ref='A',alts=(alt,),format={'GT':1})
    class VCF:
        def __init__(self,records):self.records=records;self.fetches=[]
        def fetch(self,c,start,end):
            self.fetches.append((start,end))
            # Include spanning overlaps (POS < start) as a real indexed fetch can.
            return (r for r in self.records if r.pos<=end)
    unrelated=rec(2);unrelated.alts=('T','C')
    records=[rec(1),unrelated,rec(9999,'T'),rec(10000),rec(10001),rec(20000),rec(30001)]
    v=VCF(records);got=list(exact_records(v,'chr22',keys,10000))
    assert [k[1] for k,r in got]==[1,10000,10001,20000,30001]
    assert len(v.fetches)==3
    with pytest.raises(ValidationError,match='Duplicate exact'):
        list(exact_records(VCF([rec(1),rec(1)]),'chr22',keys,10000))
    multi=rec(10001);multi.alts=('T','G')
    with pytest.raises(ValidationError,match='Multiallelic'):
        list(exact_records(VCF([multi]),'chr22',keys,10000))
    for value in [0,-1,True,1.5]:
        with pytest.raises(ValidationError):list(candidate_windows({},value))


def psam(path):
    path.write_text('#IID\tSEX\tparticipant_id\tfrequency_representative\tunrelated\n'
        'S4\t0\tP4\t1\t1\nS2\t2\tP1\t0\t1\nS1\t1\tP1\t1\t0\nS3\tNA\tP3\t1\t1\n')
    return path


def test_metadata_does_not_subset_carriers_or_change_order(tmp_path):
    a=fixture(tmp_path);a.psam=str(psam(tmp_path/'samples.psam'))
    extract(a);out=Path(a.outdir)
    assert len(rows(out/'carriers.tsv.gz'))==9
    assert [r['#IID'] for r in rows(out/'sample_metadata.tsv')]==['S1','S2','S3','S4']
    assert len(rows(out/'samples.tsv'))==4
    assert any(r['sample']=='S2' for r in rows(out/'carriers.tsv.gz'))
    r=json.loads((out/'receipt.json').read_text())
    assert r['psam']['representative_flags']==3 and r['psam']['unrelated_flags']==3
    assert r['input_identities']['psam']['sha256']
    assert r['frequency_status'].startswith('deferred')


@pytest.mark.parametrize('change',['duplicate','missing','flag','header'])
def test_psam_errors(change,tmp_path):
    p=psam(tmp_path/'x.psam');s=p.read_text()
    if change=='duplicate':s+='S1\t1\tP1\t1\t0\n'
    if change=='missing':s=s.replace('S4\t0\tP4\t1\t1\n','')
    if change=='flag':s=s.replace('S1\t1\tP1\t1\t0','S1\t1\tP1\tmaybe\t0')
    if change=='header':s=s.replace('#IID','sample')
    p.write_text(s)
    with pytest.raises(ValidationError):load_psam(p,['S1','S2','S3','S4'])


@pytest.mark.parametrize('chromosome',['X','Y'])
def test_sex_chromosomes_raw_extraction_without_frequency_policy(tmp_path,chromosome):
    import pysam
    from test_filtered_carriers import mutate
    a=fixture(tmp_path)
    source=tmp_path/'original.vcf'
    s=source.read_text().replace('chr22','chr'+chromosome)
    s=s.replace('##fileformat=VCFv4.2\n','##fileformat=VCFv4.2\n##INFO=<ID=AC,Number=A,Type=Integer,Description="Source AC">\n##INFO=<ID=AN,Number=1,Type=Integer,Description="Source AN">\n##INFO=<ID=AF,Number=A,Type=Float,Description="Source AF">\n')
    lines=s.splitlines()
    for i,line in enumerate(lines):
        if line.startswith('chr'+chromosome+'\t10\t'):
            fields=line.split('\t');fields[7]='AC=3;AN=8;AF=0.25';lines[i]='\t'.join(fields)
    source.write_text('\n'.join(lines)+'\n')
    pysam.tabix_compress(str(source),a.vcf,force=True);pysam.tabix_index(a.vcf,preset='vcf',force=True)
    for p in [a.missense,a.lof_hc]:mutate(p,f"UPDATE c SET CHROM='{chromosome}'")
    meta=json.loads(Path(a.metadata).read_text());meta['chromosome']='chr'+chromosome;Path(a.metadata).write_text(json.dumps(meta))
    extract(a)
    out=Path(a.outdir)
    assert len(rows(out/'carriers.tsv.gz'))==9
    frequency=next(r for r in rows(out/'source_frequencies.tsv') if r['POS']=='10')
    assert frequency['source_info_ac']=='3' and frequency['source_info_an']=='8'
    assert float(frequency['source_info_af'])==0.25
    assert float(frequency['source_ac_an'])==0.375
    assert json.loads((out/'receipt.json').read_text())['frequency_status'].startswith('deferred')
