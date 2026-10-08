import gzip
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from argparse import Namespace
import duckdb
import pysam
ROOT=Path(__file__).resolve().parents[3]
sys.path[:0]=[str(ROOT/'benchmarks/carrier'),str(ROOT/'scripts')]
from benchmark import compare,run


def fixture(root,index='tbi',float_gq=False,empty=False):
    root.mkdir(parents=True,exist_ok=True)
    vcf=root/'source.vcf'
    header='''##fileformat=VCFv4.2
##contig=<ID=chr22,length=1000>
##FILTER=<ID=q10,Description="Synthetic">
##FILTER=<ID=Other,Description="Synthetic">
##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">
##FORMAT=<ID=GQ,Number=1,Type=GQTYPE,Description="GQ">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="DP">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="AD">
##FORMAT=<ID=FT,Number=1,Type=String,Description="FT">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4\tS5\tZERO
'''.replace('GQTYPE','Float' if float_gq else 'Integer')
    lines=[
        'chr22\t10\t.\tA\tG\t.\tPASS\t.\tGT:GQ:DP:AD:FT\t0/1:20:10:5,5:PASS\t1|1:.:11:.,11:.\t1:3:2:1,1:LowGQ\t1/.:2:3:.:.\t./1:30:20:20,.:PASS\t0/0:20:10:10,0:PASS',
        'chr22\t10\t.\tA\tT\t.\tPASS\t.\tGT\t1/1\t1/1\t1/1\t1/1\t1/1\t1/1',
        'chr22\t20\t.\tAC\t*\t.\tq10;Other\t.\tGT\t1|.\t.|1\t1\t0/.\t./.\t0/0',
        'chr22\t30\t.\tG\tC\t.\t.\t.\tGT:AD\t0/1:4,6\t0/0:10,0\t0:9,0\t.:.\t./.:.\t0/0:10,0',
        'chr22\t40\t.\tT\tA\t.\tPASS\t.\tGT\t0/0\t0/0\t0\t0/.\t./.\t0/0',
        'chr22\t50\t.\tA\tC\t.\tPASS\t.\tGT\t1/0/1\t1|1|0\t1\t1/.\t./1\t0/0',
    ]
    if float_gq:lines[0]=lines[0].replace('0/1:20:','0/1:20.123:').replace('1:3:','1:nan:').replace('1/.:2:','1/.:inf:')
    vcf.write_text(header+'\n'.join(lines)+'\n')
    gz=root/'source.vcf.gz';pysam.tabix_compress(str(vcf),str(gz),force=True)
    pysam.tabix_index(str(gz),preset='vcf',force=True,csi=index=='csi')
    c=duckdb.connect();c.execute('CREATE TABLE c(CHROM VARCHAR,POS BIGINT,REF VARCHAR,ALT VARCHAR,Gene VARCHAR,Feature VARCHAR,SYMBOL VARCHAR,Consequence VARCHAR,LoF VARCHAR,tier VARCHAR,allele_class VARCHAR,pcf_retained BOOLEAN)')
    for kind in ['missense','lof_hc']:
        c.execute('DELETE FROM c')
        for pos,ref,alt in ([] if empty else ([(10,'A','G'),(30,'G','C')] if kind=='missense' else [(10,'A','G'),(20,'AC','*'),(40,'T','A'),(50,'A','C'),(70,'C','*')])):
            c.execute('INSERT INTO c VALUES (?,?,?,?,?,?,?,?,?,?,?,?)',['22',pos,ref,alt,'GENE_'+kind,'TX_'+kind,'SYM','missense_variant' if kind=='missense' else 'frameshift_variant',None if kind=='missense' else 'HC','miss_t1' if kind=='missense' else None,'spanning_deletion' if alt=='*' else 'sequence',True])
        c.execute('COPY c TO ? (FORMAT PARQUET)',[str(root/(kind+'.parquet'))])
    c.close()
    meta=root/'metadata.json';meta.write_text(json.dumps(dict(unit_id='synthetic',chromosome='chr22')))
    return Namespace(metadata=str(meta),missense=str(root/'missense.parquet'),lof_hc=str(root/'lof_hc.parquet'),vcf=str(gz),index=str(gz)+'.'+index,
                     expected_missense=None,expected_hc=None,outdir=str(root/'benchmark'),pairs=2,container_receipt=None)


def backend(a,name,out):
    args=[sys.executable,str(ROOT/'benchmarks/carrier/run_backend.py'),'--backend',name,'--outdir',str(out)]
    for field in ['metadata','missense','lof_hc','vcf','index']:args+=['--'+field.replace('_','-'),str(getattr(a,field))]
    return subprocess.run(args,capture_output=True,text=True)


class Parity(unittest.TestCase):
    def test_formats_partial_haploid_polyploid_alias_stars_and_indexes(self):
        for index in ['tbi','csi']:
            for floats in [False,True]:
                with self.subTest(index=index,floats=floats),tempfile.TemporaryDirectory() as d:
                    root=Path(d);a=fixture(root,index,float_gq=floats)
                    for b in ['pysam','cyvcf2']:
                        p=backend(a,b,root/b);self.assertEqual(p.returncode,0,p.stderr)
                    result=compare(root/'pysam',root/'cyvcf2');self.assertTrue(result['passed'],result)
                    report=json.loads((root/'pysam/receipt.json').read_text())
                    self.assertEqual(report['overlapping_type_alleles'],1)
                    self.assertEqual(report['samples_in_source'],6)
                    self.assertEqual(report['distinct_allele_audit_by_class']['spanning_deletion']['carrier_variant_sample_records'],3)

    def test_out_of_range_gt_matches_pysam_missing_alleles(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);a=fixture(root);p=root/'source.vcf'
            p.write_text(p.read_text().replace('0/1:20:','1/2:20:'))
            pysam.tabix_compress(str(p),a.vcf,force=True);pysam.tabix_index(a.vcf,preset='vcf',force=True)
            for b in ['pysam','cyvcf2']:self.assertEqual(backend(a,b,root/b).returncode,0)
            self.assertTrue(compare(root/'pysam',root/'cyvcf2')['passed'])

    def test_mixed_polyploid_phasing(self):
        for gt in ['1|1/0','1/1|0']:
            with self.subTest(gt=gt),tempfile.TemporaryDirectory() as d:
                root=Path(d);a=fixture(root);p=root/'source.vcf'
                p.write_text(p.read_text().replace('1|1|0',gt))
                pysam.tabix_compress(str(p),a.vcf,force=True);pysam.tabix_index(a.vcf,preset='vcf',force=True)
                for b in ['pysam','cyvcf2']:self.assertEqual(backend(a,b,root/b).returncode,0)
                self.assertTrue(compare(root/'pysam',root/'cyvcf2')['passed'])

    def test_empty_candidates(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);a=fixture(root,empty=True)
            for b in ['pysam','cyvcf2']:self.assertEqual(backend(a,b,root/b).returncode,0)
            self.assertTrue(compare(root/'pysam',root/'cyvcf2')['passed'])

    def test_reject_multiallelic_and_missing_index(self):
        for case in ['multiallelic','index','duplicate','missing_gt']:
            with self.subTest(case=case),tempfile.TemporaryDirectory() as d:
                root=Path(d);a=fixture(root)
                if case=='index':Path(a.index).unlink()
                else:
                    p=root/'source.vcf';s=p.read_text()
                    if case=='multiallelic':s=s.replace('\tAC\t*\t','\tAC\t*,A\t')
                    if case=='duplicate':
                        lines=s.splitlines();i=next(i for i,x in enumerate(lines) if x.startswith('chr22\t20\t'));lines.insert(i,lines[i]);s='\n'.join(lines)+'\n'
                    if case=='missing_gt':s=s.replace('\tGT:AD\t','\tFT:AD\t')
                    p.write_text(s);pysam.tabix_compress(str(p),a.vcf,force=True);pysam.tabix_index(a.vcf,preset='vcf',force=True)
                for b in ['pysam','cyvcf2']:self.assertNotEqual(backend(a,b,root/b).returncode,0)

    def test_driver_reversed_order_and_full_record_comparison(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);a=fixture(root)
            import contextlib,io
            stdout=io.StringIO()
            with contextlib.redirect_stdout(stdout):run(a)
            self.assertNotIn('S1',stdout.getvalue())
            out=Path(a.outdir);r=json.loads((out/'benchmark.json').read_text())
            self.assertEqual(r['status'],'passed');self.assertEqual([x['backend'] for x in r['runs']],['pysam','cyvcf2','cyvcf2','pysam'])
            self.assertTrue(all(x['wall_seconds']>0 and x['cpu_seconds']>0 for x in r['runs']))
            # Equal aggregate row counts are insufficient: changing a genotype must fail.
            p=out/'02-cyvcf2/carriers.tsv.gz'
            with gzip.open(p,'rt') as f:s=f.read()
            with gzip.open(p,'wt') as f:f.write(s.replace('0/1','1/0',1))
            self.assertFalse(compare(out/'01-pysam',out/'02-cyvcf2')['passed'])
            with self.assertRaises(FileExistsError):run(a)


if __name__=='__main__':unittest.main()
