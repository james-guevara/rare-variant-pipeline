import csv,json,sqlite3,subprocess,sys,time
from pathlib import Path
import pysam,pyBigWig,pyarrow.parquet as pq
root=Path(sys.argv[1]);work=Path(sys.argv[2]);work.mkdir(exist_ok=True)
ann=root/'targeted-annotation/ensembl-115'
chroms=[str(i) for i in range(1,23)]+['X','Y']
report={'status':'running','resource_root':str(root),'chromosomes':{},'parquets':{},'synthetic_only':True}
fasta=pysam.FastaFile(str(ann/'Homo_sapiens.GRCh38.dna.primary_assembly.fa'))
assert set(chroms)<=set(fasta.references)
ancestor=pysam.FastaFile(str(root/'loftee-grch38/human_ancestor.fa.gz'))
assert set(chroms)<=set(ancestor.references)
gerp=pyBigWig.open(str(root/'loftee-grch38/gerp_conservation_scores.homo_sapiens.GRCh38.bw'))
assert set(chroms)<=set(gerp.chroms())
for name in ['ensembl-115/transcripts.sqlite','loftee-grch38/loftee.sql']:
 with sqlite3.connect('file:'+str(root/name)+'?mode=ro',uri=True) as db:
  assert db.execute('PRAGMA integrity_check').fetchall()==[('ok',)]
  if name.startswith('ensembl'):
   counts=dict(db.execute('select seqname,count(*) from transcript group by seqname'))
   assert set(chroms)<=set(counts)
   report['transcript_counts']=counts
  report[name]='integrity_check passed'
with open(ann/'vep115.transcript-priority.tsv') as f:
 report['priority_data_rows']=sum(1 for _ in f)-1
assert report['priority_data_rows']>500000
for kind in ['parquet_scores_af','parquet_expanded_mane_select']:
 for chrom in chroms:
  p=root/'dbNSFP/5.3.1a'/kind/('chr'+chrom+'.parquet')
  meta=pq.read_metadata(p)
  assert meta.num_rows>0
  report['parquets'][kind+'/chr'+chrom]={'rows':meta.num_rows,'row_groups':meta.num_row_groups}
for name in ['genomicSuperDups','rmsk','simpleRepeat']:
 seen=set()
 with open(root/'problematic-regions'/(name+'.bed')) as f:
  for line in f:
   if line.startswith(('#','track','browser')):continue
   fields=line.split('\t');assert len(fields)>=3
   assert int(fields[2])>=int(fields[1])>=0
   seen.add(fields[0].removeprefix('chr'))
 assert set(chroms)<=seen,(name,'missing chromosomes')
 report['regions_'+name]=sorted(seen)
for chrom in chroms:
 start=time.perf_counter();tag='chr'+chrom
 vcf=work/(tag+'.synthetic.vcf');picked=work/(tag+'.picked.tsv')
 positions=[];seen=set()
 with open(ann/('Homo_sapiens.GRCh38.115.'+tag+'.gff3')) as f:
  for line in f:
   if line.startswith('#'):continue
   c=line.rstrip().split('\t')
   if len(c)<9 or c[0]!=chrom or c[2]!='CDS' or int(c[4])-int(c[3])<30:continue
   pos=int(c[3])+10
   if pos in seen:continue
   seq=fasta.fetch(chrom,pos-1,pos+1).upper()
   if len(seq)!=2 or any(x not in 'ACGT' for x in seq):continue
   positions.append((pos,seq));seen.add(pos)
   if len(positions)==20:break
 assert positions
 with vcf.open('w') as f:
  f.write('##fileformat=VCFv4.2\n##contig=<ID='+chrom+'>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
  for i,(pos,seq) in enumerate(sorted(positions)):
   f.write(f'{chrom}\t{pos}\tSYNTHETIC_{tag}_{i}\t{seq}\t{seq[0]}\t.\tPASS\t.\n')
 with (work/(tag+'.private.log')).open('w') as log:
  a=subprocess.Popen(['fastvep','annotate','--input',str(vcf),'--gff3',str(ann/('Homo_sapiens.GRCh38.115.'+tag+'.gff3')),'--fasta',str(ann/'Homo_sapiens.GRCh38.dna.primary_assembly.fa'),'--transcript-cache',str(ann/('Homo_sapiens.GRCh38.115.'+tag+'.gff3.fastvep.cache')),'--hgvs','--symbol','--canonical','--output-format','vcf','--output','-'],stdout=subprocess.PIPE,stderr=log)
  b=subprocess.Popen(['fastvep-picker','--fastvep','-','--transcript-priority',str(ann/'vep115.transcript-priority.tsv'),'--consequence-ranks',str(ann/'vep115.consequence-ranks.tsv'),'--output',str(picked)],stdin=a.stdout,stdout=log,stderr=log)
  a.stdout.close();br=b.wait();ar=a.wait();assert ar==br==0
  rows=list(csv.DictReader(picked.open(),delimiter='\t'));assert len(rows)==len(positions)
  assert any('frameshift_variant' in row['Consequence'] for row in rows)
  lof=work/(tag+'.loftee.tsv')
  subprocess.run([sys.executable,'/opt/rvp/scripts/run_standalone_loftee.py','--input',str(picked),'--transcripts',str(root/'ensembl-115/transcripts.sqlite'),'--reference',str(ann/'Homo_sapiens.GRCh38.dna.primary_assembly.fa'),'--ancestor',str(root/'loftee-grch38/human_ancestor.fa.gz'),'--gerp',str(root/'loftee-grch38/gerp_conservation_scores.homo_sapiens.GRCh38.bw'),'--conservation',str(root/'loftee-grch38/loftee.sql'),'--output',str(lof)],check=True,stdout=log,stderr=log)
  lofrows=list(csv.DictReader(lof.open(),delimiter='\t'));assert len(lofrows)>0
 report['chromosomes'][tag]={'synthetic_sites':len(positions),'picked_rows':len(rows),'loftee_rows':len(lofrows),'seconds':time.perf_counter()-start}
 print('PASS '+tag+' synthetic='+str(len(positions))+' picked='+str(len(rows))+' loftee='+str(len(lofrows)),flush=True)
 (work/'validation-progress.json').write_text(json.dumps(report,indent=2)+'\n')
fasta.close();ancestor.close();gerp.close()
report['status']='passed';report['finished_unix']=time.time()
(work/'VALIDATION_RECEIPT.json').write_text(json.dumps(report,indent=2)+'\n')
print('VALIDATION COMPLETE: 24 chromosomes; FastVEP, Rust picker, LOFTEE, 48 Parquets, reference/index, SQLite integrity, region BEDs.',flush=True)
