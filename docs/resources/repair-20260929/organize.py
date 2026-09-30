import csv,hashlib,json,os,time
from pathlib import Path
root=Path('/expanse/projects/sebat1/resources/rare-variant-pipeline')
work=root/'maintenance/repair-20260929';release=root/'releases/v1'
assert json.loads((work/'validation/VALIDATION_RECEIPT.json').read_text())['status']=='passed'
assert len(json.loads((release/'RESTORE_RECEIPT.json').read_text())['files'])==113
legacy=root/'legacy/registry-before-20260929';legacy.mkdir(parents=True,exist_ok=True)
ann=release/'targeted-annotation/ensembl-115'
old=root/'annotation/ensembl-115';new=root/'annotation/.ensembl-115-repaired'
assert not new.exists() and not (legacy/'annotation-ensembl-115').exists()
new.mkdir()
report={'status':'running','release':str(release),'legacy':str(legacy),'actions':[], 'legacy_byte_differences':{}}
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
def link(dest,source):
 assert source.exists(),str(source)
 dest.symlink_to(source)
# Preserve every existing filename, including mixed pilot files under chr22.
for d in old.iterdir():
 if not d.is_dir():raise RuntimeError('Unexpected annotation entry '+str(d))
 nd=new/d.name;nd.mkdir(exist_ok=True)
 for f in d.iterdir():
  if not f.is_file():raise RuntimeError('Unexpected annotation file '+str(f))
  if (ann/f.name).exists():
   before,after=sha(f),sha(ann/f.name)
   if before!=after:report['legacy_byte_differences'][str(f)]={'old_sha256':before,'v1_sha256':after,'old_file_preserved_under':str(legacy/'annotation-ensembl-115'/d.name/f.name)}
   link(nd/f.name,ann/f.name)
  else:
   (nd/f.name).symlink_to(legacy/'annotation-ensembl-115'/d.name/f.name)
for c in [str(i) for i in range(1,23)]+['X','Y']:
 d=new/('chr'+c);d.mkdir(exist_ok=True)
 for name in ['Homo_sapiens.GRCh38.115.chr'+c+'.gff3','Homo_sapiens.GRCh38.115.chr'+c+'.gff3.fastvep.cache','Homo_sapiens.GRCh38.dna.primary_assembly.fa','Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai','vep115.transcript-priority.tsv','vep115.consequence-ranks.tsv']:
  dest=d/name;src=ann/name;previous=old/d.name/name
  if dest.is_symlink():
   before,after=sha(previous),sha(src)
   if before!=after:report['legacy_byte_differences'][str(previous)]={'old_sha256':before,'v1_sha256':after,'old_file_preserved_under':str(legacy/'annotation-ensembl-115'/d.name/name)}
   dest.unlink()
  link(dest,src)
 d.chmod(0o2775)
old.rename(legacy/'annotation-ensembl-115')
try:new.rename(old)
except Exception:
 (legacy/'annotation-ensembl-115').rename(old);raise
report['actions'].append('All 24 annotation chromosome directories now use the same six consolidated resource aliases; all prior filenames retained.')
# Replace the legacy LOFTEE root alias with a complete compatibility directory.
oldlof=root/'loftee/GRCh38';newlof=root/'loftee/.GRCh38-repaired'
assert not newlof.exists() and not (legacy/'loftee-GRCh38').exists()
newlof.mkdir();(newlof/'ensembl-115').mkdir()
link(newlof/'loftee-grch38',release/'loftee-grch38')
link(newlof/'ensembl-115/transcripts.sqlite',release/'ensembl-115/transcripts.sqlite')
for d in [oldlof/'ensembl-115',root/'loftee/ensembl-115']:
 for f in d.glob('*.sqlite'):
  dest=newlof/'ensembl-115'/f.name
  if not dest.exists():link(dest,f.resolve())
oldlof.rename(legacy/'loftee-GRCh38')
try:newlof.rename(oldlof)
except Exception:
 (legacy/'loftee-GRCh38').rename(oldlof);raise
link(root/'loftee/ensembl-115/transcripts.sqlite',release/'ensembl-115/transcripts.sqlite')
report['actions'].append('Repaired active LOFTEE support paths; added consolidated transcript database, preserved old chromosome databases.')
# Preserve the broken GeneBayes alias as evidence, then replace it atomically.
gene=root/'genebayes/GeneBayes.Supplementary_Table_1.tsv'
assert gene.is_symlink()
(legacy/'GeneBayes.old-link').symlink_to(os.readlink(gene))
tmp=gene.with_name('.GeneBayes.repaired');link(tmp,release/'targeted-annotation/GeneBayes.Supplementary_Table_1.tsv');os.replace(str(tmp),str(gene))
link(root/'current',release)
report['actions'].append('Repaired GeneBayes alias and published current -> releases/v1.')
# Environment binding contains no cohort identity or scientific changes.
binding=root/'deployments/expanse-v1';binding.mkdir(parents=True,exist_ok=True)
paths=dict(resource_root=str(release),annotation_root=str(ann),loftee_root=str(release),genebayes=str(release/'targeted-annotation/GeneBayes.Supplementary_Table_1.tsv'),postprocess_config=str(release/'postprocess/config.json'),dbnsfp_scores_af=str(release/'dbNSFP/5.3.1a/parquet_scores_af'))
(binding/'resources.json').write_text(json.dumps(paths,indent=2)+'\n')
(binding/'resources.env').write_text(''.join(k.upper()+"='"+v+"'\n" for k,v in paths.items()))
# Verify active repaired namespaces resolve completely; archived broken links remain evidence.
for base in [root/'annotation/ensembl-115',root/'loftee/GRCh38',root/'genebayes']:
 for parent,dirs,files in os.walk(str(base),followlinks=True):
  for name in dirs+files:assert Path(parent,name).exists(),str(Path(parent,name))
assert len([p for p in (root/'annotation/ensembl-115').iterdir() if p.is_dir()])==24
assert not (root/'annotation/ensembl-115/chr22').is_symlink()
report['status']='passed';report['finished_unix']=time.time()
(work/'ORGANIZATION_RECEIPT.json').write_text(json.dumps(report,indent=2)+'\n')
print('ORGANIZATION COMPLETE: 24 consistent annotation directories; repaired GeneBayes/LOFTEE; old paths preserved.',flush=True)
