import csv,hashlib,json,os,shutil,time
from pathlib import Path
root=Path('/expanse/projects/sebat1/resources/rare-variant-pipeline');release=root/'releases/v1';work=root/'maintenance/repair-20260929'
rows=list(csv.DictReader((release/'runtime-manifest.tsv').open(),delimiter='\t'))
expected={r['relative_path']:r['sha256'] for r in rows}
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
changes=[]
ann=release/'targeted-annotation/ensembl-115'
for c in [str(i) for i in range(1,23)]+['X','Y']:
 gff=ann/('Homo_sapiens.GRCh38.115.chr'+c+'.gff3');cache=Path(str(gff)+'.fastvep.cache');want=expected[str(cache.relative_to(release))]
 before=sha(cache)
 if before!=want:
  source=root/'legacy/registry-before-20260929/annotation-ensembl-115'/('chr'+c)/cache.name
  assert sha(source)==want
  tmp=cache.with_name(cache.name+'.restoring');shutil.copyfile(source,tmp);assert sha(tmp)==want;tmp.replace(cache)
 ns=gff.stat().st_mtime_ns+1000000000
 os.utime(str(cache),ns=(ns,ns))
 assert cache.stat().st_mtime_ns>gff.stat().st_mtime_ns
 changes.append(dict(chromosome='chr'+c,expected_sha256=want,before_sha256=before,restored_bytes=before!=want,cache_mtime_ns=cache.stat().st_mtime_ns,gff_mtime_ns=gff.stat().st_mtime_ns))
for p in release.rglob('*'):
 if p.is_file():p.chmod(0o444)
 elif p.is_dir():p.chmod(0o2555)
release.chmod(0o2555)
(work/'CACHE_PINNING_RECEIPT.json').write_text(json.dumps(dict(status='passed',reason='Initial validation rebuilt 3 caches because download-time mtime ordering made them stale; restored original manifest bytes, enforced cache>GFF mtime, made release read-only.',files=changes),indent=2)+'\n')
print('Restored '+str(sum(r['restored_bytes'] for r in changes))+' caches; normalized 24 cache mtimes; release made read-only.')
