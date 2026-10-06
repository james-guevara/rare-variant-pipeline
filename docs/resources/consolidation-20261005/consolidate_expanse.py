import pathlib,hashlib,json,shutil,os
src=pathlib.Path('/expanse/projects/sebat1/s3/data/sebat/resources/dbNSFP/5.3.1a/parquet_expanded')
root=pathlib.Path('/expanse/projects/sebat1/resources/rare-variant-pipeline/releases/v1')
parent=root/'dbNSFP/5.3.1a'
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
mode=parent.stat().st_mode & 0o7777
os.chmod(parent,mode|0o200)
dest=parent/'parquet_expanded'
dest.mkdir(exist_ok=True)
os.chmod(parent,mode)
rows=[]
for chrom in [str(i) for i in range(1,23)]+['X','Y']:
 s=src/f'chr{chrom}.parquet'; d=dest/s.name
 digest=sha(s)
 if not d.exists():
  t=d.with_suffix('.copying'); shutil.copy2(s,t)
  if sha(t)!=digest:raise RuntimeError('copy mismatch '+s.name)
  os.chmod(t,0o444);os.replace(t,d)
 elif sha(d)!=digest:raise RuntimeError('existing destination mismatch '+s.name)
 rows.append(dict(path=str(d.relative_to(root)),bytes=d.stat().st_size,sha256=digest,chromosome=chrom,source=str(s)))
 pathlib.Path('/tmp/rare-expanded-inventory.json').write_text(json.dumps(rows,indent=2)+'\n')
 print(s.name,d.stat().st_size,digest,flush=True)
os.chmod(dest,0o2555)
print('COMPLETE',len(rows),flush=True)
