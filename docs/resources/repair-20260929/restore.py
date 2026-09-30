import concurrent.futures,csv,hashlib,json,os,time,urllib.request,shutil
from pathlib import Path
base=Path('/expanse/projects/sebat1/resources/rare-variant-pipeline')
work=base/'maintenance/repair-20260929'
stage=base/'releases/.v1-restoring-20260929'
stage.mkdir(exist_ok=True)
rows=list(csv.DictReader((work/'runtime-manifest.tsv').open(),delimiter='\t'))
urls=json.loads((work/'private-download-urls.json').read_text())
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
def restore(row):
 name=row['relative_path'];dest=stage/name;dest.parent.mkdir(parents=True,exist_ok=True)
 expected=row['sha256'];size=int(row['size_bytes'])
 if dest.exists() and dest.stat().st_size==size and sha(dest)==expected:
  print('VERIFIED cached '+name,flush=True);return dict(path=name,bytes=size,sha256=expected)
 for attempt in range(3):
  temp=dest.with_name(dest.name+'.partial')
  try:
   h=hashlib.sha256();n=0
   with urllib.request.urlopen(urls[name],timeout=120) as src,temp.open('wb') as dst:
    while True:
     b=src.read(8*1024*1024)
     if not b:break
     dst.write(b);h.update(b);n+=len(b)
   if n!=size or h.hexdigest()!=expected:raise ValueError('Downloaded bytes/checksum mismatch')
   if sha(temp)!=expected:raise ValueError('On-disk checksum mismatch')
   temp.chmod(0o664);temp.replace(dest)
   print('VERIFIED '+name,flush=True)
   return dict(path=name,bytes=size,sha256=expected)
  except Exception as e:
   print('RETRY '+name+' '+type(e).__name__,flush=True)

   if temp.exists():temp.unlink()
   if attempt==2:raise RuntimeError('Resource download failed: '+name) from None
   time.sleep(5)
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:verified=list(pool.map(restore,rows))
assert len(verified)==113
for p in stage.rglob('*'):
 if p.is_dir():p.chmod(0o2775)
for name in ['runtime-manifest.tsv','RESOURCE_SHA256SUMS','DBNSFP_SHA256SUMS']:
 shutil.copyfile(work/name,stage/name)
(stage/'RESTORE_RECEIPT.json').write_text(json.dumps(dict(status='verified',source='s3://sebat-genomics-work/resources/rare-variant-pipeline/v1/',verified_unix=time.time(),method='Every downloaded file SHA-256 checked during transfer and independently re-read from disk',files=verified),indent=2)+'\n')
if (base/'releases/v1').exists():raise RuntimeError('Refusing to overwrite existing published release')
stage.rename(base/'releases/v1')
(work/'private-download-urls.json').unlink()
print('RESTORE COMPLETE: 113/113; '+str(base/'releases/v1'),flush=True)
