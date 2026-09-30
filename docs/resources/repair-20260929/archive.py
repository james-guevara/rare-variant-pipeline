import hashlib,json,os,urllib.request,time
from pathlib import Path
root=Path('/expanse/projects/sebat1/resources/rare-variant-pipeline')
work=root/'maintenance/repair-20260929';release=root/'releases/v1'
rows=json.loads((work/'private-archive-downloads.json').read_text());receipt=[]
for row in rows:
 p=release/row['path'];p.parent.mkdir(parents=True,exist_ok=True)
 tmp=p.with_name(p.name+'.partial');h=hashlib.sha256()
 with urllib.request.urlopen(row['url'],timeout=120) as src,tmp.open('wb') as dst:
  while True:
   b=src.read(8*1024*1024)
   if not b:break
   dst.write(b);h.update(b)
 assert tmp.stat().st_size==row['bytes'] and h.hexdigest()==row['sha256']
 h=hashlib.sha256()
 with tmp.open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 assert h.hexdigest()==row['sha256']
 tmp.chmod(0o664);tmp.replace(p)
 receipt.append({k:row[k] for k in ['path','bytes','sha256']})
 print('ARCHIVE VERIFIED '+row['path'],flush=True)
(work/'private-archive-downloads.json').unlink()
(work/'ARCHIVE_RECEIPT.json').write_text(json.dumps(dict(status='passed',files=receipt,finished_unix=time.time()),indent=2)+'\n')
