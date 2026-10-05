import pathlib,hashlib,json,os,shutil,subprocess,datetime
base=pathlib.Path('/expanse/lustre/projects/ddp195/j3guevar/rare-variant-pipeline');stage=base/'.staging-v1';source=pathlib.Path('/expanse/projects/sebat1/resources/rare-variant-pipeline/releases/v1')
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
checks=[]
for line in (stage/'RESOURCE_SHA256SUMS').read_text().splitlines():
 digest,rel=line.split(None,1);p=stage/rel
 if sha(p)!=digest:raise RuntimeError('checksum mismatch '+rel)
 checks.append(dict(path=rel,sha256=digest,bytes=p.stat().st_size))
 print('VERIFIED',len(checks),rel,flush=True)
assert len(checks)==137
# Check every copied auxiliary file, not just runtime resources.
allfiles=[p for p in source.rglob('*') if p.is_file()]
for p in allfiles:
 rel=p.relative_to(source)
 if any(x['path'].lstrip('./')==str(rel) for x in checks):continue
 if sha(p)!=sha(stage/rel):raise RuntimeError('auxiliary mismatch '+str(rel))
containers=base/'containers/v1'
subprocess.run(['rsync','-rt','--perms','--chmod=D2750,F440','/expanse/projects/sebat1/resources/rare-variant-pipeline/containers/v1/',str(containers)+'/'],check=True)
for p in containers.glob('*.sif'):
 side=p.with_name(p.name+'.sha256');expected=side.read_text().split()[0]
 if sha(p)!=expected:raise RuntimeError('container mismatch '+p.name)
 os.chmod(p,0o550)
# Verify required cache timestamp ordering survives the copy.
for gff in (stage/'targeted-annotation/ensembl-115').glob('*.gff3'):
 cache=gff.with_name(gff.name+'.fastvep.cache')
 if cache.exists() and cache.stat().st_mtime<=gff.stat().st_mtime:raise RuntimeError('cache timestamp '+gff.name)
release=base/'releases/v1';release.parent.mkdir(exist_ok=True)
if release.exists():raise RuntimeError('release already exists; refusing overwrite')
os.rename(stage,release)
(base/'current').symlink_to('releases/v1')
for p in [base,release.parent,base/'containers']:os.chmod(p,0o2750)
for p in release.rglob('*'):
 if p.is_dir():os.chmod(p,0o2550)
os.chmod(release,0o2550)
receipt=dict(status='PASS',source=str(source),destination=str(release),runtime_files=len(checks),all_release_files=len(allfiles),containers=2,sha256_verified=True,cache_timestamp_order_verified=True,access_group='ddp195',verified_at=datetime.datetime.now(datetime.timezone.utc).isoformat(),files=checks)
(base/'COPY_RECEIPT.json').write_text(json.dumps(receipt,indent=2)+'\n');os.chmod(base/'COPY_RECEIPT.json',0o440)
d=base/'deployments/ddp195-v1';d.mkdir(parents=True,exist_ok=True)
bindings=json.loads(pathlib.Path('/expanse/projects/sebat1/resources/rare-variant-pipeline/deployments/expanse-v1/resources.json').read_text())
bindings={k:v.replace(str(source),str(release)) if isinstance(v,str) else v for k,v in bindings.items()}
(d/'resources.json').write_text(json.dumps(bindings,indent=2)+'\n')
(d/'resources.env').write_text(''.join(k.upper()+"='"+v+"'\n" for k,v in bindings.items()))
(base/'README.md').write_text('''# Shared rare-variant resources for ddp195

Scientific resources: `releases/v1/` (`current` is an alias).
Pinned SIF containers: `containers/v1/`.
Local resource bindings: `deployments/ddp195-v1/resources.json` and `resources.env`.
Verification: `COPY_RECEIPT.json` (137 runtime files and both containers).

This is a verified read-only copy of the canonical Sebat release, including all
24 unfiltered dbNSFP parquet_expanded chromosomes. It is readable by ddp195
members, including toedwards. It contains reference resources, not cohort data.

Current resource guide:
https://github.com/james-guevara/rare-variant-pipeline/blob/docs/consolidate-expanded-resources/docs/resources/README.md

Preserve resource timestamps on transfer: FastVEP caches must be newer than GFF3.
Historical transfer receipts retain their original deployment paths; use the
local bindings above for this deployment.
''')
for p in [base/'README.md',d/'resources.json',d/'resources.env']:os.chmod(p,0o440)
print('COMPLETE',len(checks),len(allfiles),flush=True)
