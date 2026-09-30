import concurrent.futures,csv,hashlib,json,time
from pathlib import Path
root=Path('/expanse/projects/sebat1/resources/rare-variant-pipeline');release=root/'releases/v1';work=root/'maintenance/repair-20260929'
rows=list(csv.DictReader((release/'runtime-manifest.tsv').open(),delimiter='\t'))
def verify(row):
 p=release/row['relative_path'];h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 assert p.stat().st_size==int(row['size_bytes']) and h.hexdigest()==row['sha256'],str(p)
 assert not p.stat().st_mode & 0o222
 print('FINAL HASH PASS '+row['relative_path'],flush=True)
 return dict(path=row['relative_path'],bytes=p.stat().st_size,sha256=h.hexdigest())
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:verified=list(pool.map(verify,rows))
assert len(verified)==113
logs=list((work/'validation-readonly').glob('*.private.log'));assert len(logs)==24
for p in logs:
 text=p.read_text()
 assert 'transcripts from cache' in text,str(p)
 assert 'Cache is stale' not in text and 'Saved transcript cache' not in text and 'cache load failed' not in text,str(p)
for tag in [str(i) for i in range(1,23)]+['X','Y']:
 name='chr'+tag+'.picked.tsv'
 assert (work/'validation'/name).read_bytes()==(work/'validation-readonly'/name).read_bytes(),'Synthetic output drift '+name
receipt=dict(initial_vs_pinned_synthetic_outputs_identical_for_24_chromosomes=True,status='passed',finished_unix=time.time(),files=verified,all_24_fastvep_runs_loaded_pinned_caches=True,release_read_only=True,method='Independent SHA-256 of all 113 files AFTER read-only container validation')
(work/'POST_VALIDATION_HASHES.json').write_text(json.dumps(receipt,indent=2)+'\n')
print('FINAL VERIFICATION COMPLETE: 113/113 hashes unchanged after validation; all 24 pinned caches loaded directly.')
