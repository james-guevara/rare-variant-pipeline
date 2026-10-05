import pathlib,json,os,shutil,hashlib
root=pathlib.Path('/expanse/projects/sebat1/resources/rare-variant-pipeline')
release=root/'releases/v1';stage=pathlib.Path('/tmp/rare-consolidation-install');audit=root/'maintenance/consolidation-20261005'
audit.mkdir(parents=True,exist_ok=True)
def replace(p,data):
 mode=p.stat().st_mode&0o7777 if p.exists() else 0o444
 pmode=p.parent.stat().st_mode&0o7777
 os.chmod(p.parent,pmode|0o200)
 t=p.with_name(p.name+'.consolidating');t.write_bytes(data);os.chmod(t,mode);os.replace(t,p);os.chmod(p.parent,pmode)
for name in ['RESOURCE_SHA256SUMS','DBNSFP_SHA256SUMS']:
 current=release/name;old=audit/('original-'+name)
 if not old.exists():shutil.copy2(current,old)
 expected=(stage/('original-'+name)).read_bytes()
 if current.read_bytes() not in [expected,(stage/name).read_bytes()]:raise RuntimeError('unexpected current manifest '+name)
 replace(current,(stage/name).read_bytes())
for name in ['verified.json','inventory.json','consolidate_expanse.py','upload.py','transfer.py','install_registry.py']:
 if (stage/name).exists():shutil.copy2(stage/name,audit/name)
replace(release/'EXPANDED_CONSOLIDATION_20261005.json',(stage/'verified.json').read_bytes())
p=root/'deployments/expanse-v1/resources.json';d=json.loads(p.read_text());d['dbnsfp_expanded']=str(release/'dbNSFP/5.3.1a/parquet_expanded');replace(p,(json.dumps(d,indent=2)+'\n').encode())
p=root/'deployments/expanse-v1/resources.env';s=p.read_text()
if 'DBNSFP_EXPANDED=' not in s:s+="DBNSFP_EXPANDED='"+str(release/'dbNSFP/5.3.1a/parquet_expanded')+"'\n"
replace(p,s.encode())
p=root/'README.md';s=p.read_text();marker='## Unfiltered expanded dbNSFP consolidation (2026-10-05)'
if marker not in s:
 s+='\n'+marker+'\n\nAll 24 chr1–22/X/Y files are now installed as read-only regular files at\n`releases/v1/dbNSFP/5.3.1a/parquet_expanded/` and uploaded under the same\nrelative path in the S3 v1 release. Source/destination SHA-256 and S3\nchecksum validation are recorded in\n`releases/v1/EXPANDED_CONSOLIDATION_20261005.json`.\nThe runtime inventory now contains 137 files: the original 113 verified in\nSeptember plus these 24 verified additions. Existing MANE-filtered files\nremain for reproducibility. No scientific configuration was changed.\n\nThe shared resource bindings now expose `dbnsfp_expanded` / `DBNSFP_EXPANDED`.\nCurrent GitHub guide: https://github.com/james-guevara/rare-variant-pipeline/blob/docs/consolidate-expanded-resources/docs/resources/README.md\n'
replace(p,s.encode())
print('INSTALLED manifests, receipt, bindings, and registry README')
