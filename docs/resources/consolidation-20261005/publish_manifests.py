import boto3,pathlib,json,hashlib,base64,csv
s=boto3.client('s3');bucket='sebat-genomics-work';prefix='resources/rare-variant-pipeline/v1/'
root=pathlib.Path('/Users/jamesguevara/work/rare-resource-expanded/docs/resources')
rows=json.load(open('/tmp/rare-consolidation/verified.json'));assert len(rows)==24
for name in ['RESOURCE_SHA256SUMS','DBNSFP_SHA256SUMS']:
 key=prefix+name;o=s.get_object(Bucket=bucket,Key=key);old=o['Body'].read();expected=pathlib.Path('/tmp/rare-consolidation/original-'+name).read_bytes();new=(root/('aws-'+name)).read_bytes()
 assert old in (expected,new)
 if old==expected:
  archive=prefix+'transfer-receipts/consolidation-20261005/original-'+name
  s.put_object(Bucket=bucket,Key=archive,Body=old,ChecksumSHA256=base64.b64encode(hashlib.sha256(old).digest()).decode(),IfNoneMatch='*')
  s.put_object(Bucket=bucket,Key=key,Body=new,IfMatch=o['ETag'],ChecksumSHA256=base64.b64encode(hashlib.sha256(new).digest()).decode())
 assert s.get_object(Bucket=bucket,Key=key)['Body'].read()==new
 print('Manifest verified',name)
receipt=pathlib.Path('/tmp/rare-consolidation/verified.json').read_bytes();s.put_object(Bucket=bucket,Key=prefix+'EXPANDED_CONSOLIDATION_20261005.json',Body=receipt,ChecksumSHA256=base64.b64encode(hashlib.sha256(receipt).digest()).decode(),IfNoneMatch='*')
catalog=root/'aws-file-catalog.tsv'
with catalog.open() as f:old=list(csv.DictReader(f,delimiter='\t'));fields=list(old[0])
newkeys={r['s3_uri'] for r in rows};assert not any(r['s3_uri'] in newkeys for r in old)
old.extend(dict(s3_uri=r['s3_uri'],bytes=r['bytes'],manifest_sha256=r['sha256'],last_modified=r['s3_last_modified']) for r in rows)
for r in old:
 if r['s3_uri'].split('/')[-1] in ('RESOURCE_SHA256SUMS','DBNSFP_SHA256SUMS'):
  h=s.head_object(Bucket=bucket,Key=r['s3_uri'].split(bucket+'/')[1]);r['bytes']=h['ContentLength'];r['last_modified']=h['LastModified'].isoformat()
with catalog.open('w',newline='') as f:w=csv.DictWriter(f,fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(old)
print('S3 manifests, original snapshots and receipt published')
