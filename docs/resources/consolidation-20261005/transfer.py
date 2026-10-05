import boto3,subprocess,json,base64,time,pathlib,concurrent.futures
from botocore.config import Config
s3=boto3.client('s3',region_name='us-east-1',config=Config(signature_version='s3v4'));bucket='sebat-genomics-work';prefix='resources/rare-variant-pipeline/v1/'
root='/expanse/projects/sebat1/resources/rare-variant-pipeline/releases/v1/'
def transfer(row):
 key=prefix+row['path'];checksum=base64.b64encode(bytes.fromhex(row['sha256'])).decode()
 try:
  h=s3.head_object(Bucket=bucket,Key=key,ChecksumMode='ENABLED')
 except s3.exceptions.ClientError as e:
  if e.response['Error']['Code'] not in ('404','NoSuchKey'):raise
  url=s3.generate_presigned_url('put_object',Params=dict(Bucket=bucket,Key=key,ChecksumSHA256=checksum,IfNoneMatch='*'),ExpiresIn=21600)
  data=dict(url=url,file=root+row['path'],headers={'Content-Length':str(row['bytes']),'x-amz-checksum-sha256':checksum,'If-None-Match':'*'})
  p=subprocess.run(['ssh','expanse','python3 /tmp/rare-expanded-upload.py'],input=json.dumps(data),text=True,capture_output=True)
  if p.returncode:raise RuntimeError('Upload failed for '+row['path']+' '+p.stderr[-500:])
  h=s3.head_object(Bucket=bucket,Key=key,ChecksumMode='ENABLED')
 if h['ContentLength']!=row['bytes'] or h.get('ChecksumSHA256')!=checksum:raise RuntimeError('S3 checksum mismatch '+key)
 return dict(**row,s3_uri='s3://'+bucket+'/'+key,s3_checksum_sha256=h['ChecksumSHA256'],s3_last_modified=h['LastModified'].isoformat(),verification='source and Expanse SHA256; S3 checksum-validated PUT and HEAD')
seen=set();done=[];pending={}
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
 while len(done)<24:
  raw=subprocess.check_output(['ssh','expanse','cat /tmp/rare-expanded-inventory.json'],text=True)
  rows=json.loads(raw)
  for row in rows:
   if row['path'] not in seen:pending[pool.submit(transfer,row)]=row;seen.add(row['path'])
  for future in list(pending):
   if future.done():
    result=future.result();done.append(result);del pending[future]
    pathlib.Path('/tmp/rare-consolidation/verified.json').write_text(json.dumps(done,indent=2)+'\n')
    print('VERIFIED',len(done),result['path'],flush=True)
  if len(done)<24:time.sleep(10)
print('COMPLETE',sum(x['bytes'] for x in done),flush=True)
