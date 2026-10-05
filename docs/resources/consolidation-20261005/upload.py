import json,http.client,urllib.parse,pathlib,sys
item=json.load(sys.stdin);u=urllib.parse.urlsplit(item['url'])
c=http.client.HTTPSConnection(u.hostname,timeout=1800)
c.putrequest('PUT',u.path+'?'+u.query)
for k,v in item['headers'].items():c.putheader(k,v)
c.endheaders()
try:
 with open(item['file'],'rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):c.send(b)
except (BrokenPipeError,ConnectionResetError):
 r=c.getresponse();print('HTTP',r.status,r.read().decode()[:500]);sys.exit(1)
r=c.getresponse();r.read()
if r.status not in (200,201):raise RuntimeError('S3 upload HTTP '+str(r.status))
print('uploaded',pathlib.Path(item['file']).name)
