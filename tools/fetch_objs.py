import json,urllib.request,zlib,struct,os,time,concurrent.futures as cf
URL='https://dbarchive.biosciencedbc.jp/data/bodyparts3d/LATEST/isa_BP3D_4.0_obj_99.zip'
os.makedirs('obj', exist_ok=True)
idx=json.load(open('zip_index.json')); g=json.load(open('groups.json'))
ids=sorted({f for v in g.values() for f in v})
def rng(a,b):
    for k in range(6):
        try: return urllib.request.urlopen(urllib.request.Request(URL,headers={'Range':f'bytes={a}-{b}'}),timeout=300).read()
        except Exception as e: print('retry',a,e,flush=True); time.sleep(3)
    raise RuntimeError
def get(fid):
    out=f'obj/{fid}.obj'
    if os.path.exists(out): return fid,0
    name,method,csz,usz,off=idx[fid+'.obj']
    h=rng(off,off+29+400+csz)
    nl,el=struct.unpack('<HH',h[26:30]); d=h[30+nl+el:30+nl+el+csz]
    data=zlib.decompress(d,-15) if method==8 else d
    assert len(data)==usz
    open(out,'wb').write(data); return fid,csz
t=time.time(); tot=0
with cf.ThreadPoolExecutor(12) as ex:
    for fid,c in ex.map(get,ids): tot+=c
print(len(ids),'files',tot/1e6,'MB compressed',round(time.time()-t),'s')
