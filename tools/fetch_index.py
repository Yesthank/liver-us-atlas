import urllib.request, zipfile, io, sys, os, json, concurrent.futures as cf, zlib, struct, time
URL='https://dbarchive.biosciencedbc.jp/data/bodyparts3d/LATEST/isa_BP3D_4.0_obj_99.zip'
SIZE=142903898
def rng(a,b):
    for k in range(5):
        try:
            r=urllib.request.Request(URL,headers={'Range':f'bytes={a}-{b}'})
            return urllib.request.urlopen(r,timeout=120).read()
        except Exception as e: print('retry',a,e); time.sleep(2)
    raise
# EOCD
tail=rng(SIZE-70000,SIZE-1)
i=tail.rfind(b'PK\x05\x06')
cd_size,cd_off=struct.unpack('<II',tail[i+12:i+20])
cd=rng(cd_off,cd_off+cd_size-1) if cd_off < SIZE-70000 else tail[cd_off-(SIZE-70000):]
entries={}
p=0
while p<len(cd) and cd[p:p+4]==b'PK\x01\x02':
    method,=struct.unpack('<H',cd[p+10:p+12]); csz,usz=struct.unpack('<II',cd[p+20:p+28])
    nl,el,cl=struct.unpack('<HHH',cd[p+28:p+34]); off,=struct.unpack('<I',cd[p+42:p+46])
    name=cd[p+46:p+46+nl].decode('latin1'); entries[os.path.basename(name)]=(name,method,csz,usz,off)
    p+=46+nl+el+cl
json.dump(entries,open('zip_index.json','w'))
print(len(entries),'entries'); 
