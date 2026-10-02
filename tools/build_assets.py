"""Pack meshes + label volume into a classic <script> data file (works from file://).

Usage: python build_assets.py <scratch_dir> <out_js>
Reads <scratch>/obj/*.obj, <scratch>/groups.json, <scratch>/volume.npz, <scratch>/labels.json
Positions are stored in BodyParts3D mm relative to CENTER, quantized x20 to int16.
"""
import sys, json, zlib, base64, numpy as np, trimesh, fast_simplification
from scipy import ndimage as ndi
from skimage import measure

SC = sys.argv[1]
OUT = sys.argv[2]
CENTER = np.array([-10.0, -120.0, 1100.0])
Q = 20.0

g = json.load(open(f'{SC}/groups.json'))
vol = np.load(f'{SC}/volume.npz')
lab, LO, P = vol['lab'], vol['lo'], float(vol['pitch'])
L = json.load(open(f'{SC}/labels.json'))


def b64(a):
    return base64.b64encode(zlib.compress(a.tobytes(), 9)).decode()


def simplify(m, target):
    if len(m.faces) <= target:
        return m
    v, f = fast_simplification.simplify(m.vertices.astype(np.float32), m.faces.astype(np.int32),
                                        target_reduction=1 - target / len(m.faces))
    return trimesh.Trimesh(v, f, process=True)


def load(ids):
    return trimesh.util.concatenate([trimesh.load(f'{SC}/obj/{i}.obj', process=True) for i in ids])


def pack(key, m):
    v = np.round((m.vertices - CENTER) * Q)
    assert np.abs(v).max() < 32767, key
    f = m.faces.astype(np.uint32 if len(m.vertices) > 65535 else np.uint16)
    return dict(key=key, nv=len(m.vertices), nf=len(m.faces), i32=bool(f.dtype == np.uint32),
                pos=b64(v.astype(np.int16)), idx=b64(f))


meshes = []
# --- Couinaud segments via marching cubes on label volume
seg_codes = {1: [L['s1']], 2: [L['s2']], 3: [L['s3']], 4: [L['s4'], L['s4b']], 5: [L['s5']],
             6: [L['s6']], 7: [L['s7']], 8: [L['s8']]}
liver_codes = sum(seg_codes.values(), [])
vessel_codes = [L[k] for k in ('mpv', 'rpv', 'lpv', 'rhv', 'mhv', 'lhv', 'ivc', 'gb', 'bile')]
for s, codes in seg_codes.items():
    m = np.isin(lab, codes)
    # let vessels inside the liver count as parenchyma so segment surfaces stay closed
    near = ndi.binary_dilation(m, iterations=2) & np.isin(lab, vessel_codes)
    m = m | near
    sl = ndi.find_objects(m.astype(np.uint8))[0]
    sl = tuple(slice(max(0, a.start - 3), a.stop + 3) for a in sl)
    f = ndi.gaussian_filter(m[sl].astype(np.float32), 1.0)
    v, fc, _, _ = measure.marching_cubes(f, 0.5)
    v = (v + [a.start for a in sl] + 0.5) * P + LO
    tm = trimesh.Trimesh(v, fc[:, ::-1], process=True)
    trimesh.smoothing.filter_taubin(tm, iterations=12)
    tm = simplify(tm, 5000)
    meshes.append(pack(f'seg{s}', tm))
    print('seg', s, len(tm.faces))

targets = dict(mpv=3000, rpv=5000, lpv=5000, splv=1200, smv=800, rhv=4000, mhv=4000, lhv=2500, ivc=2500,
               renv=2000, aorta=3000, celiac=800, sma=1800, spla=2500, rena=500, gb=2400, bile=1800,
               pduct=1500, pancreas=4000, spleen=1400, rkid=2400, lkid=2400, adrenal=1500, stomach=1900,
               duod=1900, colon=4000, diaph=9000, ribs=24000, cart=5000, spine=8000, psoas=2500)
for k, t in targets.items():
    m = simplify(load(g[k]), t)
    meshes.append(pack(k, m))
    print(k, len(m.faces))

# skin: torso crop
skin = load(g['skin'])
fv = skin.vertices[skin.faces]
keep = (fv[:, :, 2].min(1) > 905) & (fv[:, :, 2].max(1) < 1290) & (np.abs(fv[:, :, 0]).max(1) < 190)
skin.update_faces(keep); skin.remove_unreferenced_vertices()
skin = simplify(skin, 14000)
meshes.append(pack('skin', skin))
print('skin', len(skin.faces))

# --- volume (x fastest in JS: store as [z][y][x])
volT = np.ascontiguousarray(np.transpose(lab, (2, 1, 0)))
data = dict(center=CENTER.tolist(), q=Q, meshes=meshes,
            volume=dict(dims=list(lab.shape), lo=LO.tolist(), pitch=P, data=b64(volT)),
            labels=L)
js = '/* BodyParts3D, (c) The Database Center for Life Science, CC BY-SA 2.1 JP. Derived data. */\n'
js += 'window.ANATOMY=' + json.dumps(data, separators=(',', ':')) + ';\n'
open(OUT, 'w').write(js)
print('wrote', OUT, round(len(js) / 1e6, 2), 'MB')
