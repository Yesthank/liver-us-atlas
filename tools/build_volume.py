"""BodyParts3D OBJ -> label volume (1.5 mm) + Couinaud segmentation by portal territory.

Usage: python build_volume.py <obj_dir> <out_dir>
Outputs: volume.npz (labels uint8, origin, pitch), segment masks for meshing.
Coordinates are BodyParts3D mm: +X = patient left, +Y = posterior, +Z = superior.
"""
import sys, json, numpy as np, trimesh
from scipy import ndimage as ndi

OBJ = sys.argv[1] if len(sys.argv) > 1 else 'obj'
OUT = sys.argv[2] if len(sys.argv) > 2 else '.'
P = 1.5
LO = np.array([-178.0, -240.0, 880.0])
HI = np.array([178.0, 52.0, 1262.0])
SHAPE = np.ceil((HI - LO) / P).astype(int)
print('grid', SHAPE, SHAPE.prod() / 1e6, 'M')

L = dict(
    out=0, soft=1, fat=2, skin=3, muscle=4,
    s1=11, s2=12, s3=13, s4=14, s4b=19, s5=15, s6=16, s7=17, s8=18,
    mpv=20, rpv=21, lpv=22, splv=23, smv=24, rhv=25, mhv=26, lhv=27, ivc=28, renv=29,
    aorta=30, celiac=31, sma=32, spla=33, rena=34,
    gb=40, bile=41, pduct=42, pancreas=45, spleen=46, rkid=47, lkid=48, sinus=49, adrenal=50,
    stomach=55, duod=56, colon=57, diaph=60, psoas=61,
    ribs=70, cart=71, spine=72, lung=80, heart=81,
)


def load(ids):
    ms = [trimesh.load(f'{OBJ}/{i}.obj', process=True) for i in ids]
    return trimesh.util.concatenate(ms)


def surf_mask(mesh, step=0.5):
    """Mark grid voxels hit by mesh surface (dense point sampling on triangles)."""
    m = np.zeros(SHAPE, bool)
    area = mesh.area_faces
    n = np.maximum(1, (area / (step * step)).astype(int))
    n = np.minimum(n, 400)
    tri = mesh.triangles
    rep = np.repeat(np.arange(len(tri)), n)
    r1 = np.random.rand(len(rep)); r2 = np.random.rand(len(rep))
    s = np.sqrt(r1)
    a, b, c = tri[rep, 0], tri[rep, 1], tri[rep, 2]
    pts = (1 - s)[:, None] * a + (s * (1 - r2))[:, None] * b + (s * r2)[:, None] * c
    pts = np.vstack([pts, mesh.vertices])
    ijk = np.floor((pts - LO) / P).astype(int)
    ok = np.all((ijk >= 0) & (ijk < SHAPE), 1)
    ijk = ijk[ok]
    m[ijk[:, 0], ijk[:, 1], ijk[:, 2]] = True
    return m


def fill2d(m, ax):
    out = np.zeros_like(m)
    sl = [slice(None)] * 3
    for i in range(m.shape[ax]):
        sl[ax] = i
        out[tuple(sl)] = ndi.binary_fill_holes(m[tuple(sl)])
    return out


def solid(ids, close=1):
    """Surface voxels -> solid. 3D fill, plus the union of per-axis 2D fills (closed slice contours only) for leaky meshes."""
    m = surf_mask(load(ids))
    if close:
        m = ndi.binary_closing(m, iterations=close)
    bb = ndi.find_objects(m.astype(np.uint8))[0]
    bb = tuple(slice(max(0, b.start - 2), b.stop + 2) for b in bb)
    sub = m[bb]
    f3 = ndi.binary_fill_holes(sub)
    union = fill2d(sub, 0) | fill2d(sub, 1) | fill2d(sub, 2)
    out = np.zeros_like(m)
    out[bb] = f3 | union
    return out


g = json.load(open(f'{OBJ}/../groups.json'))
lab = np.zeros(SHAPE, np.uint8)

# ---- body (skin) : per axial slice 2D fill
skin = load(g['skin'])
keep = (skin.vertices[skin.faces].min(1)[:, 2] > LO[2] - 20) & (skin.vertices[skin.faces].max(1)[:, 2] < HI[2] + 20)
skin.update_faces(keep); skin.remove_unreferenced_vertices()
sm = ndi.binary_dilation(surf_mask(skin), iterations=1)
# body = span between first and last skin hit along Y (anterior->posterior) for each (x, z)
yi = np.arange(SHAPE[1])[None, :, None]
has = sm.any(1)
first = np.argmax(sm, axis=1)
last = SHAPE[1] - 1 - np.argmax(sm[:, ::-1, :], axis=1)
body = (yi >= first[:, None, :]) & (yi <= last[:, None, :]) & has[:, None, :] & ((last - first) > 40)[:, None, :]
# drop arms: keep largest connected component per slice
for k in range(SHAPE[2]):
    lb, n = ndi.label(body[:, :, k])
    if n > 1:
        sizes = ndi.sum(body[:, :, k], lb, range(1, n + 1))
        body[:, :, k] = lb == (1 + np.argmax(sizes))
print('body voxels', body.sum())
depth = ndi.distance_transform_edt(body) * P
lab[body] = L['soft']
lab[body & (depth < 16)] = L['fat']
lab[body & (depth >= 16) & (depth < 26)] = L['muscle']
lab[body & (depth < 2.2)] = L['skin']

# ---- diaphragm (sheet) & thorax contents above it
dia = ndi.binary_dilation(surf_mask(load(g['diaph'])), iterations=1)
zi = np.arange(SHAPE[2])[None, None, :]
dtop = np.where(dia.any(2), SHAPE[2] - 1 - np.argmax(dia[:, :, ::-1], axis=2), -1)
above = (zi > dtop[:, :, None]) & (dtop[:, :, None] >= 0) & body & (depth > 22)
xs = LO[0] + (np.arange(SHAPE[0]) + .5) * P
ys = LO[1] + (np.arange(SHAPE[1]) + .5) * P
X, Y = np.meshgrid(xs, ys, indexing='ij')
Zs = LO[2] + (np.arange(SHAPE[2]) + .5) * P
heart = (((X[:, :, None] - 22) / 62) ** 2 + ((Y[:, :, None] + 140) / 50) ** 2 + ((Zs[None, None, :] - 1228) / 60) ** 2) < 1
lab[above] = L['lung']
lab[above & heart] = L['heart']
lab[dia & body] = L['diaph']

order = ['psoas', 'colon', 'stomach', 'duod', 'rkid', 'lkid', 'adrenal', 'spleen', 'pancreas']
for k in order:
    m = solid(g[k], close=3 if k in ('stomach', 'colon', 'duod') else 1)
    lab[m] = L[k]
    if k in ('rkid', 'lkid'):
        d = ndi.distance_transform_edt(m) * P
        lab[m & (d > 9)] = L['sinus']

# ---- liver parenchyma + Couinaud by portal territory
segfiles = sum([g[f'seg{i}'] for i in range(1, 9)], [])
liver = solid(segfiles, close=2)
liver = ndi.binary_opening(liver, iterations=1)
print('liver vol cm3', liver.sum() * P ** 3 / 1000)
# Couinaud (Bismuth) by hepatic-vein scissurae around the IVC axis + portal plane.
# theta: angle about IVC centre in axial plane, 0 = anterior, negative = patient right.
ivcV = load(g['ivc']).vertices
def ivc_center(z):
    m = np.abs(ivcV[:, 2] - z) < 4
    return ivcV[m, :2].mean(0) if m.sum() > 10 else np.array([-15.0, -116.0])
def hv_profile(fid):
    V = load([fid]).vertices
    zz, tt = [], []
    for z in np.arange(1080, 1200, 5.0):
        m = np.abs(V[:, 2] - z) < 4
        if m.sum() < 5:
            continue
        c = ivc_center(z)
        tt.append(np.median(np.degrees(np.arctan2(V[m, 0] - c[0], -(V[m, 1] - c[1])))))
        zz.append(z)
    tt = ndi.median_filter(np.array(tt), 3)
    return np.array(zz), tt
zR, tR = hv_profile('FJ2416'); tR = np.clip(tR, -125, -95); zM, tM = hv_profile('FJ2414'); zLh, tLh = hv_profile('FJ2415')
print('RHV', tR.round(0), 'MHV', tM.round(0), 'LHV', tLh.round(0))
TH_UF = 30.0
Z_PORTAL_R, Z_PORTAL_IV, Z_II = 1110.0, 1124.0, 1118.0
zs = LO[2] + (np.arange(SHAPE[2]) + .5) * P
seg = np.zeros(SHAPE, np.uint8)
for k, z in enumerate(zs):
    sl = liver[:, :, k]
    if not sl.any():
        continue
    c = ivc_center(min(max(z, 1096), 1185))
    th = np.degrees(np.arctan2(X - c[0], -(Y - c[1])))
    tr, tm, tl = np.interp(z, zR, tR), np.interp(z, zM, tM), np.interp(z, zLh, tLh)
    s_ = np.zeros(sl.shape, np.uint8)
    post = (th < tr) | (th > 150)
    ant = (th >= tr) & (th < tm)
    med = (th >= tm) & (th < TH_UF)
    lat = (th >= TH_UF) & (th <= 150)
    s_[post] = 7 if z >= Z_PORTAL_R else 6
    s_[ant] = 8 if z >= Z_PORTAL_R else 5
    s_[med] = 4
    s_[lat] = 3
    s_[lat & (th > tl) & (z > Z_II)] = 2
    seg[:, :, k] = np.where(sl, s_, 0)
caud = solid(['FJ2816'])
seg[caud & liver] = 1
seg[~liver] = 0
for s in range(1, 9):
    lab[seg == s] = L[f's{s}']
lab[(seg == 4) & (zi < np.searchsorted(zs, Z_PORTAL_IV))] = L['s4b']
np.save(f'{OUT}/seg_tmp.npy', seg)

# ---- vessels / ducts (thin structures win)
for k in ['ivc', 'renv', 'smv', 'splv', 'mpv', 'rpv', 'lpv', 'rhv', 'mhv', 'lhv',
          'aorta', 'celiac', 'sma', 'spla', 'rena', 'gb', 'bile', 'pduct']:
    m = solid(g[k], close=1)
    lab[m] = L[k]

for k in ['cart', 'spine', 'ribs']:
    m = solid(g[k], close=1)
    lab[m] = L[k]

np.savez_compressed(f'{OUT}/volume.npz', lab=lab, lo=LO, pitch=P)
json.dump(L, open(f'{OUT}/labels.json', 'w'))
u, c = np.unique(lab, return_counts=True)
inv = {v: k for k, v in L.items()}
print({inv.get(int(a), a): int(b) for a, b in zip(u, c)})
