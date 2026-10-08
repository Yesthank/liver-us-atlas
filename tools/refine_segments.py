"""Refine the S4 boundaries (Cantlie line and umbilical fissure) in an existing data/anatomy.js.

Usage: python refine_segments.py data/anatomy.js
Run after build_assets.py. Safe to run again (the result does not change).

build_volume.py splits S4 from S5/S8 with the MHV trunk mesh angle, held constant below the
mesh. That left the gallbladder fossa entirely in S5. Here the boundary angle (about the IVC
axis) runs from the gallbladder fossa at the GB top to the MHV trunk just below the
confluence (Cantlie line), and is linear in z between them.
Portal branches labelled LPV that lie in the right liver are relabelled RPV.
build_volume.py also puts the umbilical fissure (falciform ligament line, S4 vs S2/S3) at a fixed
30 deg, about 3 cm left of the LPV umbilical portion. Here it is the angle of that portion, which
lies just in front of the IVC; S4 beyond it joins the left lateral section.
Changed segment meshes are rebuilt with the same marching-cubes settings as build_assets.py.
"""
import sys, json, zlib, base64, numpy as np, trimesh, fast_simplification
from scipy import ndimage as ndi
from skimage import measure

PATH = sys.argv[1] if len(sys.argv) > 1 else 'data/anatomy.js'
Z_PORTAL_R, Z_PORTAL_IV, Z_II = 1110.0, 1124.0, 1118.0   # same planes as build_volume.py
Z_TRUNK = 1150.0                            # MHV is a single trunk above this level
LPV_MARGIN = 5.0                            # degrees right of the boundary before an LPV branch counts as right liver

src = open(PATH).read()
head = src[:src.index('window.ANATOMY=')]
A = json.loads(src[src.index('{'):src.rindex('}') + 1])
L, CENTER, Q = A['labels'], np.array(A['center']), A['q']
V = A['volume']
NX, NY, NZ = V['dims']
LO, P = np.array(V['lo']), V['pitch']
lab = np.frombuffer(zlib.decompress(base64.b64decode(V['data'])), np.uint8).reshape(NZ, NY, NX).transpose(2, 1, 0).copy()
xs = LO[0] + (np.arange(NX) + .5) * P
ys = LO[1] + (np.arange(NY) + .5) * P
zs = LO[2] + (np.arange(NZ) + .5) * P
X, Y = np.meshgrid(xs, ys, indexing='ij')

# IVC centre per axial slice (fallback outside the IVC)
ivc = lab == L['ivc']
cx, cy = np.full(NZ, -15.0), np.full(NZ, -115.0)
for k in range(NZ):
    m = ivc[:, :, k]
    if m.sum() > 20:
        cx[k], cy[k] = X[m].mean(), Y[m].mean()


def angle(x, y, k):
    """Angle about the IVC centre in the axial plane: 0 = anterior, negative = patient right."""
    return np.degrees(np.arctan2(x - cx[k], -(y - cy[k])))


# gallbladder fossa: mean angle of the GB, at its top level
gb = lab == L['gb']
gk = np.nonzero(gb.any((0, 1)))[0]
th_gb = np.mean(np.concatenate([angle(X[gb[:, :, k]], Y[gb[:, :, k]], k) for k in gk]))
z_gb = zs[gk[-1]]

# MHV trunk below the confluence: angle of the outer half of the thickest MHV section per slice
mhv = lab == L['mhv']
edt = ndi.distance_transform_edt(mhv)
top = []
for k in np.nonzero(zs >= Z_TRUNK)[0]:
    lb, n = ndi.label(mhv[:, :, k])
    if n == 0:
        continue
    b = 1 + int(np.argmax(ndi.maximum(edt[:, :, k], lb, range(1, n + 1))))
    sel = lb == b
    r = np.hypot(X[sel] - cx[k], Y[sel] - cy[k])
    top.append(angle(X[sel], Y[sel], k)[r >= np.median(r)].mean())
th_top = float(np.median(top))
tM = np.interp(zs, [z_gb, Z_TRUNK], [th_gb, th_top])
print(f'GB fossa {th_gb:.1f} deg at z {z_gb:.0f}, MHV trunk {th_top:.1f} deg from z {Z_TRUNK:.0f}')

# relabel S4/S4b/S5/S8 by the new boundary (S1, S2/S3, S6/S7 and vessels untouched)
before = lab.copy()
codes = [L['s4'], L['s4b'], L['s5'], L['s8']]
for k in range(NZ):
    sl = lab[:, :, k]
    m = np.isin(sl, codes)
    if not m.any():
        continue
    right = angle(X, Y, k) < tM[k]
    sl[m & right] = L['s8'] if zs[k] >= Z_PORTAL_R else L['s5']
    sl[m & ~right] = L['s4'] if zs[k] >= Z_PORTAL_IV else L['s4b']


def right_liver(x, y, k):
    return (angle(x, y, k) < tM[k] - LPV_MARGIN) & (np.hypot(x - cx[k], y - cy[k]) > 20)


i, j, k = np.nonzero(lab == L['lpv'])
mv = right_liver(xs[i], ys[j], k)
lab[i[mv], j[mv], k[mv]] = L['rpv']

# umbilical fissure: median angle of the thick LPV umbilical portion in front of the hilum
i, j, k = np.nonzero(ndi.distance_transform_edt(lab == L['lpv']) * P >= 2.5)
far = np.hypot(xs[i] - cx[k], ys[j] - cy[k]) >= 25
th_uf = float(np.median(angle(xs[i], ys[j], k)[far]))
print(f'umbilical fissure {th_uf:.1f} deg')
# S4 left of it joins the left lateral section, taking S2 or S3 from the nearest lateral voxel in its slice
for k in range(NZ):
    sl = lab[:, :, k]
    a = angle(X, Y, k)
    m = np.isin(sl, [L['s4'], L['s4b']]) & (a >= th_uf) & (a <= 150)
    if not m.any():
        continue
    lat = np.isin(sl, [L['s2'], L['s3']])
    if lat.any():
        _, (ii, jj) = ndi.distance_transform_edt(~lat, return_indices=True)
        sl[m] = sl[ii[m], jj[m]]
    else:
        sl[m] = L['s2'] if zs[k] > Z_II else L['s3']
print('voxels changed', int((lab != before).sum()), '- LPV to RPV', int(mv.sum()))


# ---- meshes
def b64(a):
    return base64.b64encode(zlib.compress(a.tobytes(), 9)).decode()


def unpack(m):
    v = np.frombuffer(zlib.decompress(base64.b64decode(m['pos'])), np.int16).reshape(-1, 3) / Q + CENTER
    f = np.frombuffer(zlib.decompress(base64.b64decode(m['idx'])), np.uint32 if m['i32'] else np.uint16).reshape(-1, 3)
    return trimesh.Trimesh(v, f.astype(np.int64), process=False)


def pack(key, m):
    v = np.round((m.vertices - CENTER) * Q)
    assert np.abs(v).max() < 32767, key
    f = m.faces.astype(np.uint32 if len(m.vertices) > 65535 else np.uint16)
    return dict(key=key, nv=len(m.vertices), nf=len(m.faces), i32=bool(f.dtype == np.uint32),
                pos=b64(v.astype(np.int16)), idx=b64(f))


def simplify(m, target):
    if len(m.faces) <= target:
        return m
    v, f = fast_simplification.simplify(m.vertices.astype(np.float32), m.faces.astype(np.int32),
                                        target_reduction=1 - target / len(m.faces))
    return trimesh.Trimesh(v, f, process=True)


vessel_codes = [L[k] for k in ('mpv', 'rpv', 'lpv', 'rhv', 'mhv', 'lhv', 'ivc', 'gb', 'bile')]


def segment_mesh(codes):  # as in build_assets.py
    m = np.isin(lab, codes)
    m = m | (ndi.binary_dilation(m, iterations=2) & np.isin(lab, vessel_codes))
    sl = ndi.find_objects(m.astype(np.uint8))[0]
    sl = tuple(slice(max(0, a.start - 3), a.stop + 3) for a in sl)
    f = ndi.gaussian_filter(m[sl].astype(np.float32), 1.0)
    v, fc, _, _ = measure.marching_cubes(f, 0.5)
    v = (v + [a.start for a in sl] + 0.5) * P + LO
    tm = trimesh.Trimesh(v, fc[:, ::-1], process=True)
    trimesh.smoothing.filter_taubin(tm, iterations=12)
    return simplify(tm, 5000)


idx = {m['key']: n for n, m in enumerate(A['meshes'])}
changed = {s for s in range(1, 9)
           if any(((lab == L[x]) != (before == L[x])).any() for x in (['s4', 's4b'] if s == 4 else [f's{s}']))}
for s in sorted(changed):
    codes = [L['s4'], L['s4b']] if s == 4 else [L[f's{s}']]
    A['meshes'][idx[f'seg{s}']] = pack(f'seg{s}', segment_mesh(codes))
    print('rebuilt seg', s)

# LPV mesh faces in the right liver move to the RPV mesh
lpv, rpv = unpack(A['meshes'][idx['lpv']]), unpack(A['meshes'][idx['rpv']])
c = lpv.triangles_center
kc = np.clip(np.round((c[:, 2] - LO[2]) / P - .5).astype(int), 0, NZ - 1)
mv = right_liver(c[:, 0], c[:, 1], kc)
if mv.any():
    moved, kept = lpv.submesh([np.nonzero(mv)[0]], append=True), lpv.submesh([np.nonzero(~mv)[0]], append=True)
    A['meshes'][idx['lpv']] = pack('lpv', kept)
    A['meshes'][idx['rpv']] = pack('rpv', trimesh.util.concatenate([rpv, moved]))
print('LPV mesh faces to RPV', int(mv.sum()))

V['data'] = b64(np.ascontiguousarray(np.transpose(lab, (2, 1, 0))))
open(PATH, 'w').write(head + 'window.ANATOMY=' + json.dumps(A, separators=(',', ':')) + ';\n')
names = ['s1', 's2', 's3', 's4', 's4b', 's5', 's6', 's7', 's8']
cnt = {n: int((lab == L[n]).sum()) for n in names}
tot = sum(cnt.values())
print('segment volume %:', ' '.join(f'{n.upper()} {cnt[n] / tot * 100:.1f}' for n in names))
