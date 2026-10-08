"""Couinaud segments: sections from hepatic veins and fissures, segments from portal territories.

Usage: python portal_segments.py data/anatomy.js
Run after build_assets.py. It replaces the angle-based split from build_volume.py and is safe to run again.

Couinaud (1957) separates sectors by the scissurae that carry the hepatic veins and defines segments as the
territories of the segmental portal branches. Section names follow Brisbane 2000 (Strasberg et al.).
1. Sections. Angles are taken about the IVC axis in each axial slice.
   - Midplane (right vs left hemiliver): from the gallbladder fossa at the GB top to the MHV trunk below
     the confluence (Cantlie line).
   - Right intersectional plane (anterior vs posterior right section): the RHV trunk in each slice.
     Where the slice has no RHV, the nearest portal branch decides.
   - Umbilical fissure (S4 vs left lateral section): the umbilical portion of the LPV.
2. Segments. The portal labels are skeletonized into one tree rooted at the lower end of the MPV. Each
   SEEDS point sits on a segmental branch, and it and everything downstream belong to that segment. Inside
   its section, a liver voxel takes the segment of the nearest segmental-branch lumen voxel (nearest-branch
   territories as in Soler et al. 2001 and Selle et al. 2002). This splits S5/S8, S6/S7, S2/S3 and S4a/S4b.
3. Portal voxels are renamed by the hemiliver they supply: RPV for S5-S8, LPV for S2-S4, MPV for the shared trunk.
S1 keeps the BodyParts3D caudate lobe. Segment and portal meshes are rebuilt with the same settings as build_assets.py.
"""
import sys, json, zlib, base64, numpy as np, trimesh, fast_simplification
from scipy import ndimage as ndi
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import dijkstra
from skimage import measure
from skimage.morphology import skeletonize

PATH = sys.argv[1] if len(sys.argv) > 1 else 'data/anatomy.js'
# Points on the segmental branches of this BodyParts3D portal tree (mm), each about 4 mm past its
# branching point. Right portal vein: anterior branch -> S5 (inferior), S8 (superior); posterior
# branch -> S6 (inferior), S7 (superior). Umbilical portion of the LPV -> S2 (posterior-superior),
# S3 (anterior-inferior), S4a/S4b (rightward, above/below). The second S8 point is a branch that
# leaves the LPV and runs right into the anterior-superior right liver (variant origin, S8 territory).
SEEDS = {
    's2': [(-3.25, -153.75, 1126.75), (-6.25, -155.25, 1134.25)],
    's3': [(-3.25, -165.75, 1131.25), (-4.75, -170.25, 1119.25)],
    's4': [(-15.25, -171.75, 1129.75)],
    's4b': [(-15.25, -176.25, 1122.25)],
    's5': [(-49.75, -146.25, 1113.25)],
    's6': [(-64.75, -114.75, 1105.75)],
    's7': [(-61.75, -111.75, 1111.75)],
    's8': [(-48.25, -146.25, 1120.75), (-25.75, -152.25, 1125.25)],
}
RIGHT = ('s5', 's6', 's7', 's8')
SECTIONS = {'ant': ('s5', 's8'), 'post': ('s6', 's7'), 'med': ('s4', 's4b'), 'lat': ('s2', 's3')}
Z_TRUNK = 1150.0     # MHV is a single trunk above this level
TRUNK_SHARE = 0.75   # a trunk voxel is RPV/LPV when this share of its downstream territory is on one side
SMOOTH = 1.5         # voxels, Gaussian vote that removes single-voxel jaggies at segment borders

src = open(PATH).read()
head = src[:src.index('window.ANATOMY=')]
A = json.loads(src[src.index('{'):src.rindex('}') + 1])
L, CENTER, Q = A['labels'], np.array(A['center']), A['q']
V = A['volume']
NX, NY, NZ = V['dims']
LO, P = np.array(V['lo']), V['pitch']
lab = np.frombuffer(zlib.decompress(base64.b64decode(V['data'])), np.uint8).reshape(NZ, NY, NX).transpose(2, 1, 0).copy()
before = lab.copy()
SEG = ['s1', 's2', 's3', 's4', 's4b', 's5', 's6', 's7', 's8']
PORTAL = [L['mpv'], L['rpv'], L['lpv']]

# ---- portal skeleton tree
portal = np.isin(lab, PORTAL)
sk = skeletonize(portal) > 0
node = np.argwhere(sk)
xyz = LO + (node + .5) * P
nid = np.full(lab.shape, -1, np.int32)
nid[tuple(node.T)] = np.arange(len(node))
rows, cols, w = [], [], []
for d in np.argwhere(np.ones((3, 3, 3), bool)) - 1:
    if not d.any():
        continue
    nb = node + d
    ok = np.all((nb >= 0) & (nb < lab.shape), 1)
    j = np.full(len(node), -1)
    j[ok] = nid[tuple(nb[ok].T)]
    m = j >= 0
    rows += list(np.nonzero(m)[0]); cols += list(j[m]); w += [float(np.linalg.norm(d)) * P] * int(m.sum())
graph = coo_matrix((w, (rows, cols)), shape=(len(node), len(node))).tocsr()
mpv_nodes = np.nonzero(lab[tuple(node.T)] == L['mpv'])[0]
root = int(mpv_nodes[np.argmin(xyz[mpv_nodes, 2])])
dist, parent = dijkstra(graph, indices=root, return_predecessors=True)
assert np.isfinite(dist).all(), 'portal skeleton is not one connected tree'
children = [[] for _ in node]
for n, p in enumerate(parent):
    if p >= 0:
        children[p].append(n)


def downstream(n):
    out, stack = [], [n]
    while stack:
        a = stack.pop()
        out.append(a)
        stack += children[a]
    return out


seg_of = np.full(len(node), '', object)
for s, pts in SEEDS.items():
    for p in pts:
        d = np.linalg.norm(xyz - p, axis=1)
        n = int(np.argmin(d))
        assert d[n] < 2.5, f'seed {s} {p} is {d[n]:.1f} mm from the portal skeleton; update SEEDS'
        sub = downstream(n)
        assert not seg_of[sub].any(), f'seed {s} {p} overlaps another segment'
        seg_of[sub] = s
for s in SEEDS:
    print(f'{s.upper():4s} branch nodes {int((seg_of == s).sum())}')

# ---- portal lumen voxels inherit the segment of their nearest skeleton node
_, near = ndi.distance_transform_edt(~sk, return_indices=True)
pv = np.argwhere(portal)
pnode = nid[tuple(near[:, pv[:, 0], pv[:, 1], pv[:, 2]])]
pseg = seg_of[pnode]

# ---- section boundaries (angles about the IVC centre; 0 = anterior, negative = patient right)
xs = LO[0] + (np.arange(NX) + .5) * P
ys = LO[1] + (np.arange(NY) + .5) * P
zs = LO[2] + (np.arange(NZ) + .5) * P
X, Y = np.meshgrid(xs, ys, indexing='ij')
ivc = lab == L['ivc']
cx, cy = np.full(NZ, -15.0), np.full(NZ, -115.0)
for k in range(NZ):
    m = ivc[:, :, k]
    if m.sum() > 20:
        cx[k], cy[k] = X[m].mean(), Y[m].mean()


def angle(x, y, k):
    return np.degrees(np.arctan2(x - cx[k], -(y - cy[k])))


def trunk_angle(code, k, rmin=20.0):
    """Angle of the outer half of the thickest section of a vein in slice k (None if absent)."""
    m = lab[:, :, k] == code
    lb, n = ndi.label(m)
    if n == 0:
        return None
    e = ndi.distance_transform_edt(m)
    sel = lb == 1 + int(np.argmax(ndi.maximum(e, lb, range(1, n + 1))))
    r = np.hypot(X[sel] - cx[k], Y[sel] - cy[k])
    if r.max() < rmin:
        return None
    return float(angle(X[sel], Y[sel], k)[r >= np.median(r)].mean())


gb = lab == L['gb']
gk = np.nonzero(gb.any((0, 1)))[0]
th_gb = float(np.mean(np.concatenate([angle(X[gb[:, :, k]], Y[gb[:, :, k]], k) for k in gk])))
th_top = float(np.median([a for k in np.nonzero(zs >= Z_TRUNK)[0] if (a := trunk_angle(L['mhv'], k)) is not None]))
tM = np.interp(zs, [zs[gk[-1]], Z_TRUNK], [th_gb, th_top])
tR = np.full(NZ, np.nan)
for k in range(NZ):
    a = trunk_angle(L['rhv'], k)
    if a is not None:
        tR[k] = a
ok = np.isfinite(tR)
tR[ok] = ndi.median_filter(tR[ok], 5, mode='nearest')

# ---- portal names: segmental branches by hemiliver; trunk by the side it mostly supplies
r_sum = np.isin(seg_of, RIGHT).astype(float)
l_sum = (seg_of != '').astype(float) - r_sum
for n in np.argsort(-dist):                # leaves first
    if parent[n] >= 0:
        r_sum[parent[n]] += r_sum[n]; l_sum[parent[n]] += l_sum[n]
name = np.full(len(node), 'mpv', object)
for n in np.argsort(dist):                 # root first, so a twig with no segment inherits its parent
    tot = r_sum[n] + l_sum[n]
    if tot == 0:
        name[n] = name[parent[n]] if parent[n] >= 0 else 'mpv'
    elif r_sum[n] >= TRUNK_SHARE * tot:
        name[n] = 'rpv'
    elif l_sum[n] >= TRUNK_SHARE * tot:
        name[n] = 'lpv'
lab[tuple(pv.T)] = [L[s] for s in name[pnode]]

# umbilical fissure: thick LPV (left tree) in front of the hilum, i.e. the umbilical portion
pr = ndi.distance_transform_edt(portal)[tuple(pv.T)] * P
core = pv[(name[pnode] == 'lpv') & (pr >= 2.5)]
ck = core[:, 2]
ux, uy = xs[core[:, 0]], ys[core[:, 1]]
far = np.hypot(ux - cx[ck], uy - cy[ck]) >= 25
tU = float(np.median(angle(ux[far], uy[far], ck[far])))
print(f'midplane {th_gb:.1f} -> {th_top:.1f} deg, umbilical fissure {tU:.1f} deg, RHV found in {int(ok.sum())} slices')

# ---- liver voxels outside S1: section from the planes, then segment from the nearest segmental branch
liver = np.isin(lab, [L[s] for s in SEG])
free = liver & (lab != L['s1'])
box = ndi.find_objects(liver.astype(np.uint8))[0]
box = tuple(slice(max(0, b.start - 2), b.stop + 2) for b in box)
seed = np.zeros(lab.shape, np.uint8)
has = pseg != ''
seed[tuple(pv[has].T)] = [L[s] for s in pseg[has]]
seed = seed[box]


def nearest(segs):
    _, ind = ndi.distance_transform_edt(~np.isin(seed, [L[s] for s in segs]), return_indices=True)
    return seed[tuple(ind)]


ks = np.arange(box[2].start, box[2].stop)
th = np.stack([angle(X[box[:2]], Y[box[:2]], k) for k in ks], -1)
tMb, tRb = tM[ks][None, None, :], tR[ks][None, None, :]
right = th < tMb
near_right = nearest(RIGHT)
ant = np.where(np.isfinite(tRb), th >= np.nan_to_num(tRb), np.isin(near_right, [L['s5'], L['s8']]))
section = {'ant': right & ant, 'post': right & ~ant, 'med': ~right & (th < tU), 'lat': ~right & (th >= tU)}
assigned = np.zeros(seed.shape, np.uint8)
for sec, segs in SECTIONS.items():
    assigned[section[sec]] = nearest(segs)[section[sec]]
fb = free[box]
votes = np.stack([ndi.gaussian_filter((assigned == L[s]).astype(np.float32), SMOOTH) for s in SEEDS])
best = np.array([L[s] for s in SEEDS], np.uint8)[np.argmax(votes, 0)]
lab[box][fb] = best[fb]

changed = int((lab != before).sum())
print('voxels changed', changed)


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
for s in range(1, 9):
    codes = [L['s4'], L['s4b']] if s == 4 else [L[f's{s}']]
    if any(((lab == c) != (before == c)).any() for c in codes):
        A['meshes'][idx[f'seg{s}']] = pack(f'seg{s}', segment_mesh(codes))
        print('rebuilt seg', s)

# portal meshes: each face goes to the name of the nearest portal voxel
if any(((lab == c) != (before == c)).any() for c in PORTAL):
    pm = trimesh.util.concatenate([unpack(A['meshes'][idx[k]]) for k in ('mpv', 'rpv', 'lpv')])
    _, pind = ndi.distance_transform_edt(~portal, return_indices=True)
    c = np.clip(np.floor((pm.triangles_center - LO) / P).astype(int), 0, np.array(lab.shape) - 1)
    face_lab = lab[tuple(pind[:, c[:, 0], c[:, 1], c[:, 2]])]
    for k in ('mpv', 'rpv', 'lpv'):
        A['meshes'][idx[k]] = pack(k, pm.submesh([np.nonzero(face_lab == L[k])[0]], append=True))
        print(f'{k} mesh faces', int((face_lab == L[k]).sum()))

V['data'] = b64(np.ascontiguousarray(np.transpose(lab, (2, 1, 0))))
open(PATH, 'w').write(head + 'window.ANATOMY=' + json.dumps(A, separators=(',', ':')) + ';\n')

cnt = {s: int((lab == L[s]).sum()) for s in SEG}
tot = sum(cnt.values())
print('segment volume %:', ' '.join(f'{s.upper()} {cnt[s] / tot * 100:.1f}' for s in SEG))
