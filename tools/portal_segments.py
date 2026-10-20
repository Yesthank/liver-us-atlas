"""Couinaud segments: sections from hepatic veins and fissures, segments from portal territories.

Usage: python portal_segments.py data/anatomy.js [out.js] [orig.js]
  out.js defaults to the input; orig.js (the volume before tools/anatomy_fixes.py) makes every segment mesh whose voxels
  differ from it be rebuilt, so liver voxels added upstream reach the 3D meshes.
Run after build_assets.py. It replaces the angle-based split from build_volume.py and is safe to run again.

Couinaud (1957) separates sectors by the scissurae that carry the hepatic veins and defines segments as the
territories of the segmental portal branches. Section names follow Brisbane 2000 (Strasberg et al.).
1. Sectors. Angles are taken about the IVC axis in each axial slice.
   - Midplane (right vs left hemiliver): from the gallbladder fossa at the GB top to the MHV trunk below
     the confluence (Cantlie line).
   - Right portal scissura (S5/S8 vs S6/S7): the RHV, leaving out tributaries within 10 mm of the S7/S8
     pedicles, and kept between the anterior and posterior pedicles.
     Left portal scissura (S2 vs S3): the LHV. Each is the median angle of the vein in the slice,
     smoothed over z and held beyond the vein's ends.
   - Umbilical fissure (S4 vs S3): the umbilical portion of the LPV.
2. Segments. The portal labels are skeletonized into one tree rooted at the lower end of the MPV. Each
   SEEDS point sits on a segmental branch, and it and everything downstream belong to that segment. Inside
   its sector, a liver voxel takes the segment of the nearest segmental-branch lumen voxel (nearest-branch
   territories as in Soler et al. 2001 and Selle et al. 2002). This splits S5/S8, S6/S7 and S4a/S4b.
3. Portal voxels are renamed by the hemiliver they supply: RPV for S5-S8, LPV for S2-S4, MPV for the shared trunk.
S1 keeps the BodyParts3D caudate lobe (plus stray boundary voxels that touch only S1). Segment meshes are rebuilt with the same settings as build_assets.py,
and the existing portal meshes are re-split by the new portal names.
"""
import sys, os, json, zlib, base64, numpy as np
from scipy import ndimage as ndi
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import dijkstra
from skimage import measure
from skimage.morphology import skeletonize
try:
    import trimesh, fast_simplification
    HAVE_TRIMESH = True
except ImportError:                      # fallback only: numpy surface-nets mesher from anatomy_fixes.py
    HAVE_TRIMESH = False
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import anatomy_fixes as AFX

# ---- audit 2026-10 (part D2) switches
LAT_BY_PORTAL = True     # S2/S3 = nearest P2/P3 territory inside the lateral section (was: LHV half-plane)
POST_KZ = 2.0            # cranio-caudal distance weight for the S6/S7 nearest-branch split only (was 1)
RHV_TRUNK_ANCHOR = True  # keep tR within +-5 deg of the thick RHV trunk (EDT >= 2 mm); pedicle clip applied last
S8B_SEED_MOVE = True     # second S8 seed moved along its own branch to just right of the Cantlie plane
FORCE_SEG_MESH = os.environ.get('FORCE_SEG_MESH', '1') == '1'   # rebuild every segment mesh whose voxels differ from ORIG
if os.environ.get('D2_OFF') == '1':      # audit baseline: behave exactly like the repo script
    LAT_BY_PORTAL, POST_KZ, RHV_TRUNK_ANCHOR, S8B_SEED_MOVE = False, 1.0, False, False

PATH = sys.argv[1] if len(sys.argv) > 1 else 'data/anatomy.js'
OUT = sys.argv[2] if len(sys.argv) > 2 else PATH
ORIG = sys.argv[3] if len(sys.argv) > 3 else None
# Points on the segmental branches of this BodyParts3D portal tree (mm), each on an unbranched run a few
# mm past its branching point. Right portal vein: anterior branch -> S8 (ascending) and S5 (the inferior
# branch plus a short branch that leaves the S8 stem and descends to fuse with it); posterior branch ->
# S6 (inferior), S7 (superior). Umbilical portion of the LPV -> S2 (posterior-superior), S3 (anterior-
# inferior), S4a/S4b (rightward, above/below). The second S8 point is an anterior-superior branch that
# BodyParts3D labels RPV; it touches the LPV in the voxel model, so the skeleton tree hangs it there.
SEEDS = {
    's2': [(-3.25, -153.75, 1126.75), (-6.25, -155.25, 1134.25)],
    's3': [(-3.25, -165.75, 1131.25), (-4.75, -170.25, 1119.25)],
    's4': [(-15.25, -171.75, 1129.75)],
    's4b': [(-15.25, -174.75, 1123.75)],
    's5': [(-49.75, -146.25, 1113.25), (-51.25, -152.25, 1120.75)],
    's6': [(-64.75, -114.75, 1105.75)],
    's7': [(-61.75, -111.75, 1111.75)],
    's8': [(-52.75, -146.25, 1125.25), (-25.75, -152.25, 1125.25)],
}
RIGHT = ('s5', 's6', 's7', 's8')
SECTORS = {'ant': ('s5', 's8'), 'post': ('s6', 's7'), 'med': ('s4', 's4b'), 's3': ('s3',), 's2': ('s2',)}
if LAT_BY_PORTAL:
    SECTORS = {'ant': ('s5', 's8'), 'post': ('s6', 's7'), 'med': ('s4', 's4b'), 'lat': ('s2', 's3')}
Z_TRUNK = 1150.0     # MHV is a single trunk above this level
TRUNK_SHARE = 0.75   # a trunk voxel is RPV/LPV when this share of its downstream segmental-branch nodes is on one side
SMOOTH = 1.5         # voxels, Gaussian vote that removes single-voxel jaggies at segment borders

src = open(PATH, encoding='utf-8').read()
head = src[:src.index('window.ANATOMY=')]
A = json.loads(src[src.index('{'):src.rindex('}') + 1])
L, CENTER, Q = A['labels'], np.array(A['center']), A['q']
V = A['volume']
NX, NY, NZ = V['dims']
LO, P = np.array(V['lo']), V['pitch']
lab = np.frombuffer(zlib.decompress(base64.b64decode(V['data'])), np.uint8).reshape(NZ, NY, NX).transpose(2, 1, 0).copy()
before = lab.copy()
SEG = ['s1', 's2', 's3', 's4', 's4b', 's5', 's6', 's7', 's8']
L.setdefault('rapv', 37); L.setdefault('rppv', 38)   # right anterior (S5/S8) / posterior (S6/S7) sectional portal branches
PORTAL = [L['mpv'], L['rpv'], L['lpv'], L['rapv'], L['rppv']]

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


if S8B_SEED_MOVE:
    # The second S8 seed sits on an anterior-superior S8 branch that the voxel tree hangs on the LPV; its first
    # ~10 mm run through S4a. Move the seed along the same branch (largest-subtree child each step) to the first
    # unbranched node right of x = -36 mm (just past the Cantlie plane). The stem left behind seeds nothing and
    # keeps the LPV name of its parent, so no 'RPV' crosses S4a.
    n = int(np.argmin(np.linalg.norm(xyz - np.array(SEEDS['s8'][1]), axis=1)))
    size = np.zeros(len(node))
    for a in np.argsort(-dist):
        size[a] += 1
        if parent[a] >= 0:
            size[parent[a]] += size[a]
    while not (xyz[n, 0] <= -36 and len(children[n]) == 1 and len(children[parent[n]]) == 1):
        assert children[n], 'S8b branch ended before x = -36'
        n = max(children[n], key=lambda c: size[c])
    print('S8b seed moved', SEEDS['s8'][1], '->', tuple(np.round(xyz[n], 2)))
    SEEDS['s8'][1] = tuple(float(v) for v in xyz[n])
seg_of = np.full(len(node), '', object)
for s, pts in SEEDS.items():
    for p in pts:
        d = np.linalg.norm(xyz - p, axis=1)
        n = int(np.argmin(d))
        assert d[n] < 2.5, f'seed {s} {p} is {d[n]:.1f} mm from the portal skeleton; update SEEDS'
        assert len(children[n]) == 1 and len(children[parent[n]]) == 1, f'seed {s} {p} is at a branching point'
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


def vein_plane(code, keep=None, rmin=20.0):
    """Scissura angle per slice: median angle of the vein beyond rmin of the IVC centre (and inside keep),
    smoothed over z and held constant beyond the vein's ends."""
    t = np.full(NZ, np.nan)
    for k in range(NZ):
        m = (lab[:, :, k] == code) & (True if keep is None else keep[:, :, k])
        a = angle(X[m], Y[m], k)[np.hypot(X[m] - cx[k], Y[m] - cy[k]) >= rmin]
        if len(a) >= 5:
            t[k] = np.median(a)
    ok = np.isfinite(t)
    t[ok] = ndi.median_filter(t[ok], 5, mode='nearest')
    return ndi.gaussian_filter1d(np.interp(np.arange(NZ), np.nonzero(ok)[0], t[ok]), 3, mode='nearest'), ok


gb = lab == L['gb']
gk = np.nonzero(gb.any((0, 1)))[0]
th_gb = float(np.mean(np.concatenate([angle(X[gb[:, :, k]], Y[gb[:, :, k]], k) for k in gk])))
th_top = float(np.median([a for k in np.nonzero(zs >= Z_TRUNK)[0] if (a := trunk_angle(L['mhv'], k)) is not None]))
tM = np.interp(zs, [zs[gk[-1]], Z_TRUNK], [th_gb, th_top])
# RHV tributaries running alongside the S7 or S8 pedicle drain that segment; leave them out so the plane
# follows the scissural course between the sectors
upper_ped = np.zeros(lab.shape, bool)
upper_ped[tuple(pv[np.isin(pseg, ('s7', 's8'))].T)] = True
rhv_thick = ndi.distance_transform_edt(lab == L['rhv']) * P >= 2.0          # RHV trunk (radius >= 2 mm)
keepR = ndi.distance_transform_edt(~upper_ped) * P >= 10
if RHV_TRUNK_ANCHOR:
    keepR = keepR | rhv_thick                    # the pedicle-proximity exclusion never removes trunk voxels
tR, okR = vein_plane(L['rhv'], keep=keepR)
tT, okT = vein_plane(L['rhv'], keep=rhv_thick)
tL, okL = vein_plane(L['lhv'])
# the right portal scissura runs between the posterior (S6/S7) and anterior (S5/S8) pedicles
pk = pv[:, 2]
pa = angle(xs[pv[:, 0]], ys[pv[:, 1]], pk)
beyond = np.hypot(xs[pv[:, 0]] - cx[pk], ys[pv[:, 1]] - cy[pk]) >= 30
lo, hi = np.full(NZ, -np.inf), np.full(NZ, np.inf)
for k in np.unique(pk):
    a_post = pa[(pk == k) & beyond & np.isin(pseg, ['s6', 's7'])]
    a_ant = pa[(pk == k) & beyond & np.isin(pseg, ['s5', 's8'])]
    if len(a_post):
        lo[k] = a_post.max() + 3
    if len(a_ant):
        hi[k] = a_ant.min() - 3
tR = ndi.gaussian_filter1d(np.where(lo < hi, np.clip(tR, lo, hi), tR), 1.5, mode='nearest')
if RHV_TRUNK_ANCHOR:
    tR = np.where(okT, np.clip(tR, tT - 5, tT + 5), tR)      # anchor to the trunk where the trunk exists
tR = np.where(lo < hi, np.clip(tR, lo, hi), tR)              # pedicle gap rule stays the last word

# ---- portal names: segmental branches by hemiliver; trunk by the side it mostly supplies
r_sum = np.isin(seg_of, RIGHT).astype(float)
l_sum = (seg_of != '').astype(float) - r_sum
a_sum = np.isin(seg_of, ('s5', 's8')).astype(float)    # S5/S8 branch nodes downstream
p_sum = np.isin(seg_of, ('s6', 's7')).astype(float)    # S6/S7 branch nodes downstream
for n in np.argsort(-dist):                # leaves first
    if parent[n] >= 0:
        r_sum[parent[n]] += r_sum[n]; l_sum[parent[n]] += l_sum[n]
        a_sum[parent[n]] += a_sum[n]; p_sum[parent[n]] += p_sum[n]
name = np.full(len(node), 'mpv', object)
for n in np.argsort(dist):                 # root first, so a twig with no segment inherits its parent
    tot = r_sum[n] + l_sum[n]
    if tot == 0:
        name[n] = name[parent[n]] if parent[n] >= 0 else 'mpv'
    elif r_sum[n] >= TRUNK_SHARE * tot:
        # right side: sectional branches get their own names; only the short common trunk stays RPV
        name[n] = 'rapv' if p_sum[n] == 0 else 'rppv' if a_sum[n] == 0 else 'rpv'
    elif l_sum[n] >= TRUNK_SHARE * tot:
        name[n] = 'lpv'
name_core = name.copy()      # umbilical-portion detection below uses the names before the stem rename
if S8B_SEED_MOVE:
    # the stem between the LPV and the moved S8 seed lies in S4a: name it after the umbilical-portion side it
    # hangs from (LPV), so no 'RPV' label runs through S4a. Only this one stem is renamed.
    n = int(np.argmin(np.linalg.norm(xyz - np.array(SEEDS['s8'][1]), axis=1)))
    stem = []
    a = parent[n]
    while a >= 0 and name[a] in ('rpv', 'rapv', 'rppv'):
        stem.append(a); a = parent[a]
    if a >= 0 and name[a] == 'lpv':
        name[stem] = 'lpv'
        print(f'S8b stem renamed (right portal)->lpv: {len(stem)} skeleton nodes, from {np.round(xyz[stem[-1]], 1)} to {np.round(xyz[stem[0]], 1)}')
lab[tuple(pv.T)] = [L[s] for s in name[pnode]]

# umbilical fissure: thick LPV (left tree) in front of the hilum, i.e. the umbilical portion
pr = ndi.distance_transform_edt(portal)[tuple(pv.T)] * P
core = pv[(name_core[pnode] == 'lpv') & (pr >= 2.5)]
ck = core[:, 2]
ux, uy = xs[core[:, 0]], ys[core[:, 1]]
far = np.hypot(ux - cx[ck], uy - cy[ck]) >= 25
tU = float(np.median(angle(ux[far], uy[far], ck[far])))
print(f'midplane {th_gb:.1f} -> {th_top:.1f} deg, umbilical fissure {tU:.1f} deg, '
      f'RHV {int(okR.sum())} slices, LHV {int(okL.sum())} slices')

# ---- liver voxels outside S1: sector from the planes, then segment from the nearest segmental branch
liver = np.isin(lab, [L[s] for s in SEG])
free = liver & (lab != L['s1'])
box = ndi.find_objects(liver.astype(np.uint8))[0]
box = tuple(slice(max(0, b.start - 2), b.stop + 2) for b in box)
seed = np.zeros(lab.shape, np.uint8)
has = pseg != ''
seed[tuple(pv[has].T)] = [L[s] for s in pseg[has]]
seed = seed[box]


def nearest(segs, sampling=None):
    _, ind = ndi.distance_transform_edt(~np.isin(seed, [L[s] for s in segs]), sampling=sampling, return_indices=True)
    return seed[tuple(ind)]


ks = np.arange(box[2].start, box[2].stop)
th = np.stack([angle(X[box[:2]], Y[box[:2]], k) for k in ks], -1)
tMb, tRb, tLb = (t[ks][None, None, :] for t in (tM, tR, tL))
right = th < tMb
ant = th >= tRb
lat = ~right & (th >= tU)
sector = {'ant': right & ant, 'post': right & ~ant, 'med': ~right & ~lat, 's3': lat & (th < tLb), 's2': lat & (th >= tLb),
          'lat': lat}
assigned = np.zeros(seed.shape, np.uint8)
for sec, segs in SECTORS.items():
    smp = (1, 1, POST_KZ) if (sec == 'post' and POST_KZ != 1) else None
    assigned[sector[sec]] = L[segs[0]] if len(segs) == 1 else nearest(segs, smp)[sector[sec]]
# ---- audit diagnostics (D2): right scissura vs RHV trunk, pedicle-in-wrong-sector counts
for z in (1110, 1115, 1120, 1128, 1133, 1137, 1140):
    k = int(np.argmin(np.abs(zs - z)))
    print(f'  z{z}: tR {tR[k]:7.1f}  RHV trunk {tT[k] if okT[k] else float("nan"):7.1f}  offset {tR[k] - tT[k] if okT[k] else float("nan"):+5.1f}')
pvb = pv - np.array([b.start for b in box])
inb = np.all((pvb >= 0) & (pvb < np.array(seed.shape)), 1)
a_sec = ant[tuple(pvb[inb].T)] & right[tuple(pvb[inb].T)]
p_sec = (~ant[tuple(pvb[inb].T)]) & right[tuple(pvb[inb].T)]
ps = pseg[inb]
print('  pedicle voxels in wrong right sector: S6/S7 in ant', int((np.isin(ps, ('s6', 's7')) & a_sec).sum()),
      ' S5/S8 in post', int((np.isin(ps, ('s5', 's8')) & p_sec).sum()))
fb = free[box]
votes = np.stack([ndi.gaussian_filter((assigned == L[s]).astype(np.float32), SMOOTH) for s in SEEDS])
best = np.array([L[s] for s in SEEDS], np.uint8)[np.argmax(votes, 0)]
lab[box][fb] = best[fb]
# drop specks: segment pieces under 20 voxels join the most common adjacent segment; stray voxels on the
# caudate surface whose only liver neighbour is S1 join S1
for s in SEG[1:]:
    pieces, n = ndi.label(lab == L[s])
    if n < 2:
        continue
    size = ndi.sum(np.ones_like(pieces), pieces, range(1, n + 1))
    for i in np.nonzero(size < 20)[0] + 1:
        m = pieces == i
        ring = lab[ndi.binary_dilation(m) & ~m]
        ring = ring[np.isin(ring, [L[t] for t in SEG if t != s])]
        if len(ring):
            lab[m] = np.bincount(ring).argmax()

changed = int((lab != before).sum())
print('voxels changed', changed)




# ---- meshes
def b64(a):
    return base64.b64encode(zlib.compress(a.tobytes(), 9)).decode()


idx = {m['key']: n for n, m in enumerate(A['meshes'])}
ref = before
if ORIG:                                   # compare with the original repo volume so voxels added upstream count
    _s = open(ORIG, encoding='utf-8').read(); _A = json.loads(_s[_s.index('{'):_s.rindex('}') + 1])
    ref = np.frombuffer(zlib.decompress(base64.b64decode(_A['volume']['data'])), np.uint8).reshape(NZ, NY, NX).transpose(2, 1, 0)
vessel_codes = [L[k] for k in ('mpv', 'rpv', 'lpv', 'rapv', 'rppv', 'rhv', 'mhv', 'lhv', 'ivc', 'gb', 'bile')]

if HAVE_TRIMESH:
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

    def seg_pack(key, codes):
        return pack(key, segment_mesh(codes))
else:
    def seg_pack(key, codes):           # prototype fallback: surface nets on a 3 mm grid (about 3-6k faces)
        m = np.isin(lab, codes)
        m = m | (ndi.binary_dilation(m, iterations=2) & np.isin(lab, vessel_codes))
        v, f = AFX.surface_nets(m, LO, P, sigma=1.0, step=2, smooth=12)
        return AFX.pack(A, key, v, f)

for s in range(1, 9):
    codes = [L['s4'], L['s4b']] if s == 4 else [L[f's{s}']]
    if any(((lab == c) != ((ref if FORCE_SEG_MESH else before) == c)).any() for c in codes):
        A['meshes'][idx[f'seg{s}']] = seg_pack(f'seg{s}', codes)
        print('rebuilt seg', s, A['meshes'][idx[f'seg{s}']]['nf'], 'faces')

# portal meshes: each face goes to the name of the nearest portal voxel
if any(((lab == c) != (before == c)).any() for c in PORTAL):
    if HAVE_TRIMESH:
        pm = trimesh.util.concatenate([unpack(A['meshes'][idx[k]]) for k in ('mpv', 'rpv', 'lpv')])
        _, pind = ndi.distance_transform_edt(~portal, return_indices=True)
        c = np.clip(np.floor((pm.triangles_center - LO) / P).astype(int), 0, np.array(lab.shape) - 1)
        face_lab = lab[tuple(pind[:, c[:, 0], c[:, 1], c[:, 2]])]
        for k, codes in (('mpv', [L['mpv']]), ('rpv', [L['rpv'], L['rapv'], L['rppv']]), ('lpv', [L['lpv']])):
            A['meshes'][idx[k]] = pack(k, pm.submesh([np.nonzero(np.isin(face_lab, codes))[0]], append=True))
            print(f'{k} mesh faces', int(np.isin(face_lab, codes).sum()))
    else:
        vs, fs, off = [], [], 0
        for k in ('mpv', 'rpv', 'lpv'):
            v, f = AFX.unpack(A, A['meshes'][idx[k]]); vs.append(v); fs.append(f + off); off += len(v)
        parts = AFX.split_faces(np.concatenate(vs), np.concatenate(fs), lab, LO, P, PORTAL)
        for k, codes in (('mpv', ['mpv']), ('rpv', ['rpv', 'rapv', 'rppv']), ('lpv', ['lpv'])):
            vv, ff, o = [], [], 0
            for c in codes:
                v, f = parts[L[c]]; vv.append(v); ff.append(f + o); o += len(v)
            A['meshes'][idx[k]] = AFX.pack(A, k, np.concatenate(vv), np.concatenate(ff))
            print(f'{k} mesh faces', sum(len(f) for f in ff))

V['data'] = b64(np.ascontiguousarray(np.transpose(lab, (2, 1, 0))))
open(OUT, 'w', encoding='utf-8').write(head + 'window.ANATOMY=' + json.dumps(A, separators=(',', ':')) + ';\n')

cnt = {s: int((lab == L[s]).sum()) for s in SEG}
tot = sum(cnt.values())
print('segment volume %:', ' '.join(f'{s.upper()} {cnt[s] / tot * 100:.1f}' for s in SEG))
print('segment volume mL:', ' '.join(f'{s.upper()} {cnt[s] * P ** 3 / 1000:.0f}' for s in SEG))
