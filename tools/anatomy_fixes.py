"""Deterministic in-place anatomy fixes for the liver US atlas label volume and meshes (audit 2026-10, part D2).

Usage: python tools/anatomy_fixes.py <in anatomy.js> <out anatomy.js>
  Run on the anatomy.js produced by build_assets.py (it stops if the new labels gda/ha/lhd already exist), then run
  tools/portal_segments.py on its output:
    python tools/anatomy_fixes.py data/anatomy.js stage1.js
    python tools/portal_segments.py stage1.js data/anatomy.js data/anatomy.js.orig

Pattern follows tools/portal_segments.py: decode the 1.5 mm label volume, edit, re-encode, rewrite meshes that
the edits invalidate. No randomness; numpy + scipy only (meshes are rebuilt with a numpy surface-nets mesher so
trimesh is not needed). Run tools/portal_segments.py (patched) AFTER this script: it re-segments the liver
(Couinaud) and rebuilds the segment meshes, which also absorbs the liver voxels added in steps 8 and 9.

Steps (applied in this order; each prints a one-line log):
  1 renv->ivc    renal-vein voxels inside the per-slice convex hull of the IVC become IVC (paint-order artefact)
  2 ha           'celiac' voxels with x < -6 mm (CHA/PHA) get the new label 'ha' (36)
  3 gda          synthetic gastroduodenal artery tube (r 2.0 mm, new label 'gda' 35), painted over soft/fat only
  4 cbd          supraduodenal CBD translated 12 mm anteriorly per axial slice (smoothstep plateau z 1080-1088, ramps
                 1068-1080 and 1088-1098) so it lies anterolateral to the MPV with the HA anteromedial; painted over
                 soft/fat only, vacated voxels -> soft, calibre restored in-slice away from vessels/duodenum
  5 lhd          bile voxels z >= 1108 and x > -24 (left hepatic duct inside S3/S4b) get new label 'lhd' (43)
  6 sv-route     splenic vein body (x < 50) cleared and repainted as a 3 mm-radius tube along the back of the pancreas
                 body (flat, z 1052-1075) that crosses in front of the SMA at z ~1052 and joins the SMV/MPV confluence;
                 painted over pancreas/soft/fat only (duct, SMA, aorta, spla never overwritten)
  7 calibre      SV (tail/hilar part, x >= 50) and SMV grown by one voxel (EDT <= 1.5 mm) over soft/fat/pancreas only
  8 diaph        diaphragm voxels anterior (y < -140) and caudal to the costal margin (per-x margin from rib/cart,
                 capped by the cartilage-9/10 table) are removed and filled from the nearest muscle/soft/fat/liver
  8b diaph-slip  anterior diaphragm just above the margin (z < margin+25) with no lung within 6 mm -> 'muscle' echo
                 (the dome and costophrenic recess keep 'diaph', so the deep bright curve and mirror stay)
  9 hepatorenal  constrained growth of the liver (label s6) into soft/fat between S6 and the anterolateral upper
                 right kidney (z 1024-1068, x < -60, gap <= 26 mm, 1-voxel rim kept, never posterior to the kidney,
                 >= 4.5 mm from duodenum/GB/colon)
Meshes: celiac split into celiac/ha; new gda, lhd meshes; bile mesh shifted with the step-4 field and split
bile/lhd; splv and smv rebuilt from voxels; diaph faces in the step-8 region dropped. Segment meshes: see above.
"""
import sys, os, json, zlib, base64
import numpy as np
from scipy import ndimage as ndi
from scipy.spatial import Delaunay
from scipy.sparse import coo_matrix

NEW_CODES = {'gda': 35, 'ha': 36, 'lhd': 43}


# ------------------------------------------------------------------ io
def b64(a):
    return base64.b64encode(zlib.compress(np.ascontiguousarray(a).tobytes(), 9)).decode()


def read_js(path):
    src = open(path, encoding='utf-8').read()
    head = src[:src.index('window.ANATOMY=')]
    A = json.loads(src[src.index('{'):src.rindex('}') + 1])
    V = A['volume']
    NX, NY, NZ = V['dims']
    lab = np.frombuffer(zlib.decompress(base64.b64decode(V['data'])), np.uint8).reshape(NZ, NY, NX).transpose(2, 1, 0).copy()
    return head, A, lab


def write_js(path, head, A, lab):
    A['volume']['data'] = b64(np.transpose(lab, (2, 1, 0)))
    with open(path, 'w', encoding='utf-8') as f:
        f.write(head + 'window.ANATOMY=' + json.dumps(A, separators=(',', ':')) + ';\n')


# ------------------------------------------------------------------ meshes (numpy only)
def unpack(A, m):
    C, Q = np.array(A['center']), A['q']
    v = np.frombuffer(zlib.decompress(base64.b64decode(m['pos'])), np.int16).reshape(-1, 3) / Q + C
    f = np.frombuffer(zlib.decompress(base64.b64decode(m['idx'])), np.uint32 if m['i32'] else np.uint16).reshape(-1, 3)
    return v.astype(np.float64), f.astype(np.int64)


def pack(A, key, v, f):
    C, Q = np.array(A['center']), A['q']
    q = np.round((v - C) * Q)
    assert np.abs(q).max() < 32767, key
    f = f.astype(np.uint32 if len(v) > 65535 else np.uint16)
    return dict(key=key, nv=int(len(v)), nf=int(len(f)), i32=bool(f.dtype == np.uint32), pos=b64(q.astype(np.int16)), idx=b64(f))


def compact(v, f):
    used = np.unique(f)
    remap = np.full(len(v), -1, np.int64); remap[used] = np.arange(len(used))
    return v[used], remap[f]


def taubin(v, f, iters=10, lam=.5, mu=-.53):
    e = np.concatenate([f[:, [0, 1]], f[:, [1, 2]], f[:, [2, 0]]])
    e = np.concatenate([e, e[:, ::-1]])
    W = coo_matrix((np.ones(len(e)), (e[:, 0], e[:, 1])), shape=(len(v), len(v))).tocsr()
    W.data[:] = 1.0
    deg = np.asarray(W.sum(1)).ravel(); deg[deg == 0] = 1
    for _ in range(iters):
        for s in (lam, mu):
            v = v + s * (W @ v / deg[:, None] - v)
    return v


def surface_nets(mask, LO, P, sigma=1.0, step=1, smooth=10):
    """Closed surface of a voxel mask: Gaussian-smoothed field, naive surface nets at iso 0.5, Taubin smoothing.
    Faces are oriented outward (checked by signed volume)."""
    sl = ndi.find_objects(mask.astype(np.uint8))[0]
    sl = tuple(slice(max(0, s.start - 4), min(n, s.stop + 4)) for s, n in zip(sl, mask.shape))
    fld = ndi.gaussian_filter(mask[sl].astype(np.float32), sigma)[::step, ::step, ::step]
    fld = np.pad(fld, 1)
    ins = fld > .5
    nx, ny, nz = fld.shape
    cs = (nx - 1, ny - 1, nz - 1)
    corner = [(a, b, c) for a in (0, 1) for b in (0, 1) for c in (0, 1)]
    cv = np.stack([fld[a:a + cs[0], b:b + cs[1], c:c + cs[2]] for a, b, c in corner], -1)
    ci = cv > .5
    act = ci.any(-1) & ~ci.all(-1)
    cid = np.full(cs, -1, np.int64)
    idx = np.argwhere(act)
    cid[tuple(idx.T)] = np.arange(len(idx))
    vals = cv[act]
    acc = np.zeros((len(idx), 3)); cnt = np.zeros(len(idx))
    cpos = np.array(corner, float)
    for p in range(8):
        for q in range(p + 1, 8):
            if np.abs(cpos[p] - cpos[q]).sum() != 1:
                continue
            fa, fb = vals[:, p], vals[:, q]
            m = (fa > .5) != (fb > .5)
            t = (.5 - fa[m]) / (fb[m] - fa[m])
            acc[m] += cpos[p] + t[:, None] * (cpos[q] - cpos[p])
            cnt[m] += 1
    vg = idx + acc / cnt[:, None]
    faces = []
    # x edges
    a = ins[:-1, 1:-1, 1:-1]; b = ins[1:, 1:-1, 1:-1]
    for i, j, k in np.argwhere(a != b):
        j += 1; k += 1
        q = [cid[i, j - 1, k - 1], cid[i, j, k - 1], cid[i, j, k], cid[i, j - 1, k]]
        faces.append(q if ins[i, j, k] else q[::-1])
    a = ins[1:-1, :-1, 1:-1]; b = ins[1:-1, 1:, 1:-1]
    for i, j, k in np.argwhere(a != b):
        i += 1; k += 1
        q = [cid[i - 1, j, k - 1], cid[i - 1, j, k], cid[i, j, k], cid[i, j, k - 1]]
        faces.append(q if ins[i, j, k] else q[::-1])
    a = ins[1:-1, 1:-1, :-1]; b = ins[1:-1, 1:-1, 1:]
    for i, j, k in np.argwhere(a != b):
        i += 1; j += 1
        q = [cid[i - 1, j - 1, k], cid[i, j - 1, k], cid[i, j, k], cid[i - 1, j, k]]
        faces.append(q if ins[i, j, k] else q[::-1])
    qd = np.array(faces, np.int64)
    f = np.concatenate([qd[:, [0, 1, 2]], qd[:, [0, 2, 3]]])
    origin = np.array([s.start for s in sl], float)
    v = LO + ((vg + .5 - 1) * step + origin) * P            # grid point (padded, decimated) -> mm (voxel centres)
    v = taubin(v, f, smooth)
    vol = np.einsum('ij,ij->i', v[f[:, 0]], np.cross(v[f[:, 1]], v[f[:, 2]])).sum() / 6
    if vol < 0:
        f = f[:, ::-1]
    return v, f


def split_faces(v, f, lab, LO, P, codes):
    """Assign each face to the nearest voxel among codes (as portal_segments.py does for MPV/RPV/LPV)."""
    m = np.isin(lab, codes)
    _, ind = ndi.distance_transform_edt(~m, return_indices=True)
    c = np.clip(np.floor((v[f].mean(1) - LO) / P).astype(int), 0, np.array(lab.shape) - 1)
    face_lab = lab[tuple(ind[:, c[:, 0], c[:, 1], c[:, 2]])]
    return {code: compact(v, f[face_lab == code]) for code in codes if (face_lab == code).any()}


# ------------------------------------------------------------------ geometry helpers
class Grid:
    def __init__(self, lab, LO, P):
        self.LO, self.P, self.shape = LO, P, lab.shape
        self.x, self.y, self.z = [LO[i] + (np.arange(n) + .5) * P for i, n in enumerate(lab.shape)]

    def i(self, v, ax):
        return int(np.floor((v - self.LO[ax]) / self.P))

    def tube(self, pts, r):
        """Voxels whose centre is within r mm of the polyline pts (mm)."""
        m = np.zeros(self.shape, bool)
        pts = np.asarray(pts, float)
        for a, b in zip(pts[:-1], pts[1:]):
            lo = np.minimum(a, b) - r - 2; hi = np.maximum(a, b) + r + 2
            ix = [np.nonzero((g >= l) & (g <= h))[0] for g, l, h in zip((self.x, self.y, self.z), lo, hi)]
            if any(len(t) == 0 for t in ix):
                continue
            X, Y, Z = np.meshgrid(self.x[ix[0]], self.y[ix[1]], self.z[ix[2]], indexing='ij')
            p = np.stack([X, Y, Z], -1)
            ab = b - a
            t = np.clip(((p - a) @ ab) / (ab @ ab), 0, 1)
            d = np.linalg.norm(p - (a + t[..., None] * ab), axis=-1)
            sub = m[ix[0][0]:ix[0][-1] + 1, ix[1][0]:ix[1][-1] + 1, ix[2][0]:ix[2][-1] + 1]
            sub |= d <= r
        return m


def smoothstep(x, a, b):
    t = np.clip((x - a) / (b - a), 0, 1)
    return t * t * (3 - 2 * t)


CBD_S0 = (0.0, -12.0)                 # mm: +x = patient left (medial for the CBD), -y = anterior
CBD_RAMP = (1068.0, 1080.0, 1088.0, 1098.0)   # w: 0 below z0, smooth up to 1 at z1, 1 until z2, smooth down to 0 at z3


def cbd_shift(lab, L, G, S0, ramp, SOFT):
    """Move each axial slice of 'bile' by w(z)*S0 (rounded to voxels). Old voxels -> soft; new voxels painted
    only over soft/fat. Returns (z, dx, dy) per slice, painted count, blocked-by-label counts, components."""
    P = G.P
    ks = np.nonzero((G.z >= ramp[0] - 2) & (G.z <= ramp[3] + 2))[0]
    zk = G.z[ks]
    w = smoothstep(zk, ramp[0], ramp[1]) * (1 - smoothstep(zk, ramp[2], ramp[3]))
    dx, dy = w * S0[0], w * S0[1]
    new = []
    for n, k in enumerate(ks):
        di, dj = int(np.round(dx[n] / P)), int(np.round(dy[n] / P))
        if di == 0 and dj == 0:
            continue
        wv = np.argwhere(lab[:, :, k] == L['bile'])
        lab[wv[:, 0], wv[:, 1], k] = L['soft']
        new.append((k, wv[:, 0] + di, wv[:, 1] + dj))
    moved, blk = 0, {}
    inv = {c: a for a, c in L.items()}
    keepout = [L[k] for k in ('mpv', 'rpv', 'lpv', 'ha', 'gda', 'celiac', 'smv', 'ivc', 'duod') if k in L]
    for k, ii, jj in new:
        cur = lab[ii, jj, k]
        ok = np.isin(cur, SOFT + [L['bile']])
        for c, n in zip(*np.unique(cur[~ok], return_counts=True)):
            blk[inv[int(c)]] = blk.get(inv[int(c)], 0) + int(n)
        lab[ii[ok], jj[ok], k] = L['bile']; moved += int(ok.sum())
        # restore the calibre where neighbours clipped the moved section: grow in-slice into soft/fat that does
        # not touch a vessel or the duodenum (keeps a >= 1 voxel wall of soft tissue around the duct)
        s = lab[:, :, k]
        for _ in range(2):
            have = int((s == L['bile']).sum())
            if have >= .9 * len(ii):
                break
            near = ndi.binary_dilation(np.isin(s, keepout), np.ones((3, 3), bool))
            cand = ndi.binary_dilation(s == L['bile'], np.ones((3, 3), bool)) & np.isin(s, SOFT) & ~near
            c = np.argwhere(cand)[:len(ii) - have]
            s[c[:, 0], c[:, 1]] = L['bile']; moved += len(c)
    _, nb = ndi.label(lab == L['bile'], np.ones((3, 3, 3)))
    return zk, dx, dy, moved, blk, nb


# ------------------------------------------------------------------ main
def main(src, dst):
    head, A, lab = read_js(src)
    L = A['labels']
    for k, c in NEW_CODES.items():
        assert k not in L and c not in L.values(), f'label {k}/{c} already present: run on the original anatomy.js'
        L[k] = c
    V = A['volume']; LO = np.array(V['lo'], float); P = float(V['pitch'])
    G = Grid(lab, LO, P)
    orig = lab.copy()
    LOG = []

    def log(s):
        print(s); LOG.append(s)

    SEG = [L[s] for s in ('s1', 's2', 's3', 's4', 's4b', 's5', 's6', 's7', 's8')]
    SOFT = [L['soft'], L['fat']]

    # ---- 1 renal vein inside the IVC -> IVC
    n1 = 0
    for k in np.nonzero((G.z > 1012) & (G.z < 1046))[0]:
        s = lab[:, :, k]
        pts = np.argwhere(s == L['ivc'])
        cand = np.argwhere(s == L['renv'])
        if len(pts) < 8 or len(cand) == 0:
            continue
        inside = Delaunay(pts + 0.0).find_simplex(cand + 0.0) >= 0
        sel = cand[inside]
        s[sel[:, 0], sel[:, 1]] = L['ivc']; n1 += len(sel)
    log(f'step1 renv->ivc: {n1} voxels')

    # ---- 2 hepatic artery split
    m2 = (lab == L['celiac']) & (G.x[:, None, None] < -6)
    lab[m2] = L['ha']
    log(f'step2 celiac->ha (x<-6): {int(m2.sum())} voxels; celiac left {int((lab == L["celiac"]).sum())}')

    # ---- 3 GDA: from the lower end of the CHA, in the soft-tissue corridor between duodenal bulb and pancreatic
    #      neck (anterior to MPV), then along the anterior surface of the head to its lateral edge (view-13 level)
    GDA = [(-8.0, -140.0, 1077.0), (-9.0, -141.0, 1067.0), (-9.0, -142.5, 1060.0), (-13.0, -144.0, 1054.0),
           (-20.0, -142.5, 1048.0), (-26.0, -141.0, 1042.0), (-30.0, -142.0, 1035.0), (-31.0, -142.5, 1029.0)]
    t3 = G.tube(GDA, 2.0)
    ok3 = np.isin(lab, SOFT)
    blocked = {k: int(n) for k, n in zip(*np.unique(lab[t3 & ~ok3], return_counts=True)) if k != L['ha']}
    lab[t3 & ok3] = L['gda']
    log(f'step3 gda tube: painted {int((t3 & ok3).sum())}, blocked by ' + str({[n for n, c in L.items() if c == k][0]: v for k, v in blocked.items()}))

    # ---- 4 CBD anterolateral to the MPV: translate the duct per axial slice by w(z)*S0 (S0 in mm, x/y), where
    #      w is a smooth plateau over the supraduodenal hepatoduodenal-ligament segment. The retroduodenal /
    #      intrapancreatic CBD below and the hilar confluence above stay where they are.
    zk, dx, dy, moved, blk, nb = cbd_shift(lab, L, G, CBD_S0, CBD_RAMP, SOFT)
    field4 = np.stack([zk, dx, dy], 1)          # displacement field, reused for the bile mesh below
    log(f'step4 cbd shift S0={CBD_S0} ramp={CBD_RAMP}: painted {moved}, blocked {blk}, bile components {nb}')

    # ---- 5 left hepatic duct label
    m5 = (lab == L['bile']) & (G.z[None, None, :] >= 1108) & (G.x[:, None, None] > -24)
    lab[m5] = L['lhd']
    lb, nb = ndi.label(np.isin(lab, [L['bile'], L['lhd']]), np.ones((3, 3, 3)))
    log(f'step5 bile->lhd: {int(m5.sum())} voxels; bile+lhd components {nb}')

    # ---- 6 splenic vein re-route in front of the SMA (keep SMA and aorta as they are)
    X3 = G.x[:, None, None]
    # The tail/hilar SV (x >= 50) is kept. The new body segment runs along the back of the pancreas body with a
    # flat course (z 1052-1057 up to x 32, then rising to the tail; was 1052-1078 at x 0-40), keeps >= 1.5 mm from the main pancreatic
    # duct, crosses the midline at z ~1052 in the 5-6 mm gap between the duct (y <= -134) and the SMA
    # (y >= -128), and joins the SMV/MPV confluence at the start of the MPV.
    clr = (lab == L['splv']) & (X3 < 50)
    lab[clr] = L['soft']
    SV = [(52.0, -100.0, 1075.0), (46.0, -101.0, 1067.0), (40.0, -102.5, 1060.0), (32.0, -105.5, 1057.0),
          (24.0, -114.0, 1056.0), (14.0, -127.0, 1055.0), (6.0, -131.0, 1053.0), (-2.0, -131.5, 1052.0),
          (-8.0, -128.5, 1057.0), (-11.0, -126.5, 1062.0)]
    t6 = G.tube(SV, 3.0)
    ok6 = np.isin(lab, SOFT + [L['pancreas']])
    b6 = {[n for n, c in L.items() if c == k][0]: int(n) for k, n in zip(*np.unique(lab[t6 & ~ok6], return_counts=True)) if k != L['splv']}
    lab[t6 & ok6] = L['splv']
    # drop stray splv specks (< 10 voxels) that are not connected to the main vein
    lb, n6 = ndi.label(lab == L['splv'], np.ones((3, 3, 3)))
    if n6 > 1:
        sz = ndi.sum(np.ones_like(lb), lb, range(1, n6 + 1))
        for c in np.nonzero(sz < 10)[0] + 1:
            lab[lb == c] = L['soft']
    lb, n6b = ndi.label(lab == L['splv'], np.ones((3, 3, 3)))
    log(f'step6 sv re-route: cleared {int(clr.sum())}, painted {int((t6 & ok6).sum())}, blocked {b6}, splv components {n6b}')

    # ---- 7 calibre: one-voxel isotropic growth of SV (distal) and SMV
    for nm, extra in (('splv', X3 >= 50), ('smv', np.ones((1, 1, 1), bool))):
        m = lab == L[nm]
        d = ndi.distance_transform_edt(~m) * P
        g = (d > 0) & (d <= 1.5) & np.isin(lab, SOFT + [L['pancreas']]) & extra
        lab[g] = L[nm]
        log(f'step7 calibre {nm}: +{int(g.sum())} voxels ({g.sum() * P ** 3 / 1000:.2f} mL)')

    # ---- 8 diaphragm below the anterior costal margin
    bone = np.isin(lab, [L['ribs'], L['cart']]) & (G.y[None, :, None] < -140)
    zmin = np.where(bone.any((1, 2)), np.where(bone, G.z[None, None, :], 9e9).min((1, 2)), 9e9)
    zdata = ndi.minimum_filter1d(zmin, 3)                       # +-1 voxel column
    tx = np.array([-122, -115, -101, -60, -35, 0, 35, 60, 101, 115, 122.])
    tz = np.array([1055, 1080, 1103, 1109, 1133, 1164, 1133, 1109, 1103, 1080, 1055.])
    margin = np.minimum(zdata, np.interp(G.x, tx, tz))
    TOL = 3.0
    XLAT = 150.0     # lateral limit; beyond |x| 122 the margin is the rib-10 tip level (table clamps at 1055)
    m8 = (lab == L['diaph']) & (G.y[None, :, None] < -140) & (G.z[None, None, :] < margin[:, None, None] - TOL) & (np.abs(X3) <= XLAT)
    src_codes = [L['muscle'], L['soft'], L['fat']] + SEG
    srcm = np.isin(lab, src_codes)
    bb = ndi.find_objects(m8.astype(np.uint8))[0]
    bb = tuple(slice(max(0, s.start - 12), min(n, s.stop + 12)) for s, n in zip(bb, lab.shape))
    _, ind = ndi.distance_transform_edt(~srcm[bb], return_indices=True)
    sub = lab[bb]
    mm8 = m8[bb]
    fill = sub[tuple(ind[:, mm8])] if False else sub[ind[0][mm8], ind[1][mm8], ind[2][mm8]]
    sub[mm8] = fill
    lab[bb] = sub
    fc = {[n for n, c in L.items() if c == k][0]: int(n) for k, n in zip(*np.unique(fill, return_counts=True))}
    log(f'step8 diaph below costal margin: removed {int(m8.sum())} voxels ({m8.sum() * P ** 3 / 1000:.1f} mL), filled with {fc}')
    # 8b: the anterior costal slips that remain just above the margin (z < margin + 25 mm) are real diaphragm, but
    #     with no lung behind them they must not render as the bright diaphragm/lung line. Relabel the ones more
    #     than 6 mm from lung as 'muscle' (hypoechoic, like the slips of transversus they interdigitate with).
    #     The dome and the costophrenic recess (lung within 6 mm) keep 'diaph', so the deep bright curve and the
    #     mirror artefact are unchanged. The 3D diaphragm mesh keeps these faces (anatomically present).
    dl8 = ndi.distance_transform_edt(lab != L['lung']) * P
    dh8 = ndi.distance_transform_edt(lab != L['heart']) * P     # substernal/xiphoid diaphragm (views 1-2) is left alone
    m8b = ((lab == L['diaph']) & (G.y[None, :, None] < -140) & (G.z[None, None, :] < margin[:, None, None] + 25) & (dl8 > 6)
           & (dh8 > 10) & (np.abs(X3) <= XLAT) & (np.abs(X3) >= 30))
    del dh8
    lab[m8b] = L['muscle']
    del dl8
    log(f'step8b anterior costal-slip diaph without lung behind -> muscle: {int(m8b.sum())} voxels ({m8b.sum() * P ** 3 / 1000:.1f} mL)')

    # ---- 9 hepatorenal contact: grow the liver (as s6) into soft/fat between S6 and the upper right kidney
    liver = np.isin(lab, SEG)
    kid = np.isin(lab, [L['rkid'], L['sinus']]) & (X3 < 0)
    dk, kidx = ndi.distance_transform_edt(~kid, return_indices=True); dk = dk * P
    post_ok = G.y[None, :, None] <= G.y[kidx[1]] + 3.0
    del kidx
    dl = ndi.distance_transform_edt(~liver) * P
    Z3 = G.z[None, None, :]
    GAP = 26.0
    dd = ndi.distance_transform_edt(~np.isin(lab, [L['duod'], L['gb'], L['colon']])) * P   # keep off duodenum/GB/colon
    allow = (np.isin(lab, SOFT) & (Z3 >= 1024) & (Z3 <= 1068) & (dk >= 1.6) & (dk <= GAP) & (dl + dk <= GAP + 2) & (X3 < -60)
             & post_ok & (dd > 4.5))
    del dd
    grown = liver.copy()
    for _ in range(40):
        nb_ = ndi.binary_dilation(grown) & allow & ~grown
        if not nb_.any():
            break
        grown |= nb_
    add = ndi.binary_opening(grown & ~liver, iterations=1) & allow
    lab[add] = L['s6']
    log(f'step9 hepatorenal growth: +{int(add.sum())} voxels ({add.sum() * P ** 3 / 1000:.1f} mL) labelled s6')

    # ------------------------------------------------------------------ meshes
    idx = {m['key']: n for n, m in enumerate(A['meshes'])}
    def put(key, v, f):
        mesh = pack(A, key, v, f)
        if key in idx:
            A['meshes'][idx[key]] = mesh
        else:
            idx[key] = len(A['meshes']); A['meshes'].append(mesh)
        log(f'mesh {key}: {len(v)} v / {len(f)} f')

    # celiac -> celiac + ha (face split by nearest label)
    v, f = unpack(A, A['meshes'][idx['celiac']])
    parts = split_faces(v, f, lab, LO, P, [L['celiac'], L['ha']])
    put('celiac', *parts[L['celiac']]); put('ha', *parts[L['ha']])
    # gda from voxels
    put('gda', *surface_nets(lab == L['gda'], LO, P, sigma=.6, step=1, smooth=8))
    # bile: same displacement field as step 4 (continuous in z), then split bile / lhd
    v, f = unpack(A, A['meshes'][idx['bile']])
    fld = field4
    v = v + np.stack([np.interp(v[:, 2], fld[:, 0], fld[:, 1]), np.interp(v[:, 2], fld[:, 0], fld[:, 2]), np.zeros(len(v))], 1)
    parts = split_faces(v, f, lab, LO, P, [L['bile'], L['lhd']])
    put('bile', *parts[L['bile']]); put('lhd', *parts[L['lhd']])
    # splv, smv rebuilt from the edited voxels
    put('splv', *surface_nets(lab == L['splv'], LO, P, sigma=.7, step=1, smooth=10))
    put('smv', *surface_nets(lab == L['smv'], LO, P, sigma=.7, step=1, smooth=10))
    # diaphragm: drop faces whose three vertices all lie in the step-8 deletion region
    v, f = unpack(A, A['meshes'][idx['diaph']])
    mx = np.interp(v[:, 0], G.x, margin)
    inreg = (v[:, 1] < -140) & (v[:, 2] < mx - TOL) & (np.abs(v[:, 0]) <= XLAT)
    keep = ~inreg[f].all(1)
    put('diaph', *compact(v, f[keep]))
    log(f'diaph mesh faces dropped: {int((~keep).sum())}')

    write_js(dst, head, A, lab)
    changed = {}
    for c in np.unique(np.concatenate([np.unique(orig), np.unique(lab)])):
        a, b = int((orig == c).sum()), int((lab == c).sum())
        if a != b:
            changed[[n for n, cc in L.items() if cc == c][0]] = (a, b)
    log('label voxel counts (before, after): ' + json.dumps(changed))
    log(f'voxels changed: {int((orig != lab).sum())}')


if __name__ == '__main__':
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2])

