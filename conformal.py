"""
conformal.py - conforming RF coil winding designer (phantom -> shell -> field of view -> coil).

This is coil_conformal_phantom_v2.py reorganised so the website can call it:
  * settings come in as a dictionary (the form on the website)
  * the phantom STL is the file the visitor uploads (ASCII or binary); no 'trimesh' needed
  * results go out as numbers, tables, text instructions and PNG graphs
  * NO former STL, wire-path CSV or other export (instructions for building the former are given as text)
  * the integer program uses scipy's MILP solver instead of PuLP (PuLP cannot run in a browser)
"""
import io, re, json, base64, traceback
import numpy as np
from scipy.constants import mu_0
from scipy.ndimage import uniform_filter
from scipy.spatial.distance import cdist
from scipy.optimize import milp, LinearConstraint, Bounds
import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt

C_BLUE, C_ORANGE, C_VERM, C_GREEN, C_BLACK, C_SKY = "#0072B2", "#E69F00", "#D55E00", "#009E73", "#000000", "#56B4E9"

# Numerics (same as the original script)
N_THETA = 180      # points around each loop
Z_STEP = 1e-3      # slice spacing for the phantom radius map (m)
SLAB = 3e-3        # half-thickness of each slice (m)
FOV_GRID = (7, 7, 9)
FOV_INSIDE_FRACTION = 0.9
TH = np.linspace(0, 2 * np.pi, N_THETA, endpoint=False)


def png(fig):
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=110, bbox_inches="tight")
    plt.close(fig)
    return base64.b64encode(buf.getvalue()).decode()


# =============================================================================
# Phantom loading (STL only)
# =============================================================================
def read_stl_triangles(path):
    data = open(path, "rb").read()
    if len(data) >= 84:
        n = int(np.frombuffer(data[80:84], dtype="<u4")[0])
        if len(data) == 84 + 50 * n and n > 0:       # binary STL
            rec = np.frombuffer(data, dtype=[("n", "<f4", 3), ("v", "<f4", (3, 3)), ("a", "<u2")], count=n, offset=84)
            return rec["v"].astype(float)
    text = data.decode("utf-8", errors="ignore")      # ASCII STL
    v = re.findall(r"vertex\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)", text)
    if len(v) < 3 or len(v) % 3:
        raise ValueError("Could not read that file as an STL. Export it again as a binary or ASCII STL.")
    return np.array(v, dtype=float).reshape(-1, 3, 3)


def sample_surface(tri, n=300000):
    cross = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    area = 0.5 * np.linalg.norm(cross, axis=1)
    if area.sum() <= 0:
        raise ValueError("The STL has no surface area (it looks empty or flat).")
    rng = np.random.default_rng(0)
    idx = rng.choice(len(tri), size=n, p=area / area.sum())
    u, v = rng.random(n), rng.random(n)
    r1 = np.sqrt(u)
    a, b, c = 1 - r1, r1 * (1 - v), r1 * v
    t = tri[idx]
    return a[:, None] * t[:, 0] + b[:, None] * t[:, 1] + c[:, None] * t[:, 2]


def synthetic_head(n=300000, semi=(0.075, 0.095, 0.115)):
    v = np.random.default_rng(0).normal(size=(n, 3))
    v /= np.linalg.norm(v, axis=1)[:, None]
    return v * np.array(semi)


def load_phantom(path, units_to_m, long_axis):
    p = sample_surface(read_stl_triangles(path)) * units_to_m
    ax = long_axis
    if ax == "auto":
        ax = "xyz"[int(np.argmax(np.ptp(p, axis=0)))]
    return p[:, {"z": [0, 1, 2], "y": [0, 2, 1], "x": [1, 2, 0]}[ax]], ax


# =============================================================================
# Geometry helpers
# =============================================================================
def get_groove_positions(length, M, spacing):
    if spacing == "uniform":
        return np.linspace(-length / 2, length / 2, M)
    if spacing == "sqrt":
        t = np.linspace(-1, 1, M)
        return np.sign(t) * (length / 2) * np.sqrt(np.abs(t))
    if spacing == "cosine":
        return np.sort((length / 2) * np.cos(np.pi * np.arange(M) / (M - 1)))
    raise ValueError(f"Unknown groove spacing: {spacing}")


def radial_profile(points, z_lo, z_hi, slide, strict=True):
    """Outer radius of the phantom per z slice and angle, plus the slide-on envelope."""
    ZF = np.arange(points[:, 2].min(), points[:, 2].max() + Z_STEP / 2, Z_STEP)
    ang = np.arctan2(points[:, 1], points[:, 0]) % (2 * np.pi)
    rad = np.hypot(points[:, 0], points[:, 1])
    kb = (ang / (2 * np.pi / N_THETA)).astype(int) % N_THETA
    order = np.argsort(points[:, 2])
    zs, kb, rad = points[order, 2], kb[order], rad[order]
    R = np.zeros((len(ZF), N_THETA))
    for i, zg in enumerate(ZF):
        a, b = np.searchsorted(zs, [zg - SLAB, zg + SLAB])
        if b > a:
            np.maximum.at(R[i], kb[a:b], rad[a:b])
    idx = np.arange(N_THETA)
    for i in range(len(R)):
        ok = R[i] > 0
        if ok.sum() >= N_THETA // 4 and not ok.all():
            R[i] = np.interp(idx, idx[ok], R[i][ok], period=N_THETA)
    R = uniform_filter(R, size=(5, 5), mode=("nearest", "wrap"))
    env = R.copy()
    if slide == "top":
        env = np.maximum.accumulate(R[::-1], axis=0)[::-1]
    elif slide == "bottom":
        env = np.maximum.accumulate(R, axis=0)
    keep = (ZF >= z_lo - 1e-9) & (ZF <= z_hi + 1e-9)
    R, env = R[keep], env[keep]
    for i in range(len(R)):
        if (R[i] > 0).sum() < N_THETA // 4:
            raise ValueError(f"The phantom has too little surface at z = {ZF[keep][i]*1000:.0f} mm. "
                             "Move the field of view or shrink the coil margin so the coil stays on the head.")
    return ZF[keep], R, env


def loop_path(z, r):
    return np.c_[r * np.cos(TH), r * np.sin(TH), np.full(N_THETA, z)]


def bz_from_loop(path, P):
    nxt = np.roll(path, -1, axis=0)
    mid, dl = (path + nxt) / 2, nxt - path
    r = P[:, None, :] - mid[None, :, :]
    d3 = np.linalg.norm(r, axis=2) ** 3
    cz = dl[None, :, 0] * r[:, :, 1] - dl[None, :, 1] * r[:, :, 0]
    return mu_0 / (4 * np.pi) * (cz / d3).sum(axis=1)


def optimize(B, target, n_max, n_min, allow_bucking, time_limit):
    """Minimise the worst relative error over the field of view (scipy MILP instead of PuLP)."""
    npts, mc = B.shape
    s = 1e6
    Bs, ts = B * s, target * s
    c = np.r_[np.zeros(mc), 1.0]
    A = np.vstack([np.c_[Bs, -ts], np.c_[-Bs, -ts]])
    lb_n = n_min if allow_bucking else 0
    res = milp(c, constraints=LinearConstraint(A, -np.inf, np.r_[ts, -ts]),
               integrality=np.r_[np.ones(mc), 0],
               bounds=Bounds(np.r_[np.full(mc, lb_n), 0], np.r_[np.full(mc, n_max), np.inf]),
               options={"time_limit": time_limit, "mip_rel_gap": 0.02})
    if res.x is None:
        raise ValueError("The optimiser could not find a winding pattern. Try more grooves or a larger turn limit.")
    return np.round(res.x[:mc]).astype(int), res.status == 1


# =============================================================================
# Main entry point called by the website
# =============================================================================
def _run(P, stl_path):
    mm = lambda k: float(P[k]) * 1e-3
    opt = lambda k: None if P.get(k) in (None, "", "auto") else float(P[k]) * 1e-3

    STANDOFF, former_thickness = mm("standoff_mm"), mm("former_thickness_mm")
    GROOVE_DEPTH, GROOVE_CLEARANCE, FORMER_END_MARGIN = mm("groove_depth_mm"), mm("groove_clearance_mm"), mm("end_margin_mm")
    wire_bare, wire_ins = mm("wire_bare_mm"), mm("wire_ins_mm")
    M, n_max, n_min = int(P["M"]), int(P["n_max"]), int(P["n_min"])
    ALLOW_BUCKING = bool(P["bucking"])
    spacing = P["spacing"]
    TARGET_FRACTION = float(P["target_fraction"])
    f_op = float(P["freq_mhz"]) * 1e6
    SLIDE_FROM = None if P["slide_from"] == "none" else P["slide_from"]
    FOV_MODE = P["fov_mode"]
    FOV_CENTER = (mm("fov_cx_mm"), mm("fov_cy_mm"), mm("fov_cz_mm"))
    fx, fy, fz = opt("fov_x_mm"), opt("fov_y_mm"), opt("fov_z_mm")
    ROI_SHRINK, ROI_SPACING = float(P["roi_shrink"]), mm("roi_spacing_mm")
    COIL_MARGIN = opt("coil_margin_mm")
    time_limit = max(float(P.get("solver_s", 20)), 2)

    if M < 2:
        raise ValueError("Use at least 2 grooves.")
    if n_max < 1:
        raise ValueError("Max forward windings per groove must be at least 1.")
    if ALLOW_BUCKING and n_min > 0:
        n_min = 0
    if GROOVE_DEPTH >= former_thickness:
        raise ValueError("Groove depth must be smaller than the former thickness.")
    if wire_ins < wire_bare:
        raise ValueError("The insulated wire diameter must be at least the bare diameter.")
    if FOV_MODE not in ("anatomy", "box"):
        raise ValueError("Field-of-view mode must be 'anatomy' or 'box'.")

    # --- phantom
    notes = []
    if stl_path:
        pts, axis_used = load_phantom(stl_path, float(P["units_to_m"]), P["axis"])
        source = f"uploaded STL (file's {axis_used.upper()} axis used as the coil axis)"
    else:
        pts, axis_used = synthetic_head(), "z"
        source = "built-in ellipsoid test head (no STL uploaded)"
    lo, hi = pts.min(0), pts.max(0)
    pts = pts - (lo + hi) / 2
    ext = hi - lo
    size_txt = f"{ext[0]*1e3:.0f} x {ext[1]*1e3:.0f} x {ext[2]*1e3:.0f} mm"
    if ext.max() > 1.0 or ext.max() < 0.02:
        raise ValueError(f"The phantom measures {size_txt}, which does not look like a head. Check the file units setting.")

    # --- coil extent and shell
    FOV_SIZE = (fx or 0.6 * ext[0], fy or 0.6 * ext[1], fz or 0.4 * ext[2])
    if COIL_MARGIN is None:
        COIL_MARGIN = 0.1 * ext[2]
    coil_length = FOV_SIZE[2] + 2 * COIL_MARGIN
    z_loops = FOV_CENTER[2] + get_groove_positions(coil_length, M, spacing)
    z_lo, z_hi = z_loops.min() - FORMER_END_MARGIN, z_loops.max() + FORMER_END_MARGIN
    if z_lo < pts[:, 2].min() or z_hi > pts[:, 2].max():
        raise ValueError("The coil is longer than the phantom along the coil axis. Check the long-axis and units settings, "
                         "or make the field-of-view height, coil margin or former end margin smaller.")
    ZG, R_ph, R_env = radial_profile(pts, z_lo, z_hi, SLIDE_FROM)

    def _interp(R, z):
        i = np.clip(np.searchsorted(ZG, z) - 1, 0, len(ZG) - 2)
        t = np.clip((z - ZG[i]) / (ZG[i + 1] - ZG[i]), 0, 1)
        return (1 - t) * R[i] + t * R[i + 1]

    radial_at = lambda z: _interp(R_env, z)
    phantom_at = lambda z: _interp(R_ph, z)
    r_inner = lambda z: radial_at(z) + STANDOFF
    r_outer = lambda z: r_inner(z) + former_thickness
    r_wire = lambda z: r_outer(z) + wire_ins / 2

    R_cyl = max(r_wire(z).max() for z in z_loops)
    if min(phantom_at(z).min() for z in z_loops) < 0.3 * R_ph.max():
        notes.append("The coil reaches a narrow part of the phantom (under 30% of its widest radius), where the radial offset is less "
                     "accurate. Consider moving the field of view.")
    gaps = np.array([r_inner(z) - phantom_at(z) for z in z_loops]) * 1000

    # --- field of view points
    def inside_phantom(Pt, frac):
        r = np.hypot(Pt[:, 0], Pt[:, 1])
        k = np.round((np.arctan2(Pt[:, 1], Pt[:, 0]) % (2 * np.pi)) / (2 * np.pi / N_THETA)).astype(int) % N_THETA
        i = np.clip(np.round((Pt[:, 2] - ZG[0]) / Z_STEP).astype(int), 0, len(ZG) - 1)
        return r < frac * R_ph[i, k]

    zc0 = FOV_CENTER[2]
    if FOV_MODE == "box":
        g = [np.linspace(-sz / 2, sz / 2, n) + c for sz, n, c in zip(FOV_SIZE, FOV_GRID, FOV_CENTER)]
        F = np.array(np.meshgrid(*g, indexing="ij")).reshape(3, -1).T
        F = F[inside_phantom(F, FOV_INSIDE_FRACTION)]
    else:
        xr = R_ph.max()
        gx = np.arange(-xr, xr + ROI_SPACING / 2, ROI_SPACING)
        gz = np.arange(zc0 - FOV_SIZE[2] / 2, zc0 + FOV_SIZE[2] / 2 + 1e-9, ROI_SPACING)
        F = np.array(np.meshgrid(gx, gx, gz, indexing="ij")).reshape(3, -1).T
        F = F[inside_phantom(F, ROI_SHRINK)]
        if len(F) > 800:
            F = F[::int(np.ceil(len(F) / 800))]
    if len(F) < 20:
        raise ValueError("Fewer than 20 sample points in the imaging region. Check the field-of-view size and centre, "
                         "or reduce the sample spacing.")

    # --- fields, optimisation
    paths_con = [loop_path(z, r_wire(z)) for z in z_loops]
    paths_cyl = [loop_path(z, np.full(N_THETA, R_cyl)) for z in z_loops]
    B_con = np.column_stack([bz_from_loop(p, F) for p in paths_con])
    B_cyl = np.column_stack([bz_from_loop(p, F) for p in paths_cyl])
    target = np.full(len(F), TARGET_FRACTION * B_cyl.sum(axis=1).mean())
    n_con, hit1 = optimize(B_con, target, n_max, n_min, ALLOW_BUCKING, time_limit)
    n_cyl, hit2 = optimize(B_cyl, target, n_max, n_min, ALLOW_BUCKING, time_limit)
    if hit1 or hit2:
        notes.append(f"The optimiser stopped at its {time_limit:.0f} s time limit for at least one coil and kept the best pattern found.")

    # --- comparison metrics
    def inductance(paths, n):
        segs = [((np.roll(p, -1, 0) + p) / 2, np.roll(p, -1, 0) - p) for p in paths]
        Lm = np.zeros((len(paths), len(paths)))
        for i, (mi, di) in enumerate(segs):
            Req = np.linalg.norm(di, axis=1).sum() / (2 * np.pi)
            Lm[i, i] = mu_0 * Req * (np.log(8 * Req / (wire_bare / 2)) - 2)
            for j in range(i + 1, len(paths)):
                mj, dj = segs[j]
                Lm[i, j] = Lm[j, i] = mu_0 / (4 * np.pi) * np.sum((di @ dj.T) / cdist(mi, mj))
        return float(n @ Lm @ n)

    def metrics(paths, B, n):
        b = B @ n
        per = [np.linalg.norm(np.roll(p, -1, 0) - p, axis=1).sum() for p in paths]
        L = inductance(paths, n)
        omega = 2 * np.pi * f_op
        return {"err": float(np.max(np.abs(b - target) / target) * 100), "b1": float(b.mean() * 1e6),
                "turns": int(np.abs(n).sum()), "wire": float(np.sum(np.abs(n) * per)), "L": L * 1e6,
                "C": float(1e12 / (omega ** 2 * L)) if L > 0 else None}

    mc, my = metrics(paths_con, B_con, n_con), metrics(paths_cyl, B_cyl, n_cyl)
    fmt = lambda v, f: f.format(v) if v is not None else "n/a"
    comparison = [
        ["Max field error (%)", fmt(mc["err"], "{:.2f}"), fmt(my["err"], "{:.2f}")],
        ["Mean B1 in region (uT/A)", fmt(mc["b1"], "{:.3f}"), fmt(my["b1"], "{:.3f}")],
        ["Total turns", mc["turns"], my["turns"]],
        ["Wire length (m)", fmt(mc["wire"], "{:.1f}"), fmt(my["wire"], "{:.1f}")],
        ["Inductance (uH)", fmt(mc["L"], "{:.2f}"), fmt(my["L"], "{:.2f}")],
        ["Tuning capacitance (pF)", fmt(mc["C"], "{:.0f}"), fmt(my["C"], "{:.0f}")],
    ]

    # --- winding table and former instructions
    gw = wire_ins + GROOVE_CLEARANCE
    used = [(m + 1, z, int(n)) for m, (z, n) in enumerate(zip(z_loops, n_con)) if n != 0]
    for (_, a, _), (_, b, _) in zip(used[:-1], used[1:]):
        if (b - a) - gw < 0.5e-3:
            notes.append(f"Grooves at {a*1e3:+.1f} and {b*1e3:+.1f} mm leave less than 0.5 mm of wall between them.")
    windings = [{"groove": g, "from_end_mm": round((z - z_lo) * 1e3, 2), "conformal": int(n_con[g - 1]),
                 "cylinder": int(n_cyl[g - 1])} for g, z, _ in used]
    cyl_only = [{"groove": m + 1, "from_end_mm": round((z - z_lo) * 1e3, 2), "conformal": 0, "cylinder": int(n_cyl[m])}
                for m, z in enumerate(z_loops) if n_con[m] == 0 and n_cyl[m] != 0]
    shell_in = [r_inner(z) for z in z_loops]
    n_used = len(used)
    inst = [
        f"Place the phantom with its {axis_used.upper()} axis (the long axis of the head) as the coil axis. Height z = 0 is the middle of the phantom's bounding box.",
        f"Make the shell: offset the phantom surface outward by {STANDOFF*1e3:.1f} mm for the inner surface, then thicken outward by "
        f"{former_thickness*1e3:.1f} mm. Along the coil the inner radius runs from {min(r.min() for r in shell_in)*1e3:.1f} to "
        f"{max(r.max() for r in shell_in)*1e3:.1f} mm.",
    ]
    if SLIDE_FROM:
        inst.append(f"So the head can slide in, the shell must be at least as wide as every part of the head that passes through it. "
                    f"The coil slides on from the {SLIDE_FROM} of the head. On average the gap to the head is {gaps.mean():.1f} mm "
                    f"(largest {gaps.max():.1f} mm). The graph 'Former profile' shows the shape to follow.")
    else:
        inst.append("The shell follows the head exactly, so it has undercuts. Build it as a split (clamshell) former.")
    inst += [
        f"Trim the shell to run from z = {z_lo*1e3:+.1f} mm to z = {z_hi*1e3:+.1f} mm ({(z_hi-z_lo)*1e3:.1f} mm long).",
        f"Cut {n_used} grooves around the outside of the shell, each {gw*1e3:.2f} mm wide and {GROOVE_DEPTH*1e3:.2f} mm deep, "
        "centred at the 'from former end' distances in the table. Each groove follows the outline of the shell at that height (it is not a circle).",
        f"Wind the wire ({wire_bare*1e3:.3f} mm bare, {wire_ins*1e3:.3f} mm insulated) along each groove using the turns in the table. "
        "Extra turns sit on top of the first one." + (" Bucking (negative) turns run the opposite way." if ALLOW_BUCKING else ""),
        "Connect all turns in series in the order of the table.",
    ]
    if mc["C"]:
        inst.append(f"To resonate at {f_op/1e6:.3f} MHz, start with about {mc['C']:.0f} pF of tuning capacitance and fine-tune on the bench.")

    # --- figures
    figs = []
    zc = zc0
    fig = plt.figure(figsize=(16, 5.5))
    ax = fig.add_subplot(1, 3, 1, projection="3d")
    ax.scatter(*(pts[::60][:5000] * 1e3).T, s=1, c="#BBBBBB", alpha=.4)
    for p, n in zip(paths_con, n_con):
        if n:
            q = np.vstack([p, p[:1]]) * 1e3
            ax.plot(q[:, 0], q[:, 1], q[:, 2], color=C_BLUE if n > 0 else C_ORANGE, lw=1.2)
    lim = max(R_cyl, np.abs(pts).max()) * 1e3
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_zlim(-lim, lim)
    ax.set_xlabel("x (mm)"); ax.set_ylabel("y (mm)"); ax.set_zlabel("z (mm)")
    ax.set_title("Conformal loops around the phantom")
    ax = fig.add_subplot(1, 3, 2)

    def contour(r):
        return np.r_[r * np.cos(TH), r[0] * np.cos(TH[0])] * 1e3, np.r_[r * np.sin(TH), r[0] * np.sin(TH[0])] * 1e3

    ax.plot(*contour(phantom_at(zc)), color=C_BLACK, lw=2, label="Phantom")
    ax.plot(*contour(r_outer(zc)), color=C_ORANGE, lw=2, label="Shell outer surface")
    ax.plot(*contour(r_wire(zc)), color=C_BLUE, lw=1.5, ls="-.", label="Conformal wire")
    ax.plot(*contour(np.full(N_THETA, R_cyl)), color=C_VERM, lw=2, ls="--", label="Cylinder wire")
    zlev = np.unique(F[:, 2]); zlev = zlev[np.argmin(np.abs(zlev - zc))]
    sl = F[np.isclose(F[:, 2], zlev)]
    ax.scatter(sl[:, 0] * 1e3, sl[:, 1] * 1e3, s=14, c=C_GREEN, marker="s", label="Imaging region")
    ax.set_aspect("equal"); ax.set_xlabel("x (mm)"); ax.set_ylabel("y (mm)")
    ax.set_title(f"Axial slice at z = {zc*1e3:.0f} mm"); ax.legend(fontsize=8, loc="upper right")
    ax = fig.add_subplot(1, 3, 3)
    bw = 0.35 * np.min(np.diff(np.sort(z_loops))) * 1e3
    ax.bar(z_loops * 1e3 - bw / 2, n_con, width=bw, color=C_BLUE, edgecolor="black", lw=.5, label="Conformal")
    ax.bar(z_loops * 1e3 + bw / 2, n_cyl, width=bw, color=C_VERM, edgecolor="black", lw=.5, hatch="//", label="Cylinder")
    ax.axhline(0, color="black", lw=.8)
    ax.set_xlabel("Groove position z (mm)"); ax.set_ylabel("Turns"); ax.set_title("Winding pattern")
    ax.legend(); ax.grid(alpha=.3, axis="y")
    fig.tight_layout()
    figs.append({"title": "Overview", "caption": "Left: loops around the head. Middle: one slice through phantom, shell and wire. Right: turns in each groove.",
                 "png": png(fig)})

    # former profile (the graph that explains how to build the shell)
    zs = np.linspace(z_lo, z_hi, 120)
    ph = np.array([phantom_at(z) for z in zs]); ro = np.array([r_outer(z) for z in zs]); ri = np.array([r_inner(z) for z in zs])
    fig, axs = plt.subplots(1, 2, figsize=(11, 6), sharey=True)
    for a, (name, kr, kl) in zip(axs, [("x-z section", 0, N_THETA // 2), ("y-z section", N_THETA // 4, 3 * N_THETA // 4)]):
        for arr, col, lab, ls in [(ph, C_BLACK, "Phantom", "-"), (ri, C_SKY, "Shell inner", "--"), (ro, C_ORANGE, "Shell outer", "-")]:
            a.plot(arr[:, kr] * 1e3, zs * 1e3, color=col, ls=ls, lw=2, label=lab)
            a.plot(-arr[:, kl] * 1e3, zs * 1e3, color=col, ls=ls, lw=2)
        for _, z, _ in used:
            a.axhline(z * 1e3, color=C_BLUE, lw=.6, alpha=.6)
        a.set_aspect("equal"); a.set_xlabel("radius (mm)"); a.set_title(f"Former profile, {name}"); a.grid(alpha=.2)
    axs[0].set_ylabel("z (mm)"); axs[0].legend(fontsize=8, loc="lower center")
    fig.tight_layout()
    figs.append({"title": "Former profile", "caption": "Phantom, shell inner and outer surfaces. Blue lines mark the groove heights.",
                 "png": png(fig)})

    # field error maps
    xr2 = max(R_cyl, R_ph.max()) * 1.05
    gx2 = np.linspace(-xr2, xr2, 61)
    gzm = np.linspace(zc - coil_length / 2, zc + coil_length / 2, 61)
    GX, GZ = np.meshgrid(gx2, gzm)
    G = np.c_[GX.ravel(), np.zeros(GX.size), GZ.ravel()]
    mask = ~inside_phantom(G, 0.95)

    def dev_map(paths, n):
        Bm = sum(int(nn) * bz_from_loop(p, G) for p, nn in zip(paths, n) if nn)
        d = np.abs(Bm - target[0]) / target[0] * 100
        d[mask] = np.nan
        return np.clip(d.reshape(GX.shape), 0, 40)

    near = F[np.abs(F[:, 1]) <= np.abs(F[:, 1]).min() + 1e-6]
    fig, axs = plt.subplots(1, 2, figsize=(12, 5.5), sharey=True)
    for a, (name, paths, n) in zip(axs, [("Conformal", paths_con, n_con), ("Cylinder", paths_cyl, n_cyl)]):
        im = a.imshow(dev_map(paths, n), extent=[gx2[0] * 1e3, gx2[-1] * 1e3, gzm[0] * 1e3, gzm[-1] * 1e3],
                      origin="lower", cmap="cividis", vmin=0, vmax=40, aspect="equal")
        a.scatter(near[:, 0] * 1e3, near[:, 2] * 1e3, s=8, c="white", marker="s", edgecolors="black", linewidths=.4)
        a.set_title(f"{name}: field error inside the head (y = 0)"); a.set_xlabel("x (mm)")
    axs[0].set_ylabel("z (mm)")
    fig.colorbar(im, ax=axs, label="Deviation from target (%), capped at 40", shrink=.85)
    figs.append({"title": "Field error maps", "caption": "How far the field is from uniform. White squares are the sample points used for the design.",
                 "png": png(fig)})

    summary = [
        ["Phantom", source], ["Phantom size (x, y, z)", size_txt],
        ["Coil length", f"{coil_length*1e3:.1f} mm"], ["Imaging region", f"{FOV_MODE}, {len(F)} sample points"],
        ["Mean gap, head to shell", f"{gaps.mean():.1f} mm (largest {gaps.max():.1f} mm)"],
        ["Baseline cylinder radius", f"{R_cyl*1e3:.1f} mm"],
        ["Field error, conformal vs cylinder", f"{mc['err']:.2f} % vs {my['err']:.2f} %"],
    ]
    return {"summary": summary, "comparison": comparison, "windings": windings + cyl_only,
            "instructions": inst, "notes": notes, "figures": figs}


def run(params_json, stl_path=""):
    try:
        return json.dumps(_run(json.loads(params_json), stl_path))
    except Exception as e:
        traceback.print_exc()
        return json.dumps({"error": str(e)})
