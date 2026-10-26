"""
solenoid.py - regular (cylindrical) RF coil winding designer.

This is coil_optimisation_v1_5.py reorganised so the website can call it:
  * all coil settings come in as a dictionary (the form on the website)
  * results go out as numbers, tables, text instructions and PNG graphs
  * NO former / STL / Fusion 360 export (instructions for building the former are given as text)
  * the integer program uses scipy's MILP solver instead of PuLP (PuLP cannot run in a browser)
"""
import io, json, base64, traceback
import numpy as np
from scipy.constants import mu_0
from scipy.special import ellipk, ellipe
from scipy.optimize import milp, LinearConstraint, Bounds
import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Rectangle

# Okabe-Ito colours (safe for most types of colour blindness)
BLUE, VERM, ORANGE, SKY, GREEN = "#0072B2", "#D55E00", "#E69F00", "#56B4E9", "#009E73"
DARK_BLUE, DARK_VERM = "#003F63", "#7A3300"


def png(fig):
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=110, bbox_inches="tight")
    plt.close(fig)
    return base64.b64encode(buf.getvalue()).decode()


# =============================================================================
# Groove positions and target field
# =============================================================================
def get_groove_positions(coil_length, M, spacing):
    if spacing == "uniform":
        return np.linspace(-coil_length / 2, coil_length / 2, M)
    if spacing == "sqrt":
        t = np.linspace(-1, 1, M)
        return np.sign(t) * (coil_length / 2) * np.sqrt(np.abs(t))
    if spacing == "cosine":
        k = np.arange(M)
        return np.sort((coil_length / 2) * np.cos(np.pi * k / (M - 1)))
    raise ValueError(f"Unknown groove spacing: {spacing}")


def get_target_field(x_eval, profile, sigma=0.25):
    x_norm = (x_eval - x_eval.min()) / (x_eval.max() - x_eval.min())
    profiles = {
        "sqrt": (np.sqrt(x_norm), "sqrt(x) (Bloch-Siegert linear gradient)"),
        "uniform": (np.ones_like(x_eval), "Uniform (flat)"),
        "linear": (x_norm, "Linear gradient"),
        "quadratic": (x_norm ** 2, "Quadratic (x^2)"),
        "cosine": (0.5 * (1 + np.cos(np.pi * (2 * x_norm - 1))), "Cosine bell"),
        "gaussian": (np.exp(-0.5 * ((x_norm - 0.5) / sigma) ** 2), f"Gaussian (sigma={sigma})"),
        "linear_symmetric": (np.abs(2 * x_norm - 1), "Linear symmetric |x|"),
    }
    if profile not in profiles:
        raise ValueError(f"Unknown target profile: {profile}")
    B, desc = profiles[profile]
    return np.maximum(B, 0), desc


# =============================================================================
# Biot-Savart
# =============================================================================
def b1_single_loop_on_axis(loop_x, loop_radius, eval_x):
    dx = eval_x - loop_x
    return mu_0 * loop_radius ** 2 / (2.0 * (loop_radius ** 2 + dx ** 2) ** 1.5)


def b1_field_map(groove_positions, n_windings, coil_radius, former_thickness, wire_ins, coil_length, res):
    """2D |B1| map in the x-y plane. Same method as the original script, but vectorised row by row."""
    fov_x, fov_y = coil_length * 0.85, coil_radius * 1.5
    x_map = np.arange(-fov_x / 2, fov_x / 2 + res, res)
    y_map = np.arange(-fov_y / 2, fov_y / 2 + res, res)
    B1 = np.zeros((len(y_map), len(x_map)))
    n_seg = 360
    theta = np.linspace(0, 2 * np.pi, n_seg, endpoint=False)
    dtheta = 2 * np.pi / n_seg
    for pos, n_w in zip(groove_positions, n_windings):
        if n_w == 0:
            continue
        for layer in range(abs(int(n_w))):
            R = coil_radius + former_thickness + wire_ins * (layer + 0.5)
            w_y, w_z = R * np.cos(theta), R * np.sin(theta)
            dl_y, dl_z = -R * np.sin(theta) * dtheta, R * np.cos(theta) * dtheta
            rx = (x_map - pos)[:, None]
            rz = -w_z[None, :]
            sign = np.sign(n_w)
            for iy, yp in enumerate(y_map):
                ry = yp - w_y[None, :]
                r3 = np.sqrt(rx ** 2 + ry ** 2 + rz ** 2) ** 3
                Bx = np.sum((dl_y[None, :] * rz - dl_z[None, :] * ry) / r3, axis=1)
                By = np.sum((dl_z[None, :] * rx) / r3, axis=1)
                Bz = np.sum((-dl_y[None, :] * rx) / r3, axis=1)
                B1[iy] += sign * mu_0 / (4 * np.pi) * np.sqrt(Bx ** 2 + By ** 2 + Bz ** 2)
    return x_map, y_map, B1


# =============================================================================
# Integer optimisation (scipy MILP instead of PuLP)
# =============================================================================
def optimize_windings(B, target, n_max, n_min, allow_bucking, time_limit=20):
    """Minimise the worst relative error: |B n - target| <= eps * target at every point."""
    npts, mc = B.shape
    s = 1e6                                   # work in microtesla so the solver sees sensible numbers
    Bs, ts = B * s, target * s
    c = np.r_[np.zeros(mc), 1.0]
    A = np.vstack([np.c_[Bs, -ts], np.c_[-Bs, -ts]])
    ub = np.r_[ts, -ts]
    lb_n = n_min if allow_bucking else 0
    res = milp(c, constraints=LinearConstraint(A, -np.inf, ub),
               integrality=np.r_[np.ones(mc), 0],
               bounds=Bounds(np.r_[np.full(mc, lb_n), 0], np.r_[np.full(mc, n_max), np.inf]),
               options={"time_limit": time_limit, "mip_rel_gap": 0.02})
    if res.x is None:
        raise ValueError("The optimiser could not find a winding pattern. Try more grooves or a larger turn limit.")
    n = np.round(res.x[:mc]).astype(int)
    err = float(np.max(np.abs(B @ n - target) / target))
    return n, err, res.status == 1


# =============================================================================
# Inductance
# =============================================================================
def self_inductance_loop(R, a_wire):
    return mu_0 * R * (np.log(8.0 * R / a_wire) - 2.0)


def mutual_inductance_coaxial(R1, R2, d):
    k_sq = min(4.0 * R1 * R2 / ((R1 + R2) ** 2 + d ** 2), 0.999999)
    k = np.sqrt(k_sq)
    return mu_0 * np.sqrt(R1 * R2) * ((2.0 / k - k) * ellipk(k_sq) - (2.0 / k) * ellipe(k_sq))


def calculate_inductance(n_windings, groove_positions, coil_radius, wire_radius, wire_ins, former_thickness):
    turns = []
    for m, n_w in enumerate(n_windings):
        if n_w == 0:
            continue
        for t in range(abs(int(n_w))):
            turns.append({"x": groove_positions[m],
                          "R": coil_radius + former_thickness + wire_ins * (t + 0.5),
                          "sign": 1 if n_w > 0 else -1})
    nt = len(turns)
    if nt == 0:
        return 0.0, 0.0, 0.0, None, []
    L = np.zeros((nt, nt))
    for i in range(nt):
        L[i, i] = self_inductance_loop(turns[i]["R"], wire_radius)
        for j in range(i + 1, nt):
            d = abs(turns[i]["x"] - turns[j]["x"])
            R1, R2 = turns[i]["R"], turns[j]["R"]
            if d < 1e-10 and abs(R1 - R2) < 1e-10:
                Mij = L[i, i] * 0.98
            elif d < 1e-10:
                Mij = mutual_inductance_coaxial(R1, R2, wire_ins * 0.1)
            else:
                Mij = mutual_inductance_coaxial(R1, R2, d)
            L[i, j] = L[j, i] = Mij
    s = np.array([t["sign"] for t in turns], dtype=float)
    total = float(s @ L @ s)
    self_total = float(np.trace(L))
    return total, self_total, total - self_total, L, turns


# =============================================================================
# Main entry point called by the website
# =============================================================================
def _run(P):
    mm = lambda k: float(P[k]) * 1e-3
    coil_length, coil_radius, former_thickness = mm("coil_length_mm"), mm("coil_radius_mm"), mm("former_thickness_mm")
    wire_bare, wire_ins = mm("wire_bare_mm"), mm("wire_ins_mm")
    margin, clearance, groove_depth = mm("end_margin_mm"), mm("groove_clearance_mm"), mm("groove_depth_mm")
    res_map = mm("map_res_mm")
    M, N = int(P["M"]), int(P["N"])
    n_max, n_min = int(P["n_max"]), int(P["n_min"])
    allow_bucking = bool(P["bucking"])
    spacing, profile = P["spacing"], P["profile"]
    sigma = float(P.get("sigma", 0.25))
    f_op = float(P["freq_mhz"]) * 1e6
    time_limit = max(float(P.get("solver_s", 20)), 2)

    if M < 2 or N < 5:
        raise ValueError("Use at least 2 grooves and 5 field evaluation points.")
    if n_max < 1:
        raise ValueError("Max forward windings per groove must be at least 1.")
    if allow_bucking and n_min > 0:
        n_min = 0
    if groove_depth >= former_thickness:
        raise ValueError("Groove depth must be smaller than the former thickness.")
    if min(coil_length, coil_radius, former_thickness, wire_bare, wire_ins, res_map) <= 0 or wire_ins < wire_bare:
        raise ValueError("Lengths must be positive, and the insulated wire diameter must be at least the bare diameter.")

    wire_radius = wire_bare / 2
    groove_positions = get_groove_positions(coil_length, M, spacing)
    x_eval = np.linspace(-coil_length / 2 * 0.8, coil_length / 2 * 0.8, N)
    B_target_norm, profile_desc = get_target_field(x_eval, profile, sigma)

    # Field matrix (on axis) and target scaling, as in the original script
    R_wire = coil_radius + former_thickness + wire_ins * 0.5
    B_matrix = np.column_stack([b1_single_loop_on_axis(p, R_wire, x_eval) for p in groove_positions])
    if profile == "uniform":
        scale = np.mean(B_matrix.sum(axis=1)) * 0.5
    else:
        nz = B_target_norm > 0.01
        scale = (n_max * np.mean(B_matrix.sum(axis=1))) / np.mean(B_target_norm[nz]) if np.any(nz) else 1.0
        scale *= 0.5
    B_target = B_target_norm * scale
    valid = B_target > 1e-12
    n_windings, max_error, hit_limit = optimize_windings(B_matrix[valid], B_target[valid], n_max, n_min, allow_bucking, time_limit)
    n_fwd, n_buck = np.maximum(n_windings, 0), np.minimum(n_windings, 0)

    # Inductance and RF numbers
    L_total, L_self, L_mut, L_mat, turns = calculate_inductance(
        n_windings, groove_positions, coil_radius, wire_radius, wire_ins, former_thickness)
    omega = 2 * np.pi * f_op
    rho = 1.68e-8
    skin = float(np.sqrt(2 * rho / (omega * mu_0)))
    wire_len = float(sum(2 * np.pi * t["R"] for t in turns))
    R_dc = rho * wire_len / (np.pi * wire_radius ** 2) if wire_len > 0 else 0.0
    R_ac = rho * wire_len / (2 * np.pi * wire_radius * skin) if wire_radius > skin else R_dc
    Q = omega * L_total / R_ac if R_ac > 0 else 0.0
    C_res = 1.0 / (omega ** 2 * L_total) if L_total > 0 else None

    # ---------------------------------------------------------------- former instructions
    groove_w = wire_ins + clearance
    R_in, R_out = coil_radius, coil_radius + former_thickness
    former_len = coil_length + 2 * margin
    used = [(m, groove_positions[m], int(n_windings[m])) for m in range(M) if n_windings[m] != 0]
    notes = []
    if hit_limit:
        notes.append(f"The optimiser stopped at its {time_limit:.0f} s time limit and kept the best pattern found. "
                     "Raise the time limit for a possibly slightly better pattern.")
    xs = [u[1] for u in used]
    for a, b in zip(xs[:-1], xs[1:]):
        wall = (b - a) - groove_w
        if wall < 0.5e-3:
            notes.append(f"Grooves at {a*1e3:+.2f} mm and {b*1e3:+.2f} mm leave only {wall*1e3:.2f} mm of wall between them. "
                         "Use fewer grooves or a longer coil.")
    rows = [{"groove": m + 1,
             "from_end_mm": round((margin + p + coil_length / 2) * 1e3, 2),
             "from_center_mm": round(p * 1e3, 2),
             "turns": n, "direction": "Forward" if n > 0 else "Bucking"} for m, p, n in used]
    instructions = [
        f"Make a cylindrical tube (the former): {2*R_in*1e3:.1f} mm inner diameter, {2*R_out*1e3:.1f} mm outer diameter "
        f"({former_thickness*1e3:.1f} mm wall) and {former_len*1e3:.1f} mm long.",
        f"Cut {len(used)} grooves around the outside of the tube, each {groove_w*1e3:.2f} mm wide and {groove_depth*1e3:.2f} mm deep. "
        "Centre them at the 'from former end' distances in the table. Grooves with no turns are not cut.",
        f"Wind the wire ({wire_bare*1e3:.3f} mm bare, {wire_ins*1e3:.3f} mm insulated) into the grooves, using the number of turns in the table. "
        "Extra turns in a groove sit on top of the first turn, one layer at a time.",
        "Forward turns carry current in one direction. Bucking turns must carry current the opposite way, "
        "so wind them in the opposite sense (or reverse the connection).",
        "Connect all turns in series in the order of the table.",
    ]
    if C_res is not None:
        instructions.append(f"To resonate at {f_op/1e6:.3f} MHz, start with about {C_res*1e12:.0f} pF of tuning capacitance "
                            "and fine-tune on the bench.")

    # ---------------------------------------------------------------- figures
    figs = []
    B_ach = B_matrix @ n_windings

    # 1. field performance
    fig, axs = plt.subplots(3, 1, figsize=(10, 9))
    ax = axs[0]
    ax.plot(x_eval * 1e3, B_target * 1e6, color=BLUE, lw=2, label=f"Target: {profile_desc}")
    ax.plot(x_eval * 1e3, B_ach * 1e6, color=VERM, ls="--", lw=2, label="Achieved")
    ax.set_xlabel("Position (mm)"); ax.set_ylabel("B1 (uT per amp)")
    ax.set_title("Target vs achieved B1 field (on axis)"); ax.legend(); ax.grid(alpha=.3)
    ax = axs[1]
    rel = np.zeros(N)
    rel[valid] = (B_ach[valid] - B_target[valid]) / B_target[valid]
    ax.plot(x_eval[valid] * 1e3, rel[valid] * 100, color="black", lw=1.5)
    ax.axhline(0, color="gray", ls="--")
    ax.set_xlabel("Position (mm)"); ax.set_ylabel("Relative error (%)")
    ax.set_title(f"Relative error (max = {max_error*100:.2f}%)"); ax.grid(alpha=.3)
    ax = axs[2]
    cols = [BLUE if n > 0 else (VERM if n < 0 else "#CCCCCC") for n in n_windings]
    bw = 0.7 * np.min(np.abs(np.diff(np.sort(groove_positions)))) * 1e3
    ax.bar(groove_positions * 1e3, n_windings, width=bw, color=cols, edgecolor="black", lw=.5)
    ax.axhline(0, color="black", lw=.8)
    ax.set_xlabel("Groove position (mm)"); ax.set_ylabel("Windings")
    ax.set_title(f"Winding pattern (forward: {n_fwd.sum()}, bucking: {abs(n_buck.sum())})")
    ax.set_ylim(min(n_min, n_windings.min()) - 1, max(n_max, n_windings.max()) + 1)
    ax.grid(alpha=.3, axis="y")
    ax.legend(handles=[mpatches.Patch(color=BLUE, label="Forward"), mpatches.Patch(color=VERM, label="Bucking")])
    fig.tight_layout()
    figs.append({"title": "Field performance and winding pattern",
                 "caption": "Top: target and achieved on-axis field. Middle: error. Bottom: turns in each groove.",
                 "png": png(fig)})

    # 2. side view (former drawn as a guide)
    Rm, Lm, ftm, wdm = coil_radius * 1e3, coil_length * 1e3, former_thickness * 1e3, wire_ins * 1e3
    max_layers = max(int(np.max(np.abs(n_windings))), 1)
    r_outer = Rm + ftm + max_layers * wdm
    wire_vis = max(wdm / 2, Rm * 0.025)
    fig, ax = plt.subplots(figsize=(12, 6))
    for y0 in (Rm, -Rm - ftm):
        ax.add_patch(Rectangle((-Lm / 2, y0), Lm, ftm, ec="saddlebrown", fc="wheat", lw=1.5, zorder=1))
    ax.axhline(0, color="gray", ls="-.", lw=.5)
    for pos, n_w in zip(groove_positions, n_windings):
        if n_w == 0:
            continue
        c, ec = (BLUE, DARK_BLUE) if n_w > 0 else (VERM, DARK_VERM)
        for w in range(abs(int(n_w))):
            r = Rm + ftm + wdm * (w + .5)
            for sgn in (1, -1):
                ax.add_patch(plt.Circle((pos * 1e3, sgn * r), wire_vis, fc=c, ec=ec, lw=.4, zorder=3))
    dim_y = -(r_outer + Rm * .25)
    ax.annotate("", xy=(Lm / 2, dim_y), xytext=(-Lm / 2, dim_y), arrowprops=dict(arrowstyle="<->", lw=1.5))
    ax.text(0, dim_y - Rm * .12, f"L = {Lm:.1f} mm", ha="center", fontweight="bold")
    ax.set_aspect("equal")
    ax.set_xlim(-Lm / 2 - max(Lm * .08, 3), Lm / 2 + max(Lm * .08, 3))
    ax.set_ylim(dim_y - Rm * .3, r_outer + Rm * .2)
    ax.set_xlabel("Axial position (mm)"); ax.set_ylabel("Radial position (mm)")
    ax.set_title("Side view: wire cross-sections around the former\n"
                 f"ID {2*Rm:.1f} mm | up to {max_layers} layer(s) | L = {L_total*1e6:.1f} uH")
    ax.grid(alpha=.2)
    ax.legend(handles=[mpatches.Patch(fc=BLUE, ec=DARK_BLUE, label="Forward"),
                       mpatches.Patch(fc=VERM, ec=DARK_VERM, label="Bucking"),
                       mpatches.Patch(fc="wheat", ec="saddlebrown", label="Former")], loc="lower right", fontsize=9)
    fig.tight_layout()
    figs.append({"title": "Side view", "caption": "Wire cross-sections above and below the axis, with the former wall in tan.",
                 "png": png(fig)})

    # 3. developed view
    layers_top = max(n_max, int(np.max(np.abs(n_windings))))
    fig, ax = plt.subplots(figsize=(12, 4.2))
    ax.fill_between([-Lm / 2, Lm / 2], -.3, 0, color="wheat", ec="saddlebrown", lw=1.5)
    for pos, n_w in zip(groove_positions, n_windings):
        if n_w == 0:
            continue
        c, mk = (BLUE, "o") if n_w > 0 else (VERM, "s")
        for w in range(abs(int(n_w))):
            ax.plot(pos * 1e3, w + 1, mk, color=c, ms=9, mec="black", mew=.6)
        ax.text(pos * 1e3, abs(n_w) + .3, f"{int(n_w):+d}", ha="center", fontsize=7, color=c, fontweight="bold")
    ax.set_yticks(range(0, layers_top + 1))
    ax.set_yticklabels(["Former"] + [f"Layer {i}" for i in range(1, layers_top + 1)])
    ax.set_xlabel("Axial position (mm)"); ax.set_title("Developed view: layers in each groove")
    ax.set_xlim(-Lm / 2 - 10, Lm / 2 + 10); ax.set_ylim(-.6, layers_top + 1.5); ax.grid(alpha=.3)
    ax.legend(handles=[plt.Line2D([0], [0], marker="o", color="w", mfc=BLUE, ms=10, mec="k", label="Forward"),
                       plt.Line2D([0], [0], marker="s", color="w", mfc=VERM, ms=10, mec="k", label="Bucking")],
              loc="upper right")
    fig.tight_layout()
    figs.append({"title": "Developed view", "caption": "Each groove unrolled: how many turns stack up and in which direction.",
                 "png": png(fig)})

    # 4. B1 field map
    x_map, y_map, B1 = b1_field_map(groove_positions, n_windings, coil_radius, former_thickness, wire_ins, coil_length, res_map)
    fig, axs = plt.subplots(2, 1, figsize=(11, 8))
    ax = axs[0]
    im = ax.imshow(B1 * 1e6, extent=[x_map[0] * 1e3, x_map[-1] * 1e3, y_map[0] * 1e3, y_map[-1] * 1e3],
                   aspect="auto", cmap="cividis", origin="lower", vmin=0, vmax=np.percentile(np.abs(B1), 99) * 1e6)
    plt.colorbar(im, ax=ax, label="|B1| (uT/A)")
    ax.set_xlabel("x (mm)"); ax.set_ylabel("y (mm)"); ax.set_title(f"B1 field map: {profile_desc}")
    ax = axs[1]
    ny = len(y_map)
    for iy in [0, ny // 4, ny // 2, 3 * ny // 4, ny - 1]:
        ax.plot(x_map * 1e3, B1[iy] * 1e6, lw=1.4, label=f"y = {y_map[iy]*1e3:.1f} mm")
    ax.set_xlabel("x (mm)"); ax.set_ylabel("|B1| (uT/A)"); ax.set_title("B1 profiles at different y offsets")
    ax.legend(fontsize=8, ncol=2); ax.grid(alpha=.3)
    fig.tight_layout()
    figs.append({"title": "B1 field map", "caption": "Field inside and around the coil, computed loop by loop (Biot-Savart).",
                 "png": png(fig)})

    # 5. inductance matrices
    if L_mat is not None:
        fig, ax = plt.subplots(figsize=(6, 5))
        im = ax.imshow(L_mat * 1e6, cmap="viridis", origin="lower")
        plt.colorbar(im, ax=ax, label="L or M (uH)")
        ax.set_xlabel("Turn index"); ax.set_ylabel("Turn index")
        ax.set_title(f"Turn-by-turn inductance ({len(turns)} turns)")
        fig.tight_layout()
        figs.append({"title": "Inductance matrix", "caption": "Self inductance on the diagonal, mutual inductance between turns elsewhere.",
                     "png": png(fig)})

    summary = [
        ["Target profile", profile_desc],
        ["Max field error", f"{max_error*100:.2f} %"],
        ["Forward turns", int(n_fwd.sum())],
        ["Bucking turns", int(abs(n_buck.sum()))],
        ["Total inductance", f"{L_total*1e6:.3f} uH"],
        ["Reactance", f"{omega*L_total:.2f} ohm"],
        ["Tuning capacitance", f"{C_res*1e12:.1f} pF" if C_res else "n/a"],
        ["Wire length", f"{wire_len:.3f} m"],
        ["Skin depth", f"{skin*1e6:.1f} um"],
        ["DC / AC resistance", f"{R_dc*1e3:.2f} / {R_ac*1e3:.2f} mohm"],
        ["Q factor (unloaded)", f"{Q:.0f}"],
    ]
    return {"summary": summary, "windings": rows, "instructions": instructions, "notes": notes, "figures": figs}


def run(params_json, stl_path=""):
    try:
        return json.dumps(_run(json.loads(params_json)))
    except Exception as e:
        traceback.print_exc()
        return json.dumps({"error": str(e)})
