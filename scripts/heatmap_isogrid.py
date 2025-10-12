#!/usr/bin/env python3
"""
Generate 3D heatmap data for a skinless isogrid cylinder.

Defaults:
- R = 550 mm (outer radius)
- Material = AISI 304
- 3D axes: b [mm], t [mm], n_theta [-]
- Color metric: mass-specific global buckling capacity N_cr_global / m

If Plotly is installed, produces an interactive 3D scatter; otherwise writes CSV.
"""
import argparse
import csv
import math
import sys
from pathlib import Path

# Allow running this script from any working directory by adding the repo root to sys.path
_ROOT = Path(__file__).resolve().parents[1]
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

try:
    import plotly.graph_objects as go  # type: ignore
    import plotly.io as pio  # type: ignore
    _HAS_PLOTLY = True
    # Prefer opening in browser on macOS/desktop
    try:
        if pio.renderers.default == "notebook":
            pio.renderers.default = "browser"
    except Exception:
        pass
except Exception:
    _HAS_PLOTLY = False
from isogrid.materials import STEEL_304
from isogrid.geometry import a_from_n_theta, n_theta_from_a
from isogrid.mechanics import (
    GlobalAssumptions,
    metric_mass_specific_capacity,
    Ncr_global,
    total_mass,
    areal_mass,
    local_safety_factors,
    bending_global_sf,
    torsion_sf,
)


def mm(val_mm: float) -> float:
    return val_mm / 1000.0


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    p.add_argument("--R_mm", type=float, default=550.0, help="Cylinder outer radius [mm] (default: 550)")
    p.add_argument("--L_m", type=float, default=0.8, help="Segment length [m] (default: 0.8)")
    p.add_argument("--b_mm", type=float, nargs=3, metavar=("start","stop","count"), default=[6.0, 25.0, 15], help="Rib width sweep [mm]: start stop count")
    p.add_argument("--t_mm", type=float, nargs=3, metavar=("start","stop","count"), default=[1.0, 4.0, 12], help="Rib thickness sweep [mm]: start stop count")
    p.add_argument("--n_theta", type=int, nargs=3, metavar=("start","stop","count"), default=[40, 140, 21], help="Circumferential cell count sweep: start stop count")
    p.add_argument("--a_mm", type=float, nargs=3, metavar=("start","stop","count"), default=[25.0, 90.0, 20], help="Cell edge length a sweep [mm]: start stop count (used if axis is 'a')")
    p.add_argument("--L_m_sweep", type=float, nargs=3, metavar=("start","stop","count"), default=[0.3, 1.2, 10], help="Segment length L sweep [m]: start stop count (used if axis is 'L')")
    p.add_argument("--K", type=float, default=1.0, help="Euler K-factor for global buckling (default: 1.0)")
    p.add_argument("--N_req", type=float, default=None, help="Optional axial force [N] to compute safety factors")
    p.add_argument("--M_req", type=float, default=None, help="Optional bending moment [N·m] for bending safety factor")
    p.add_argument("--T_req", type=float, default=None, help="Optional torque [N·m] for torsion safety factor")
    p.add_argument("--out", type=Path, default=Path("out/isogrid_heatmap.csv"), help="Output CSV path (written always)")
    p.add_argument("--show", action="store_true", help="Open interactive plot output")
    p.add_argument("--no_save_html", action="store_true", help="Do not write HTML files; only show plots if --show is set")
    p.add_argument("--renderer", type=str, default="browser", help="Plotly renderer to use for fig.show(); e.g., browser, vega, png, notebok")
    # Plot mode: scatter over 3 params (current behavior) or surface (2 params -> z)
    p.add_argument("--plot", choices=["scatter", "surface"], default="surface", help="Plot type: scatter (3 params) or surface (2D sweep)")
    # Interactive slider across a third parameter for surface plots
    p.add_argument("--interactive", action="store_true", help="Create an interactive surface with a slider over a third parameter")
    p.add_argument("--slider_axis", choices=["b", "t", "n_theta", "a", "L"], default="n_theta", help="Which parameter is controlled by the slider in interactive mode")
    p.add_argument(
        "--metric",
        choices=["ncr_over_m", "ncr_global", "mass", "areal_mass", "sf_min"],
        default="ncr_global",
        help="Quantity to visualize: mass-specific Ncr, absolute Ncr, total mass, areal mass density, or min safety factor (requires --N_req)",
    )
    p.add_argument("--x_axis", choices=["b", "t", "n_theta", "a", "L"], default="b", help="Surface X axis parameter")
    p.add_argument("--y_axis", choices=["b", "t", "n_theta", "a", "L"], default="t", help="Surface Y axis parameter")
    p.add_argument("--fix_b_mm", type=float, default=12.0, help="Fixed b [mm] if b is not an axis")
    p.add_argument("--fix_t_mm", type=float, default=2.0, help="Fixed t [mm] if t is not an axis")
    p.add_argument("--fix_n_theta", type=int, default=100, help="Fixed n_theta if n_theta is not an axis")
    p.add_argument("--fix_a_mm", type=float, default=60.0, help="Fixed a [mm] if a is not an axis (overrides fix_n_theta)")
    p.add_argument("--fix_L_m", type=float, default=0.8, help="Fixed L [m] if L is not an axis")
    p.add_argument("--batch", action="store_true", help="Generate multiple preset surface plots as separate HTML files")
    return p.parse_args()


def frange(start: float, stop: float, count: int):
    if count <= 1:
        yield start
        return
    step = (stop - start) / (count - 1)
    for i in range(count):
        yield start + i * step


def main():
    args = parse_args()

    R = mm(args.R_mm)
    L = args.L_m
    K = args.K
    asum = GlobalAssumptions()

    # Configure renderer preference if Plotly present
    if _HAS_PLOTLY:
        try:
            pio.renderers.default = args.renderer
        except Exception:
            pass

    b_vals_mm = list(frange(args.b_mm[0], args.b_mm[1], int(args.b_mm[2])))
    t_vals_mm = list(frange(args.t_mm[0], args.t_mm[1], int(args.t_mm[2])))
    ntheta_vals = list(range(int(args.n_theta[0]), int(args.n_theta[1]) + 1, max(1, int((args.n_theta[1]-args.n_theta[0])/(args.n_theta[2]-1) if args.n_theta[2] > 1 else 1))))
    a_vals_mm = list(frange(args.a_mm[0], args.a_mm[1], int(args.a_mm[2])))
    L_vals_m = list(frange(args.L_m_sweep[0], args.L_m_sweep[1], int(args.L_m_sweep[2])))

    rows = []
    for bmm in b_vals_mm:
        for tmm in t_vals_mm:
            for n_th in ntheta_vals:
                b = mm(bmm)
                t = mm(tmm)
                a = a_from_n_theta(R, n_th)

                try:
                    ncr = Ncr_global(STEEL_304, R, L, b, t, a, K)
                    mass_total = total_mass(STEEL_304, R, L, b, t, a)
                    if args.metric == "ncr_over_m":
                        metric_val = ncr / max(mass_total, 1e-12)
                    elif args.metric == "ncr_global":
                        metric_val = ncr
                    elif args.metric == "mass":
                        metric_val = mass_total
                    elif args.metric == "sf_min":
                        if args.N_req is None or args.N_req <= 0:
                            raise ValueError("--metric sf_min requires a positive --N_req")
                        sfs = local_safety_factors(
                            STEEL_304, b, t, a, L, R, args.N_req, asum
                        )
                        sf_local = sfs["sf_min"]
                        sf_global = ncr / args.N_req
                        metric_val = min(sf_local, sf_global)
                    elif args.metric == "sf_bending":
                        metric_val = bending_global_sf(STEEL_304, R, b, t, a, args.M_req, asum)
                    elif args.metric == "sf_torsion":
                        metric_val = torsion_sf(STEEL_304, R, b, t, a, args.T_req, asum)
                    elif args.metric == "areal_mass":
                        metric_val = areal_mass(STEEL_304, b, t, a)
                    else:
                        metric_val = float("nan")
                except Exception:
                    metric_val = float("nan")
                    ncr = float("nan")
                    mass_total = float("nan")

                rows.append({
                    "b_mm": bmm,
                    "t_mm": tmm,
                    "n_theta": n_th,
                    "a_mm": a * 1000.0,
                    "metric_value": metric_val,
                    "Ncr_global_N": ncr,
                    "mass_total_kg": mass_total,
                })

    # Ensure output directory exists
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote CSV: {args.out}")

    if _HAS_PLOTLY:
        # Helper to build a surface figure for arbitrary axes/metric and fixed parameters
        def build_surface_figure(x_axis: str, y_axis: str, metric: str,
                                 fix_b_mm: float, fix_t_mm: float, fix_n_theta: int,
                                 fix_a_mm: float, fix_L_m: float):
            if x_axis == y_axis:
                raise ValueError("x_axis and y_axis must differ")

            def get_axis_values(axis_name: str):
                if axis_name == "b":
                    return b_vals_mm
                if axis_name == "t":
                    return t_vals_mm
                if axis_name == "n_theta":
                    return ntheta_vals
                if axis_name == "a":
                    return a_vals_mm
                if axis_name == "L":
                    return L_vals_m
                raise ValueError(axis_name)

            X_vals = list(get_axis_values(x_axis))
            Y_vals = list(get_axis_values(y_axis))

            Z = []
            for yv in Y_vals:
                row = []
                for xv in X_vals:
                    bmm = fix_b_mm if x_axis != "b" and y_axis != "b" else (xv if x_axis == "b" else yv)
                    tmm = fix_t_mm if x_axis != "t" and y_axis != "t" else (xv if x_axis == "t" else yv)
                    L_here = fix_L_m if x_axis != "L" and y_axis != "L" else (xv if x_axis == "L" else yv)

                    if x_axis == "a" or y_axis == "a":
                        a_mm_here = float(xv) if x_axis == "a" else float(yv)
                        a_here = mm(a_mm_here)
                        nth = n_theta_from_a(R, a_here)
                    else:
                        if x_axis == "n_theta" or y_axis == "n_theta":
                            nth = int(xv) if x_axis == "n_theta" else int(yv)
                            a_here = a_from_n_theta(R, nth)
                        elif fix_a_mm is not None:
                            a_here = mm(float(fix_a_mm))
                            nth = n_theta_from_a(R, a_here)
                        else:
                            nth = int(fix_n_theta)
                            a_here = a_from_n_theta(R, nth)

                    b_here = mm(float(bmm))
                    t_here = mm(float(tmm))
                    try:
                        ncr_val = Ncr_global(STEEL_304, R, float(L_here), b_here, t_here, a_here, K)
                        mass_val = total_mass(STEEL_304, R, float(L_here), b_here, t_here, a_here)
                        if metric == "ncr_over_m":
                            zval = ncr_val / max(mass_val, 1e-12)
                        elif metric == "ncr_global":
                            zval = ncr_val
                        elif metric == "mass":
                            zval = mass_val
                        elif metric == "areal_mass":
                            zval = areal_mass(STEEL_304, b_here, t_here, a_here)
                        elif metric == "sf_min":
                            if args.N_req is None or args.N_req <= 0:
                                zval = float("nan")
                            else:
                                sfs = local_safety_factors(
                                    STEEL_304, b_here, t_here, a_here, float(L_here), R, args.N_req, asum
                                )
                                sf_local = sfs["sf_min"]
                                sf_global = ncr_val / args.N_req
                                zval = min(sf_local, sf_global)
                        elif metric == "sf_bending":
                            zval = bending_global_sf(STEEL_304, R, b_here, t_here, a_here, args.M_req, asum)
                        elif metric == "sf_torsion":
                            zval = torsion_sf(STEEL_304, R, b_here, t_here, a_here, args.T_req, asum)
                        else:
                            zval = float("nan")
                    except Exception:
                        zval = float("nan")
                    row.append(zval)
                Z.append(row)

            axis_title = {
                "b": "Rippenbreite b [mm]",
                "t": "Rippendicke t [mm]",
                "n_theta": "Umfangszellzahl n_theta [-]",
                "a": "Zellkante a [mm]",
                "L": "Segmentlänge L [m]",
            }
            ctitle = {
                "ncr_over_m": "Ncr/m [m/s^2]",
                "ncr_global": "Ncr [N]",
                "mass": "Masse [kg]",
                "areal_mass": "Arealmasse m̄ [kg/m²]",
                "sf_min": "Sicherheitsfaktor min [-]",
                "sf_bending": "Sicherheitsfaktor Biegung [-]",
                "sf_torsion": "Sicherheitsfaktor Torsion [-]",
            }[args.metric if metric is None else metric]
            ctick = {
                "ncr_over_m": ".2e",
                "ncr_global": ".2e",
                "mass": ".2e",
                "areal_mass": ".2e",
                "sf_min": ".2f",
                "sf_bending": ".2f",
                "sf_torsion": ".2f",
            }[args.metric if metric is None else metric]

            fixed_info = []
            if "b" not in (x_axis, y_axis):
                fixed_info.append(f"b={fix_b_mm} mm")
            if "t" not in (x_axis, y_axis):
                fixed_info.append(f"t={fix_t_mm} mm")
            if "n_theta" not in (x_axis, y_axis) and ("a" in (x_axis, y_axis) or fix_a_mm is None):
                fixed_info.append(f"n_theta={fix_n_theta}")
            if "a" not in (x_axis, y_axis):
                fixed_info.append(f"a={fix_a_mm} mm")
            if "L" not in (x_axis, y_axis):
                fixed_info.append(f"L={fix_L_m} m")
            title = f"Isogrid Surface — {ctitle}; fixed: " + ", ".join(fixed_info)

            fig = go.Figure(
                data=[
                    go.Surface(
                        x=X_vals,
                        y=Y_vals,
                        z=Z,
                        colorscale="RdBu",
                        colorbar=dict(title=ctitle, tickformat=ctick),
                    )
                ]
            )
            fig.update_layout(
                title=title,
                scene=dict(
                    xaxis_title=axis_title[x_axis],
                    yaxis_title=axis_title[y_axis],
                    zaxis_title=ctitle,
                ),
                margin=dict(l=0, r=0, b=0, t=40),
            )
            return fig
        # Decide filename for the primary plot
        if args.plot == "scatter":
            # 3-parameter scatter
            x = [r["b_mm"] for r in rows]
            y = [r["t_mm"] for r in rows]
            z = [r["n_theta"] for r in rows]
            c = [r["metric_value"] for r in rows]
            ctitle = {
                "ncr_over_m": "Ncr/m [m/s^2]",
                "ncr_global": "Ncr [N]",
                "mass": "Masse [kg]",
                "areal_mass": "Arealmasse m̄ [kg/m²]",
                "sf_min": "Sicherheitsfaktor min [-]",
                "sf_bending": "Sicherheitsfaktor Biegung [-]",
                "sf_torsion": "Sicherheitsfaktor Torsion [-]",
            }[args.metric]
            ctick = {
                "ncr_over_m": ".2e",
                "ncr_global": ".2e",
                "mass": ".2e",
                "areal_mass": ".2e",
                "sf_min": ".2f",
                "sf_bending": ".2f",
                "sf_torsion": ".2f",
            }[args.metric]
            fig = go.Figure(
                data=[
                    go.Scatter3d(
                        x=x, y=y, z=z, mode="markers",
                        marker=dict(
                            size=3,
                            color=c,
                            colorscale="Viridis",
                            colorbar=dict(title=ctitle, tickformat=ctick),
                            showscale=True,
                        ),
                    )
                ]
            )
            fig.update_layout(
                title=f"Isogrid 3D Scatter — color: {ctitle}",
                scene=dict(
                    xaxis_title="Rippenbreite b [mm]",
                    yaxis_title="Rippendicke t [mm]",
                    zaxis_title="Umfangszellzahl n_theta [-]",
                ),
                margin=dict(l=0, r=0, b=0, t=40),
            )
            html_name = f"isogrid_scatter_{args.metric}.html"
        else:
            # True surface: choose two axes; fix the third parameter
            if args.x_axis == args.y_axis:
                raise SystemExit("x_axis and y_axis must differ for surface mode")

            # Build 2D grids
            def get_axis_values(axis_name: str):
                if axis_name == "b":
                    return b_vals_mm
                if axis_name == "t":
                    return t_vals_mm
                if axis_name == "n_theta":
                    return ntheta_vals
                if axis_name == "a":
                    return a_vals_mm
                if axis_name == "L":
                    return L_vals_m
                raise ValueError(axis_name)

            X_vals = list(get_axis_values(args.x_axis))
            Y_vals = list(get_axis_values(args.y_axis))

            # Map numeric values per axis to metric
            Z = []  # rows over Y, cols over X (Plotly surface expects 2D array)
            for yv in Y_vals:
                row = []
                for xv in X_vals:
                    # Determine b, t, n_theta for this cell
                    bmm = args.fix_b_mm if args.x_axis != "b" and args.y_axis != "b" else (xv if args.x_axis == "b" else yv)
                    tmm = args.fix_t_mm if args.x_axis != "t" and args.y_axis != "t" else (xv if args.x_axis == "t" else yv)
                    # Determine L
                    L_here = args.fix_L_m if args.x_axis != "L" and args.y_axis != "L" else (xv if args.x_axis == "L" else yv)
                    # Determine a (mm) and n_theta
                    if args.x_axis == "a" or args.y_axis == "a":
                        a_mm_here = float(xv) if args.x_axis == "a" else float(yv)
                        a = mm(a_mm_here)
                        nth = n_theta_from_a(R, a)
                    else:
                        # If n_theta is an axis, take it and derive a
                        if args.x_axis == "n_theta" or args.y_axis == "n_theta":
                            nth = int(xv) if args.x_axis == "n_theta" else int(yv)
                            a = a_from_n_theta(R, nth)
                        # Else prefer fixed a if provided; derive n_theta from it
                        elif args.fix_a_mm is not None:
                            a = mm(float(args.fix_a_mm))
                            nth = n_theta_from_a(R, a)
                        # Else use fixed n_theta and derive a
                        else:
                            nth = int(args.fix_n_theta)
                            a = a_from_n_theta(R, nth)

                    b = mm(float(bmm))
                    t = mm(float(tmm))
                    try:
                        ncr = Ncr_global(STEEL_304, R, float(L_here), b, t, a, K)
                        mass_total = total_mass(STEEL_304, R, float(L_here), b, t, a)
                        if args.metric == "ncr_over_m":
                            zval = ncr / max(mass_total, 1e-12)
                        elif args.metric == "ncr_global":
                            zval = ncr
                        elif args.metric == "mass":
                            zval = mass_total
                        elif args.metric == "areal_mass":
                            zval = areal_mass(STEEL_304, b, t, a)
                        elif args.metric == "sf_min":
                            if args.N_req is None or args.N_req <= 0:
                                zval = float("nan")
                            else:
                                sfs = local_safety_factors(
                                    STEEL_304, b, t, a, float(L_here), R, args.N_req, asum
                                )
                                sf_local = sfs["sf_min"]
                                sf_global = ncr / args.N_req
                                zval = min(sf_local, sf_global)
                        else:
                            zval = float("nan")
                    except Exception:
                        zval = float("nan")
                    row.append(zval)
                Z.append(row)

            # Axes labels
            axis_title = {
                "b": "Rippenbreite b [mm]",
                "t": "Rippendicke t [mm]",
                "n_theta": "Umfangszellzahl n_theta [-]",
                "a": "Zellkante a [mm]",
                "L": "Segmentlänge L [m]",
            }
            fixed_info = []
            if "b" not in (args.x_axis, args.y_axis):
                fixed_info.append(f"b={args.fix_b_mm} mm")
            if "t" not in (args.x_axis, args.y_axis):
                fixed_info.append(f"t={args.fix_t_mm} mm")
            if "n_theta" not in (args.x_axis, args.y_axis) and ("a" in (args.x_axis, args.y_axis) or args.fix_a_mm is None):
                fixed_info.append(f"n_theta={args.fix_n_theta}")
            if "a" not in (args.x_axis, args.y_axis):
                fixed_info.append(f"a={args.fix_a_mm} mm")
            if "L" not in (args.x_axis, args.y_axis):
                fixed_info.append(f"L={args.fix_L_m} m")
            ctitle = {
                "ncr_over_m": "Ncr/m [m/s^2]",
                "ncr_global": "Ncr [N]",
                "mass": "Masse [kg]",
                "areal_mass": "Arealmasse m̄ [kg/m²]",
                "sf_min": "Sicherheitsfaktor min [-]",
            }[args.metric]
            ctick = {
                "ncr_over_m": ".2e",
                "ncr_global": ".2e",
                "mass": ".2e",
                "areal_mass": ".2e",
                "sf_min": ".2f",
            }[args.metric]
            title = f"Isogrid Surface — {ctitle}; fixed: " + ", ".join(fixed_info)

            fig = go.Figure(
                data=[
                    go.Surface(
                        x=X_vals,
                        y=Y_vals,
                        z=Z,
                        colorscale="RdBu",
                        colorbar=dict(title=ctitle, tickformat=ctick),
                    )
                ]
            )
            fig.update_layout(
                title=title,
                scene=dict(
                    xaxis_title=axis_title[args.x_axis],
                    yaxis_title=axis_title[args.y_axis],
                    zaxis_title=ctitle,
                ),
                margin=dict(l=0, r=0, b=0, t=40),
            )
            html_name = f"isogrid_surface_{args.metric}_{args.x_axis}x{args.y_axis}.html"

            # Optional interactive slider over a third parameter
            if args.interactive:
                if args.slider_axis in (args.x_axis, args.y_axis):
                    print("Slider axis must differ from x_axis and y_axis — ignoring --interactive.")
                else:
                    def get_axis_values(axis_name: str):
                        if axis_name == "b":
                            return b_vals_mm
                        if axis_name == "t":
                            return t_vals_mm
                        if axis_name == "n_theta":
                            return ntheta_vals
                        if axis_name == "a":
                            return a_vals_mm
                        if axis_name == "L":
                            return L_vals_m
                        raise ValueError(axis_name)

                    slider_vals = list(get_axis_values(args.slider_axis))
                    frames = []
                    step_labels = []

                    for sval in slider_vals:
                        # Recompute Z for this slider value
                        Zs = []
                        for yv in Y_vals:
                            row = []
                            for xv in X_vals:
                                # Determine current parameters
                                bmm = args.fix_b_mm if args.x_axis != "b" and args.y_axis != "b" else (xv if args.x_axis == "b" else yv)
                                tmm = args.fix_t_mm if args.x_axis != "t" and args.y_axis != "t" else (xv if args.x_axis == "t" else yv)
                                L_here = args.fix_L_m if args.x_axis != "L" and args.y_axis != "L" else (xv if args.x_axis == "L" else yv)

                                # Apply slider override
                                if args.slider_axis == "b":
                                    bmm = sval
                                elif args.slider_axis == "t":
                                    tmm = sval
                                elif args.slider_axis == "L":
                                    L_here = sval

                                # a / n_theta handling with slider
                                if args.x_axis == "a" or args.y_axis == "a":
                                    a_mm_here = float(xv) if args.x_axis == "a" else float(yv)
                                    a_here = mm(a_mm_here)
                                    nth = n_theta_from_a(R, a_here)
                                else:
                                    if args.x_axis == "n_theta" or args.y_axis == "n_theta":
                                        nth = int(xv) if args.x_axis == "n_theta" else int(yv)
                                        a_here = a_from_n_theta(R, nth)
                                    else:
                                        # slider may control 'a' or 'n_theta', otherwise use fixeds
                                        if args.slider_axis == "a":
                                            a_here = mm(float(sval))
                                            nth = n_theta_from_a(R, a_here)
                                        elif args.slider_axis == "n_theta":
                                            nth = int(sval)
                                            a_here = a_from_n_theta(R, nth)
                                        elif args.fix_a_mm is not None:
                                            a_here = mm(float(args.fix_a_mm))
                                            nth = n_theta_from_a(R, a_here)
                                        else:
                                            nth = int(args.fix_n_theta)
                                            a_here = a_from_n_theta(R, nth)

                                b_here = mm(float(bmm))
                                t_here = mm(float(tmm))

                                try:
                                    ncr_val = Ncr_global(STEEL_304, R, float(L_here), b_here, t_here, a_here, K)
                                    mass_val = total_mass(STEEL_304, R, float(L_here), b_here, t_here, a_here)
                                    if args.metric == "ncr_over_m":
                                        zval = ncr_val / max(mass_val, 1e-12)
                                    elif args.metric == "ncr_global":
                                        zval = ncr_val
                                    elif args.metric == "mass":
                                        zval = mass_val
                                    elif args.metric == "areal_mass":
                                        zval = areal_mass(STEEL_304, b_here, t_here, a_here)
                                    elif args.metric == "sf_min":
                                        if args.N_req is None or args.N_req <= 0:
                                            zval = float("nan")
                                        else:
                                            sfs = local_safety_factors(
                                                STEEL_304, b_here, t_here, a_here, float(L_here), R, args.N_req, asum
                                            )
                                            sf_local = sfs["sf_min"]
                                            sf_global = ncr_val / args.N_req
                                            zval = min(sf_local, sf_global)
                                    elif args.metric == "sf_bending":
                                        zval = bending_global_sf(STEEL_304, R, b_here, t_here, a_here, args.M_req, asum)
                                    elif args.metric == "sf_torsion":
                                        zval = torsion_sf(STEEL_304, R, b_here, t_here, a_here, args.T_req, asum)
                                    else:
                                        zval = float("nan")
                                except Exception:
                                    zval = float("nan")
                                row.append(zval)
                            Zs.append(row)

                        frames.append(go.Frame(data=[go.Surface(x=X_vals, y=Y_vals, z=Zs, colorscale="RdBu")], name=str(sval)))
                        step_labels.append(str(sval))

                    if frames:
                        fig.frames = frames
                        steps = [dict(method="animate", args=[[lab], {"mode": "immediate", "frame": {"duration": 0, "redraw": True}, "transition": {"duration": 0}}], label=str(lab)) for lab in step_labels]
                        fig.update_layout(
                            sliders=[{
                                "active": 0,
                                "currentvalue": {"prefix": f"{args.slider_axis} = "},
                                "pad": {"t": 30},
                                "steps": steps,
                            }]
                        )

        html_path = args.out.parent / html_name
        if args.no_save_html:
            if args.show:
                try:
                    fig.show()
                    print("Opened interactive plot (no HTML saved).")
                except Exception as e:
                    print(f"Plot show failed: {e}")
        else:
            try:
                fig.write_html(str(html_path), auto_open=bool(args.show))
                print(f"Wrote HTML: {html_path}")
            except Exception as e:
                print(f"Plotly write_html failed: {e}")

        # Batch preset plots
        if args.plot == "surface" and args.batch:
            presets = [
                ("b", "t", "ncr_global"),
                ("b", "t", "mass"),
                ("a", "t", "areal_mass"),
                ("a", "L", "ncr_global"),
                ("t", "n_theta", "ncr_global"),
            ]
            if args.N_req and args.N_req > 0:
                presets.append(("b", "t", "sf_min"))

            for xax, yax, met in presets:
                try:
                    fig2 = build_surface_figure(
                        xax, yax, met,
                        args.fix_b_mm, args.fix_t_mm, int(args.fix_n_theta),
                        args.fix_a_mm, args.fix_L_m,
                    )
                    name = f"isogrid_surface_{met}_{xax}x{yax}.html"
                    path = args.out.parent / name
                    if args.no_save_html:
                        if args.show:
                            try:
                                fig2.show()
                            except Exception as e:
                                print(f"Plot show failed for {name}: {e}")
                        else:
                            print(f"Skipped saving {name} due to --no_save_html (use --show to open).")
                    else:
                        fig2.write_html(str(path), auto_open=bool(args.show))
                        print(f"Wrote HTML: {path}")
                except Exception as e:
                    print(f"Preset plot {met} {xax}x{yax} failed: {e}")
    elif args.show and not _HAS_PLOTLY:
        print("Plotly not available. Install plotly to visualize or open the CSV in your tool.")


if __name__ == "__main__":
    sys.exit(main())
