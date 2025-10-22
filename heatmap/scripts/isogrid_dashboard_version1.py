#!/usr/bin/env python3
from pathlib import Path
import argparse
import json
import sys

_ROOT = Path(__file__).resolve().parents[2]
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

# This is a frozen copy of heatmap/scripts/isogrid_dashboard.py at "version1".
# No further changes should be applied here; new features go into version2.

from heatmap.scripts.isogrid_dashboard import HTML_TEMPLATE as _T

def main():
    ap = argparse.ArgumentParser(description="Frozen version1 of the Isogrid Dashboard generator")
    default_out = Path(__file__).resolve().parents[1] / "out" / "isogrid_dashboard_version1.html"
    ap.add_argument("--out", type=Path, default=default_out)
    ap.add_argument("--R_mm", type=float, default=550.0)
    ap.add_argument("--K", type=float, default=1.0)
    ap.add_argument("--N_req", type=float, default=150000.0)
    ap.add_argument("--M_req", type=float, default=500.0)
    ap.add_argument("--T_req", type=float, default=800.0)
    # default ranges
    ap.add_argument("--b_mm", nargs=2, type=float, default=[6.0, 25.0])
    ap.add_argument("--t_mm", nargs=2, type=float, default=[1.0, 4.0])
    ap.add_argument("--n_theta", nargs=2, type=int, default=[40, 140])
    ap.add_argument("--a_mm", nargs=2, type=float, default=[25.0, 90.0])
    ap.add_argument("--L", nargs=2, type=float, default=[0.3, 1.2])
    ap.add_argument("--nx", type=int, default=15)
    ap.add_argument("--ny", type=int, default=12)
    ap.add_argument("--nz", type=int, default=12)
    args = ap.parse_args()

    state = {
        "xAxis": "b", "yAxis": "t", "zAxis": "n_theta",
        "plotType": "surface",
        "material": "AISI 304",
        "metric": "ncr_global",
        "KDF": 0.85,
        "nLayers": 6,
        "layerAlpha": 0.85,
        "nx": args.nx, "ny": args.ny, "nz": args.nz,
        "b_min_mm": args.b_mm[0], "b_max_mm": args.b_mm[1],
        "t_min_mm": args.t_mm[0], "t_max_mm": args.t_mm[1],
        "n_min": args.n_theta[0], "n_max": args.n_theta[1],
        "a_min_mm": args.a_mm[0], "a_max_mm": args.a_mm[1],
        "L_min": args.L[0], "L_max": args.L[1],
        "b_fix_mm": (args.b_mm[0] + args.b_mm[1]) / 2,
        "t_fix_mm": (args.t_mm[0] + args.t_mm[1]) / 2,
        "n_fix": (args.n_theta[0] + args.n_theta[1]) / 2,
        "a_fix_mm": (args.a_mm[0] + args.a_mm[1]) / 2,
        "L_fix": (args.L[0] + args.L[1]) / 2,
        "R_mm": args.R_mm, "K": args.K,
        "N_req": args.N_req, "M_req": args.M_req, "T_req": args.T_req,
    }

    html = _T % {
        "nx": args.nx, "ny": args.ny, "nz": args.nz,
        "n_layers": 6, "layer_alpha": 0.85,
        "b_min_mm": args.b_mm[0], "b_max_mm": args.b_mm[1],
        "t_min_mm": args.t_mm[0], "t_max_mm": args.t_mm[1],
        "n_min": args.n_theta[0], "n_max": args.n_theta[1],
        "a_min_mm": args.a_mm[0], "a_max_mm": args.a_mm[1],
        "L_min": args.L[0], "L_max": args.L[1],
        "R_mm": args.R_mm, "K": args.K, "KDF": 0.85, "N_req": args.N_req, "M_req": args.M_req, "T_req": args.T_req,
        "state_json": json.dumps(state),
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(html, encoding="utf-8")
    print(f"Wrote version1 dashboard: {args.out}")


if __name__ == "__main__":
    main()
