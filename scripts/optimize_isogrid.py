#!/usr/bin/env python3
"""
Simple derivative-free optimizer for the skinless isogrid model.

Default objective: minimize mass subject to constraints:
  - SF_min (global vs. local) >= sf_target
  - optional bending/torsion SF >= 1 if M_req/T_req provided

Variables: by default (b, t, n_theta). You can switch to (b, t, a) with --use_a.

Writes best result to stdout and a CSV of evaluated points.
"""
import argparse
import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Tuple, Dict, Any
import sys

# Allow running this script from any working directory by adding the repo root to sys.path
_ROOT = Path(__file__).resolve().parents[1]
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from isogrid.materials import STEEL_304, IsotropicMaterial
from isogrid.geometry import a_from_n_theta
from isogrid.mechanics import (
    GlobalAssumptions,
    Ncr_global,
    total_mass,
    areal_mass,
    local_safety_factors,
    bending_global_sf,
    torsion_sf,
)


def mm(x: float) -> float:
    return x / 1000.0


MATERIALS: Dict[str, IsotropicMaterial] = {
    "AISI 304": STEEL_304,
    "Al 6061-T6": IsotropicMaterial(name="Al 6061-T6", E=69e9, nu=0.33, rho=2700.0, sigma_y=276e6, sigma_ult=310e6),
    "AISI 316L": IsotropicMaterial(name="AISI 316L", E=193e9, nu=0.30, rho=8000.0, sigma_y=170e6, sigma_ult=485e6),
}


@dataclass
class Bounds:
    b_mm: Tuple[float, float]
    t_mm: Tuple[float, float]
    n_theta: Tuple[int, int]
    a_mm: Tuple[float, float]


def eval_point(material: IsotropicMaterial, R: float, L: float, K: float, KDF: float,
               use_a: bool, b_mm: float, t_mm: float, n_theta: int, a_mm: float,
               N_req: float, M_req: float, T_req: float, sf_target: float,
               assumptions: GlobalAssumptions) -> Dict[str, Any]:
    b = mm(b_mm)
    t = mm(t_mm)
    a = mm(a_mm) if use_a else a_from_n_theta(R, n_theta)
    mat = material

    ncr = KDF * Ncr_global(mat, R, L, b, t, a, K)
    mass = total_mass(mat, R, L, b, t, a)

    # Local/global SF
    sfs_local = local_safety_factors(mat, b, t, a, L, R, N_req or 0.0, assumptions)
    sf_local_min = sfs_local["sf_min"] if N_req and N_req > 0 else float('inf')
    sf_global = (ncr / N_req) if N_req and N_req > 0 else float('inf')
    sf_min = min(sf_local_min, sf_global)

    # Bending / torsion
    sf_b = bending_global_sf(mat, R, b, t, a, M_req or 0.0, assumptions) if M_req else float('inf')
    sf_t = torsion_sf(mat, R, b, t, a, T_req or 0.0, assumptions) if T_req else float('inf')

    feasible = (
        (sf_min >= sf_target if sf_target and sf_target > 0 else True)
        and (sf_b >= 1.0 if M_req else True)
        and (sf_t >= 1.0 if T_req else True)
    )

    return {
        "b_mm": b_mm, "t_mm": t_mm, "n_theta": n_theta, "a_mm": a * 1000.0,
        "Ncr_N": ncr, "mass_kg": mass,
        "SF_min": sf_min, "SF_b": sf_b, "SF_t": sf_t,
        "feasible": feasible,
    }


def pattern_search(material: IsotropicMaterial, R: float, L: float, K: float, KDF: float,
                   use_a: bool, bounds: Bounds, N_req: float, M_req: float, T_req: float,
                   sf_target: float, objective: str, iters: int = 200) -> Tuple[Dict[str, Any], list]:
    """Coordinate pattern search with shrinking steps and feasibility penalty."""
    asum = GlobalAssumptions()

    # Initial guess: mid of bounds
    b = 0.5 * (bounds.b_mm[0] + bounds.b_mm[1])
    t = 0.5 * (bounds.t_mm[0] + bounds.t_mm[1])
    if use_a:
        a = 0.5 * (bounds.a_mm[0] + bounds.a_mm[1])
        nth = None
    else:
        nth = int(0.5 * (bounds.n_theta[0] + bounds.n_theta[1]))
        a = 0.0

    sb, st, sa, sn = (b * 0.2 or 1.0), (t * 0.2 or 0.2), (10.0), (10)

    def score(pt):
        res = eval_point(material, R, L, K, KDF, use_a, pt['b'], pt['t'], pt.get('n_theta', 0), pt.get('a', 0.0), N_req, M_req, T_req, sf_target, asum)
        # Objective
        if objective == 'min_mass':
            obj = res['mass_kg']
        elif objective == 'max_ncr_over_m':
            obj = - (res['Ncr_N'] / max(res['mass_kg'], 1e-12))
        elif objective == 'max_ncr':
            obj = - res['Ncr_N']
        else:
            obj = res['mass_kg']
        # Penalty for infeasible
        if not res['feasible']:
            obj += 1e9
        return obj, res

    cur = {'b': b, 't': t}
    if use_a:
        cur['a'] = a
    else:
        cur['n_theta'] = nth

    history = []
    best_obj, best_res = score(cur)
    history.append(best_res)

    for _ in range(iters):
        improved = False
        # Try coordinate moves
        for k, step in [('b', sb), ('t', st), ('a' if use_a else 'n_theta', sa if use_a else sn)]:
            for sgn in (+1, -1):
                trial = cur.copy()
                if k == 'b':
                    trial['b'] = min(max(cur['b'] + sgn * step, bounds.b_mm[0]), bounds.b_mm[1])
                elif k == 't':
                    trial['t'] = min(max(cur['t'] + sgn * step, bounds.t_mm[0]), bounds.t_mm[1])
                elif k == 'a':
                    trial['a'] = min(max(cur['a'] + sgn * step, bounds.a_mm[0]), bounds.a_mm[1])
                else:
                    trial['n_theta'] = int(min(max(cur['n_theta'] + sgn * step, bounds.n_theta[0]), bounds.n_theta[1]))
                obj, res = score(trial)
                history.append(res)
                if obj < best_obj:
                    best_obj, best_res = obj, res
                    cur = trial
                    improved = True
        if not improved:
            # shrink steps
            sb *= 0.5; st *= 0.5
            if use_a:
                sa *= 0.5
            else:
                sn = max(1, int(sn * 0.5))
            if sb < 0.05 and st < 0.05 and ((use_a and sa < 0.5) or (not use_a and sn <= 1)):
                break

    return best_res, history


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument('--R_mm', type=float, default=550.0)
    ap.add_argument('--L_m', type=float, default=0.8)
    ap.add_argument('--K', type=float, default=1.0)
    ap.add_argument('--KDF', type=float, default=0.85)
    ap.add_argument('--material', type=str, default='AISI 304', choices=list(MATERIALS.keys()))
    ap.add_argument('--use_a', action='store_true', help='Use (b,t,a) as design variables instead of (b,t,n_theta)')
    ap.add_argument('--b_mm', nargs=2, type=float, default=[6.0, 25.0])
    ap.add_argument('--t_mm', nargs=2, type=float, default=[1.0, 4.0])
    ap.add_argument('--n_theta', nargs=2, type=int, default=[40, 140])
    ap.add_argument('--a_mm', nargs=2, type=float, default=[25.0, 90.0])
    ap.add_argument('--N_req', type=float, default=150000.0)
    ap.add_argument('--M_req', type=float, default=0.0)
    ap.add_argument('--T_req', type=float, default=0.0)
    ap.add_argument('--sf_target', type=float, default=1.2)
    ap.add_argument('--objective', type=str, default='min_mass', choices=['min_mass', 'max_ncr_over_m', 'max_ncr'])
    ap.add_argument('--iters', type=int, default=200)
    ap.add_argument('--out', type=Path, default=Path('out/optimize_trace.csv'))
    args = ap.parse_args()

    mat = MATERIALS[args.material]
    bounds = Bounds(b_mm=tuple(args.b_mm), t_mm=tuple(args.t_mm), n_theta=tuple(args.n_theta), a_mm=tuple(args.a_mm))

    best, hist = pattern_search(
        material=mat,
        R=mm(args.R_mm), L=args.L_m, K=args.K, KDF=args.KDF,
        use_a=args.use_a, bounds=bounds,
        N_req=args.N_req, M_req=args.M_req, T_req=args.T_req,
        sf_target=args.sf_target, objective=args.objective, iters=args.iters,
    )

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(hist[0].keys()))
        writer.writeheader(); writer.writerows(hist)

    print('Best design:')
    for k in ['b_mm','t_mm','n_theta','a_mm','Ncr_N','mass_kg','SF_min','SF_b','SF_t','feasible']:
        print(f'  {k}: {best.get(k)}')
    print(f'Trace written to {args.out}')


if __name__ == '__main__':
    main()
