import argparse
import json
import os
import sys
from pathlib import Path
from typing import Any, Dict

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from mapdl.utils import resolve_material, derive_params, ensure_run_dir, material_code


TEMPLATE_FILE = Path(__file__).parent / "templates" / "isogrid_0deg_tri_template.inp"


def load_cases(path: Path) -> Dict[str, Any]:
    data = json.loads(path.read_text())
    if "cases" not in data or not isinstance(data["cases"], list):
        raise ValueError("Cases file must contain a 'cases' list")
    return data


def render_template(
    template: str,
    mat: Dict[str, float],
    p: Dict[str, float],
    extras: Dict[str, str],
) -> str:
    # Flatten material and params into a single dict for format()
    ctx = {
        "MAT_NAME": mat["name"],
        "E": mat["E"],
        "NU": mat["nu"],
        "RHO": mat["rho"],
        "SIGMA_Y": mat["sigma_y"],
        "SIGMA_ULT": mat["sigma_ult"],
        "R": p["R"],
        "L": p["L"],
        "A": p["a"],
        "NTHETA": p["n_theta"],
        "B": p["b"],
        "T": p["t"],
        "DZ": p["dz"],
    }
    ctx.update(extras)
    return template.format(**ctx)


def try_launch_mapdl(input_file: Path, workdir: Path) -> None:
    try:
        from ansys.mapdl.core import launch_mapdl  # type: ignore
    except Exception:
        print("ansys.mapdl.core not available; skipping MAPDL run. Generated:", input_file)
        return

    mapdl = launch_mapdl(run_location=str(workdir))
    try:
        mapdl.input(str(input_file))
        mapdl.finish()
        print("MAPDL run completed for", workdir.name)
    finally:
        mapdl.exit()


def main() -> None:
    ap = argparse.ArgumentParser(description="Run parametric MAPDL cases for isogrid shells")
    ap.add_argument("--cases", type=Path, required=True, help="Path to JSON cases file")
    ap.add_argument("--dry-run", action="store_true", help="Only write input files; do not run MAPDL")
    args = ap.parse_args()

    data = load_cases(args.cases)
    template = TEMPLATE_FILE.read_text()
    root = Path(__file__).parent

    description = data.get("description")
    if description:
        print(f"Case library: {description}")
    print(f"Total cases: {len(data['cases'])}")

    for case in data["cases"]:
        case_id = case["id"]
        mat_name = case["material"]
        params = derive_params(case["params"])
        material = resolve_material(mat_name)
        mat_dict = {
            "name": material.name,
            "E": material.E,
            "nu": material.nu,
            "rho": material.rho,
            "sigma_y": material.sigma_y,
            "sigma_ult": material.sigma_ult,
        }
        mat_code = material_code(mat_name)

        run_dir = ensure_run_dir(root, case_id)
        extras = {
            "CASE_ID": case_id,
            "MAT_CODE": mat_code,
        }
        inp_text = render_template(template, mat_dict, params, extras)
        inp_path = run_dir / f"{case_id}.inp"
        inp_path.write_text(inp_text)

        print(f"Prepared case '{case_id}' in {run_dir}")

        if not args.dry_run:
            try_launch_mapdl(inp_path, run_dir)


if __name__ == "__main__":
    main()
