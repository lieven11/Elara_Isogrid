import argparse
import json
import os
import sys
from pathlib import Path
from typing import Any, Dict

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from isogrid.mechanics import total_mass  # type: ignore
from isogrid.geometry import annulus_area

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
    for key, value in p.items():
        upper = key.upper()
        if upper not in ctx:
            ctx[upper] = value
    return template.format(**ctx)


ERROR_MARKERS = (
    "*** ERROR ***",
    "*** FATAL ERROR ***",
    "***FATAL ERROR***",
    "*** SEVERE ERROR ***",
    "SOLUTION ABORTED",
)


def try_launch_mapdl(input_file: Path, workdir: Path) -> bool:
    try:
        from ansys.mapdl.core import launch_mapdl  # type: ignore
    except Exception:
        print("ansys.mapdl.core not available; skipping MAPDL run. Generated:", input_file)
        return False

    mapdl = launch_mapdl(run_location=str(workdir))
    try:
        mapdl.input(str(input_file))
        mapdl.finish()
        return True
    except Exception as exc:  # pragma: no cover - MAPDL specific
        print(f"MAPDL run raised an exception for {workdir.name}: {exc}")
        return False
    finally:
        mapdl.exit()


def _logs_contain_errors(run_dir: Path) -> bool:
    for path in sorted(run_dir.glob("file*.err")):
        try:
            text = path.read_text(encoding="utf-8", errors="ignore")
        except OSError:
            continue
        for marker in ERROR_MARKERS:
            if marker in text:
                return True
    return False


def main() -> None:
    ap = argparse.ArgumentParser(description="Run parametric MAPDL cases for isogrid shells")
    ap.add_argument("--cases", type=Path, required=True, help="Path to JSON cases file")
    ap.add_argument("--dry-run", action="store_true", help="Only write input files; do not run MAPDL")
    args = ap.parse_args()

    data = load_cases(args.cases)
    template = TEMPLATE_FILE.read_text()
    root = Path(__file__).parent

    metadata = data.get("metadata", {})
    default_params = {}
    if isinstance(metadata, dict):
        for key in ("fixed", "defaults"):
            cfg = metadata.get(key)
            if isinstance(cfg, dict):
                default_params.update(cfg)

    description = data.get("description")
    if description:
        print(f"Case library: {description}")
    print(f"Total cases: {len(data['cases'])}")

    for case in data["cases"]:
        case_id = case["id"]
        mat_name = case["material"]
        case_params = case["params"]
        raw_params = {**default_params, **case_params}
        params = derive_params(raw_params)
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

        load_face = float(
            raw_params.get(
                "axial_load_N",
                raw_params.get("load_total", raw_params.get("load", 100.0)),
            )
        )
        face_area = annulus_area(params["R"], params["t"])
        if face_area <= 0.0:
            face_area = 1e-12
        face_pressure = load_face / face_area
        mass_total = total_mass(
            material,
            params["R"],
            params["L"],
            params["b"],
            params["t"],
            params["a"],
        )

        extras = {
            "CASE_ID": case_id,
            "MAT_CODE": mat_code,
            "AXIAL_LOAD_N": load_face,
            "FACE_AREA": face_area,
            "FACE_PRESSURE": face_pressure,
            "TOTAL_MASS": mass_total,
        }
        inp_text = render_template(template, mat_dict, params, extras)
        inp_path = run_dir / f"{case_id}.inp"
        inp_path.write_text(inp_text)

        print(f"Prepared case '{case_id}' in {run_dir}")

        if not args.dry_run:
            launched = try_launch_mapdl(inp_path, run_dir)
            if launched:
                if _logs_contain_errors(run_dir):
                    print(f"MAPDL run completed for {run_dir.name} [logs report errors]")
                else:
                    print(f"MAPDL run completed for {run_dir.name} [logs clean]")


if __name__ == "__main__":
    main()
