import math
from pathlib import Path
from typing import Dict, Any

# Reuse materials and geometry helpers from existing package
from isogrid.materials import (
    IsotropicMaterial,
    MATERIAL_LIBRARY,
)
from isogrid.geometry import n_theta_from_a, a_from_n_theta, delta_z_from_a


MATERIAL_REGISTRY: Dict[str, IsotropicMaterial] = {
    key: mat for key, mat in MATERIAL_LIBRARY.items()
}

# Provide some human-readable aliases matching dashboard/heatmap usage
MATERIAL_REGISTRY.update({
    "AISI 304": MATERIAL_LIBRARY["AISI_304"],
    "AISI 316L": MATERIAL_LIBRARY["AISI_316L"],
    "Al 6061-T6": MATERIAL_LIBRARY["AL_6061_T6"],
    "Al 2024-T3": MATERIAL_LIBRARY["AL_2024_T3"],
    "Al 7075-T6": MATERIAL_LIBRARY["AL_7075_T6"],
})


MATERIAL_CODE_MAP = {
    "AISI_304": "304",
    "STEEL_304": "304",
    "AISI_316L": "316L",
    "AL_6061_T6": "AL6061",
    "AL_2024_T3": "AL2024",
    "AL_7075_T6": "AL7075",
}


def resolve_material(name: str) -> IsotropicMaterial:
    if name in MATERIAL_REGISTRY:
        return MATERIAL_REGISTRY[name]

    normalized = name.replace(" ", "_").upper()
    if normalized in MATERIAL_REGISTRY:
        return MATERIAL_REGISTRY[normalized]

    raise KeyError(f"Unknown material '{name}'. Known: {sorted(MATERIAL_REGISTRY)}")


def canonical_material_key(name: str) -> str:
    mat = resolve_material(name)
    for key, candidate in MATERIAL_LIBRARY.items():
        if candidate == mat:
            return key
    normalized = mat.name.replace(" ", "_").upper()
    MATERIAL_REGISTRY[normalized] = mat
    return normalized


def material_code(name: str) -> str:
    key = canonical_material_key(name)
    if key in MATERIAL_CODE_MAP:
        return MATERIAL_CODE_MAP[key]

    mat = resolve_material(name)
    alnum = ''.join(ch for ch in mat.name if ch.isalnum())
    return alnum.upper()


def derive_params(params: Dict[str, Any]) -> Dict[str, float]:
    """Return a dict with required lattice params.

    Required inputs include: R, L, (a or n_theta), b, t
    Outputs include: R, L, a, n_theta, b, t, dz
    """
    R = float(params["R"])  # m
    L = float(params["L"])  # m
    b = float(params["b"])  # m
    t = float(params["t"])  # m

    a = params.get("a")
    n_theta = params.get("n_theta")

    if a is None and n_theta is None:
        raise ValueError("Provide either 'a' or 'n_theta' in params")

    if a is None:
        a = a_from_n_theta(R, int(n_theta))
    else:
        a = float(a)

    if n_theta is None:
        n_theta = n_theta_from_a(R, a)
    else:
        n_theta = int(n_theta)

    dz = delta_z_from_a(a)

    return {
        "R": R,
        "L": L,
        "a": float(a),
        "n_theta": int(n_theta),
        "b": b,
        "t": t,
        "dz": dz,
    }


def ensure_run_dir(root: Path, case_id: str) -> Path:
    run_dir = root / "runs" / case_id
    run_dir.mkdir(parents=True, exist_ok=True)
    return run_dir
