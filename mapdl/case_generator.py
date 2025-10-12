#!/usr/bin/env python3
"""Generate MAPDL case libraries from high-level sweep specifications."""

import argparse
import json
import sys
from decimal import Decimal, getcontext
from itertools import product
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from mapdl.utils import resolve_material, material_code, canonical_material_key  # noqa: E402


getcontext().prec = 16


def decimal_sequence(cfg: Dict[str, float]) -> List[float]:
    if "values" in cfg:
        return [float(v) for v in cfg["values"]]
    if "value" in cfg:
        return [float(cfg["value"])]

    try:
        start = Decimal(str(cfg["start"]))
        stop = Decimal(str(cfg["stop"]))
        step = Decimal(str(cfg["step"]))
    except KeyError as exc:  # pragma: no cover - config error surfaces to user
        raise ValueError(f"Range config missing key: {exc!s}") from exc

    if step <= 0:
        raise ValueError("Range 'step' must be positive")
    if stop < start:
        raise ValueError("Range 'stop' must be >= 'start'")

    values: List[float] = []
    current = start
    epsilon = step * Decimal("1e-6")
    while current <= stop + epsilon:
        values.append(float(current))
        current += step
    return values


def format_param_snippet(name: str, value: float) -> str:
    if name in {"t", "b"}:
        mm = value * 1000.0
        return f"{name}{mm:05.1f}".replace(".", "p")
    if name == "a":
        cm = value * 100.0
        return f"{name}{cm:03.0f}"
    if name == "L" or name == "R":
        cm = value * 100.0
        return f"{name}{cm:03.0f}"
    if name == "n_theta":
        return f"nth{int(round(value))}"
    if abs(value) >= 1.0:
        return f"{name}{value:04.1f}".replace(".", "p")
    return f"{name}{value:05.3f}".replace(".", "p")


def build_cases(spec: Dict[str, object]) -> Tuple[List[Dict[str, object]], Dict[str, object]]:
    prefix = spec.get("id_prefix", "case")
    materials = spec.get("materials")
    if not isinstance(materials, list) or not materials:
        raise ValueError("'materials' must be a non-empty list in spec")

    parameters_cfg = spec.get("parameters")
    if not isinstance(parameters_cfg, dict) or not parameters_cfg:
        raise ValueError("'parameters' must be a dict in spec")

    param_items: List[Tuple[str, List[float]]] = []
    for name, cfg in parameters_cfg.items():
        if not isinstance(cfg, dict):
            raise ValueError(f"Parameter '{name}' must map to a dict config")
        values = decimal_sequence(cfg)
        if not values:
            raise ValueError(f"Parameter '{name}' produced no values")
        param_items.append((name, values))

    varying_params = [name for name, values in param_items if len(values) > 1]
    id_param_order = spec.get("id_parameters") or varying_params
    id_param_order = [p for p in id_param_order if p in {name for name, _ in param_items}]

    total_per_material = 1
    for _, values in param_items:
        total_per_material *= len(values)

    cases: List[Dict[str, object]] = []
    metadata = {
        "id_prefix": prefix,
        "materials": [],
        "parameters": {},
        "total_per_material": total_per_material,
    }

    for name, values in param_items:
        metadata["parameters"][name] = {
            "count": len(values),
            "values": values,
        }

    for mat_name in materials:
        mat_obj = resolve_material(mat_name)
        canonical_key = canonical_material_key(mat_name)
        mat_code = material_code(mat_name)
        metadata["materials"].append({
            "input": mat_name,
            "canonical_key": canonical_key,
            "display_name": mat_obj.name,
            "code": mat_code,
        })

        for combo_index, combo in enumerate(product(*(values for _, values in param_items)), start=1):
            params = {name: value for (name, _), value in zip(param_items, combo)}
            snippets = [format_param_snippet(name, params[name]) for name in id_param_order]
            case_id_parts = [prefix, mat_code]
            if snippets:
                case_id_parts.extend(snippets)
            else:
                case_id_parts.append(f"{combo_index:05d}")
            case_id = "_".join(case_id_parts)

            cases.append({
                "id": case_id,
                "material": canonical_key,
                "params": params,
            })

    metadata["total_cases"] = len(cases)
    return cases, metadata


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--spec", type=Path, required=True, help="JSON spec describing sweeps")
    ap.add_argument("--out", type=Path, required=True, help="Output JSON file for generated cases")
    ap.add_argument("--max-cases", type=int, default=200000, help="Safety cap on total cases")
    ap.add_argument("--dry-run", action="store_true", help="Parse spec and report counts without writing")
    args = ap.parse_args()

    spec = json.loads(args.spec.read_text())
    cases, metadata = build_cases(spec)

    if metadata["total_cases"] > args.max_cases:
        raise ValueError(
            f"Spec expands to {metadata['total_cases']} cases which exceeds max {args.max_cases}. "
            "Adjust --max-cases or coarsen ranges."
        )

    print(f"Materials: {[m['canonical_key'] for m in metadata['materials']]}")
    print(
        " per material combinations:",
        metadata["total_per_material"],
        " total cases:",
        metadata["total_cases"],
    )

    if args.dry_run:
        return

    payload = {
        "description": spec.get("description", f"Generated from {args.spec.name}"),
        "spec_source": str(args.spec),
        "metadata": metadata,
        "cases": cases,
    }

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2))
    print(f"Wrote {metadata['total_cases']} cases to {args.out}")


if __name__ == "__main__":  # pragma: no cover - CLI entry point
    main()

