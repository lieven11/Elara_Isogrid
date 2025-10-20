#!/usr/bin/env python3
"""Generate a readable summary of MAPDL sweep results.

Reads the consolidated ``summary.csv`` written by ``summarize_runs.py`` and
derives a ranked view of the best performing cases together with a few quick
analytics. By default, a Markdown report is emitted next to the input CSV.
"""

from __future__ import annotations

import argparse
import csv
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from isogrid.mechanics import total_mass  # type: ignore
from mapdl.utils import resolve_material  # type: ignore


DEFAULT_SUMMARY = Path(__file__).parent / "runs" / "summary.csv"
DEFAULT_REPORT = Path(__file__).parent / "runs" / "summary_report.md"

# Metrics we consider for ranking. The weight is only used if the caller does
# not override it via --weight tip_defl=0.5 etc.
METRICS = {
    "tip_per_load": {
        "label": "Tip Deflection per Load [m/N]",
        "higher_is_better": False,
        "default_weight": 0.35,
    },
    "buckling_per_mass": {
        "label": "Buckling Capacity per Mass [N/kg]",
        "higher_is_better": True,
        "default_weight": 0.40,
    },
    "sigma_per_load": {
        "label": "Stress per Load [Pa/N]",
        "higher_is_better": False,
        "default_weight": 0.25,
    },
}

# Columns we try to coerce to floats for nicer formatting later on.
NUMERIC_COLUMNS = {
    "R_m",
    "L_m",
    "a_m",
    "b_m",
    "t_m",
    "n_theta",
    "tip_defl",
    "tip_defl_m",
    "buckling",
    "buckling_factor",
    "sigma_ma",
    "sigma_max_pa",
    "axial_load_face_N",
    "face_pressure_pa",
    "load_pair_total_N",
    "buckling_face_N",
    "buckling_pair_N",
    "total_mass_kg",
    "tip_per_load_m_per_N",
    "sigma_per_load_pa_per_N",
    *METRICS.keys(),
}


@dataclass
class Row:
    data: Dict[str, Optional[float]]
    raw: Dict[str, str]
    score: Optional[float] = None

    def value(self, key: str) -> Optional[float]:
        return self.data.get(key)

    def text(self, key: str) -> str:
        return self.raw.get(key, "")


@dataclass
class MetricStats:
    name: str
    higher_is_better: bool
    minimum: float
    maximum: float

    @property
    def span(self) -> float:
        return self.maximum - self.minimum

    def normalised(self, value: float) -> float:
        if self.span <= 0:
            return 1.0
        if self.higher_is_better:
            return (value - self.minimum) / self.span
        return (self.maximum - value) / self.span


def _first_present(data: Dict[str, Optional[float]], keys: Sequence[str]) -> Optional[float]:
    for key in keys:
        value = data.get(key)
        if value is not None:
            return value
    return None


def augment_row(row: "Row") -> None:
    data = row.data

    tip_defl = _first_present(data, ("tip_defl_m", "tip_defl"))
    sigma_max = _first_present(data, ("sigma_max_pa", "sigma_ma"))
    buckling_factor = _first_present(data, ("buckling_factor", "buckling"))

    load_face = _first_present(data, ("axial_load_face_N",))
    load_pair = _first_present(data, ("load_pair_total_N",))

    if load_face is None and load_pair is not None:
        load_face = 0.5 * load_pair
    if load_pair is None and load_face is not None:
        load_pair = 2.0 * load_face

    data["axial_load_face_N"] = load_face
    data["load_pair_total_N"] = load_pair

    buckling_pair = _first_present(data, ("buckling_pair_N",))
    if buckling_pair is None and buckling_factor is not None and load_pair is not None:
        buckling_pair = buckling_factor * load_pair
    buckling_face = _first_present(data, ("buckling_face_N",))
    if buckling_face is None and buckling_pair is not None:
        buckling_face = 0.5 * buckling_pair
    data["buckling_pair_N"] = buckling_pair
    data["buckling_face_N"] = buckling_face

    mass = _first_present(data, ("total_mass_kg",))
    R = data.get("R_m")
    L = data.get("L_m")
    a = data.get("a_m")
    b = data.get("b_m")
    t = data.get("t_m")

    if mass is None and None not in (R, L, a, b, t):
        material_key = row.text("material") or row.text("mat_code")
        try:
            material = resolve_material(material_key)
        except Exception:
            material = None
        if material is not None:
            mass = total_mass(material, R, L, b, t, a)
    data["total_mass_kg"] = mass

    if buckling_pair is not None and mass not in (None, 0.0):
        data["buckling_per_mass"] = buckling_pair / mass
    else:
        data.setdefault("buckling_per_mass", None)

    load_ref = None
    if load_pair is not None and abs(load_pair) > 0.0:
        load_ref = abs(load_pair)

    tip_per_load = _first_present(data, ("tip_per_load", "tip_per_load_m_per_N"))
    sigma_per_load = _first_present(data, ("sigma_per_load", "sigma_per_load_pa_per_N"))

    if tip_per_load is None and tip_defl is not None and load_ref:
        tip_per_load = tip_defl / load_ref
    if sigma_per_load is None and sigma_max is not None and load_ref:
        sigma_per_load = sigma_max / load_ref

    data["tip_per_load"] = tip_per_load
    data["sigma_per_load"] = sigma_per_load


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--summary",
        type=Path,
        default=DEFAULT_SUMMARY,
        help="Path to the consolidated summary CSV.",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=DEFAULT_REPORT,
        help="Destination path for the Markdown report.",
    )
    parser.add_argument(
        "--top",
        type=int,
        default=10,
        help="Number of high-scoring rows to include in the main table.",
    )
    parser.add_argument(
        "--weight",
        action="append",
        default=[],
        metavar="METRIC=VALUE",
        help="Override default metric weights, e.g. --weight tip_defl=0.5",
    )
    return parser.parse_args(argv)


def parse_weight_overrides(weight_flags: Iterable[str]) -> Dict[str, float]:
    weights = {name: cfg["default_weight"] for name, cfg in METRICS.items()}
    if not weight_flags:
        return weights
    for flag in weight_flags:
        if "=" not in flag:
            raise ValueError(f"Invalid weight override '{flag}'. Use name=value.")
        name, raw_value = flag.split("=", 1)
        name = name.strip()
        if name not in METRICS:
            known = ", ".join(sorted(METRICS))
            raise ValueError(f"Unknown metric '{name}'. Known metrics: {known}.")
        try:
            weight = float(raw_value)
        except ValueError as exc:
            raise ValueError(f"Weight for '{name}' must be numeric") from exc
        if weight < 0:
            raise ValueError(f"Weight for '{name}' must be >= 0.")
        weights[name] = weight
    total = sum(weights.values())
    if total <= 0:
        raise ValueError("Sum of metric weights must be > 0.")
    return weights


def coerce_float(value: str) -> Optional[float]:
    value = value.strip()
    if not value or set(value) == {"*"}:
        return None
    try:
        return float(value)
    except ValueError:
        return None


def load_rows(path: Path) -> List[Row]:
    if not path.exists():
        raise FileNotFoundError(f"Summary CSV '{path}' does not exist.")
    rows: List[Row] = []
    with path.open(newline="") as fp:
        data_lines = (line for line in fp if not line.lstrip().startswith("#"))
        reader = csv.DictReader(data_lines)
        for raw_row in reader:
            raw = {key: (value.strip() if value is not None else "") for key, value in raw_row.items()}
            coerced: Dict[str, Optional[float]] = {}
            for key, value in raw.items():
                if key in NUMERIC_COLUMNS:
                    coerced[key] = coerce_float(value)
            row = Row(data=coerced, raw=raw)
            augment_row(row)
            rows.append(row)
    return rows


def compute_metric_stats(rows: Iterable[Row]) -> Dict[str, MetricStats]:
    stats: Dict[str, MetricStats] = {}
    for name, cfg in METRICS.items():
        values = [row.value(name) for row in rows if row.value(name) is not None]
        if not values:
            continue
        stats[name] = MetricStats(
            name=name,
            higher_is_better=cfg["higher_is_better"],
            minimum=min(values),
            maximum=max(values),
        )
    return stats


def assign_scores(rows: Iterable[Row], weights: Dict[str, float], stats: Dict[str, MetricStats]) -> None:
    for row in rows:
        weighted_sum = 0.0
        weight_total = 0.0
        for name, weight in weights.items():
            metric_stats = stats.get(name)
            value = row.value(name)
            if metric_stats is None or value is None:
                continue
            weighted_sum += weight * metric_stats.normalised(value)
            weight_total += weight
        if weight_total > 0:
            row.score = weighted_sum / weight_total
        else:
            row.score = None


def sorted_rows(rows: Iterable[Row]) -> List[Row]:
    eligible = [row for row in rows if row.score is not None]
    eligible.sort(key=lambda row: row.score if row.score is not None else -1.0, reverse=True)
    return eligible


def format_number(value: Optional[float], precision: int = 3) -> str:
    if value is None:
        return "–"
    if abs(value) >= 1000 or (0 < abs(value) < 10 ** -(precision - 1)):
        return f"{value:.{precision}e}"
    return f"{value:.{precision}f}"


def render_table(rows: Sequence[Tuple[int, Row]]) -> str:
    columns = [
        ("Rank", lambda idx, row: str(idx)),
        ("Case", lambda idx, row: row.text("case_id")),
        ("Material", lambda idx, row: row.text("material") or row.text("mat_code")),
        ("t [m]", lambda idx, row: format_number(row.value("t_m"), precision=4)),
        ("a [m]", lambda idx, row: format_number(row.value("a_m"), precision=4)),
        ("Mass [kg]", lambda idx, row: format_number(row.value("total_mass_kg"), precision=4)),
        ("Load [N]", lambda idx, row: format_number(row.value("load_pair_total_N"), precision=4)),
        ("Buckling/Mass", lambda idx, row: format_number(row.value("buckling_per_mass"), precision=4)),
        ("Tip/Load", lambda idx, row: format_number(row.value("tip_per_load"), precision=4)),
        ("Sigma/Load", lambda idx, row: format_number(row.value("sigma_per_load"), precision=4)),
        ("Score", lambda idx, row: format_number(row.score, precision=4)),
    ]
    header = " | ".join(name for name, _ in columns)
    divider = " | ".join("---" for _ in columns)
    lines = [f"| {header} |", f"| {divider} |"]
    for rank, row in rows:
        cells = [func(rank, row) for _, func in columns]
        lines.append("| " + " | ".join(cells) + " |")
    return "\n".join(lines)


def render_metric_highlights(rows: Sequence[Row], stats: Dict[str, MetricStats]) -> str:
    lines: List[str] = []
    for name, cfg in METRICS.items():
        metric_stats = stats.get(name)
        if metric_stats is None:
            lines.append(f"- **{cfg['label']}**: no valid data available.")
            continue
        if cfg["higher_is_better"]:
            ordered = sorted(
                (row for row in rows if row.value(name) is not None),
                key=lambda row: row.value(name),
                reverse=True,
            )
        else:
            ordered = sorted(
                (row for row in rows if row.value(name) is not None),
                key=lambda row: row.value(name),
            )
        top = ordered[:5]
        if not top:
            lines.append(f"- **{cfg['label']}**: no valid data available.")
            continue
        best_lines = []
        for row in top:
            best_lines.append(
                f"{row.text('case_id')} ({row.text('material')}) -> {format_number(row.value(name), precision=4)}"
            )
        lines.append(f"- **{cfg['label']}**: " + "; ".join(best_lines))
    return "\n".join(lines)


def render_stats_section(stats: Dict[str, MetricStats]) -> str:
    lines = []
    for name, cfg in METRICS.items():
        metric_stats = stats.get(name)
        if metric_stats is None:
            lines.append(f"- **{cfg['label']}**: no valid samples.")
            continue
        lines.append(
            f"- **{cfg['label']}**: min {format_number(metric_stats.minimum, 4)}, "
            f"max {format_number(metric_stats.maximum, 4)}"
        )
    return "\n".join(lines)


def generate_report(
    summary_path: Path,
    rows: Sequence[Row],
    ranked: Sequence[Row],
    top_n: int,
    stats: Dict[str, MetricStats],
    weights: Dict[str, float],
) -> str:
    total_rows = len(rows)
    ranked_rows = len(ranked)
    selected = list(enumerate(ranked, start=1))[:top_n]

    lines: List[str] = []
    lines.append("# MAPDL Run Summary\n")
    lines.append(f"- Source file: `{summary_path}`")
    lines.append(f"- Total rows parsed: {total_rows}")
    lines.append(f"- Rows with complete score: {ranked_rows}")
    weight_text = ", ".join(f"{name}={weights[name]:.2f}" for name in sorted(weights))
    lines.append(f"- Weights: {weight_text}\n")

    lines.append("## Top Results\n")
    if not selected:
        lines.append("No rows contained enough data to compute a score.")
    else:
        lines.append(render_table(selected))
    lines.append("")

    lines.append("## Metric Highlights\n")
    lines.append(render_metric_highlights(ranked, stats))
    lines.append("")

    lines.append("## Metric Ranges\n")
    lines.append(render_stats_section(stats))
    lines.append("")

    return "\n".join(lines).rstrip() + "\n"


def write_report(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")


def main(argv: Optional[Sequence[str]] = None) -> None:
    args = parse_args(argv)
    weights = parse_weight_overrides(args.weight)
    rows = load_rows(args.summary)
    stats = compute_metric_stats(rows)
    if not stats:
        raise SystemExit("No usable metric data found. Did the runs finish correctly?")
    assign_scores(rows, weights, stats)
    ranked = sorted_rows(rows)
    report = generate_report(args.summary, rows, ranked, args.top, stats, weights)
    write_report(args.out, report)
    print(f"Wrote report with {len(ranked)} scored rows to {args.out}")


if __name__ == "__main__":  # pragma: no cover
    main()
