#!/usr/bin/env python3
"""Static 3D preview for the 0° triangular isogrid geometry."""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import matplotlib.pyplot as plt
import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from isogrid.geometry import END_RING_HEIGHT, delta_z_from_a, n_theta_from_a


Vector3 = Tuple[float, float, float]
Segment = Tuple[int, int]


def _default_from_param(param: Dict[str, float]) -> float:
    if "value" in param:
        return float(param["value"])
    if "values" in param and param["values"]:
        return float(param["values"][0])
    if "start" in param:
        return float(param["start"])
    raise ValueError("Parameter dictionary lacks a default-like entry.")


def load_spec_defaults(spec_path: Path) -> Dict[str, float]:
    data = json.loads(spec_path.read_text(encoding="utf-8"))
    params = data.get("parameters", {})
    defaults: Dict[str, float] = {}
    for key in ("R", "b", "t", "a", "n_theta"):
        if key in params:
            defaults[key] = _default_from_param(params[key])
    return defaults


def build_triangular_lattice(
    R: float,
    a: float,
    n_theta: int,
    axial_length: float,
) -> Tuple[List[Vector3], Dict[str, List[Segment]]]:
    width = n_theta * a
    dz = delta_z_from_a(a)

    rows_z: List[float] = []
    i = 0
    while True:
        z = i * dz
        if z > axial_length + 1e-9:
            break
        rows_z.append(z)
        i += 1
    if not rows_z or rows_z[-1] < axial_length - 1e-9:
        rows_z.append(axial_length)

    nodes: List[Vector3] = []
    node_index: Dict[Tuple[int, int], int] = {}

    def to_xyz(x_unwrapped: float, z_val: float) -> Vector3:
        theta = (x_unwrapped % width) / max(R, 1e-12)
        return (
            R * math.cos(theta),
            R * math.sin(theta),
            z_val,
        )

    for r, z_val in enumerate(rows_z):
        offset = 0.0 if r % 2 == 0 else 0.5
        for c in range(n_theta):
            x_unwrapped = (c + offset) * a
            idx = len(nodes)
            nodes.append(to_xyz(x_unwrapped, z_val))
            node_index[(r, c)] = idx

    segments: Dict[str, List[Segment]] = {"0": [], "+60": [], "-60": []}

    for r, z_val in enumerate(rows_z):
        for c in range(n_theta):
            here = node_index[(r, c)]
            right = node_index[(r, (c + 1) % n_theta)]
            segments["0"].append((here, right))

            if r + 1 >= len(rows_z):
                continue

            if r % 2 == 0:
                up_right = node_index[(r + 1, c)]
                up_left = node_index[(r + 1, (c - 1) % n_theta)]
                segments["+60"].append((here, up_right))
                segments["-60"].append((here, up_left))
            else:
                up_left = node_index[(r + 1, c)]
                up_right = node_index[(r + 1, (c + 1) % n_theta)]
                segments["-60"].append((here, up_left))
                segments["+60"].append((here, up_right))

    return nodes, segments


def draw_ring(
    ax,
    radius: float,
    z_values: Iterable[float],
    *,
    color: str,
    linewidth: float = 1.2,
    inner_radius: float | None = None,
) -> None:
    theta = np.linspace(0.0, 2.0 * math.pi, 240)
    for z in z_values:
        x_outer = radius * np.cos(theta)
        y_outer = radius * np.sin(theta)
        z_line = np.full_like(theta, z)
        ax.plot(x_outer, y_outer, z_line, color=color, linewidth=linewidth)
        if inner_radius is not None and inner_radius > 0.0:
            x_inner = inner_radius * np.cos(theta)
            y_inner = inner_radius * np.sin(theta)
            ax.plot(
                x_inner,
                y_inner,
                z_line,
                color=color,
                linewidth=linewidth * 0.8,
                linestyle="--",
            )


def set_equal_axes(ax, nodes: Iterable[Vector3], extra: float = 0.05) -> None:
    xs, ys, zs = zip(*nodes)
    x_min, x_max = min(xs), max(xs)
    y_min, y_max = min(ys), max(ys)
    z_min, z_max = min(zs), max(zs)

    span_x = x_max - x_min
    span_y = y_max - y_min
    span_z = z_max - z_min
    span = max(span_x, span_y, span_z)
    pad = span * extra

    x_mid = (x_min + x_max) / 2.0
    y_mid = (y_min + y_max) / 2.0
    z_mid = (z_min + z_max) / 2.0

    ax.set_xlim(x_mid - span / 2.0 - pad, x_mid + span / 2.0 + pad)
    ax.set_ylim(y_mid - span / 2.0 - pad, y_mid + span / 2.0 + pad)
    ax.set_zlim(z_mid - span / 2.0 - pad, z_mid + span / 2.0 + pad)


def plot_preview(
    nodes: List[Vector3],
    segments: Dict[str, List[Segment]],
    *,
    ring_bottom: float,
    ring_height: float,
    radius: float,
    wall_thickness: float,
    save_path: Path | None,
) -> None:
    fig = plt.figure(figsize=(9, 7))
    ax = fig.add_subplot(111, projection="3d")

    family_colors = {"0": "#1f77b4", "+60": "#ff7f0e", "-60": "#2ca02c"}
    family_labels = {"0": "0° hoops", "+60": "+60° ribs", "-60": "-60° ribs"}

    for family, segs in segments.items():
        color = family_colors.get(family, "#444444")
        label = family_labels.get(family)
        used_label = label
        for idx0, idx1 in segs:
            x0, y0, z0 = nodes[idx0]
            x1, y1, z1 = nodes[idx1]
            ax.plot(
                [x0, x1],
                [y0, y1],
                [z0, z1],
                color=color,
                linewidth=1.0,
                alpha=0.85,
                label=used_label,
            )
            used_label = None

    ring_color = "#444444"
    ring_top = ring_bottom + ring_height
    inner_radius = max(radius - wall_thickness, 0.0)
    draw_ring(ax, radius, [ring_bottom, ring_top], color=ring_color, linewidth=1.6, inner_radius=inner_radius)

    for angle in np.linspace(0.0, 2.0 * math.pi, 8, endpoint=False):
        x_outer = radius * math.cos(angle)
        y_outer = radius * math.sin(angle)
        ax.plot(
            [x_outer, x_outer],
            [y_outer, y_outer],
            [ring_bottom, ring_top],
            color=ring_color,
            linewidth=1.0,
        )
        if inner_radius > 0.0:
            x_inner = inner_radius * math.cos(angle)
            y_inner = inner_radius * math.sin(angle)
            ax.plot(
                [x_inner, x_inner],
                [y_inner, y_inner],
                [ring_bottom, ring_top],
                color=ring_color,
                linewidth=0.9,
                linestyle="--",
            )

    axis_extent = radius * 1.2
    ax.plot([0.0, axis_extent], [0.0, 0.0], [0.0, 0.0], color="k", linewidth=1.0)
    ax.plot([0.0, 0.0], [0.0, axis_extent], [0.0, 0.0], color="k", linewidth=1.0)
    ax.plot([0.0, 0.0], [0.0, 0.0], [0.0, ring_top * 1.1], color="k", linewidth=1.0)
    ax.text(axis_extent, 0.0, 0.0, "X", fontsize=10)
    ax.text(0.0, axis_extent, 0.0, "Y", fontsize=10)
    ax.text(0.0, 0.0, ring_top * 1.1, "Z", fontsize=10)

    ax.set_xlabel("X [m]")
    ax.set_ylabel("Y [m]")
    ax.set_zlabel("Z [m]")
    ax.set_title("Triangular Isogrid (0°) Preview")

    set_equal_axes(ax, nodes + [
        (radius, 0.0, ring_top),
        (-radius, 0.0, ring_top),
        (0.0, radius, ring_top),
        (0.0, -radius, ring_top),
    ])

    handles, labels = ax.get_legend_handles_labels()
    if labels:
        ax.legend(loc="upper right")

    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches="tight")
    else:
        plt.show()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    default_spec = Path(__file__).resolve().parents[1] / "specs" / "tri0_dense_template.json"
    parser.add_argument(
        "--spec",
        type=Path,
        default=default_spec,
        help="Path to a tri0 spec JSON used to grab representative parameters.",
    )
    parser.add_argument(
        "--length",
        type=float,
        default=1.0,
        help="Axial preview length in meters (default: 1.0 m).",
    )
    parser.add_argument(
        "--save",
        type=Path,
        help="Optional path to save the rendered preview instead of showing a window.",
    )
    args = parser.parse_args()

    defaults = load_spec_defaults(args.spec)
    R = defaults.get("R", 0.55)
    a = defaults.get("a")
    n_theta = defaults.get("n_theta")

    if a is None and n_theta is None:
        raise ValueError("Spec must provide either 'a' or 'n_theta'.")
    if a is None:
        n_theta = int(max(1, round(n_theta)))
        a = 2.0 * math.pi * R / n_theta
    else:
        n_theta = n_theta_from_a(R, a)

    wall_thickness = defaults.get("t", 0.002)
    axial_length = max(args.length, 1e-6)
    nodes, segments = build_triangular_lattice(R, a, n_theta, axial_length)

    plot_preview(
        nodes,
        segments,
        ring_bottom=axial_length,
        ring_height=END_RING_HEIGHT,
        radius=R,
        wall_thickness=wall_thickness,
        save_path=args.save,
    )


if __name__ == "__main__":
    main()
