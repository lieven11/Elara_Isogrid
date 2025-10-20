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
    *,
    theta_fraction: float = 1.0,
    height_fraction: float = 1.0,
) -> Tuple[List[Vector3], Dict[str, List[Segment]]]:
    theta_fraction = max(1e-6, min(theta_fraction, 1.0))
    height_fraction = max(1e-6, min(height_fraction, 1.0))
    effective_length = axial_length * height_fraction
    dz = delta_z_from_a(a)

    rows_z: List[float] = []
    i = 0
    while True:
        z = i * dz
        if z > effective_length + 1e-9:
            break
        rows_z.append(z)
        i += 1
    if not rows_z or rows_z[-1] < effective_length - 1e-9:
        rows_z.append(effective_length)

    nodes: List[Vector3] = []
    node_index: Dict[Tuple[int, int], int] = {}

    theta_step = a / max(R, 1e-12)
    n_cols = max(2, int(math.ceil(n_theta * theta_fraction)))
    theta_min = float("inf")
    theta_max = float("-inf")

    for r, z_val in enumerate(rows_z):
        offset = 0.0 if r % 2 == 0 else 0.5
        for c in range(n_cols):
            col_pos = c + offset
            idx = len(nodes)
            theta = col_pos * theta_step
            nodes.append(
                (
                    R * math.cos(theta),
                    R * math.sin(theta),
                    z_val,
                )
            )
            node_index[(r, c)] = idx
            theta_min = min(theta_min, theta)
            theta_max = max(theta_max, theta)

    segments: Dict[str, List[Segment]] = {"0": [], "+60": [], "-60": []}

    for r, z_val in enumerate(rows_z):
        for c in range(n_cols):
            here = node_index[(r, c)]

            if r + 1 >= len(rows_z):
                upper_exists = False
            else:
                upper_exists = True

            if c + 1 < n_cols:
                right = node_index[(r, c + 1)]
                segments["0"].append((here, right))

            if not upper_exists:
                continue

            if r % 2 == 0:
                up_right = node_index[(r + 1, c)]
                segments["+60"].append((here, up_right))
                if c - 1 >= 0:
                    up_left = node_index[(r + 1, c - 1)]
                    segments["-60"].append((here, up_left))
            else:
                up_left = node_index[(r + 1, c)]
                segments["-60"].append((here, up_left))
                if c + 1 < n_cols:
                    up_right = node_index[(r + 1, c + 1)]
                    segments["+60"].append((here, up_right))

    if not math.isfinite(theta_min) or not math.isfinite(theta_max):
        theta_min, theta_max = 0.0, 2.0 * math.pi

    return nodes, segments, (theta_min, theta_max)


def draw_ring(
    ax,
    radius: float,
    z_values: Iterable[float],
    *,
    color: str,
    linewidth: float = 1.2,
    inner_radius: float | None = None,
    theta_limits: Tuple[float, float] | None = None,
) -> None:
    if theta_limits is None:
        theta_limits = (0.0, 2.0 * math.pi)
    theta = np.linspace(theta_limits[0], theta_limits[1], 240)
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
    theta_limits: Tuple[float, float] | None = None,
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

    if theta_limits is not None:
        theta_min, theta_max = theta_limits
    else:
        theta_values = np.array([math.atan2(y, x) for x, y, _ in nodes])
        if theta_values.size == 0:
            theta_min, theta_max = 0.0, 2.0 * math.pi
        else:
            theta_unwrapped = np.unwrap(theta_values)
            theta_min = float(theta_unwrapped.min())
            theta_max = float(theta_unwrapped.max())
    if theta_max < theta_min:
        theta_min, theta_max = theta_max, theta_min

    ring_color = "#444444"
    ring_top = ring_bottom + ring_height
    inner_radius = max(radius - wall_thickness, 0.0)
    bottom_ring_top = 0.0
    bottom_ring_bottom = -ring_height
    draw_ring(
        ax,
        radius,
        [ring_bottom, ring_top],
        color=ring_color,
        linewidth=1.6,
        inner_radius=inner_radius,
        theta_limits=(theta_min, theta_max),
    )
    draw_ring(
        ax,
        radius,
        [bottom_ring_bottom, bottom_ring_top],
        color=ring_color,
        linewidth=1.6,
        inner_radius=inner_radius,
        theta_limits=(theta_min, theta_max),
    )

    support_angles = np.linspace(theta_min, theta_max, 4)
    for angle in support_angles:
        x_outer = radius * math.cos(angle)
        y_outer = radius * math.sin(angle)
        ax.plot(
            [x_outer, x_outer],
            [y_outer, y_outer],
            [bottom_ring_top, ring_top],
            color=ring_color,
            linewidth=1.0,
        )
        if inner_radius > 0.0:
            x_inner = inner_radius * math.cos(angle)
            y_inner = inner_radius * math.sin(angle)
            ax.plot(
                [x_inner, x_inner],
                [y_inner, y_inner],
                [bottom_ring_top, ring_top],
                color=ring_color,
                linewidth=0.9,
                linestyle="--",
            )

    axis_extent = radius * 1.2
    z_axis_top = ring_top * 1.1
    z_axis_bottom = bottom_ring_bottom * 1.1
    ax.plot([0.0, axis_extent], [0.0, 0.0], [0.0, 0.0], color="k", linewidth=1.0)
    ax.plot([0.0, 0.0], [0.0, axis_extent], [0.0, 0.0], color="k", linewidth=1.0)
    ax.plot([0.0, 0.0], [0.0, 0.0], [z_axis_bottom, z_axis_top], color="k", linewidth=1.0)
    ax.text(axis_extent, 0.0, 0.0, "X", fontsize=10)
    ax.text(0.0, axis_extent, 0.0, "Y", fontsize=10)
    ax.text(0.0, 0.0, z_axis_top, "Z", fontsize=10)
    ax.text(0.0, 0.0, z_axis_bottom, "-Z", fontsize=10)

    arrow_length = radius * 0.25
    inner_radius = max(radius - wall_thickness, 0.0)
    theta_samples = np.linspace(theta_min, theta_max, 120)
    if np.isclose(inner_radius, radius):
        radii = np.array([0.0, radius])
    else:
        radii = np.array([inner_radius, radius])
    theta_grid, radius_grid = np.meshgrid(theta_samples, radii)
    x_grid = radius_grid * np.cos(theta_grid)
    y_grid = radius_grid * np.sin(theta_grid)
    z_top = np.full_like(x_grid, ring_top)
    if radius > 0.0:
        ax.plot_surface(
            x_grid,
            y_grid,
            z_top,
            color="#d62728",
            alpha=0.18,
            linewidth=0.0,
            antialiased=False,
        )

    load_theta = 0.5 * (theta_min + theta_max)
    ring_mid_radius = (inner_radius + radius) * 0.5 if radius > 0.0 else 0.0
    load_origin = (
        ring_mid_radius * math.cos(load_theta),
        ring_mid_radius * math.sin(load_theta),
        ring_top + arrow_length * 0.6,
    )
    load_direction = (0.0, 0.0, -arrow_length)
    arrow_angles = np.linspace(theta_min, theta_max, 8)
    for angle in arrow_angles:
        origin = (
            ring_mid_radius * math.cos(angle),
            ring_mid_radius * math.sin(angle),
            ring_top + arrow_length * 0.4,
        )
        ax.quiver(
            origin[0],
            origin[1],
            origin[2],
            load_direction[0],
            load_direction[1],
            load_direction[2],
            color="#d62728",
            linewidth=1.6,
            arrow_length_ratio=0.2,
        )
    text_offset = arrow_length * 0.2
    ax.text(
        load_origin[0],
        load_origin[1],
        load_origin[2] + load_direction[2] - text_offset,
        "Axial load on annulus",
        color="#d62728",
        fontsize=10,
        ha="center",
    )

    z_bottom = np.full_like(x_grid, bottom_ring_bottom)
    if radius > 0.0:
        ax.plot_surface(
            x_grid,
            y_grid,
            z_bottom,
            color="#333333",
            alpha=0.12,
            linewidth=0.0,
            antialiased=False,
        )
    base_label_radius = ring_mid_radius
    ax.text(
        base_label_radius * math.cos(load_theta),
        base_label_radius * math.sin(load_theta),
        bottom_ring_bottom - arrow_length * 0.3,
        "Fixed base annulus",
        color="#333333",
        fontsize=9,
        ha="center",
    )

    ax.set_xlabel("X [m]")
    ax.set_ylabel("Y [m]")
    ax.set_zlabel("Z [m]")
    ax.set_title("Triangular Isogrid (0°) Preview")

    set_equal_axes(
        ax,
        nodes
        + [
            (radius * math.cos(theta_min), radius * math.sin(theta_min), ring_top),
            (radius * math.cos(theta_max), radius * math.sin(theta_max), ring_top),
            (radius * math.cos(theta_min), radius * math.sin(theta_min), bottom_ring_bottom),
            (radius * math.cos(theta_max), radius * math.sin(theta_max), bottom_ring_bottom),
        ],
    )

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
    theta_fraction = 0.25
    height_fraction = 0.5
    nodes, segments, theta_limits = build_triangular_lattice(
        R,
        a,
        n_theta,
        axial_length,
        theta_fraction=theta_fraction,
        height_fraction=height_fraction,
    )
    isogrid_height = axial_length * height_fraction

    plot_preview(
        nodes,
        segments,
        ring_bottom=isogrid_height,
        ring_height=END_RING_HEIGHT,
        radius=R,
        wall_thickness=wall_thickness,
        save_path=args.save,
        theta_limits=theta_limits,
    )


if __name__ == "__main__":
    main()
