import math
from dataclasses import dataclass


SQRT3 = math.sqrt(3.0)


@dataclass(frozen=True)
class LatticeParams:
    R: float          # Cylinder outer radius [m]
    L: float          # Segment length [m]
    b: float          # Rib width (tangential) [m]
    t: float          # Rib thickness (radial) [m]
    n_theta: int      # Number of cells around circumference [-]


def a_from_n_theta(R: float, n_theta: int) -> float:
    """Equilateral cell edge length a from circumferential cell count.

    a = 2*pi*R / n_theta (on the unwrapped cylinder)
    """
    return 2.0 * math.pi * R / max(1, n_theta)


def n_theta_from_a(R: float, a: float) -> int:
    """Circumferential cell count from desired edge length a (rounded to >=1)."""
    if a <= 0:
        return 1
    return max(1, int(round(2.0 * math.pi * R / a)))


def delta_z_from_a(a: float) -> float:
    """Axial node spacing Δz for 60° triangular lattice."""
    return (SQRT3 / 2.0) * a


def line_length_density_total(a: float) -> float:
    """Total bar length per unit surface area of the triangular lattice.

    For an equilateral triangular lattice: L_total / Area = 2*sqrt(3) / a
    Proof via edge sharing over equilateral triangles.
    """
    return 2.0 * SQRT3 / a


def line_length_density_per_family(a: float) -> float:
    """Bar length per unit area for one orientation family (0°, +60°, -60°)."""
    return line_length_density_total(a) / 3.0


END_RING_HEIGHT = 0.04  # m, fixed axial height of continuous end-ring segments


def annulus_area(outer_radius: float, wall_thickness: float) -> float:
    """Cross-sectional area of a thin cylindrical ring (outer - inner circle).

    Args:
        outer_radius: Outer radius of the ring [m].
        wall_thickness: Radial wall thickness [m].
    """
    inner_radius = max(outer_radius - wall_thickness, 0.0)
    return math.pi * max(outer_radius ** 2 - inner_radius ** 2, 0.0)


def end_ring_volume(outer_radius: float, wall_thickness: float, height: float = END_RING_HEIGHT) -> float:
    """Volume of a continuous cylindrical ring segment at one end of the shell."""
    return annulus_area(outer_radius, wall_thickness) * max(height, 0.0)


def surface_area(R: float, L: float) -> float:
    return 2.0 * math.pi * R * L
