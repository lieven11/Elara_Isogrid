import math
from dataclasses import dataclass

from .materials import IsotropicMaterial
from .geometry import (
    a_from_n_theta,
    delta_z_from_a,
    line_length_density_total,
    line_length_density_per_family,
    surface_area,
)


SQRT3 = math.sqrt(3.0)


@dataclass(frozen=True)
class GlobalAssumptions:
    K_euler: float = 1.0  # Effective length factor (pinned-pinned default)
    k_crippling: float = 4.0  # Crippling coefficient (conservative sheet-like)
    fos_yield: float = 1.25   # Factor of safety for yield allowables
    k_torsion: float = 0.2    # Reduction factor for torsional rigidity of lattice (empirical)


def areal_mass(material: IsotropicMaterial, b: float, t: float, a: float) -> float:
    """Areal mass density m̄ [kg/m^2] for the lattice on the cylindrical surface."""
    l_tot = line_length_density_total(a)
    return material.rho * (b * t) * l_tot


def total_mass(material: IsotropicMaterial, R: float, L: float, b: float, t: float, a: float) -> float:
    return areal_mass(material, b, t, a) * surface_area(R, L)


def EI_equivalent(material: IsotropicMaterial, R: float, b: float, t: float, a: float) -> float:
    """Equivalent bending stiffness EI about a lateral axis for global Euler buckling.

    Derived via energy equivalence assuming uniform lattice and pin-jointed bars.
    Result: EI_eq = π * E * (b t) * (√3/a) * R^3
    """
    return math.pi * material.E * (b * t) * (SQRT3 / a) * (R ** 3)


def EA_ring(material: IsotropicMaterial, R: float, b: float, t: float, a: float) -> float:
    """Equivalent axial stiffness of the ring (per unit length in z):

    EA_ring = π * E * (b t) * (√3/a) * R
    (Consistent with EI_eq = EA_ring * R^2)
    """
    return math.pi * material.E * (b * t) * (SQRT3 / a) * R


def Ncr_global(material: IsotropicMaterial, R: float, L: float, b: float, t: float, a: float, K: float = 1.0) -> float:
    EI = EI_equivalent(material, R, b, t, a)
    return (math.pi ** 2) * EI / (K * L) ** 2


def GJ_equivalent(material: IsotropicMaterial, R: float, b: float, t: float, a: float, k_torsion: float = 0.2) -> float:
    """Equivalent torsional rigidity GJ for the lattice.

    Approximation: start from EA_ring scaling and convert to torsion via R^2 lever arm,
    reduced by k_torsion to reflect poor torsional continuity of a skinless lattice.
    Using G = E/(2(1+nu)).
    GJ_eq ≈ k_torsion * π * G * (b t) * (√3/a) * R^3
    """
    G = material.E / (2.0 * (1.0 + material.nu))
    return k_torsion * math.pi * G * (b * t) * (SQRT3 / a) * (R ** 3)


def bending_global_sf(material: IsotropicMaterial, R: float, b: float, t: float, a: float, M_req: float, assumptions: GlobalAssumptions) -> float:
    """Global bending safety factor from yield strain at the outer fiber.

    curvature kappa = M / EI_eq; strain at outer radius ε = kappa * R;
    stress E*ε must be <= σ_allow.
    M_cap = σ_allow * EI_eq / R; SF = M_cap / M_req.
    """
    if M_req is None or M_req <= 0:
        return float("nan")
    EI = EI_equivalent(material, R, b, t, a)
    sigma_allow = material.sigma_y / assumptions.fos_yield
    M_cap = sigma_allow * EI / max(R, 1e-12)
    return M_cap / M_req


def torsion_sf(material: IsotropicMaterial, R: float, b: float, t: float, a: float, T_req: float, assumptions: GlobalAssumptions) -> float:
    """Torsional safety factor based on shear yield proxy.

    tau_eq ≈ T * R / J_eq; tau_allow ≈ 0.58 * σ_allow.
    SF = tau_allow / tau_eq.
    """
    if T_req is None or T_req <= 0:
        return float("nan")
    J = GJ_equivalent(material, R, b, t, a, assumptions.k_torsion)
    tau_eq = T_req * R / max(J, 1e-18)
    sigma_allow = material.sigma_y / assumptions.fos_yield
    tau_allow = 0.58 * sigma_allow
    return tau_allow / max(tau_eq, 1e-18)


def rectangular_I_min(b: float, t: float) -> float:
    """Weak-axis second moment for a rectangular rib (buckling about thickness)."""
    return b * (t ** 3) / 12.0


def euler_crit_stress(material: IsotropicMaterial, length: float, r_g: float, K: float = 1.0) -> float:
    """Euler critical stress for a slender member: σ_cr = π^2 E / (λ^2),
    where λ = (K * L) / r_g is slenderness.
    """
    lam = (K * length) / max(r_g, 1e-12)
    if lam <= 0.0:
        return float("inf")
    return (math.pi ** 2) * material.E / (lam ** 2)


def crippling_stress(material: IsotropicMaterial, b: float, t: float, k: float = 4.0) -> float:
    """Conservative plate-like crippling stress for a slender rectangular strip.
    σ_cr ≈ k * π^2 E / (12(1-ν^2)) * (t/b)^2
    """
    denom = 12.0 * (1.0 - material.nu ** 2)
    return k * (math.pi ** 2) * material.E / denom * (t / max(b, 1e-12)) ** 2


def local_safety_factors(
    material: IsotropicMaterial,
    b: float,
    t: float,
    a: float,
    L: float,
    R: float,
    N_req: float,
    assumptions: GlobalAssumptions,
) -> dict:
    """Compute local safety factors for 0°, ±60° members against an applied axial load N_req.

    Approach: infer macroscopic axial strain ε from N_req using EA_ring, then
    bar axial stress = E * ε * cos^2(alpha). Compare to min of Euler/Crippling/Yield.
    """
    # Global axial strain from requested axial force
    EA = EA_ring(material, R, b, t, a)
    eps = N_req / max(EA, 1e-12)

    # Bar axial stress demand by orientation
    cos2 = {
        0.0: 1.0,
        60.0: 0.25,
        -60.0: 0.25,
    }
    sigma_demand = {ang: material.E * eps * c2 for ang, c2 in cos2.items()}

    # Critical stresses
    Imin = rectangular_I_min(b, t)
    A = b * t
    r_g = math.sqrt(max(Imin / max(A, 1e-12), 1e-18))

    dz = delta_z_from_a(a)
    lengths = {0.0: dz, 60.0: a, -60.0: a}

    sigma_allow = material.sigma_y / assumptions.fos_yield
    sigma_crip = crippling_stress(material, b, t, assumptions.k_crippling)

    sf = {}
    for ang, Lseg in lengths.items():
        sigma_euler = euler_crit_stress(material, Lseg, r_g, assumptions.K_euler)
        sigma_crit = min(sigma_euler, sigma_crip, sigma_allow)
        dem = sigma_demand[ang]
        sf[ang] = float("inf") if dem <= 0 else sigma_crit / dem

    sf_min = min(sf.values()) if sf else float("inf")
    return {"sf_per_angle": sf, "sf_min": sf_min}


def metric_mass_specific_capacity(material: IsotropicMaterial, R: float, L: float, b: float, t: float, a: float, K: float = 1.0) -> float:
    """N_cr,global / total_mass — mass-specific global buckling capacity."""
    ncr = Ncr_global(material, R, L, b, t, a, K)
    m = total_mass(material, R, L, b, t, a)
    return ncr / max(m, 1e-12)
