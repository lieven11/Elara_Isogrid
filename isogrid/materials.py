from dataclasses import dataclass


@dataclass(frozen=True)
class IsotropicMaterial:
    name: str
    E: float        # Young's modulus [Pa]
    nu: float       # Poisson's ratio [-]
    rho: float      # Density [kg/m^3]
    sigma_y: float  # Yield strength [Pa]
    sigma_ult: float  # Ultimate strength [Pa]


STEEL_304 = IsotropicMaterial(
    name="AISI 304",
    E=193e9,
    nu=0.29,
    rho=8000.0,
    sigma_y=215e6,
    sigma_ult=505e6,
)

