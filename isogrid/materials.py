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


# Alias for consistency with naming elsewhere
AISI_304 = STEEL_304


AISI_316L = IsotropicMaterial(
    name="AISI 316L",
    E=193e9,
    nu=0.30,
    rho=8000.0,
    sigma_y=170e6,
    sigma_ult=485e6,
)


AL_6061_T6 = IsotropicMaterial(
    name="Al 6061-T6",
    E=69e9,
    nu=0.33,
    rho=2700.0,
    sigma_y=276e6,
    sigma_ult=310e6,
)


AL_2024_T3 = IsotropicMaterial(
    name="Al 2024-T3",
    E=73e9,
    nu=0.33,
    rho=2780.0,
    sigma_y=324e6,
    sigma_ult=469e6,
)


AL_7075_T6 = IsotropicMaterial(
    name="Al 7075-T6",
    E=72e9,
    nu=0.33,
    rho=2810.0,
    sigma_y=503e6,
    sigma_ult=572e6,
)


MATERIAL_LIBRARY = {
    "AISI_304": AISI_304,
    "STEEL_304": STEEL_304,
    "AISI_316L": AISI_316L,
    "AL_6061_T6": AL_6061_T6,
    "AL_2024_T3": AL_2024_T3,
    "AL_7075_T6": AL_7075_T6,
}

