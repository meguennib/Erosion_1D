from dataclasses import dataclass, field
from typing import List, Optional, Union


@dataclass
class Params:
    # ---- Geometry & mesh ----
    L: float = 0.1          # conduit length (m)
    Nx: int = 600           # number of cells
    R0: float = 3.2e-3      # initial radius (m)
    R_min: float = 1e-6     # minimum allowed radius (m)
    R_break: float = 0.5    # emergency stop if max(R) > R_break (m)

    # ---- Boundary conditions ----
    Pin: float = 2300.0     # inlet pressure (Pa)
    Pout: float = 0.0       # outlet pressure (Pa)
    K_out: float = 0.0      # outlet minor-loss coefficient
    K_in: float = 0.0       # inlet minor-loss coefficient (0.0 = historical p(0)=Pin)
    inlet_kinetic: bool = False  # if True, count the rho*u0^2/2 head needed at the inlet,
                                 # i.e. solve with K_in_total = K_in + 1

    # ---- Fluid / particle properties ----
    rho_w: float = 1000.0   # water density (kg/m3)
    rho_p: float = 2650.0   # particle density (kg/m3)
    rho_s: float = 1600.0   # dry bulk density of intact soil (kg/m3)
    mu_w: float = 1e-3      # dynamic viscosity water (Pa.s)
    phi_soil: float = 0.62  # solid volume fraction of intact soil
    phi_in: float = 0.0     # inlet concentration

    # ---- Erosion parameters ----
    k_er: float = 0.005     # erosion coefficient (m3/(N.s) or m/s)
    tau_c: float = 12.0     # critical shear stress (Pa)

    # ---- Clear-water friction (Barenblatt / smooth pipe) ----
    Re_min: float = 1.0     # floor for Reynolds number

    # ---- Mixture friction multiplier (Julien) ----
    cB: float = 0.2         # dispersive coefficient (Reviewer 1 — C5)
    dp: float = 75e-6       # characteristic particle diameter (m) — silt: 75 µm
    cl: float = 0.07        # mixing length coefficient
    julien_lambda_power: float = 2.0   # exponent for lambda(phi)
    fm_max: float = 5.0     # maximum cap for fm (Julien, 2012, hyperconcentrated flow)
    fm_cap_mode: str = "hard"  # "hard" (published min[fm_max,.]), "smooth", "none" 
    lm_mode: str = "clR0"   # mixing length mode: "clR0" or other

    # ---- Transport scheme ----
    phi_scheme: str = "muscl"        # "muscl" or "upwind"
    phi_form: str = "conservative"   # "conservative" or "nonconservative"
    phi_limiter: str = "vanleer"     # "minmod", "vanleer", or "mc"

    # ---- Time stepping ----
    dt_mode: str = "adaptive"   # "adaptive" or "fixed"
    CFL: float = 0.35
    dt_min: float = 1e-6
    dt_max: float = 10.0
    t_end_factor: float = 15.7   # stop when t = t_end_factor * t_er

    # ---- Pressure solver ----
    Q_tol_abs: float = 1e-3        # absolute tolerance for outlet pressure residual (Pa)
    Q_bisect_max_iter: int = 80    # maximum bisection iterations
    Q_bracket_max_expand: int = 60 # maximum expansions when searching for bracket
    Q_bracket_growth: float = 2.0  # growth factor for Q during bracketing

    # ---- Snapshot times (in units of t_er) ----
    snap_t_over_ter: List[float] = field(default_factory=lambda: [0.0, 2.2, 4.5, 6.8, 15.7])

    # ---- I/O control ----
    log_every_steps: int = 100
    save_ts_every: int = 1    # 1 = save every step; N = save every N-th step

    def to_dict(self):
        return {
            k: v for k, v in self.__dict__.items()
            if not k.startswith('_')
        }


def params_from_dict(d: dict) -> Params:
    """
    Creates a Params object from a dictionary, using recognized fields.
    """
    known = {f.name for f in Params.__dataclass_fields__.values()}
    filtered = {k: v for k, v in d.items() if k in known}
    return Params(**filtered)
