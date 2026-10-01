from dataclasses import dataclass, field
from typing import List


@dataclass
class Params:
    # ---- Geometry & mesh ----
    L: float = 0.1
    Nx: int = 600
    R0: float = 3.2e-3
    R_min: float = 1e-6
    R_break: float = 0.5

    # ---- Boundary conditions ----
    Pin: float = 2300.0
    Pout: float = 0.0
    K_out: float = 0.0

    # ---- Fluid / particle properties ----
    rho_w: float = 1000.0
    rho_p: float = 2650.0
    mu_w: float = 1e-3
    phi_soil: float = 0.62
    phi_init: float = 0.0
    phi_in: float = 0.0

    # ---- Erosion parameters ----
    k_er: float = 0.005       # [s/m] for mdot = k_er (|tau_b|-tau_c)
    tau_c: float = 12.0       # [Pa]

    # ---- Clear-water friction: Barenblatt ----
    Re_min: float = 1.0       # numerical floor only; not a validity criterion

    # ---- Mixture-friction closure (model-specific phenomenological law) ----
    cB: float = 0.2
    dp: float = 75e-6
    cl: float = 0.07
    julien_lambda_power: float = 2.0
    fm_max: float = 5.0       # post-processing/diagnostic cap only; not used by solver
    lm_mode: str = "clR0"

    # ---- Transport scheme ----
    # "eulerian" reproduces the v1 reference transport operator.
    # "slfv" activates the experimental conservative semi-Lagrangian operator.
    transport_method: str = "eulerian"
    phi_scheme: str = "muscl"
    phi_form: str = "conservative"
    phi_limiter: str = "vanleer"

    # ---- Time stepping ----
    dt_mode: str = "adaptive"
    CFL: float = 0.35
    dt_min: float = 1e-6
    dt_max: float = 10.0
    t_end_factor: float = 15.7

    # ---- Pressure solver ----
    Q_tol_abs: float = 1e-3       # [Pa]
    Q_tol_rel: float = 1e-8       # relative pressure residual tolerance
    Q_bisect_max_iter: int = 80
    Q_bracket_max_expand: int = 60
    Q_bracket_growth: float = 2.0

    # ---- Snapshot times (in units of t_er) ----
    snap_t_over_ter: List[float] = field(
        default_factory=lambda: [0.0, 2.2, 4.5, 6.8, 15.7]
    )

    # ---- I/O control ----
    log_every_steps: int = 100
    save_ts_every: int = 1

    def to_dict(self):
        return {k: v for k, v in self.__dict__.items() if not k.startswith("_")}


def params_from_dict(d: dict) -> Params:
    known = {f.name for f in Params.__dataclass_fields__.values()}
    filtered = {k: v for k, v in d.items() if k in known}
    return Params(**filtered)
