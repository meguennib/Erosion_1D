from dataclasses import dataclass, field
from typing import List


@dataclass
class Params:
    """Validated parameters for the coupled erosion simulation.

    The scientific reference formulation transports the conservative solid
    volume variable ``S = A * phi``.  ``phi`` is reconstructed from ``S`` and
    the local conduit area when it is needed by the constitutive laws.
    """

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
    phi_in: float = 0.0

    # ---- Erosion parameters ----
    # mdot is an erosion mass flux [kg m^-2 s^-1].  The radial growth law is
    # R_t = mdot / rho_soil_sat, so k_er has the corresponding calibrated
    # units required to convert excess shear stress into a mass flux.
    k_er: float = 0.005
    tau_c: float = 12.0

    # ---- Clear-water friction (Barenblatt / smooth pipe) ----
    Re_min: float = 1.0

    # ---- Mixture friction multiplier (Julien) ----
    cB: float = 0.2
    dp: float = 75e-6
    cl: float = 0.07
    julien_lambda_power: float = 2.0
    fm_max: float = 5.0
    lm_mode: str = "clR0"

    # ---- Conservative transport scheme ----
    phi_scheme: str = "muscl"
    # ``conservative`` is retained as a compatibility spelling and is mapped
    # to the scientific reference form ``conservative_area``.
    phi_form: str = "conservative_area"
    phi_limiter: str = "minmod"

    # ---- Time stepping ----
    dt_mode: str = "adaptive"
    CFL: float = 0.35
    dt_min: float = 1e-6
    dt_max: float = 10.0
    t_end_factor: float = 15.7

    # ---- Pressure solver ----
    Q_tol_abs: float = 1e-3
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

    def __post_init__(self) -> None:
        # Normalize enum-like settings before validation and execution.
        self.phi_scheme = self.phi_scheme.lower()
        self.phi_limiter = self.phi_limiter.lower()
        self.dt_mode = self.dt_mode.lower()
        self.lm_mode = self.lm_mode.strip()
        # Older scenario files used this name for the scalar conservative
        # form.  The implementation now always transports S=A*phi.
        if self.phi_form == "conservative":
            self.phi_form = "conservative_area"
        self.validate()

    def validate(self) -> None:
        """Raise ``ValueError`` for an invalid scientific/numerical setup."""
        errors = []

        if self.L <= 0.0:
            errors.append("L must be strictly positive")
        if self.Nx < 2:
            errors.append("Nx must be at least 2")
        if self.R0 <= 0.0:
            errors.append("R0 must be strictly positive")
        if not (0.0 < self.R_min < self.R0):
            errors.append("R_min must satisfy 0 < R_min < R0")
        if self.R_break <= self.R_min:
            errors.append("R_break must be greater than R_min")

        if self.Pin < self.Pout:
            errors.append("Pin must be greater than or equal to Pout")
        if self.rho_w <= 0.0 or self.rho_p <= 0.0:
            errors.append("densities must be strictly positive")
        if self.mu_w <= 0.0:
            errors.append("mu_w must be strictly positive")
        if not (0.0 < self.phi_soil < 1.0):
            errors.append("phi_soil must lie strictly between 0 and 1")
        if not (0.0 <= self.phi_in <= self.phi_soil):
            errors.append("phi_in must lie between 0 and phi_soil")

        if self.k_er < 0.0:
            errors.append("k_er must be non-negative")
        if self.tau_c < 0.0:
            errors.append("tau_c must be non-negative")
        if self.Re_min < 1.0:
            errors.append("Re_min must be at least 1")
        if self.dp <= 0.0 or self.cl <= 0.0:
            errors.append("dp and cl must be strictly positive")
        if self.fm_max <= 1.0:
            errors.append("fm_max must be strictly greater than 1")
        if self.julien_lambda_power <= 0.0:
            errors.append("julien_lambda_power must be strictly positive")
        if self.lm_mode != "clR0":
            errors.append("only lm_mode='clR0' is currently supported")

        if self.phi_scheme.lower() not in {"muscl", "upwind"}:
            errors.append("phi_scheme must be 'muscl' or 'upwind'")
        if self.phi_form != "conservative_area":
            errors.append("the reference transport form is conservative_area")
        if self.phi_limiter.lower() not in {"minmod", "vanleer", "mc"}:
            errors.append("phi_limiter must be 'minmod', 'vanleer', or 'mc'")

        if self.dt_mode.lower() not in {"adaptive", "fixed"}:
            errors.append("dt_mode must be 'adaptive' or 'fixed'")
        if not (0.0 < self.CFL <= 1.0):
            errors.append("CFL must satisfy 0 < CFL <= 1")
        if self.dt_min <= 0.0 or self.dt_max < self.dt_min:
            errors.append("dt_max must be greater than or equal to dt_min > 0")
        if self.t_end_factor <= 0.0:
            errors.append("t_end_factor must be strictly positive")

        if self.Q_tol_abs <= 0.0:
            errors.append("Q_tol_abs must be strictly positive")
        if self.Q_bisect_max_iter < 1:
            errors.append("Q_bisect_max_iter must be at least 1")
        if self.Q_bracket_max_expand < 0:
            errors.append("Q_bracket_max_expand must be non-negative")
        if self.Q_bracket_growth <= 1.0:
            errors.append("Q_bracket_growth must be greater than 1")

        if any(float(v) < 0.0 for v in self.snap_t_over_ter):
            errors.append("snapshot times must be non-negative")
        if list(self.snap_t_over_ter) != sorted(self.snap_t_over_ter):
            errors.append("snapshot times must be sorted increasingly")
        if self.log_every_steps < 1 or self.save_ts_every < 1:
            errors.append("log_every_steps and save_ts_every must be at least 1")

        if errors:
            raise ValueError("Invalid simulation parameters: " + "; ".join(errors))

    def to_dict(self):
        """Return a JSON-serialisable dictionary of effective parameters."""
        return {k: v for k, v in self.__dict__.items() if not k.startswith("_")}


def params_from_dict(d: dict) -> Params:
    """Create validated parameters from recognised dictionary fields."""
    known = {f.name for f in Params.__dataclass_fields__.values()}
    filtered = {k: v for k, v in d.items() if k in known}
    return Params(**filtered)
