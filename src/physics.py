from __future__ import annotations

from typing import Any

import numpy as np


def rho_mix(phi, rho_w: float, rho_p: float):
    """Mixture density for a local suspended solid volume fraction ``phi``."""
    return rho_w + phi * (rho_p - rho_w)


def rho_soil_sat(rho_w: float, rho_p: float, phi_soil: float):
    """Bulk density of intact saturated soil at ``phi_soil``."""
    return rho_w + phi_soil * (rho_p - rho_w)


def barenblatt_fw(Re):
    """Return the full Darcy-Weisbach friction factor ``f_D``.

    The wall-stress law is

        tau_b = -(1/8) * f_D * rho * f_m * u**2.

    ``Re`` is floored before evaluating the Barenblatt exponent and the
    exponent is clipped to keep the constitutive law finite near ``Re=1``.
    """
    Re = np.maximum(np.asarray(Re, dtype=float), 1.0000001)
    alpha = 3.0 / (2.0 * np.log(Re))
    alpha = np.clip(alpha, 1e-3, 0.5)

    num = (2.0**alpha) * alpha * (1.0 + alpha) * (2.0 + alpha)
    den = np.exp(1.5) * (np.sqrt(3.0) + 5.0 * alpha)
    return 8.0 * (num / den) ** (2.0 / (1.0 + alpha))


def beta_barenblatt(Re):
    """Return the Barenblatt momentum/transport coefficient ``beta(Re)``."""
    Re = np.maximum(np.asarray(Re, dtype=float), 1.0000001)
    alpha = 3.0 / (2.0 * np.log(Re))
    alpha = np.clip(alpha, 1e-3, 0.5)
    return ((1.0 + alpha) * (2.0 + alpha) ** 2) / (4.0 * (1.0 + 2.0 * alpha))


def julien_lambda(phi, phi_soil: float, eps: float = 1e-12):
    """Return Julien's concentration parameter with singularity protection.

    The clipping only protects the endpoints of the admissible concentration
    interval.  The C1 rheological smoothing is applied later to ``fm_raw``;
    this function itself is not a smoothing function.
    """
    phi = np.asarray(phi, dtype=float)
    phi = np.clip(phi, eps, phi_soil * (1.0 - 1e-6))
    return 1.0 / ((phi_soil / phi) ** (1.0 / 3.0) - 1.0)


def fm_julien_raw(phi, p: Any):
    """Return the uncapped Julien mixture-friction multiplier.

    The raw law is

        fm_raw = 1 + cB * (rho_p/rho) * (dp/lm)**2 * lambda(phi)**n.

    It is retained for diagnostics only.  The hydrodynamic solver must use
    :func:`fm_julien`, which applies the active smooth ``fm_max`` cap.
    """
    rho = rho_mix(phi, p.rho_w, p.rho_p)

    if p.lm_mode != "clR0":
        raise ValueError(f"Unsupported Julien mixing-length mode: {p.lm_mode}")
    lm = p.cl * p.R0
    if lm <= 0.0:
        raise ValueError("Julien mixing length must be strictly positive")

    lam = julien_lambda(phi, p.phi_soil)
    coefficient = p.cB * (p.rho_p / rho) * (p.dp / lm) ** 2
    return 1.0 + coefficient * (lam ** p.julien_lambda_power)


def fm_julien_from_raw(fm_raw, p: Any):
    """Apply Julien's smooth active cap to a previously computed raw value."""
    fm_max = float(p.fm_max)
    if fm_max <= 1.0:
        raise ValueError("fm_max must be strictly greater than 1")

    fm_raw = np.asarray(fm_raw, dtype=float)
    span = fm_max - 1.0
    z = (fm_raw - 1.0) / span
    return 1.0 + span * np.tanh(z)


def fm_julien(phi, p: Any):
    """Return Julien's C1-smoothed, physically capped multiplier.

    The active constitutive law is

        fm = 1 + (fm_max - 1) * tanh((fm_raw - 1)/(fm_max - 1)).

    For the admissible positive parameters, ``fm_raw >= 1`` and therefore
    ``1 <= fm < fm_max``.  The hyperbolic tangent is smooth; endpoint clips
    on ``phi`` remain separate numerical domain guards.
    """
    return fm_julien_from_raw(fm_julien_raw(phi, p), p)


def shear_tau_b(rho, fw, fm, u):
    """Return wall shear stress using the full Darcy factor ``fw``."""
    return -(1.0 / 8.0) * rho * fw * fm * (u**2)


def mdot_erosion(tau_b, tau_c: float, k_er: float):
    """Return erosion mass flux ``mdot`` [kg m^-2 s^-1].

    ``mdot`` is the positive mass flux of intact saturated soil removed from
    the wall.  The corresponding radial wall-growth rate is

        R_t = mdot / rho_soil_sat.

    The coefficient ``k_er`` carries the calibrated units needed to transform
    excess shear stress into this mass flux.
    """
    return k_er * np.maximum(np.abs(tau_b) - tau_c, 0.0)


def radius_growth_rate(mdot, rho_soil_sat_value):
    """Convert erosion mass flux into radial growth rate ``R_t`` [m/s]."""
    if rho_soil_sat_value <= 0.0:
        raise ValueError("rho_soil_sat must be strictly positive")
    return np.asarray(mdot, dtype=float) / rho_soil_sat_value


def conservative_solid_source(area_old, area_new, phi_soil: float):
    """Return the solid-volume source caused by a wall-area increment.

    The newly opened saturated-soil volume carries the intact solid fraction,
    so the source in ``S=A*phi`` is ``phi_soil * (A_new-A_old)``.  The helper
    keeps this conservation statement in the physical layer rather than
    hiding it inside the time integrator.
    """
    area_old = np.asarray(area_old, dtype=float)
    area_new = np.asarray(area_new, dtype=float)
    if area_old.shape != area_new.shape:
        raise ValueError("area_old and area_new must have identical shapes")
    if not (0.0 <= phi_soil <= 1.0):
        raise ValueError("phi_soil must lie between 0 and 1")
    if np.any(area_new < area_old):
        raise ValueError("the erosion source cannot decrease the conduit area")
    return phi_soil * (area_new - area_old)
