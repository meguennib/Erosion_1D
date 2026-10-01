import numpy as np


def rho_mix(phi, rho_w, rho_p):
    """Effective mixture density: rho_m = rho_w + phi (rho_p-rho_w)."""
    return rho_w + phi * (rho_p - rho_w)


def rho_soil_sat(rho_w, rho_p, phi_soil):
    """Saturated bulk density associated with the intact-soil solid fraction."""
    return rho_w + phi_soil * (rho_p - rho_w)


def barenblatt_fw(Re):
    """
    Full Darcy-Weisbach friction factor from the Barenblatt scaling law.

    The numerical floor on Re prevents division by zero/log singularities.
    It is a numerical safeguard, not a statement that the closure is valid
    at arbitrary Reynolds number.
    """
    Re = np.maximum(np.asarray(Re, dtype=float), 1.0000001)
    alpha = 3.0 / (2.0 * np.log(Re))
    alpha = np.clip(alpha, 1e-3, 0.5)

    num = (2.0 ** alpha) * alpha * (1.0 + alpha) * (2.0 + alpha)
    den = np.exp(1.5) * (np.sqrt(3.0) + 5.0 * alpha)
    return 8.0 * (num / den) ** (2.0 / (1.0 + alpha))


def beta_barenblatt(Re):
    """
    Legacy auxiliary relation retained for exploratory studies.

    It is NOT used by the thesis reference solver.
    """
    Re = np.maximum(np.asarray(Re, dtype=float), 1.0000001)
    alpha = 3.0 / (2.0 * np.log(Re))
    alpha = np.clip(alpha, 1e-3, 0.5)
    return ((1.0 + alpha) * (2.0 + alpha) ** 2) / (
        4.0 * (1.0 + 2.0 * alpha)
    )


def julien_lambda(phi, phi_soil, eps=1e-12):
    """
    Dimensionless concentration-dependent length-scale parameter.

    lambda(phi) = 1 / ((phi_soil/phi)^(1/3) - 1)

    The clipping keeps the numerical state strictly inside 0 < phi < phi_soil.
    """
    phi = np.clip(phi, eps, phi_soil * (1.0 - 1e-6))
    return 1.0 / ((phi_soil / phi) ** (1.0 / 3.0) - 1.0)


def fm_mixture(phi, p):
    """
    Model-specific phenomenological multiplier for mixture hydraulic resistance.

    fm = 1 + cB (rho_p/rho_m) (dp/lm)^2 lambda(phi)^n

    with lm = cl*R0 for the reference closure.

    This relation is retained as a closure specific to the present model.
    It is not presented as an exact equation taken from Julien.
    No cap is applied during the numerical solve.
    """
    rho = rho_mix(phi, p.rho_w, p.rho_p)

    if p.lm_mode == "clR0":
        lm = p.cl * p.R0
    else:
        lm = p.cl * p.R0

    lam = julien_lambda(phi, p.phi_soil)
    coeff = p.cB * (p.rho_p / rho) * (p.dp / lm) ** 2
    return 1.0 + coeff * lam ** p.julien_lambda_power


# Backward-compatible alias for existing post-processing code.
fm_julien = fm_mixture


def shear_tau_b(rho, fw, fm, u):
    """Signed wall shear stress using the full Darcy friction factor."""
    return -(1.0 / 8.0) * rho * fw * fm * u * np.abs(u)


def mdot_erosion(tau_b, tau_c, k_er):
    """Surface erosion mass flux: mdot = k_er max(|tau_b|-tau_c, 0)."""
    return k_er * np.maximum(np.abs(tau_b) - tau_c, 0.0)
