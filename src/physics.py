import numpy as np

def rho_mix(phi, rho_w, rho_p):
    return rho_w + phi * (rho_p - rho_w)

def rho_soil_sat(rho_w, rho_p, phi_soil):
    return rho_w + phi_soil * (rho_p - rho_w)

def barenblatt_fw(Re):
    """
    Complete Darcy-Weisbach friction factor f_D.

    Based on Barenblatt's power law with alpha = 3/(2 ln Re).
    Returns f_D such that: tau_b = (1/8) * f_D * rho * u^2

    Correction:
      The returned value is now the FULL Darcy factor f_D.
      The 1/8 coefficient is applied explicitly in shear_tau_b,
      following the strict Darcy-Weisbach formulation.
    """
    Re = np.maximum(Re, 1.0000001)
    alpha = 3.0 / (2.0 * np.log(Re))

    # robust clamping (avoid alpha->inf when Re->1)
    alpha = np.clip(alpha, 1e-3, 0.5)

    num = (2.0 ** alpha) * alpha * (1.0 + alpha) * (2.0 + alpha)
    den = (np.exp(1.5)) * (np.sqrt(3.0) + 5.0 * alpha)

    # f_D = full Darcy factor
    fD = 8.0 * (num / den) ** (2.0 / (1.0 + alpha))
    return fD

def beta_barenblatt(Re):
    """
    beta(Re) from the same alpha parameterization.
    """
    Re = np.maximum(Re, 1.0000001)
    alpha = 3.0 / (2.0 * np.log(Re))
    alpha = np.clip(alpha, 1e-3, 0.5)
    beta = ((1.0 + alpha) * (2.0 + alpha) ** 2) / (4.0 * (1.0 + 2.0 * alpha))
    return beta

def julien_lambda(phi, phi_soil, eps=1e-12):
    """
    lambda(phi) = 1 / ((phi_soil/phi)^(1/3) - 1)
    diverges as phi -> phi_soil.
    """
    phi = np.clip(phi, eps, phi_soil * (1.0 - 1e-6))
    return 1.0 / ((phi_soil / phi) ** (1.0 / 3.0) - 1.0)

def fm_julien(phi, p):
    """
    Julien mixture friction multiplier with:
      lm_mode = "clR0"  => lm = cl*R0 constant (Option A)
      power = julien_lambda_power = 2 (requested)

    Classical scaling:
      a = cB * (rho_p/rho) * (dp/lm)^2
      fm = 1 + a * lambda(phi)^power

    NOTE — NO CAP APPLIED HERE:
      Any cap (hard min or smooth tanh) on fm during the time-stepping
      introduces a mathematical discontinuity (or near-discontinuity)
      that MUSCL/TVD slope limiters interpret as a genuine extremum,
      triggering excessive limiting and sawtooth oscillations at the
      advective front.

      The solver must operate on the raw, unbounded fm so that the
      spatial gradient of fm is smooth everywhere. The cap fm <= 5.0
      is enforced ONLY as a post-processing step inside the plotting
      scripts (np.minimum(fm_raw, 5.0)), never during the solve.

      fm_raw = 1 + a * lambda^power is C-infinity in phi because:
        - julien_lambda is C-inf away from phi=0 (phi is clipped to eps)
        - the only singularity is at phi -> phi_soil, which is also
          clipped (phi <= phi_soil*(1-1e-6)), so the formula is bounded
          by construction given the physical state constraints.
    """
    rho = rho_mix(phi, p.rho_w, p.rho_p)

    if p.lm_mode == "clR0":
        lm = p.cl * p.R0
    else:
        # fallback safe default
        lm = p.cl * p.R0

    lam = julien_lambda(phi, p.phi_soil)
    a = p.cB * (p.rho_p / rho) * (p.dp / lm) ** 2

    # Raw, unbounded fm — NO cap, NO tanh blending, NO np.minimum
    fm = 1.0 + a * (lam ** p.julien_lambda_power)
    return fm


def shear_tau_b(rho, fw, fm, u):
    """
    Wall shear stress — strict Darcy-Weisbach formulation:
        tau_b = -(1/8) * f_D * rho * fm * u^2

    Correction:
      The 1/8 coefficient is now EXPLICIT.
      fw must be the FULL Darcy factor f_D returned by barenblatt_fw.
    """
    return -(1.0 / 8.0) * rho * fw * fm * (u ** 2)

def mdot_erosion(tau_b, tau_c, k_er):
    """
    m_dot = k_er * max(|tau_b| - tau_c, 0)
    """
    return k_er * np.maximum(np.abs(tau_b) - tau_c, 0.0)
