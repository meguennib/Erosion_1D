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

def fm_julien_raw(phi, p):
    """
    Julien mixture friction multiplier, UNBOUNDED form.

        a   = cB * (rho_p / rho(phi)) * (dp/lm)^2
        fm  = 1 + a * lambda(phi)^julien_lambda_power
        lm  = cl * R0        (lm_mode = "clR0")

    This is the diagnostic form. It is NOT what the solver applies: it
    diverges as phi -> phi_soil (lambda -> +infinity). Measured consequence
    for a field-scale coarse material (dp = 1 mm): fm ~ 1.5e12 in a narrow
    zone, producing an unphysical LOCAL radius spike (R_max = 179 R0 while
    the outlet radius is only 3.5 R0) and a spurious stop on the safety
    criterion R_break.

    Use fm_julien(phi, p), which applies the cap.
    """
    rho = rho_mix(phi, p.rho_w, p.rho_p)

    if p.lm_mode == "clR0":
        lm = p.cl * p.R0
    else:
        # fallback safe default
        lm = p.cl * p.R0

    lam = julien_lambda(phi, p.phi_soil)
    a = p.cB * (p.rho_p / rho) * (p.dp / lm) ** 2
    return 1.0 + a * (lam ** p.julien_lambda_power)


def fm_julien(phi, p):
    """
    Mixture friction multiplier AS APPLIED BY THE SOLVER:
    the bounded closure actually published, i.e. with the cap enforced
    inside the time-stepping (not only in post-processing).

        fm(phi) = min[ fm_max , fm_julien_raw(phi) ]

    Cap treatment selected by p.fm_cap_mode:

      "hard"   : fm = min(fm_max, fm_raw)                      <- default,
                 this is exactly the published equation.
      "smooth" : fm = 1 + (fm_max-1)*tanh((fm_raw-1)/(fm_max-1)),
                 a C1-continuous cap, provided as a regularisation variant.
      "none"   : fm = fm_raw. Historically the only mode used. Kept for
                 reproducibility investigations only; NOT recommended.

    HISTORY / WHY THIS REPLACED THE UNBOUNDED-ONLY VERSION
      The original code applied the cap only in the plotting scripts, on the
      argument that a discontinuity in fm would be read by the TVD limiter
      as a genuine extremum and cause sawtooth oscillations at the front.
      That argument was tested: for the laboratory campaign the multiplier
      stays at fm = 1.0001 (cap inactive, all three modes give an identical
      doubling time), and for the field campaign the capped field is SMOOTHER
      than the uncapped one, not rougher -- the uncapped version is what
      creates a spurious isolated spike. Capping inside the solver is what
      makes the field campaign reproduce the published failure times
      (better than 0.5% on four configurations).

    CONCLUSION FOR THESIS REPORTING
      fm_max is an IDENTIFIED physical parameter, not a numerical guard-rail:
      with fm_max = 5 the field configuration arrests instead of running
      away, with fm_max = 2000 it reproduces the published values. It must be
      reported for every published configuration.
    """
    fm_raw = fm_julien_raw(phi, p)

    cap_mode = str(getattr(p, "fm_cap_mode", "hard")).lower()
    if cap_mode == "none":
        return fm_raw

    fm_max = float(getattr(p, "fm_max", 5.0))
    if fm_max <= 1.0:
        return np.ones_like(fm_raw)

    if cap_mode == "hard":
        return np.minimum(fm_raw, fm_max)

    # "smooth" (regularisation variant): C1 cap with asymptote fm_max
    return 1.0 + (fm_max - 1.0) * np.tanh((fm_raw - 1.0) / (fm_max - 1.0))


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
