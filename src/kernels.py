"""Hot numerical kernels used by the hydraulic pressure solve.

The kernels deliberately contain only scalar arithmetic and NumPy arrays so
that they can be compiled by Numba without coupling the compiler to Params,
logging, or Python dictionaries.  A NumPy/Python fallback is retained for
environments where Numba is not installed.
"""

from __future__ import annotations

import math

import numpy as np

try:  # Numba is an optional acceleration dependency.
    from numba import njit

    NUMBA_AVAILABLE = True
except ImportError:  # pragma: no cover - exercised only without Numba
    NUMBA_AVAILABLE = False


def _barenblatt_values(reynolds: float):
    """Scalar Barenblatt friction and beta calculation."""
    alpha = 3.0 / (2.0 * math.log(reynolds))
    alpha = min(max(alpha, 1e-3), 0.5)
    num = (2.0**alpha) * alpha * (1.0 + alpha) * (2.0 + alpha)
    den = math.exp(1.5) * (math.sqrt(3.0) + 5.0 * alpha)
    fw = 8.0 * (num / den) ** (2.0 / (1.0 + alpha))
    beta = ((1.0 + alpha) * (2.0 + alpha) ** 2) / (
        4.0 * (1.0 + 2.0 * alpha)
    )
    return fw, beta


def _pressure_residual_python(
    R,
    rho,
    fm,
    mu_w,
    Re_min,
    tau_c,
    k_er,
    Pin,
    Pout,
    K_out,
    dx,
    R_min,
    Q,
):
    pL = Pin
    u_last = 0.0
    rho_last = rho[-1]
    n = R.size
    for i in range(n):
        u = Q / (math.pi * R[i] * R[i] + 1e-30)
        reynolds = 2.0 * rho[i] * abs(u) * R[i] / mu_w
        reynolds = max(reynolds, Re_min, 1.0000001)
        fw, _ = _barenblatt_values(reynolds)
        tau_b = -0.125 * rho[i] * fw * fm[i] * u * u
        pL += (2.0 / max(R[i], R_min)) * tau_b * dx
        if i == n - 1:
            u_last = u
    pL -= K_out * 0.5 * rho_last * u_last * u_last
    return pL - Pout


def _hydraulic_fields_python(
    R,
    rho,
    fm,
    mu_w,
    Re_min,
    tau_c,
    k_er,
    R_min,
    Q,
):
    n = R.size
    tau_b = np.empty(n, dtype=np.float64)
    mdot = np.empty(n, dtype=np.float64)
    fw_arr = np.empty(n, dtype=np.float64)
    beta_arr = np.empty(n, dtype=np.float64)
    u_arr = np.empty(n, dtype=np.float64)
    for i in range(n):
        u = Q / (math.pi * R[i] * R[i] + 1e-30)
        reynolds = 2.0 * rho[i] * abs(u) * R[i] / mu_w
        reynolds = max(reynolds, Re_min, 1.0000001)
        fw, beta = _barenblatt_values(reynolds)
        tau = -0.125 * rho[i] * fw * fm[i] * u * u
        tau_b[i] = tau
        mdot[i] = k_er * max(abs(tau) - tau_c, 0.0)
        fw_arr[i] = fw
        beta_arr[i] = beta
        u_arr[i] = u
    return tau_b, mdot, fw_arr, beta_arr, u_arr


if NUMBA_AVAILABLE:

    @njit(cache=True)
    def pressure_residual_kernel(
        R,
        rho,
        fm,
        mu_w,
        Re_min,
        tau_c,
        k_er,
        Pin,
        Pout,
        K_out,
        dx,
        R_min,
        Q,
    ):
        pL = Pin
        u_last = 0.0
        rho_last = rho[rho.size - 1]
        for i in range(R.size):
            u = Q / (math.pi * R[i] * R[i] + 1e-30)
            reynolds = 2.0 * rho[i] * abs(u) * R[i] / mu_w
            if reynolds < Re_min:
                reynolds = Re_min
            if reynolds < 1.0000001:
                reynolds = 1.0000001
            alpha = 3.0 / (2.0 * math.log(reynolds))
            if alpha < 1e-3:
                alpha = 1e-3
            elif alpha > 0.5:
                alpha = 0.5
            num = (2.0**alpha) * alpha * (1.0 + alpha) * (2.0 + alpha)
            den = math.exp(1.5) * (math.sqrt(3.0) + 5.0 * alpha)
            fw = 8.0 * (num / den) ** (2.0 / (1.0 + alpha))
            tau_b = -0.125 * rho[i] * fw * fm[i] * u * u
            radius = R[i] if R[i] > R_min else R_min
            pL += (2.0 / radius) * tau_b * dx
            if i == R.size - 1:
                u_last = u
        pL -= K_out * 0.5 * rho_last * u_last * u_last
        return pL - Pout

    @njit(cache=True)
    def hydraulic_fields_kernel(
        R,
        rho,
        fm,
        mu_w,
        Re_min,
        tau_c,
        k_er,
        R_min,
        Q,
    ):
        n = R.size
        tau_b = np.empty(n, dtype=np.float64)
        mdot = np.empty(n, dtype=np.float64)
        fw_arr = np.empty(n, dtype=np.float64)
        beta_arr = np.empty(n, dtype=np.float64)
        u_arr = np.empty(n, dtype=np.float64)
        for i in range(n):
            u = Q / (math.pi * R[i] * R[i] + 1e-30)
            reynolds = 2.0 * rho[i] * abs(u) * R[i] / mu_w
            if reynolds < Re_min:
                reynolds = Re_min
            if reynolds < 1.0000001:
                reynolds = 1.0000001
            alpha = 3.0 / (2.0 * math.log(reynolds))
            if alpha < 1e-3:
                alpha = 1e-3
            elif alpha > 0.5:
                alpha = 0.5
            num = (2.0**alpha) * alpha * (1.0 + alpha) * (2.0 + alpha)
            den = math.exp(1.5) * (math.sqrt(3.0) + 5.0 * alpha)
            fw = 8.0 * (num / den) ** (2.0 / (1.0 + alpha))
            beta = ((1.0 + alpha) * (2.0 + alpha) ** 2) / (
                4.0 * (1.0 + 2.0 * alpha)
            )
            tau = -0.125 * rho[i] * fw * fm[i] * u * u
            tau_b[i] = tau
            mdot[i] = k_er * max(abs(tau) - tau_c, 0.0)
            fw_arr[i] = fw
            beta_arr[i] = beta
            u_arr[i] = u
        return tau_b, mdot, fw_arr, beta_arr, u_arr

else:  # pragma: no cover - fallback is selected only without Numba
    pressure_residual_kernel = _pressure_residual_python
    hydraulic_fields_kernel = _hydraulic_fields_python
