import numpy as np

from src.config import Params
from src.kernels import hydraulic_fields_kernel, pressure_residual_kernel
from src.physics import (
    barenblatt_fw,
    beta_barenblatt,
    fm_julien,
    mdot_erosion,
    rho_mix,
    shear_tau_b,
)


def test_hydraulic_kernels_match_vectorized_reference():
    params = Params(
        L=0.117,
        Nx=8,
        R0=0.003,
        Pin=4905.0,
        Pout=0.0,
        k_er=0.01,
        tau_c=1.0,
        Re_min=1.0,
    )
    rng = np.random.default_rng(4)
    R = np.ascontiguousarray(params.R0 * (1.0 + 0.2 * rng.random(params.Nx)))
    phi = 0.4 * rng.random(params.Nx)
    rho = np.ascontiguousarray(rho_mix(phi, params.rho_w, params.rho_p))
    fm = np.ascontiguousarray(fm_julien(phi, params))
    Q = 1.2e-4
    dx = params.L / params.Nx

    residual_kernel = pressure_residual_kernel(
        R,
        rho,
        fm,
        params.mu_w,
        params.Re_min,
        params.tau_c,
        params.k_er,
        params.Pin,
        params.Pout,
        params.K_out,
        dx,
        params.R_min,
        Q,
    )

    u = Q / (np.pi * R**2 + 1e-30)
    Re = np.maximum(2.0 * rho * np.abs(u) * R / params.mu_w, params.Re_min)
    Re = np.maximum(Re, 1.0000001)
    fw = barenblatt_fw(Re)
    beta = beta_barenblatt(Re)
    tau = shear_tau_b(rho, fw, fm, u)
    mdot = mdot_erosion(tau, params.tau_c, params.k_er)
    px = 2.0 * tau / np.maximum(R, params.R_min)
    residual_reference = (
        params.Pin
        + np.sum(px) * dx
        - params.K_out * 0.5 * rho[-1] * u[-1] ** 2
        - params.Pout
    )

    fields = hydraulic_fields_kernel(
        R,
        rho,
        fm,
        params.mu_w,
        params.Re_min,
        params.tau_c,
        params.k_er,
        params.R_min,
        Q,
    )

    assert np.isclose(residual_kernel, residual_reference, rtol=1e-12, atol=1e-10)
    for actual, reference in zip(fields, (tau, mdot, fw, beta, u)):
        assert np.allclose(actual, reference, rtol=1e-12, atol=1e-12)
