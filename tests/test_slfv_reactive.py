import numpy as np

from src.slfv_reactive import slfv_reactive_phi


def test_reactive_constant_velocity_exact_structure():
    n = 200
    L = 1.0
    dx = L / n
    x = (np.arange(n) + 0.5) * dx
    u = np.ones(n)
    phi_s = 0.62
    phi_in = 0.10
    k0 = 0.05
    phi0 = 0.20 + 0.10 * np.sin(2.0 * np.pi * x)

    for dt in (0.005, 0.05, 0.2):
        out = slfv_reactive_phi(
            phi0, u, np.full(n, k0), dx, dt, phi_s, phi_in
        )
        exact = np.empty(n)
        inside = x > dt
        exact[inside] = (
            phi_s
            + (0.20 + 0.10 * np.sin(2.0 * np.pi * (x[inside] - dt)) - phi_s)
            * np.exp(-k0 * dt)
        )
        exact[~inside] = (
            phi_s + (phi_in - phi_s) * np.exp(-k0 * x[~inside])
        )
        # Second-order reconstruction with cell-average sampling is expected
        # to converge as the mesh is refined; this tolerance is intentionally
        # loose enough for a regression test rather than a formal convergence claim.
        assert np.sqrt(np.mean((out - exact) ** 2)) < 2e-3
        assert out.min() >= -1e-12
        assert out.max() <= phi_s + 1e-12


def test_reactive_operator_accepts_large_courant():
    n = 100
    dx = 1.0 / n
    phi = np.full(n, 0.1)
    u = np.ones(n)
    k = np.full(n, 0.1)
    out = slfv_reactive_phi(phi, u, k, dx, 0.5, 0.62, 0.0)
    assert np.all(np.isfinite(out))
    assert out.min() >= 0.0
    assert out.max() < 0.62
