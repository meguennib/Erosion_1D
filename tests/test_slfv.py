import numpy as np
from src.slfv import slfv_advection_phi


def test_open_boundary_mass_conservation_large_cfl():
    L = 1.0
    n = 200
    dx = L / n
    phi = np.zeros(n)
    u = np.ones(n)
    for dt in (0.29, 0.29, 0.15):
        phi = slfv_advection_phi(phi, u, dx, dt, phi_in=1.0)
    assert abs(dx * phi.sum() - 0.73) < 1e-12
    assert phi.min() >= -1e-12
    assert phi.max() <= 1.0 + 1e-12


def test_variable_positive_velocity_finite():
    n = 200
    dx = 1.0 / n
    x = (np.arange(n) + 0.5) * dx
    u = 1.0 + 0.4 * np.sin(2.0 * np.pi * x)
    phi = 0.2 + 0.05 * np.sin(2.0 * np.pi * x)
    out = slfv_advection_phi(phi, u, dx, 20.0 * dx / u.max(), phi_in=phi[0])
    assert np.all(np.isfinite(out))
