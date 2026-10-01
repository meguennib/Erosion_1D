import numpy as np

from src.numerics_slfv import (
    conservative_slfv_positive,
    strang_transport_source_positive,
)


def test_periodic_constant_velocity_large_cfl():
    L = 1.0
    N = 200
    dx = L / N
    x = (np.arange(N) + 0.5) * dx
    u = np.full(N, 1.0)

    phi0 = 0.25 + 0.15 * np.sin(2.0 * np.pi * x / L)
    T = 0.937

    dt_nom = 50.0 * dx / 1.0
    nsteps = int(np.ceil(T / dt_nom))
    phi = phi0.copy()

    for j in range(nsteps):
        dt = dt_nom if j < nsteps - 1 else T - dt_nom * (nsteps - 1)
        phi = conservative_slfv_positive(phi, u, dx, dt, periodic=True)

    exact = 0.25 + 0.15 * np.sin(2.0 * np.pi * (x - T) / L)

    assert np.sqrt(np.mean((phi - exact) ** 2)) < 3e-5
    assert abs(dx * np.sum(phi) - dx * np.sum(phi0)) < 1e-12
    assert phi.min() >= 0.0
    assert phi.max() <= 0.5


def test_open_boundary_mass_conservation():
    L = 1.0
    N = 200
    dx = L / N
    u = np.full(N, 1.0)
    phi = np.zeros(N)
    T = 0.73

    dt = 20.0 * dx
    nsteps = int(np.ceil(T / dt))

    for j in range(nsteps):
        local_dt = dt if j < nsteps - 1 else T - dt * (nsteps - 1)
        phi = conservative_slfv_positive(
            phi, u, dx, local_dt, phi_in=1.0, periodic=False
        )

    expected_mass = T  # u*T, since T<L
    numerical_mass = dx * np.sum(phi)

    assert abs(numerical_mass - expected_mass) < 1e-12
    assert phi.min() >= -1e-12
    assert phi.max() <= 1.0 + 1e-12


def test_variable_velocity_conservation_large_cfl():
    L = 1.0
    N = 400
    dx = L / N
    x = (np.arange(N) + 0.5) * dx
    u = 1.0 + 0.4 * np.sin(2.0 * np.pi * x / L)
    phi = 0.25 + 0.08 * np.sin(2.0*np.pi*x/L) + 0.03*np.cos(4.0*np.pi*x/L)
    mass0 = dx * np.sum(phi)

    # CFL based on max(u) is approximately 50.
    dt = 50.0 * dx / np.max(u)
    T = 0.731
    nsteps = int(np.ceil(T / dt))

    for j in range(nsteps):
        local_dt = dt if j < nsteps - 1 else T - dt * (nsteps - 1)
        phi = conservative_slfv_positive(
            phi, u, dx, local_dt, periodic=True
        )

    mass1 = dx * np.sum(phi)
    assert abs(mass1 - mass0) < 1e-12
    assert np.all(np.isfinite(phi))


def test_constant_velocity_with_relaxation():
    L = 1.0
    N = 200
    dx = L/N
    x = (np.arange(N)+0.5)*dx
    u = np.full(N, 1.0)
    phi_soil = 0.62
    k = 0.05
    T = 2.0

    phi0 = 0.22 + 0.08*np.sin(2*np.pi*x/L)
    dt = 50.0*dx
    nsteps = int(np.ceil(T/dt))
    phi = phi0.copy()

    for j in range(nsteps):
        local_dt = dt if j < nsteps-1 else T-dt*(nsteps-1)
        phi = strang_transport_source_positive(
            phi, u, k, phi_soil, dx, local_dt, periodic=True
        )

    exact = phi_soil + (
        0.22 + 0.08*np.sin(2*np.pi*(x-T)/L) - phi_soil
    )*np.exp(-k*T)

    assert np.sqrt(np.mean((phi-exact)**2)) < 3e-5
    assert abs(dx*np.sum(phi)-dx*np.sum(exact)) < 1e-12
