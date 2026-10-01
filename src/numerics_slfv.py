"""
Experimental conservative semi-Lagrangian finite-volume transport.

This module is deliberately isolated from src/numerics.py because it is
part of the thesis-v2 numerical development, not yet the frozen reference
solver. It targets

    dphi/dt + d(u*phi)/dx = 0

with cell-averaged phi and a quasi-steady velocity field u(x). For the
current thesis reference case, u>0 because Pin>Pout.

Core idea:
  1. reconstruct phi with a MUSCL/Van-Leer piecewise-linear profile;
  2. trace BOTH faces of every target cell backward along characteristics;
  3. integrate the old reconstruction over the resulting departure interval;
  4. divide by the target-cell width.

Because adjacent target cells share the same departure face, the remap is
conservative up to round-off (for closed/periodic domains) and naturally
handles Courant numbers larger than one.

The current implementation uses a piecewise-linear representation of u
between faces and an exact travel-time map inside each cell. This avoids
subcycling the characteristic when the displacement spans many cells.
"""

from __future__ import annotations

import numpy as np


def vanleer_slope(dm: np.ndarray, dp: np.ndarray) -> np.ndarray:
    dm = np.asarray(dm, dtype=float)
    dp = np.asarray(dp, dtype=float)
    out = np.zeros_like(dm)
    mask = dm * dp > 0.0
    out[mask] = 2.0 * dm[mask] * dp[mask] / (dm[mask] + dp[mask])
    return out


def muscl_slopes(phi: np.ndarray, periodic: bool = False) -> np.ndarray:
    phi = np.asarray(phi, dtype=float)
    if periodic:
        dm = phi - np.roll(phi, 1)
        dp = np.roll(phi, -1) - phi
    else:
        dm = np.empty_like(phi)
        dp = np.empty_like(phi)
        dm[0] = 0.0
        dm[1:] = phi[1:] - phi[:-1]
        dp[:-1] = phi[1:] - phi[:-1]
        dp[-1] = 0.0
    return vanleer_slope(dm, dp)


def build_face_velocity(u_cell: np.ndarray) -> np.ndarray:
    u_cell = np.asarray(u_cell, dtype=float)
    if u_cell.ndim != 1 or u_cell.size < 2:
        raise ValueError("u_cell must be a 1-D array with at least 2 cells.")
    uf = np.empty(u_cell.size + 1, dtype=float)
    uf[0] = u_cell[0]
    uf[-1] = u_cell[-1]
    uf[1:-1] = 0.5 * (u_cell[:-1] + u_cell[1:])
    return uf


def _cell_travel_time(u_left: float, u_right: float, dx: float) -> float:
    if u_left <= 0.0 or u_right <= 0.0:
        raise ValueError("Positive-velocity prototype requires u_face > 0.")
    if abs(u_right - u_left) <= 1e-13 * max(1.0, abs(u_left), abs(u_right)):
        return dx / u_left
    return dx * np.log(u_right / u_left) / (u_right - u_left)


def build_travel_time_faces(u_face: np.ndarray, dx: float) -> np.ndarray:
    u_face = np.asarray(u_face, dtype=float)
    if np.any(u_face <= 0.0):
        raise ValueError("Positive-velocity prototype requires u_face > 0.")
    cell_times = np.array(
        [_cell_travel_time(u_face[i], u_face[i + 1], dx)
         for i in range(u_face.size - 1)],
        dtype=float,
    )
    return np.r_[0.0, np.cumsum(cell_times)]


def _local_travel_time(xi: float, u_left: float, u_right: float, dx: float) -> float:
    xi = float(np.clip(xi, 0.0, dx))
    slope = (u_right - u_left) / dx
    if abs(slope) <= 1e-13 * max(1.0, abs(u_left), abs(u_right)):
        return xi / u_left
    u_x = u_left + slope * xi
    return np.log(u_x / u_left) / slope


def _invert_cell_travel_time(
    tau: float, u_left: float, u_right: float, dx: float
) -> float:
    tau = max(float(tau), 0.0)
    slope = (u_right - u_left) / dx
    if abs(slope) <= 1e-13 * max(1.0, abs(u_left), abs(u_right)):
        xi = u_left * tau
    else:
        xi = u_left * np.expm1(slope * tau) / slope
    return float(np.clip(xi, 0.0, dx))


def face_travel_coordinate(
    x: float, u_face: np.ndarray, tface: np.ndarray, dx: float
) -> float:
    """Return T(x)=integral_0^x ds/u(s) for x in [0,L]."""
    n = u_face.size - 1
    if x <= 0.0:
        return 0.0
    L = n * dx
    if x >= L:
        return float(tface[-1])

    i = min(int(np.floor(x / dx)), n - 1)
    xi = x - i * dx
    return float(
        tface[i]
        + _local_travel_time(xi, u_face[i], u_face[i + 1], dx)
    )


def invert_travel_coordinate(
    target_T: float, u_face: np.ndarray, tface: np.ndarray, dx: float
) -> float:
    """Invert T(x) on [0,L], assuming strictly positive u."""
    if target_T <= 0.0:
        return 0.0
    if target_T >= tface[-1]:
        return (u_face.size - 1) * dx

    i = int(np.searchsorted(tface, target_T, side="right") - 1)
    i = max(0, min(i, u_face.size - 2))
    xi = _invert_cell_travel_time(
        target_T - tface[i], u_face[i], u_face[i + 1], dx
    )
    return i * dx + xi


def backtrace_faces_positive(
    x_faces: np.ndarray,
    u_face: np.ndarray,
    dx: float,
    dt: float,
    *,
    periodic: bool = False,
) -> np.ndarray:
    """
    Backtrace target faces for a frozen positive velocity field.

    Periodic: extend the travel-time map by integer multiples of the travel
    time through one period.

    Open domain: if the departure lies upstream of x=0, a negative
    departure coordinate is returned using the inlet face velocity as the
    local boundary extrapolation. The caller then supplies phi_in there.
    """
    x_faces = np.asarray(x_faces, dtype=float)
    tface = build_travel_time_faces(u_face, dx)
    Tperiod = float(tface[-1])
    L = (u_face.size - 1) * dx
    dep = np.empty_like(x_faces)

    for j, xf in enumerate(x_faces):
        Txf = face_travel_coordinate(xf, u_face, tface, dx)
        target = Txf - dt

        if periodic:
            nper = np.floor(target / Tperiod)
            tau = target - nper * Tperiod
            xlocal = invert_travel_coordinate(tau, u_face, tface, dx)
            dep[j] = nper * L + xlocal
        elif target < 0.0:
            # Constant boundary-face extrapolation in the upstream ghost region.
            dep[j] = target * u_face[0]
        else:
            dep[j] = invert_travel_coordinate(target, u_face, tface, dx)

    return dep


def _primitive_open(
    phi: np.ndarray, slopes: np.ndarray, dx: float, x: float
) -> float:
    """Integral of MUSCL reconstruction from 0 to x for 0<=x<=L."""
    n = phi.size
    L = n * dx
    x = float(np.clip(x, 0.0, L))
    if x <= 0.0:
        return 0.0
    if x >= L:
        return float(dx * np.sum(phi))

    i = min(int(np.floor(x / dx)), n - 1)
    xi = x - i * dx
    xcenter = (i + 0.5) * dx

    cumulative = dx * float(np.sum(phi[:i]))
    local = phi[i] * xi
    local += 0.5 * slopes[i] / dx * ((x - xcenter) ** 2 - (-0.5 * dx) ** 2)
    return cumulative + local


def _primitive_periodic(
    phi: np.ndarray, slopes: np.ndarray, dx: float, x: float
) -> float:
    """Periodic extension of the exact integral of the MUSCL reconstruction."""
    n = phi.size
    L = n * dx
    mass_period = dx * float(np.sum(phi))

    cycles = np.floor(x / L)
    r = x - cycles * L
    if r < 0.0:
        r += L
        cycles -= 1.0
    if r >= L:
        r -= L
        cycles += 1.0

    return cycles * mass_period + _primitive_open(phi, slopes, dx, r)


def conservative_slfv_positive(
    phi: np.ndarray,
    u_cell: np.ndarray,
    dx: float,
    dt: float,
    phi_in: float = 0.0,
    *,
    periodic: bool = False,
) -> np.ndarray:
    """
    One conservative semi-Lagrangian FV step for positive velocity.

    The velocity is frozen during the step. This function intentionally has
    no CFL restriction. Accuracy is controlled by the reconstruction,
    characteristic tracing, and the physical/coupling timestep selected by
    the outer erosion solver.
    """
    phi = np.asarray(phi, dtype=float)
    u_cell = np.asarray(u_cell, dtype=float)
    if phi.shape != u_cell.shape:
        raise ValueError("phi and u_cell must have the same shape.")
    if dt <= 0.0 or dx <= 0.0:
        raise ValueError("Require dt>0 and dx>0.")
    if np.any(u_cell <= 0.0):
        raise ValueError("Positive-velocity prototype requires u_cell > 0.")

    n = phi.size
    L = n * dx
    x_faces = np.linspace(0.0, L, n + 1)
    u_face = build_face_velocity(u_cell)
    dep = backtrace_faces_positive(
        x_faces, u_face, dx, dt, periodic=periodic
    )
    slopes = muscl_slopes(phi, periodic=periodic)

    out = np.empty_like(phi)

    for i in range(n):
        a, b = dep[i], dep[i + 1]
        mass = 0.0

        if periodic:
            mass = (
                _primitive_periodic(phi, slopes, dx, b)
                - _primitive_periodic(phi, slopes, dx, a)
            )
        else:
            # Upstream Dirichlet contribution outside x=0.
            if a < 0.0:
                mass += phi_in * (min(b, 0.0) - a)

            aa = max(a, 0.0)
            bb = min(b, L)
            if bb > aa:
                mass += (
                    _primitive_open(phi, slopes, dx, bb)
                    - _primitive_open(phi, slopes, dx, aa)
                )

        out[i] = mass / dx

    return out


def source_relaxation_exact(
    phi: np.ndarray, k: np.ndarray | float, phi_soil: float, dt: float
) -> np.ndarray:
    """Exact frozen-k solution of phi_t=k(phi_soil-phi)."""
    return phi_soil - (phi_soil - phi) * np.exp(-np.asarray(k) * dt)


def strang_transport_source_positive(
    phi: np.ndarray,
    u_cell: np.ndarray,
    k: np.ndarray | float,
    phi_soil: float,
    dx: float,
    dt: float,
    phi_in: float = 0.0,
    *,
    periodic: bool = False,
) -> np.ndarray:
    """
    Strang splitting:
       source(dt/2) -> conservative SL transport(dt) -> source(dt/2).
    """
    q = source_relaxation_exact(phi, k, phi_soil, 0.5 * dt)
    q = conservative_slfv_positive(
        q, u_cell, dx, dt, phi_in=phi_in, periodic=periodic
    )
    q = source_relaxation_exact(q, k, phi_soil, 0.5 * dt)
    return q
