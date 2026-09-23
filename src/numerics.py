# src/numerics.py
from __future__ import annotations

import numpy as np


# -----------------------------
# Derivatives
# -----------------------------
def ddx_centered(q: np.ndarray, dx: float) -> np.ndarray:
    """Centered derivative with one-sided differences at the boundaries."""
    q = np.asarray(q, dtype=float)
    if dx <= 0.0:
        raise ValueError("dx must be strictly positive")

    N = q.size
    dq = np.zeros_like(q)
    if N < 2:
        return dq

    dq[0] = (q[1] - q[0]) / dx
    dq[-1] = (q[-1] - q[-2]) / dx
    if N > 2:
        dq[1:-1] = (q[2:] - q[:-2]) / (2.0 * dx)
    return dq


# -----------------------------
# Limiters (TVD)
# -----------------------------
def _minmod(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Pairwise minmod limiter."""
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    out = np.zeros_like(a)
    same_sign = (a * b) > 0.0
    out[same_sign] = np.sign(a[same_sign]) * np.minimum(
        np.abs(a[same_sign]), np.abs(b[same_sign])
    )
    return out


def limiter_minmod(dL: np.ndarray, dR: np.ndarray) -> np.ndarray:
    """Classical MinMod slope."""
    return _minmod(dL, dR)


def limiter_mc(dL: np.ndarray, dR: np.ndarray) -> np.ndarray:
    """Monotonized-central limiter."""
    return _minmod(
        _minmod(2.0 * dL, 0.5 * (dL + dR)),
        2.0 * dR,
    )


def limiter_vanleer(dL: np.ndarray, dR: np.ndarray) -> np.ndarray:
    """van Leer limiter, retained for controlled comparison studies."""
    dL = np.asarray(dL, dtype=float)
    dR = np.asarray(dR, dtype=float)
    out = np.zeros_like(dL)
    mask = (dL * dR) > 0.0
    r = np.zeros_like(dL)
    r[mask] = dL[mask] / (dR[mask] + 1e-30)
    out[mask] = (r[mask] + np.abs(r[mask])) / (1.0 + np.abs(r[mask])) * dR[mask]
    return out


# -----------------------------
# Ghost cells helper
# -----------------------------
def _fill_ghost_dirichlet_neumann(q: np.ndarray, q_in: float) -> np.ndarray:
    """Build a two-ghost-cell array with inlet Dirichlet and outlet Neumann."""
    q = np.asarray(q, dtype=float)
    N = q.size
    if N == 0:
        return np.empty(0, dtype=float)

    qg = np.empty(N + 4, dtype=float)
    qg[2:-2] = q
    qg[0] = q_in
    qg[1] = q_in
    qg[-2] = q[-1]
    qg[-1] = q[-1]
    return qg


def _fill_ghost_neumann(a: np.ndarray) -> np.ndarray:
    """Build a two-ghost-cell zero-gradient array for the velocity field."""
    a = np.asarray(a, dtype=float)
    N = a.size
    if N == 0:
        return np.empty(0, dtype=float)

    ag = np.empty(N + 4, dtype=float)
    ag[2:-2] = a
    ag[0] = a[0]
    ag[1] = a[0]
    ag[-2] = a[-1]
    ag[-1] = a[-1]
    return ag


# -----------------------------
# Conservative advection
# -----------------------------
def advection_conservative(
    q: np.ndarray,
    a: np.ndarray,
    dx: float,
    dt: float,
    *,
    q_in: float = 0.0,
    scheme: str = "muscl",
    limiter: str = "minmod",
) -> np.ndarray:
    """Advect a conservative scalar ``q`` with velocity ``a``.

    The solved equation is

        q_t + (a q)_x = 0.

    In the erosion model ``q`` is the solid volume per unit axial length,
    ``S = A*phi``.  Advecting ``S`` rather than ``phi`` preserves the total
    suspended solid volume when the conduit area varies.

    The inlet uses a Dirichlet value ``q_in`` and the outlet uses a zero
    gradient condition.  MUSCL reconstruction is paired with an upwind
    Godunov flux for the variable-coefficient linear advection equation.
    """
    q = np.asarray(q, dtype=float)
    a = np.asarray(a, dtype=float)
    if q.ndim != 1 or a.ndim != 1 or q.size != a.size:
        raise ValueError("q and a must be one-dimensional arrays of equal size")
    if dx <= 0.0 or dt < 0.0:
        raise ValueError("dx must be positive and dt must be non-negative")
    if q.size == 0 or dt == 0.0:
        return q.copy()

    scheme_name = scheme.lower()
    limiter_name = limiter.lower()
    if scheme_name not in {"muscl", "upwind"}:
        raise ValueError("scheme must be 'muscl' or 'upwind'")
    if limiter_name not in {"minmod", "vanleer", "mc"}:
        raise ValueError("limiter must be 'minmod', 'vanleer', or 'mc'")

    qg = _fill_ghost_dirichlet_neumann(q, float(q_in))
    ag = _fill_ghost_neumann(a)
    a_face = 0.5 * (ag[:-1] + ag[1:])

    if scheme_name == "muscl":
        dL = qg[1:-1] - qg[:-2]
        dR = qg[2:] - qg[1:-1]
        slope = np.zeros_like(qg)
        if limiter_name == "minmod":
            slope[1:-1] = limiter_minmod(dL, dR)
        elif limiter_name == "vanleer":
            slope[1:-1] = limiter_vanleer(dL, dR)
        else:
            slope[1:-1] = limiter_mc(dL, dR)
        qL = qg[:-1] + 0.5 * slope[:-1]
        qR = qg[1:] - 0.5 * slope[1:]
    else:
        qL = qg[:-1]
        qR = qg[1:]

    flux = a_face * np.where(a_face >= 0.0, qL, qR)
    qnew = qg.copy()
    c = dt / dx
    qnew[2 : q.size + 2] = qg[2 : q.size + 2] - c * (
        flux[2 : q.size + 2] - flux[1 : q.size + 1]
    )
    return qnew[2 : q.size + 2].copy()


def advection_phi(
    phi: np.ndarray,
    a: np.ndarray,
    dx: float,
    dt: float,
    scheme: str = "muscl",
    form: str = "conservative_area",
    phi_in: float = 0.0,
    limiter: str = "minmod",
) -> np.ndarray:
    """Compatibility wrapper for scalar tests.

    The production solver transports ``S=A*phi`` through
    :func:`advection_conservative`.  This wrapper remains available for
    isolated numerical tests, but the former non-conservative form is no
    longer part of the scientific reference model.
    """
    if form.lower().startswith("non"):
        raise ValueError("non-conservative phi transport is not supported")
    return advection_conservative(
        phi,
        a,
        dx,
        dt,
        q_in=phi_in,
        scheme=scheme,
        limiter=limiter,
    )
