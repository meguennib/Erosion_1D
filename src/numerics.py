# src/numerics.py
import numpy as np


# -----------------------------
# Derivatives
# -----------------------------
def ddx_centered(q: np.ndarray, dx: float) -> np.ndarray:
    """
    Centered derivative with one-sided at boundaries.
    q is cell-centered of length N.
    """
    q = np.asarray(q, dtype=float)
    N = q.size
    dq = np.zeros_like(q)

    if N < 2:
        return dq

    # one-sided
    dq[0] = (q[1] - q[0]) / dx
    dq[-1] = (q[-1] - q[-2]) / dx

    if N > 2:
        dq[1:-1] = (q[2:] - q[:-2]) / (2.0 * dx)
    return dq


# -----------------------------
# Limiters (TVD)
# -----------------------------
def _minmod(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """
    minmod limiter.
    """
    s = np.sign(a) + np.sign(b)
    out = np.zeros_like(a)
    mask = (np.abs(s) > 1.5)  # same sign
    out[mask] = np.sign(a[mask]) * np.minimum(np.abs(a[mask]), np.abs(b[mask]))
    return out


def limiter_mc(dL: np.ndarray, dR: np.ndarray) -> np.ndarray:
    """
    Monotonized Central (MC) limiter:
      slope = minmod( 2 dL, 0.5(dL+dR), 2 dR )
    Implemented via pairwise minmod.
    """
    a = 2.0 * dL
    b = 0.5 * (dL + dR)
    c = 2.0 * dR
    return _minmod(_minmod(a, b), c)


def limiter_vanleer(dL: np.ndarray, dR: np.ndarray) -> np.ndarray:
    """
    van Leer limiter: smooth, TVD, less diffusive than minmod.
    slope = (r + |r|) / (1 + |r|)  where r = dL / dR.
    Returns 0 where signs differ.
    """
    dL = np.asarray(dL, dtype=float)
    dR = np.asarray(dR, dtype=float)
    out = np.zeros_like(dL)
    # Only apply where same sign
    mask = (dL * dR > 0)
    r = np.zeros_like(dL)
    r[mask] = dL[mask] / (dR[mask] + 1e-30)
    out[mask] = (r[mask] + np.abs(r[mask])) / (1.0 + np.abs(r[mask])) * dR[mask]
    return out


# -----------------------------
# Ghost cells helper
# -----------------------------
def _fill_ghost_dirichlet_neumann(q: np.ndarray, q_in: float) -> np.ndarray:
    """
    Build 2-ghost array for a cell-centered field q (size N):
      - inlet Dirichlet: q = q_in at x=0 boundary -> set left ghosts to q_in
      - outlet Neumann: dq/dx = 0 at x=L -> set right ghosts equal to last interior
    Returns qg of size N+4 with interior at indices [2 : N+2].
    """
    q = np.asarray(q, dtype=float)
    N = q.size
    qg = np.empty(N + 4, dtype=float)

    # interior
    qg[2:-2] = q

    # inlet Dirichlet
    qg[0] = q_in
    qg[1] = q_in

    # outlet Neumann
    qg[-2] = q[-1]
    qg[-1] = q[-1]
    return qg


def _fill_ghost_neumann(a: np.ndarray) -> np.ndarray:
    """
    Ghosting for velocity a (size N). We use zero-gradient on both ends
    (copy boundary values). This is stable for flux construction.
    """
    a = np.asarray(a, dtype=float)
    N = a.size
    ag = np.empty(N + 4, dtype=float)
    ag[2:-2] = a

    # left
    ag[0] = a[0]
    ag[1] = a[0]

    # right
    ag[-2] = a[-1]
    ag[-1] = a[-1]
    return ag


# -----------------------------
# Advection
# -----------------------------
def advection_phi(
    phi: np.ndarray,
    a: np.ndarray,
    dx: float,
    dt: float,
    scheme: str = "muscl",
    form: str = "nonconservative",
    phi_in: float = 0.0,
    limiter: str = "minmod",
) -> np.ndarray:
    """
    Advect phi by a(x) using a FV scheme with explicit ghost cells:
      - phi inlet Dirichlet: phi=phi_in
      - outlet Neumann: dphi/dx=0
    Two options for PDE form:
      - form="conservative"     : phi_t + (a phi)_x = 0
      - form="nonconservative" : phi_t + a phi_x = 0
            implemented via conservative update + dt*(a_x * phi) correction because:
              a phi_x = (a phi)_x - a_x phi
              => phi_t + (a phi)_x = a_x phi
    scheme:
      - "upwind": 1st order
      - "muscl" : MUSCL-TVD with selectable limiter
    limiter: "minmod", "vanleer", or "mc"
    """
    phi = np.asarray(phi, dtype=float)
    a = np.asarray(a, dtype=float)
    N = phi.size
    if N == 0:
        return phi.copy()

    # ghosted arrays
    qg = _fill_ghost_dirichlet_neumann(phi, phi_in)
    ag = _fill_ghost_neumann(a)

    # face velocities (between cell j and j+1)
    # faces count = (N+4)-1 = N+3
    a_face = 0.5 * (ag[:-1] + ag[1:])

    # reconstruct states for MUSCL, or use piecewise-constant for upwind
    if scheme.lower() == "muscl":
        # slopes on cells (aligned with qg indices)
        # dL and dR defined on qg[1:-1] (size N+2)
        dL = qg[1:-1] - qg[0:-2]
        dR = qg[2:  ] - qg[1:-1]
        s = np.zeros_like(qg)
        # STABILIZATION (Step 3 — selectable limiter):
        # MinMod guarantees strict TVD at all CFL but is diffusive.
        # vanLeer is smoother than MC and less diffusive than minmod.
        if limiter.lower() == "vanleer":
            s[1:-1] = limiter_vanleer(dL, dR)
        elif limiter.lower() == "mc":
            s[1:-1] = limiter_mc(dL, dR)
        else:  # default: minmod
            s[1:-1] = _minmod(dL, dR)

        # left/right states at faces
        qL = qg[:-1] + 0.5 * s[:-1]
        qR = qg[1: ] - 0.5 * s[1: ]
    else:
        # upwind (1st order): left/right are just cell values
        qL = qg[:-1]
        qR = qg[1:]

    # Upwind flux (Godunov for linear advection)
    F = a_face * np.where(a_face >= 0.0, qL, qR)

    # Conservative FV update on interior cells i = 2..N+1
    c = dt / dx
    qnew = qg.copy()
    # cell i uses outflux F[i] (i+1/2) and influx F[i-1] (i-1/2)
    qnew[2:N+2] = qg[2:N+2] - c * (F[2:N+2] - F[1:N+1])

    phi_adv = qnew[2:N+2].copy()

    # Nonconservative correction: + dt * a_x * phi (evaluated at old time)
    if form.lower().startswith("non"):
        ax = ddx_centered(a, dx)
        phi_adv = phi_adv + dt * ax * phi

    # enforce bounds softly (caller may clip with phi_soil)
    # keep inlet cell non-negative in case of small undershoot
    phi_adv = np.where(np.isfinite(phi_adv), phi_adv, 0.0)
    return phi_adv
