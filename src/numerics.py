import numpy as np


def ddx_centered(q: np.ndarray, dx: float) -> np.ndarray:
    """Second-order centered derivative in the interior; one-sided at boundaries."""
    q = np.asarray(q, dtype=float)
    N = q.size
    dq = np.zeros_like(q)
    if N < 2:
        return dq
    dq[0] = (q[1] - q[0]) / dx
    dq[-1] = (q[-1] - q[-2]) / dx
    if N > 2:
        dq[1:-1] = (q[2:] - q[:-2]) / (2.0 * dx)
    return dq


def _minmod(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    s = np.sign(a) + np.sign(b)
    out = np.zeros_like(a)
    mask = np.abs(s) > 1.5
    out[mask] = np.sign(a[mask]) * np.minimum(np.abs(a[mask]), np.abs(b[mask]))
    return out


def limiter_mc(dL: np.ndarray, dR: np.ndarray) -> np.ndarray:
    return _minmod(_minmod(2.0 * dL, 0.5 * (dL + dR)), 2.0 * dR)


def limiter_vanleer(dL: np.ndarray, dR: np.ndarray) -> np.ndarray:
    dL = np.asarray(dL, dtype=float)
    dR = np.asarray(dR, dtype=float)
    out = np.zeros_like(dL)
    mask = dL * dR > 0.0
    r = np.zeros_like(dL)
    r[mask] = dL[mask] / (dR[mask] + 1e-30)
    out[mask] = ((r[mask] + np.abs(r[mask])) /
                 (1.0 + np.abs(r[mask]))) * dR[mask]
    return out


def _fill_ghost_dirichlet_neumann(q: np.ndarray, q_in: float) -> np.ndarray:
    q = np.asarray(q, dtype=float)
    N = q.size
    qg = np.empty(N + 4, dtype=float)
    qg[2:-2] = q
    qg[0] = q_in
    qg[1] = q_in
    qg[-2] = q[-1]
    qg[-1] = q[-1]
    return qg


def _fill_ghost_neumann(a: np.ndarray) -> np.ndarray:
    a = np.asarray(a, dtype=float)
    N = a.size
    ag = np.empty(N + 4, dtype=float)
    ag[2:-2] = a
    ag[0] = a[0]
    ag[1] = a[0]
    ag[-2] = a[-1]
    ag[-1] = a[-1]
    return ag


def advection_phi(
    phi: np.ndarray,
    a: np.ndarray,
    dx: float,
    dt: float,
    scheme: str = "muscl",
    form: str = "conservative",
    phi_in: float = 0.0,
    limiter: str = "vanleer",
) -> np.ndarray:
    """
    Conservative finite-volume advection of phi.

    Reference thesis configuration:
        scheme = MUSCL
        limiter = Van Leer
        form = conservative
        velocity = u = Q/A

    The non-conservative option remains available as a numerical utility,
    but is not part of the thesis reference model.
    """
    phi = np.asarray(phi, dtype=float)
    a = np.asarray(a, dtype=float)
    N = phi.size
    if N == 0:
        return phi.copy()

    qg = _fill_ghost_dirichlet_neumann(phi, phi_in)
    ag = _fill_ghost_neumann(a)
    a_face = 0.5 * (ag[:-1] + ag[1:])

    if scheme.lower() == "muscl":
        dL = qg[1:-1] - qg[:-2]
        dR = qg[2:] - qg[1:-1]
        s = np.zeros_like(qg)
        if limiter.lower() == "vanleer":
            s[1:-1] = limiter_vanleer(dL, dR)
        elif limiter.lower() == "mc":
            s[1:-1] = limiter_mc(dL, dR)
        else:
            s[1:-1] = _minmod(dL, dR)

        qL = qg[:-1] + 0.5 * s[:-1]
        qR = qg[1:] - 0.5 * s[1:]
    else:
        qL = qg[:-1]
        qR = qg[1:]

    F = a_face * np.where(a_face >= 0.0, qL, qR)

    c = dt / dx
    qnew = qg.copy()
    qnew[2:N + 2] = qg[2:N + 2] - c * (F[2:N + 2] - F[1:N + 1])
    phi_adv = qnew[2:N + 2].copy()

    if form.lower().startswith("non"):
        ax = ddx_centered(a, dx)
        phi_adv = phi_adv + dt * ax * phi

    return np.where(np.isfinite(phi_adv), phi_adv, 0.0)
