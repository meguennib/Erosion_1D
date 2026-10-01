import numpy as np


def _face_values(q: np.ndarray) -> np.ndarray:
    q = np.asarray(q, dtype=float)
    if q.size == 0:
        return q.copy()
    f = np.empty(q.size + 1, dtype=float)
    f[0] = q[0]
    f[-1] = q[-1]
    if q.size > 1:
        f[1:-1] = 0.5 * (q[:-1] + q[1:])
    return f


def _travel_time_faces(u: np.ndarray, dx: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Build the monotone travel-time coordinate y(x)=integral_0^x dxi/u(xi).

    Face velocities are linearly interpolated inside each x-cell.
    Returns y-faces, y-cell centers and y-cell widths.
    """
    uf = _face_values(u)
    if np.any(uf <= 0.0):
        raise ValueError("Characteristic reactive SLFV requires u>0 everywhere.")

    n = u.size
    dy = np.empty(n, dtype=float)

    for i in range(n):
        a = uf[i]
        b = uf[i + 1]
        if abs(b - a) <= 1e-14 * max(1.0, abs(a), abs(b)):
            dy[i] = dx / a
        else:
            dy[i] = dx * np.log(b / a) / (b - a)

    y_faces = np.r_[0.0, np.cumsum(dy)]
    y_centers = 0.5 * (y_faces[:-1] + y_faces[1:])
    return y_faces, y_centers, dy


def _vanleer_gradient(g_left: np.ndarray, g_right: np.ndarray) -> np.ndarray:
    out = np.zeros_like(g_left, dtype=float)
    mask = g_left * g_right > 0.0
    out[mask] = (
        2.0 * g_left[mask] * g_right[mask]
        / (g_left[mask] + g_right[mask])
    )
    return out


def _reconstruct_psi(psi: np.ndarray, y_centers: np.ndarray) -> np.ndarray:
    """Second-order MUSCL/Van-Leer slope on the non-uniform travel-time grid."""
    n = psi.size
    slope = np.zeros_like(psi)
    if n < 2:
        return slope

    g = np.empty(max(n - 2, 0), dtype=float)
    if n > 2:
        left = (psi[1:-1] - psi[:-2]) / (y_centers[1:-1] - y_centers[:-2])
        right = (psi[2:] - psi[1:-1]) / (y_centers[2:] - y_centers[1:-1])
        g = _vanleer_gradient(left, right)
        slope[1:-1] = g

    slope[0] = (
        (psi[1] - psi[0]) / (y_centers[1] - y_centers[0])
        if n > 1 else 0.0
    )
    slope[-1] = (
        (psi[-1] - psi[-2]) / (y_centers[-1] - y_centers[-2])
        if n > 1 else 0.0
    )
    return slope


def _prepare_source_integrals(
    y_faces: np.ndarray,
    k: np.ndarray,
    u: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Build cumulative optical depth K(y) and source primitive H(y):

        K(y) = integral_0^y k(s) ds
        H(y) = integral_0^y k(s) u(s) exp(K(s)) ds

    Cell-center values are linearly interpolated to faces and integrated by
    trapezoidal quadrature. This is the numerical representation of the
    frozen-in-time characteristic reaction coefficient.
    """
    kf = _face_values(k)
    ku_f = _face_values(k * u)
    dy = np.diff(y_faces)

    K_faces = np.zeros_like(y_faces)
    H_faces = np.zeros_like(y_faces)

    K_faces[1:] = np.cumsum(0.5 * (kf[:-1] + kf[1:]) * dy)

    # Approximate exp(K) at faces after K has been accumulated.
    integrand = ku_f * np.exp(K_faces)
    H_faces[1:] = np.cumsum(
        0.5 * (integrand[:-1] + integrand[1:]) * dy
    )

    return K_faces, H_faces


def _interp_face_field(y, y_faces, values):
    return np.interp(
        np.asarray(y, dtype=float),
        y_faces,
        values,
        left=values[0],
        right=values[-1],
    )


def _interp_psi(y, psi, slope, y_centers, y_faces, psi_in):
    y = np.asarray(y, dtype=float)
    out = np.empty_like(y)

    neg = y < 0.0
    out[neg] = psi_in

    pos = ~neg
    if np.any(pos):
        yy = y[pos]
        idx = np.searchsorted(y_faces, yy, side="right") - 1
        idx = np.clip(idx, 0, psi.size - 1)
        out[pos] = psi[idx] + slope[idx] * (yy - y_centers[idx])

    return out


def _characteristic_reaction(
    y_departure,
    dt,
    psi,
    slope,
    y_centers,
    y_faces,
    K_faces,
    H_faces,
    psi_in,
    phi_soil,
):
    """
    Evolve psi along one frozen-u characteristic exactly in the linear
    reaction coefficient representation.

    For departure y_d >= 0:
        psi_t = exp(-(K_t-K_d))*psi_d
                + phi_soil*exp(-K_t)*(H_t-H_d)

    For a characteristic entering through y=0:
        psi_t = exp(-K_t)*(psi_in + phi_soil*H_t)
    """
    yd = np.asarray(y_departure, dtype=float)
    yt = yd + dt

    Kt = _interp_face_field(yt, y_faces, K_faces)
    Ht = _interp_face_field(yt, y_faces, H_faces)

    inside = yd >= 0.0
    out = np.empty_like(yd)

    if np.any(inside):
        yd_i = yd[inside]
        Kd = _interp_face_field(yd_i, y_faces, K_faces)
        Hd = _interp_face_field(yd_i, y_faces, H_faces)
        psi_d = _interp_psi(
            yd_i, psi, slope, y_centers, y_faces, psi_in
        )
        out[inside] = (
            np.exp(-(Kt[inside] - Kd)) * psi_d
            + phi_soil * np.exp(-Kt[inside]) * (Ht[inside] - Hd)
        )

    if np.any(~inside):
        # Before entering the erodible domain there is no soil source.
        out[~inside] = np.exp(-Kt[~inside]) * (
            psi_in + phi_soil * Ht[~inside]
        )

    return out


def slfv_reactive_phi(
    phi: np.ndarray,
    u: np.ndarray,
    k: np.ndarray,
    dx: float,
    dt: float,
    phi_soil: float,
    phi_in: float = 0.0,
) -> np.ndarray:
    """
    Conservative semi-Lagrangian finite-volume advection-reaction operator.

    It solves, during one frozen-u step,

        partial_t phi + partial_x(u phi)
            = k (phi_soil - phi)

    by introducing the travel-time coordinate

        y(x) = integral_0^x dxi/u(x)

    and the conservative variable

        psi(y,t) = u(x) phi(x,t).

    The transformed equation is

        partial_t psi + partial_y psi
            = k(y) [u(y) phi_soil - psi].

    The linear reaction is integrated analytically along each characteristic
    through cumulative K and H integrals. Target-cell mass is then obtained
    by quadrature over its departure interval. Therefore the source is applied
    continuously along the characteristic rather than as a lumped post-step
    source.

    Current scope:
      - one-dimensional positive-flow reference configuration u > 0;
      - frozen hydraulic field during one time step;
      - MUSCL/Van-Leer reconstruction in travel-time coordinate;
      - prescribed Dirichlet inlet concentration;
      - outflow at the downstream boundary;
      - arbitrary Courant number.
    """
    phi = np.asarray(phi, dtype=float)
    u = np.asarray(u, dtype=float)
    k = np.asarray(k, dtype=float)

    if phi.ndim != 1 or u.ndim != 1 or k.ndim != 1:
        raise ValueError("phi, u and k must be one-dimensional.")
    if not (phi.size == u.size == k.size):
        raise ValueError("phi, u and k must have equal size.")
    if dx <= 0.0 or dt < 0.0:
        raise ValueError("Require dx > 0 and dt >= 0.")
    if np.any(k < 0.0):
        raise ValueError("Reactive source coefficient k must be non-negative.")
    if not (0.0 <= phi_in < phi_soil):
        raise ValueError("Require 0 <= phi_in < phi_soil.")

    n = phi.size
    if n == 0 or dt == 0.0:
        return phi.copy()

    y_faces, y_centers, dy = _travel_time_faces(u, dx)

    # Cell-integrated physical mass is M_i = phi_i * dx.
    # The conservative transformed variable is psi = u*phi and
    # M_i = integral_{cell in y} psi dy.
    psi = phi * dx / dy
    psi_in = float(u[0] * phi_in)

    slope = _reconstruct_psi(psi, y_centers)
    K_faces, H_faces = _prepare_source_integrals(y_faces, k, u)

    # Three-point Gauss-Legendre quadrature in each target y-cell.
    nodes = np.array(
        [-np.sqrt(3.0 / 5.0), 0.0, np.sqrt(3.0 / 5.0)]
    )
    weights = np.array([5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0])

    new_mass = np.empty(n, dtype=float)

    for i in range(n):
        a = y_faces[i] - dt
        b = y_faces[i + 1] - dt
        mid = 0.5 * (a + b)
        half = 0.5 * (b - a)
        yd = mid + half * nodes

        psi_target = _characteristic_reaction(
            yd,
            dt,
            psi,
            slope,
            y_centers,
            y_faces,
            K_faces,
            H_faces,
            psi_in,
            phi_soil,
        )
        new_mass[i] = half * np.sum(weights * psi_target)

    phi_new = new_mass / dx

    # Numerical roundoff can produce tiny negative values. The clipping here
    # is a physical admissibility safeguard, not a stability mechanism.
    phi_new = np.clip(phi_new, 0.0, phi_soil * (1.0 - 1e-12))

    if not np.all(np.isfinite(phi_new)):
        raise FloatingPointError(
            "Non-finite value generated by reactive SLFV transport."
        )

    return phi_new
