import numpy as np


def _vanleer_slope(d_left, d_right):
    """Van-Leer limited slope expressed as a cell-to-cell difference."""
    d_left = np.asarray(d_left, dtype=float)
    d_right = np.asarray(d_right, dtype=float)
    out = np.zeros_like(d_left)
    mask = d_left * d_right > 0.0
    out[mask] = (
        2.0 * d_left[mask] * d_right[mask]
        / (d_left[mask] + d_right[mask])
    )
    return out


def _slopes(phi, phi_in):
    """MUSCL Van-Leer slopes with inlet Dirichlet/outlet zero-gradient ghosts."""
    n = len(phi)
    g = np.empty(n + 4, dtype=float)
    g[2:-2] = phi
    g[:2] = float(phi_in)
    g[-2:] = phi[-1]
    return _vanleer_slope(
        g[2:n + 2] - g[1:n + 1],
        g[3:n + 3] - g[2:n + 2],
    )


def _face_velocity(u):
    """Cell-centred to face velocity: arithmetic interior, zero-gradient ends."""
    u = np.asarray(u, dtype=float)
    uf = np.empty(len(u) + 1, dtype=float)
    uf[0] = u[0]
    uf[-1] = u[-1]
    if len(u) > 1:
        uf[1:-1] = 0.5 * (u[:-1] + u[1:])
    return uf


def _cell_travel_time(u_left, u_right, dx):
    """
    Exact travel time through a cell for linearly interpolated positive u(x).
    """
    scale = max(1.0, abs(u_left), abs(u_right))
    if u_left <= 0.0 or u_right <= 0.0:
        raise ValueError("SLFV reference operator requires strictly positive velocity.")

    du = u_right - u_left
    if abs(du) <= 1e-14 * scale:
        return dx / u_left

    return dx * np.log(u_right / u_left) / du


def _travel_time_faces(u_face, dx):
    """Cumulative steady travel time from x=0 to every physical face."""
    cell_times = np.array(
        [_cell_travel_time(u_face[i], u_face[i + 1], dx)
         for i in range(len(u_face) - 1)],
        dtype=float,
    )
    return np.r_[0.0, np.cumsum(cell_times)]


def _invert_travel_time(tau, T_faces, u_face, dx):
    """Invert T(x)=tau for a non-negative travel-time coordinate."""
    if tau <= 0.0:
        return 0.0
    if tau >= T_faces[-1]:
        return dx * (len(u_face) - 1)

    i = int(np.searchsorted(T_faces, tau, side="right") - 1)
    i = min(max(i, 0), len(u_face) - 2)
    local_tau = tau - T_faces[i]

    u_left = u_face[i]
    u_right = u_face[i + 1]
    scale = max(1.0, abs(u_left), abs(u_right))

    if abs(u_right - u_left) <= 1e-14 * scale:
        local_x = u_left * local_tau
    else:
        slope_u = (u_right - u_left) / dx
        local_x = u_left * np.expm1(slope_u * local_tau) / slope_u

    return i * dx + np.clip(local_x, 0.0, dx)


def _backtrace_faces(u, dx, dt):
    """
    Trace physical target faces backward along frozen characteristics.

    The current reference case assumes u>0. A departure face outside x=0 is
    represented by a signed upstream distance; its value is only used to
    determine the amount of prescribed inlet material entering the domain.
    """
    u_face = _face_velocity(u)
    if np.any(u_face <= 0.0):
        raise ValueError(
            "SLFV currently supports the positive-flow reference case only."
        )

    T_faces = _travel_time_faces(u_face, dx)
    dep = np.empty_like(T_faces)

    for j, target_tau in enumerate(T_faces):
        departure_tau = target_tau - dt
        if departure_tau >= 0.0:
            dep[j] = _invert_travel_time(
                departure_tau, T_faces, u_face, dx
            )
        else:
            dep[j] = departure_tau * u_face[0]

    return dep


def _primitive_piecewise_linear(phi, slope, dx, prefix_mass, x):
    """
    Primitive of the MUSCL reconstruction on [0,x].

    prefix_mass[i] is the exact reconstructed mass in cells 0,...,i-1.
    The prefix array makes each query O(1), avoiding the O(N^2) cost of
    repeatedly summing all previous cells.
    """
    n = len(phi)
    L = n * dx
    x = float(np.clip(x, 0.0, L))

    if x <= 0.0:
        return 0.0
    if x >= L:
        return float(prefix_mass[-1])

    i = min(int(x // dx), n - 1)
    local = x - i * dx
    return float(
        prefix_mass[i]
        + phi[i] * local
        + 0.5 * slope[i] / dx
        * ((local - 0.5 * dx) ** 2 - (0.5 * dx) ** 2)
    )


def _departure_mass(phi, slope, dx, prefix_mass, a, b, phi_in):
    """Reconstructed mass over the departure interval [a,b]."""
    if b <= a:
        return 0.0

    if b <= 0.0:
        return float(phi_in) * (b - a)

    if a >= 0.0:
        return (
            _primitive_piecewise_linear(phi, slope, dx, prefix_mass, b)
            - _primitive_piecewise_linear(phi, slope, dx, prefix_mass, a)
        )

    return (
        float(phi_in) * (-a)
        + _primitive_piecewise_linear(phi, slope, dx, prefix_mass, b)
    )


def slfv_advection_phi(phi, u, dx, dt, phi_in=0.0):
    """
    Conservative semi-Lagrangian finite-volume transport for u>0.

    The target-cell average is obtained by integrating the reconstructed
    old concentration over the backward-mapped departure interval.
    The method is not CFL-limited in the same way as the explicit Eulerian
    upwind/Godunov update.
    """
    phi = np.asarray(phi, dtype=float)
    u = np.asarray(u, dtype=float)

    if phi.ndim != 1 or u.ndim != 1 or len(phi) != len(u):
        raise ValueError("phi and u must be 1-D arrays of equal size.")
    if dx <= 0.0 or dt < 0.0:
        raise ValueError("Require dx > 0 and dt >= 0.")
    if len(phi) == 0 or dt == 0.0:
        return phi.copy()

    dep = _backtrace_faces(u, dx, dt)
    slope = _slopes(phi, phi_in)

    prefix_mass = np.r_[0.0, np.cumsum(phi * dx)]

    out = np.empty_like(phi)
    for i in range(len(phi)):
        out[i] = (
            _departure_mass(
                phi,
                slope,
                dx,
                prefix_mass,
                float(dep[i]),
                float(dep[i + 1]),
                float(phi_in),
            )
            / dx
        )

    if not np.all(np.isfinite(out)):
        raise FloatingPointError("Non-finite value generated by SLFV advection.")

    return out
