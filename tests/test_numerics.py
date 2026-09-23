import numpy as np

from src.numerics import advection_conservative, limiter_minmod


def test_minmod_is_zero_at_extrema():
    d_left = np.array([1.0, 1.0, -1.0, -1.0])
    d_right = np.array([1.0, -1.0, 1.0, -1.0])
    slope = limiter_minmod(d_left, d_right)
    assert np.allclose(slope, [1.0, 0.0, 0.0, -1.0])


def test_constant_conservative_state_is_preserved():
    q = np.full(32, 0.4)
    a = np.full(32, 0.7)
    out = advection_conservative(
        q,
        a,
        dx=0.1,
        dt=0.02,
        q_in=0.4,
        scheme="muscl",
        limiter="minmod",
    )
    assert np.allclose(out, q)


def test_zero_velocity_preserves_conservative_state():
    q = np.linspace(0.0, 1.0, 32)
    out = advection_conservative(
        q,
        np.zeros_like(q),
        dx=0.1,
        dt=0.5,
        q_in=0.0,
        scheme="muscl",
        limiter="minmod",
    )
    assert np.allclose(out, q)


def test_inlet_boundary_uses_conservative_value():
    q = np.zeros(20)
    a = np.ones(20)
    out = advection_conservative(
        q,
        a,
        dx=0.1,
        dt=0.02,
        q_in=1.0,
        scheme="upwind",
        limiter="minmod",
    )
    assert out[0] > 0.0
    assert np.all(out >= 0.0)
