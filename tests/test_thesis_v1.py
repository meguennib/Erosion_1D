import numpy as np

from src.config import Params
from src.physics import barenblatt_fw, mdot_erosion, rho_mix, rho_soil_sat


def test_saturated_soil_density():
    p = Params()
    expected = p.rho_w + p.phi_soil * (p.rho_p - p.rho_w)
    assert np.isclose(rho_soil_sat(p.rho_w, p.rho_p, p.phi_soil), expected)


def test_mixture_density_bounds():
    p = Params()
    rho = rho_mix(np.array([0.0, p.phi_soil]), p.rho_w, p.rho_p)
    assert np.isclose(rho[0], p.rho_w)
    assert np.isclose(rho[1], rho_soil_sat(p.rho_w, p.rho_p, p.phi_soil))


def test_erosion_threshold():
    p = Params()
    tau = np.array([p.tau_c - 1.0, p.tau_c, p.tau_c + 1.0])
    mdot = mdot_erosion(tau, p.tau_c, p.k_er)
    assert mdot[0] == 0.0
    assert mdot[1] == 0.0
    assert np.isclose(mdot[2], p.k_er)


def test_barenblatt_factor_is_finite():
    Re = np.array([1.0e3, 1.0e4, 1.0e5])
    f = barenblatt_fw(Re)
    assert np.all(np.isfinite(f))
    assert np.all(f > 0.0)
