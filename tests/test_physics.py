import numpy as np

from src.config import Params
from src.physics import (
    fm_julien,
    fm_julien_raw,
    conservative_solid_source,
    mdot_erosion,
    radius_growth_rate,
    rho_soil_sat,
)


def test_reference_transport_and_limiter_defaults():
    params = Params()
    assert params.phi_form == "conservative_area"
    assert params.phi_limiter == "minmod"


def test_julien_smoothing_is_active_and_bounded():
    params = Params(fm_max=5.0)
    phi = np.array([0.0, 0.1, 0.5, params.phi_soil * (1.0 - 1e-6)])
    raw = fm_julien_raw(phi, params)
    effective = fm_julien(phi, params)

    assert np.all(np.isfinite(raw))
    assert np.all(np.isfinite(effective))
    assert np.all(effective >= 1.0)
    assert np.all(effective <= params.fm_max)
    assert effective[-1] > 0.99 * params.fm_max
    assert raw[-1] > params.fm_max


def test_erosion_flux_has_explicit_radial_rate_conversion():
    tau = np.array([-20.0, -2.0])
    mdot = mdot_erosion(tau, tau_c=10.0, k_er=0.5)
    assert np.allclose(mdot, [5.0, 0.0])

    rho = rho_soil_sat(1000.0, 2650.0, 0.62)
    assert np.allclose(radius_growth_rate(mdot, rho), mdot / rho)


def test_area_increment_source_conserves_intact_solid_fraction():
    old_area = np.array([1.0, 2.0])
    new_area = np.array([1.5, 2.25])
    source = conservative_solid_source(old_area, new_area, 0.62)
    assert np.allclose(source, 0.62 * (new_area - old_area))
