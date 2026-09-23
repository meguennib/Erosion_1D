import logging

import numpy as np

from src.config import Params
from src.simulation import run_simulation


def _null_logger():
    logger = logging.getLogger("test_simulation")
    logger.handlers.clear()
    logger.addHandler(logging.NullHandler())
    logger.propagate = False
    return logger


def test_short_area_conservative_run():
    params = Params(
        L=0.1,
        Nx=20,
        R0=0.0032,
        Pin=2300.0,
        Pout=0.0,
        phi_in=0.2,
        tau_c=1.0e9,
        t_end_factor=1.0e-4,
        dt_max=0.1,
        snap_t_over_ter=[0.0, 1.0e-4],
        log_every_steps=100000,
        save_ts_every=1,
    )
    result = run_simulation(params, _null_logger())
    manifest = result["run_manifest"]

    assert manifest["stop_reason"] == "completed"
    assert manifest["transport_variable"] == "S=A*phi"
    assert manifest["phi_limiter_effective"] == "minmod"
    assert manifest["t_final"] <= result["t_end"] + 1e-12
    assert len(result["snapshots"]) == 2

    for snapshot in result["snapshots"]:
        assert np.allclose(snapshot["S"], snapshot["A"] * snapshot["phi"])
        assert np.all(snapshot["phi"] >= 0.0)
        assert np.all(snapshot["phi"] <= params.phi_soil + 1e-12)
        assert np.all(snapshot["fm"] <= params.fm_max + 1e-12)

    for row in result["timeseries"]:
        assert abs(row["res_pL"]) <= params.Q_tol_abs
        assert row["max_fm"] <= params.fm_max + 1e-12


def test_zero_erosion_keeps_radius_fixed():
    params = Params(
        L=0.1,
        Nx=12,
        R0=0.003,
        Pin=2300.0,
        tau_c=1.0e9,
        t_end_factor=2.0e-5,
        dt_max=0.1,
        snap_t_over_ter=[0.0, 2.0e-5],
        log_every_steps=100000,
        save_ts_every=1,
    )
    result = run_simulation(params, _null_logger())
    for snapshot in result["snapshots"]:
        assert np.allclose(snapshot["R"], params.R0)
