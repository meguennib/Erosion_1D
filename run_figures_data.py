"""Re-run the figure campaigns and generate ``figures_data.npz``.

This module is side-effect free when imported.  It is a separate, expensive
campaign runner and is not the official default entry point; use
``python run_figures_data.py`` explicitly.
"""

from __future__ import annotations

import time as time_module
from pathlib import Path

import numpy as np

from src.config import Params
from src.simulation import run_simulation
from src.utils import setup_logger


SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR
FIG_DIR = REPO_ROOT / "Figures_HD"

COMMON = dict(
    L=100.0,
    Nx=600,
    R0=0.003,
    Pin=3330000.0,
    Pout=0.0,
    rho_w=1000.0,
    rho_p=2700.0,
    phi_soil=0.5,
    phi_in=0.0,
    mu_w=0.001,
    tau_c=10.0,
    k_er=0.01,
    dp=0.001,
    cl=0.1,
    cB=0.01,
    julien_lambda_power=2.0,
    fm_max=2000.0,
    lm_mode="clR0",
    Re_min=1.0,
    R_min=1e-6,
    R_break=0.5,
    phi_scheme="muscl",
    phi_form="conservative_area",
    phi_limiter="minmod",
    dt_mode="adaptive",
    dt_min=1e-6,
    dt_max=10.0,
    CFL=0.45,
    Q_tol_abs=1e-3,
    Q_bisect_max_iter=80,
    Q_bracket_max_expand=60,
    Q_bracket_growth=2.0,
    log_every_steps=500,
    save_ts_every=20,
)

SOIL_A_KOUT10 = dict(K_out=10.0)


def make_params(**kwargs):
    known = set(Params.__dataclass_fields__)
    unknown = sorted(set(kwargs) - known)
    if unknown:
        raise ValueError(f"Unknown simulation parameters: {unknown}")
    return Params(**kwargs)


def nearest_snap(snaps, target_over_ter, ter):
    best = None
    best_dist = float("inf")
    for snapshot in snaps:
        t_s = float(snapshot.get("t_actual", snapshot.get("t", 0.0)))
        dist = abs(t_s / ter - target_over_ter)
        if dist < best_dist:
            best_dist = dist
            best = snapshot
    if best is None:
        raise RuntimeError(f"No snapshot available near t/t_er={target_over_ter}")
    return best


def main() -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    print("\n" + "=" * 70)
    print("FIG 3 — Soil A, K_out=10, Nx=600  (snap at t/ter = 0, 2.2, 4.5)")
    print("=" * 70)

    cfg_f3 = dict(
        COMMON,
        **SOIL_A_KOUT10,
        t_end_factor=4.7,
        snap_t_over_ter=[0.0, 2.2, 4.5],
    )
    p_f3 = make_params(**cfg_f3)
    logger_f3 = setup_logger(FIG_DIR / "run_fig3.log")

    t0 = time_module.time()
    res_f3 = run_simulation(p_f3, logger_f3)
    print(f"\n  Completed in {time_module.time() - t0:.1f}s")

    snaps_f3 = res_f3["snapshots"]
    ter_f3 = res_f3["t_er"]
    x_f3 = res_f3["x"]

    snap_t0 = nearest_snap(snaps_f3, 0.0, ter_f3)
    snap_t22 = nearest_snap(snaps_f3, 2.2, ter_f3)
    snap_t45 = nearest_snap(snaps_f3, 4.5, ter_f3)
    x_norm_f3 = x_f3 / p_f3.L

    conv_data = {}
    for Nx in [300, 600, 1200]:
        print("\n" + "=" * 70)
        print(f"FIG 4 — Soil A, K_out=10, Nx={Nx}  (snap at t/ter = 2.2)")
        print("=" * 70)

        cfg_nx = dict(
            COMMON,
            **SOIL_A_KOUT10,
            Nx=Nx,
            t_end_factor=2.4,
            snap_t_over_ter=[2.2],
        )
        p_nx = make_params(**cfg_nx)
        logger_nx = setup_logger(FIG_DIR / f"run_fig4_Nx{Nx}.log")

        t0 = time_module.time()
        res_nx = run_simulation(p_nx, logger_nx)
        print(f"\n  Completed in {time_module.time() - t0:.1f}s")

        snap_22 = nearest_snap(res_nx["snapshots"], 2.2, res_nx["t_er"])
        x_nx = res_nx["x"]
        conv_data[Nx] = {
            "x": x_nx / p_nx.L,
            "R": np.asarray(snap_22["R"]) / p_nx.R0,
            "p": np.asarray(snap_22["p"]) / p_nx.Pin,
            "phi": np.asarray(snap_22["phi"]) / p_nx.phi_soil,
        }

    out_dict = {
        "x": x_norm_f3,
        "ref_t0_R": np.asarray(snap_t0["R"]) / p_f3.R0,
        "ref_t0_p": np.asarray(snap_t0["p"]) / p_f3.Pin,
        "ref_t0_phi": np.asarray(snap_t0["phi"]) / p_f3.phi_soil,
        "ref_t0_fm": np.asarray(snap_t0["fm"]),
        "ref_t22_R": np.asarray(snap_t22["R"]) / p_f3.R0,
        "ref_t22_p": np.asarray(snap_t22["p"]) / p_f3.Pin,
        "ref_t22_phi": np.asarray(snap_t22["phi"]) / p_f3.phi_soil,
        "ref_t22_fm": np.asarray(snap_t22["fm"]),
        "ref_t45_R": np.asarray(snap_t45["R"]) / p_f3.R0,
        "ref_t45_p": np.asarray(snap_t45["p"]) / p_f3.Pin,
        "ref_t45_phi": np.asarray(snap_t45["phi"]) / p_f3.phi_soil,
        "ref_t45_fm": np.asarray(snap_t45["fm"]),
        "conv_Nx300_x": conv_data[300]["x"],
        "conv_Nx300_R": conv_data[300]["R"],
        "conv_Nx300_p": conv_data[300]["p"],
        "conv_Nx300_phi": conv_data[300]["phi"],
        "conv_Nx600_x": conv_data[600]["x"],
        "conv_Nx600_R": conv_data[600]["R"],
        "conv_Nx600_p": conv_data[600]["p"],
        "conv_Nx600_phi": conv_data[600]["phi"],
        "conv_Nx1200_x": conv_data[1200]["x"],
        "conv_Nx1200_R": conv_data[1200]["R"],
        "conv_Nx1200_p": conv_data[1200]["p"],
        "conv_Nx1200_phi": conv_data[1200]["phi"],
        "R0": np.array([p_f3.R0]),
        "L": np.array([p_f3.L]),
        "Pin_ref": np.array([p_f3.Pin]),
        "phi_soil": np.array([p_f3.phi_soil]),
        "R_validity": np.array([1.95]),
        "R_L_0_1": np.array([0.1]),
    }

    outpath = FIG_DIR / "figures_data.npz"
    np.savez_compressed(outpath, **out_dict)
    print(f"\nDATA SAVED  →  {outpath}")


if __name__ == "__main__":
    main()
