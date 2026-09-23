"""Profile one representative simulation with cProfile.

Example:
    python tools/profile_simulation.py --nx 100 --t-end-factor 0.001
"""

from __future__ import annotations

import argparse
import cProfile
import logging
import pstats
import time
from pathlib import Path
from typing import List, Optional

from src.config import Params
from src.simulation import run_simulation


def build_params(nx: int, t_end_factor: float) -> Params:
    return Params(
        L=0.117,
        Nx=nx,
        R0=0.003,
        Pin=4905.0,
        Pout=0.0,
        rho_w=1000.0,
        rho_p=2650.0,
        mu_w=1e-3,
        phi_soil=0.62,
        phi_in=0.0,
        dp=75e-6,
        fm_max=5.0,
        cB=0.2,
        cl=0.07,
        tau_c=1.0,
        k_er=0.01,
        dt_mode="adaptive",
        dt_min=1e-6,
        dt_max=10.0,
        CFL=0.35,
        t_end_factor=t_end_factor,
        snap_t_over_ter=[0.0, t_end_factor],
        log_every_steps=100000,
        save_ts_every=100000,
    )


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--nx", type=int, default=100)
    parser.add_argument("--t-end-factor", type=float, default=1e-3)
    parser.add_argument("--stats", type=Path, default=None)
    args = parser.parse_args(argv)

    params = build_params(args.nx, args.t_end_factor)
    logger = logging.getLogger("erosion-profile")
    logger.handlers.clear()
    logger.addHandler(logging.NullHandler())
    logger.propagate = False

    profiler = cProfile.Profile()
    start = time.perf_counter()
    profiler.enable()
    result = run_simulation(params, logger)
    profiler.disable()
    elapsed = time.perf_counter() - start

    print(
        f"elapsed={elapsed:.6f}s, steps={result['run_manifest']['n_steps']}, "
        f"backend={result['run_manifest']['hydraulic_kernel_backend']}"
    )
    if args.stats is not None:
        args.stats.parent.mkdir(parents=True, exist_ok=True)
        profiler.dump_stats(str(args.stats))
        print(f"profile data written to {args.stats}")

    pstats.Stats(profiler).strip_dirs().sort_stats("cumulative").print_stats(40)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
