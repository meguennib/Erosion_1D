"""Official campaign runner for the revised 1D erosion model.

The module is intentionally side-effect free when imported.  Use ``python
main.py`` or ``python run_revision.py`` to execute the reference cases.
"""

from __future__ import annotations

import argparse
import json
import time as time_module
from pathlib import Path
from typing import Dict, Iterable, Optional

from src.config import Params
from src.simulation import run_simulation
from src.utils import (
    make_run_dir,
    save_config,
    save_json,
    save_snapshots_npz,
    save_timeseries_csv,
    setup_logger,
)
from src.validation import validate_signatures


COMMON = dict(
    L=0.117,
    R0=0.003,
    Nx=600,
    rho_w=1000.0,
    rho_p=2650.0,
    mu_w=1e-3,
    phi_soil=0.62,
    phi_in=0.0,
    dp=75e-6,
    fm_max=5.0,
    cB=0.2,
    cl=0.07,
    julien_lambda_power=2.0,
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
    CFL=0.35,
    t_end_factor=5.0,
    Q_tol_abs=1e-3,
    Q_bisect_max_iter=80,
    Q_bracket_max_expand=60,
    Q_bracket_growth=2.0,
    snap_t_over_ter=[0.0, 2.0, 3.0, 4.0, 5.0],
    log_every_steps=100,
    save_ts_every=1,
)


SCENARIOS = {
    "SolA_Kout0": dict(
        Pin=4905.0,
        Pout=0.0,
        K_out=0.0,
        k_er=0.01,
        tau_c=1.0,
        description="Soil A (erodible silt-clay), K_out=0, Delta_p=0.5m",
    ),
    "SolB_Kout0": dict(
        Pin=4905.0,
        Pout=0.0,
        K_out=0.0,
        k_er=0.001,
        tau_c=10.0,
        description="Soil B (resistant silt-clay), K_out=0, Delta_p=0.5m",
    ),
    "SolA_Kout10": dict(
        Pin=4905.0,
        Pout=0.0,
        K_out=10.0,
        k_er=0.01,
        tau_c=1.0,
        description="Soil A (erodible), K_out=10, Delta_p=0.5m",
    ),
    "SolA_Kout0_dp0.25": dict(
        Pin=2452.5,
        Pout=0.0,
        K_out=0.0,
        k_er=0.01,
        tau_c=1.0,
        description="Soil A, K_out=0, Delta_p=0.25m",
    ),
    "SolA_Kout0_dp1.00": dict(
        Pin=9810.0,
        Pout=0.0,
        K_out=0.0,
        k_er=0.01,
        tau_c=1.0,
        description="Soil A, K_out=0, Delta_p=1.0m",
    ),
}


def build_params(overrides: Dict) -> Params:
    """Build one validated parameter object from the campaign dictionaries."""
    cfg = dict(COMMON, **overrides)
    known = set(Params.__dataclass_fields__)
    return Params(**{key: value for key, value in cfg.items() if key in known})


def _extract_summary(res: dict, p: Params, elapsed: float) -> dict:
    info = res["run_manifest"]
    ts = res["timeseries"]
    ter = float(res["t_er"])
    t_end = float(res["t_end"])
    R_out_final = float(ts[-1]["R_out"]) if ts else p.R0

    t2 = None
    for row in ts:
        if float(row.get("R_out", 0.0)) >= 2.0 * p.R0:
            t2 = float(row["t"])
            break

    return {
        "ter": ter,
        "t_end": t_end,
        "t_final": float(info.get("t_final", 0.0)),
        "t_end_over_ter": float(info.get("t_final_over_ter", 0.0)),
        "R_out_final": R_out_final,
        "R_out_R0": R_out_final / p.R0,
        "t2": t2,
        "t2_over_ter": (t2 / ter if t2 is not None and ter > 0.0 else None),
        "n_steps": int(info.get("n_steps", 0)),
        "stop_reason": info.get("stop_reason", "unknown"),
        "R_out_L": R_out_final / p.L,
        "t_star_1D": info.get("t_star_1D"),
        "elapsed_s": float(elapsed),
    }


def run_case(name: str, overrides: Dict, output_dir: str = "output_revision") -> dict:
    """Run and persist one scenario; return its machine-readable summary."""
    description = str(overrides["description"])
    p = build_params(overrides)
    run_dir = make_run_dir(
        base_dir=output_dir,
        version="v2_area_conservative",
        scenario=name,
    )
    logger = setup_logger(run_dir / "run.log")
    save_config(run_dir / "config.json", p)

    print(f"\n{'=' * 70}")
    print(f"Running: {name}")
    print(f"  {description}")
    print(f"{'=' * 70}")

    t0 = time_module.time()
    try:
        res = run_simulation(p, logger)
        elapsed = time_module.time() - t0
        metrics = validate_signatures(res["snapshots"], p, logger)
        summary = _extract_summary(res, p, elapsed)

        save_timeseries_csv(run_dir / "timeseries.csv", res["timeseries"])
        # results.csv is the canonical name expected by the campaign plotting
        # tools; timeseries.csv is retained for backwards compatibility.
        save_timeseries_csv(run_dir / "results.csv", res["timeseries"])
        save_snapshots_npz(run_dir / "snapshots.npz", res["snapshots"])

        manifest = dict(res["run_manifest"])
        manifest.update(
            {
                "name": name,
                "description": description,
                "params": p.to_dict(),
                "results": summary,
                "validation": metrics,
            }
        )
        save_json(run_dir / "validation_metrics.json", metrics)
        save_json(run_dir / "run_manifest.json", manifest)

        print(f"  Done in {elapsed:.1f}s")
        print(f"  t_er={summary['ter']:.2f}s, t_final/ter={summary['t_end_over_ter']:.3f}")
        print(
            f"  R_out_final={summary['R_out_final'] * 1000:.2f}mm, "
            f"R_out/R0={summary['R_out_R0']:.2f}"
        )
        print(f"  stop_reason={summary['stop_reason']}")
        return summary

    except Exception as exc:
        elapsed = time_module.time() - t0
        logger.exception("Simulation failed")
        failure = {
            "error": str(exc),
            "elapsed_s": float(elapsed),
            "stop_reason": "exception",
        }
        save_json(
            run_dir / "run_manifest.json",
            {
                "name": name,
                "description": description,
                "params": p.to_dict(),
                "results": failure,
            },
        )
        print(f"  FAILED after {elapsed:.1f}s: {exc}")
        return failure


def main(argv: Optional[Iterable[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--case",
        action="append",
        choices=sorted(SCENARIOS),
        help="run only the selected case; may be repeated",
    )
    parser.add_argument(
        "--output-dir",
        default="output_revision",
        help="base directory for run outputs",
    )
    args = parser.parse_args(list(argv) if argv is not None else None)

    names = args.case if args.case else list(SCENARIOS)
    results = {
        name: run_case(name, SCENARIOS[name], output_dir=args.output_dir)
        for name in names
    }

    output_path = Path(args.output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    with open(output_path / "summary.json", "w", encoding="utf-8") as f:
        json.dump(results, f, indent=2, ensure_ascii=False)

    print("\n" + "=" * 100)
    print("SUMMARY: AREA-CONSERVATIVE REFERENCE SIMULATIONS")
    print("=" * 100)
    print(
        f"{'Case':<25} {'t_final/ter':>12} {'R_out/R0':>10} "
        f"{'R_out/L':>10} {'Stop':>24}"
    )
    print("-" * 100)
    for name, result in results.items():
        if "error" in result:
            print(f"{name:<25} ERROR: {result['error'][:60]}")
        else:
            print(
                f"{name:<25} {result['t_end_over_ter']:>12.3f} "
                f"{result['R_out_R0']:>10.2f} {result['R_out_L']:>10.3f} "
                f"{result['stop_reason']:>24}"
            )
    print(f"\nSummary saved to {output_path / 'summary.json'}")
    return 0 if all("error" not in result for result in results.values()) else 1


if __name__ == "__main__":
    raise SystemExit(main())
