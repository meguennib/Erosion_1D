"""
Compare Eulerian v1, lumped-source SLFV and characteristic-reactive SLFV.

Example:
    python run_slfv_reactive_audit.py --t-end-factor 0.02
    python run_slfv_reactive_audit.py --t-end-factor 15.7
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

from src.config import Params
from src.simulation import run_simulation
from src.utils import setup_logger


def run_case(method, args, root):
    p = Params(
        transport_method=method,
        t_end_factor=args.t_end_factor,
        morph_rel_change=args.morph_rel_change,
    )
    out = root / method
    out.mkdir(parents=True, exist_ok=True)
    logger = setup_logger(out / "run.log")
    result = run_simulation(p, logger)

    rows = result["timeseries"]
    if rows:
        with open(out / "timeseries.csv", "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=rows[0].keys())
            writer.writeheader()
            writer.writerows(rows)

    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--t-end-factor", type=float, default=0.02)
    parser.add_argument("--morph-rel-change", type=float, default=0.01)
    parser.add_argument("--out", default="output_slfv_reactive_audit")
    args = parser.parse_args()

    root = Path(args.out)
    root.mkdir(parents=True, exist_ok=True)

    methods = ("eulerian", "slfv", "slfv-reactive")
    results = {m: run_case(m, args, root) for m in methods}

    print("\n=== REACTIVE SLFV AUDIT ===")
    for method, result in results.items():
        m = result["run_manifest"]
        print(
            f"{method:15s} | steps={m['n_steps']:7d} | "
            f"t/ter={m['t_final_over_ter']:.6g} | "
            f"CFLmax={m['cfl_max']:.6g} | stop={m['stop_reason']}"
        )


if __name__ == "__main__":
    main()
