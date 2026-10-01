"""
Run the decisive v1-vs-SLFV audit.

Examples:
    python run_slfv_audit.py --t-end-factor 0.02 --morph-rel-change 0.01
    python run_slfv_audit.py --t-end-factor 15.7 --morph-rel-change 0.01

The script never changes the thesis-v1 branch. It runs both transport
operators from the same physical parameters on the experimental branch.
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

from src.config import Params
from src.simulation import run_simulation
from src.utils import setup_logger


def run_case(method: str, t_end_factor: float, morph_rel_change: float, out: Path):
    p = Params(
        transport_method=method,
        t_end_factor=t_end_factor,
        morph_rel_change=morph_rel_change,
    )
    case_dir = out / method
    case_dir.mkdir(parents=True, exist_ok=True)
    logger = setup_logger(case_dir / "run.log")
    result = run_simulation(p, logger)

    rows = result["timeseries"]
    if rows:
        keys = list(rows[0].keys())
        with open(case_dir / "timeseries.csv", "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=keys)
            w.writeheader()
            w.writerows(rows)

    return result


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--t-end-factor", type=float, default=0.02)
    ap.add_argument("--morph-rel-change", type=float, default=0.01)
    ap.add_argument("--out", default="output_slfv_audit")
    args = ap.parse_args()

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    results = {}
    for method in ("eulerian", "slfv"):
        results[method] = run_case(
            method,
            args.t_end_factor,
            args.morph_rel_change,
            out,
        )

    print("\\n=== SLFV AUDIT ===")
    for method, result in results.items():
        m = result["run_manifest"]
        print(
            f"{method:9s} | steps={m['n_steps']:6d} | "
            f"t/ter={m['t_final_over_ter']:.6g} | "
            f"CFLmax={m['cfl_max']:.6g} | "
            f"stop={m['stop_reason']}"
        )

    sl = results["slfv"]["timeseries"]
    if sl:
        active = []
        for r in sl:
            vals = {
                "dt_adv": r.get("dt_adv", float("nan")),
                "dt_src": r.get("dt_src", float("nan")),
                "dt_morph": r.get("dt_morph", float("nan")),
            }
            finite = {k: v for k, v in vals.items() if v == v}
            if finite:
                active.append(min(finite, key=finite.get))
        from collections import Counter
        print("SLFV limiting constraint:", Counter(active).most_common())


if __name__ == "__main__":
    main()
