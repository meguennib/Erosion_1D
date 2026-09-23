
"""
Revision runner: regenerate Tables 3, 4, and 5 with corrected parameters.
- tau_b coefficient: 1/8 (already in shear_tau_b)
- dp = 75 µm, fm_max = 5.0 (updated in config defaults)
- R/L validity check (already in simulation loop)
"""

import sys, os, json, time as time_module
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from src.config import Params
from src.simulation import run_simulation
from src.utils import make_run_dir, setup_logger, save_config, save_timeseries_csv, save_snapshots_npz, save_json

# Common parameters matching manuscript reference case
COMMON = dict(
    L=0.117, R0=0.003, Nx=600,
    rho_w=1000.0, rho_p=2650.0,
    mu_w=1e-3, phi_soil=0.62, phi_in=0.0,
    dp=75e-6, fm_max=5.0, cB=0.2, cl=0.07,
    julien_lambda_power=2.0, lm_mode="clR0",
    Re_min=1.0, R_min=1e-6, R_break=0.5,
    phi_scheme="muscl", phi_form="conservative",
    dt_mode="adaptive", dt_min=1e-6, dt_max=10.0, CFL=0.35,
    t_end_factor=5.0,  # run to 5 * t_er
    Q_tol_abs=1e-3, Q_bisect_max_iter=80,
    Q_bracket_max_expand=60, Q_bracket_growth=2.0,
    snap_t_over_ter=[0.0, 2.0, 3.0, 4.0, 5.0],
    log_every_steps=100,
    save_ts_every=1,
)

SCENARIOS = {
    # Table 3: Soil A (erodible), K_out=0
    "SolA_Kout0": dict(
        Pin=4905.0, Pout=0.0, K_out=0.0,
        k_er=0.01, tau_c=1.0,
        description="Soil A (erodible silt-clay), K_out=0, Delta_p=0.5m"
    ),
    # Table 3: Soil B (resistant), K_out=0
    "SolB_Kout0": dict(
        Pin=4905.0, Pout=0.0, K_out=0.0,
        k_er=0.001, tau_c=10.0,
        description="Soil B (resistant silt-clay), K_out=0, Delta_p=0.5m"
    ),
    # Table 4: Soil A, K_out=10
    "SolA_Kout10": dict(
        Pin=4905.0, Pout=0.0, K_out=10.0,
        k_er=0.01, tau_c=1.0,
        description="Soil A (erodible), K_out=10, Delta_p=0.5m"
    ),
    # Table 5: Soil A, K_out=0, Delta_p=0.25m
    "SolA_Kout0_dp0.25": dict(
        Pin=2452.5, Pout=0.0, K_out=0.0,
        k_er=0.01, tau_c=1.0,
        description="Soil A, K_out=0, Delta_p=0.25m"
    ),
    # Table 5: Soil A, K_out=0, Delta_p=1.0m
    "SolA_Kout0_dp1.00": dict(
        Pin=9810.0, Pout=0.0, K_out=0.0,
        k_er=0.01, tau_c=1.0,
        description="Soil A, K_out=0, Delta_p=1.0m"
    ),
}

results = {}
for name, overrides in SCENARIOS.items():
    print(f"\n{'='*60}")
    print(f"Running: {name}")
    print(f"  {overrides['description']}")
    print(f"{'='*60}")
    
    cfg = dict(COMMON, **overrides)
    p = Params(**{k: v for k, v in cfg.items() if k in Params.__dataclass_fields__})
    
    # Account for extra fields not in Params dataclass
    for k, v in cfg.items():
        if not hasattr(p, k):
            setattr(p, k, v)
    
    run_dir = make_run_dir(base_dir="output_revision", version="v1_corrected", scenario=name)
    logger = setup_logger(run_dir / "run.log")
    save_config(run_dir / "config.json", p)
    
    t0 = time_module.time()
    try:
        res = run_simulation(p, logger)
        elapsed = time_module.time() - t0
        
        # Extract key metrics
        info = res["run_manifest"]
        ts = res["timeseries"]
        snaps = res["snapshots"]
        ter = res["t_er"]
        t = res["t_end"]
        
        # We need to get the last R array from timeseries or snaps
        # Let's get R_out_final from the last timeseries row
        R_out_final = float(ts[-1]["R_out"]) if len(ts) > 0 else p.R0
        t_end = t
        R_out_R0 = R_out_final / p.R0
        
        # Find doubling time t2: time when R_out >= 2*R0
        t2 = None
        for row in ts:
            if row.get("R_out", 0) >= 2.0 * p.R0:
                t2 = row["t"]
                break
        
        results[name] = {
            "ter": ter,
            "t_end": t_end,
            "t_end_over_ter": t_end / ter if ter else None,
            "R_out_final": R_out_final,
            "R_out_R0": R_out_R0,
            "t2": t2,
            "t2_over_ter": t2 / ter if ter and t2 else None,
            "n_steps": info.get("n_steps", 0),
            "stop_reason": info.get("stop_reason", "unknown"),
            "R_out_L": R_out_final / p.L,
            "t_star_1D": info.get("t_star_1D", None),
            "elapsed_s": elapsed,
        }
        
        # Save results
        save_timeseries_csv(run_dir / "timeseries.csv", ts)
        save_snapshots_npz(run_dir / "snapshots.npz", snaps)
        
        manifest = {
            "name": name,
            "description": overrides["description"],
            "results": results[name],
            "params": cfg,
        }
        save_json(run_dir / "run_manifest.json", manifest)
        
        print(f"  Done in {elapsed:.1f}s")
        print(f"  t_er={ter:.2f}s, t_end/ter={t_end/ter:.3f}")
        print(f"  R_out_final={R_out_final*1000:.2f}mm, R_out/R0={R_out_R0:.2f}")
        print(f"  t2={t2:.2f}s" if t2 else "  t2: NOT REACHED")
        print(f"  R_out/L={R_out_final/p.L:.3f}, t*_1D={info.get('t_star_1D','N/A')}")
        
    except Exception as e:
        elapsed = time_module.time() - t0
        print(f"  FAILED after {elapsed:.1f}s: {e}")
        import traceback
        traceback.print_exc()
        results[name] = {"error": str(e)}

# Print summary table
print("\n" + "=" * 80)
print("SUMMARY: CORRECTED SIMULATION RESULTS")
print("=" * 80)
print(f"{'Case':<25} {'t_er(s)':>10} {'t2(s)':>10} {'t_end/ter':>10} {'R_out/R0':>10} {'R_out/L':>10} {'Stop':>15}")
print("-" * 80)
for name, r in results.items():
    if "error" in r:
        print(f"{name:<25} {'ERROR: ' + r['error'][:40]}")
    else:
        print(f"{name:<25} {r['ter']:>10.1f} {r['t2'] if r['t2'] else 'N/A':>10} "
              f"{r['t_end_over_ter']:>10.3f} {r['R_out_R0']:>10.2f} "
              f"{r['R_out_L']:>10.3f} {r['stop_reason']:>15}")

# Save summary
summary = {name: r for name, r in results.items()}
import json
with open("output_revision/summary.json", "w") as f:
    json.dump(summary, f, indent=2)
print("\nSummary saved to output_revision/summary.json")
