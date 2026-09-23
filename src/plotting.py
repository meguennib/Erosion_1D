from __future__ import annotations

import csv
import json
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np

RUN_DIR_RE = re.compile(
    r"^\d{8}_\d{6}_(?P<versioncore>.+?)_(?P<stage>long|matrix)_(?P<scenario>.+)_Nx(?P<nx>\d+)_CFL(?P<cfl>[0-9.]+)$"
)


@dataclass
class RunRecord:
    run_dir: Path
    stage: str
    version_core: str
    version_stage: str
    scenario: str
    nx: int
    cfl: float
    config: Dict
    manifest: Dict
    metrics: Dict
    snapshots: List[Dict]
    results_rows: List[Dict]


def _scalarise(value: np.ndarray):
    arr = np.asarray(value)
    if arr.shape == ():
        return arr.item()
    if arr.size == 1:
        item = arr.reshape(-1)[0]
        if isinstance(item, bytes):
            return item.decode("utf-8")
        if hasattr(item, "item"):
            try:
                return item.item()
            except Exception:
                return item
        return item
    return arr.copy()


def load_snapshots_npz(path: Path) -> List[Dict]:
    data = np.load(path, allow_pickle=True)
    grouped: Dict[str, Dict] = {}
    for key in data.files:
        prefix, field = key.split("_", 1)
        grouped.setdefault(prefix, {})[field] = _scalarise(data[key])
    return [grouped[k] for k in sorted(grouped.keys(), key=lambda s: int(s[1:]))]


def load_csv_rows(path: Path) -> List[Dict]:
    with open(path, "r", encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f))


def parse_run_dir_name(name: str) -> Optional[Dict[str, object]]:
    m = RUN_DIR_RE.match(name)
    if not m:
        return None
    version_core = str(m.group("versioncore"))
    stage = str(m.group("stage"))
    return {
        "version_core": version_core,
        "version_stage": f"{version_core}_{stage}",
        "stage": stage,
        "scenario": str(m.group("scenario")),
        "nx": int(m.group("nx")),
        "cfl": float(m.group("cfl")),
    }


def discover_run_records(outdir: Path) -> List[RunRecord]:
    """Discover both campaign-style and reference-run output directories."""
    out = []
    for child in sorted(outdir.iterdir()):
        if not child.is_dir():
            continue
        req = [child / n for n in ["config.json", "results.csv", "snapshots.npz"]]
        if not all(p.exists() for p in req):
            continue

        with open(child / "config.json", "r", encoding="utf-8") as f:
            config = json.load(f)
        manifest_path = child / "run_manifest.json"
        metrics_path = child / "validation_metrics.json"
        if manifest_path.exists():
            with open(manifest_path, "r", encoding="utf-8") as f:
                manifest = json.load(f)
        else:
            manifest = {}
        if metrics_path.exists():
            with open(metrics_path, "r", encoding="utf-8") as f:
                metrics = json.load(f)
        else:
            metrics = {}

        parsed = parse_run_dir_name(child.name)
        if parsed is None:
            # run_revision.py deliberately uses a human-readable scenario
            # name rather than the matrix-campaign naming convention.  Its
            # manifest contains enough metadata to expose the same record
            # interface to the post-processing tools.
            scenario = str(manifest.get("name", child.name))
            version_core = str(manifest.get("version_core", "reference"))
            parsed = {
                "version_core": version_core,
                "version_stage": f"{version_core}_long",
                "stage": "long",
                "scenario": scenario,
                "nx": int(config.get("Nx", 0)),
                "cfl": float(config.get("CFL", 0.0)),
            }

        out.append(
            RunRecord(
                run_dir=child,
                stage=str(parsed["stage"]),
                version_core=str(parsed["version_core"]),
                version_stage=str(parsed["version_stage"]),
                scenario=str(parsed["scenario"]),
                nx=int(parsed["nx"]),
                cfl=float(parsed["cfl"]),
                config=config,
                manifest=manifest,
                metrics=metrics,
                snapshots=load_snapshots_npz(child / "snapshots.npz"),
                results_rows=load_csv_rows(child / "results.csv"),
            )
        )
    return out


def choose_latest_campaign(records: Sequence[RunRecord], scenario: Optional[str] = None) -> Tuple[RunRecord, List[RunRecord]]:
    filtered = [r for r in records if scenario is None or r.scenario == scenario]
    if not filtered:
        raise FileNotFoundError("No completed run directories were found.")
    groups: Dict[Tuple[str, str], List[RunRecord]] = {}
    for rec in filtered:
        groups.setdefault((rec.version_core, rec.scenario), []).append(rec)
    key, picked = max(groups.items(), key=lambda kv: (sum(1 for r in kv[1] if r.stage == "long") > 0, sum(1 for r in kv[1] if r.stage == "matrix"), max(r.run_dir.stat().st_mtime for r in kv[1])))
    longs = [r for r in picked if r.stage == "long"]
    matrices = sorted([r for r in picked if r.stage == "matrix"], key=lambda r: (r.nx, r.cfl))
    if not longs:
        raise FileNotFoundError(f"No long run found for campaign {key}.")
    return max(longs, key=lambda r: r.run_dir.stat().st_mtime), matrices


def _float_or_nan(v) -> float:
    try:
        return float(v)
    except Exception:
        return float("nan")


def _nearest_snapshot(snaps: Sequence[Dict], target_over_ter: float, require_reached: bool = True):
    best = None
    best_score = (float("inf"), 1.0)
    for s in snaps:
        tgt = s.get("t_target_over_ter")
        reached = bool(s.get("reached_target", True))
        if tgt is not None and abs(float(tgt) - float(target_over_ter)) < 1e-9:
            if require_reached and not reached:
                return s, "filled_missing_target"
            return s, "target"
        t = _float_or_nan(s.get("t", np.nan))
        ter = max(_float_or_nan(s.get("ter", np.nan)), 1e-30)
        score = (abs(t / ter - target_over_ter), 0.0 if reached else 1.0)
        if score < best_score:
            best = s
            best_score = score
    if best is None:
        return None, "missing"
    if require_reached and not bool(best.get("reached_target", True)):
        return best, "filled_missing_target"
    return best, "nearest"


def _initial_snapshot(snaps: Sequence[Dict]) -> Dict:
    return min(snaps, key=lambda s: _float_or_nan(s.get("t_actual", s.get("t", 0.0))))


def _save_figure_all_formats(fig: plt.Figure, stem: Path) -> Dict[str, str]:
    out = {}
    for ext in ("png", "pdf", "svg"):
        path = stem.with_suffix(f".{ext}")
        fig.savefig(path, bbox_inches="tight", dpi=200 if ext == "png" else None)
        out[ext] = str(path)
    plt.close(fig)
    return out


def _format_cell(v) -> str:
    if v is None:
        return ""
    if isinstance(v, float):
        return "NA" if math.isnan(v) else f"{v:.6g}"
    return str(v)


def _write_csv(path: Path, rows: Sequence[Dict]) -> None:
    if not rows:
        return
    with open(path, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)


def _write_md(path: Path, rows: Sequence[Dict]) -> None:
    if not rows:
        return
    cols = list(rows[0].keys())
    lines = ["| " + " | ".join(cols) + " |", "| " + " | ".join(["---"] * len(cols)) + " |"]
    for row in rows:
        lines.append("| " + " | ".join(_format_cell(row.get(c)) for c in cols) + " |")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _final_rout_over_r0(run: RunRecord) -> float:
    if not run.results_rows:
        return float("nan")
    return _float_or_nan(run.results_rows[-1].get("R_out", np.nan)) / max(_float_or_nan(run.config.get("R0", np.nan)), 1e-30)


def build_table1(matrix_runs: Sequence[RunRecord], out_dir: Path) -> List[Dict]:
    rows = []
    for run in sorted(matrix_runs, key=lambda r: (r.nx, r.cfl)):
        rows.append({
            "Nx": run.nx,
            "CFL_target": run.cfl,
            "CFL_mean": _float_or_nan(run.manifest.get("cfl_mean", np.nan)),
            "CFL_max": _float_or_nan(run.manifest.get("cfl_max", np.nan)),
            "clip_dtmax_fraction": _float_or_nan(run.manifest.get("clip_dtmax_fraction", np.nan)),
            "t_final_over_ter": _float_or_nan(run.manifest.get("t_final_over_ter", np.nan)),
            "Rout_over_R0_final": _final_rout_over_r0(run),
            "ratio_A": _float_or_nan(run.metrics.get("ratio_A", np.nan)),
            "ratio_C": _float_or_nan(run.metrics.get("ratio_C", np.nan)),
            "phi_max_2p2": _float_or_nan(run.metrics.get("phi_max_2p2", np.nan)),
            "phi_end_max_15p7": _float_or_nan(run.metrics.get("phi_end_max_15p7", np.nan)),
            "target_2p2_status": run.metrics.get("target_2p2_status", ""),
            "target_15p7_status": run.metrics.get("target_15p7_status", ""),
            "stop_reason": run.manifest.get("stop_reason", ""),
        })
    _write_csv(out_dir / "table1_cfl_diagnostics.csv", rows)
    _write_md(out_dir / "table1_cfl_diagnostics.md", rows)
    return rows


def _relative_errors(y: np.ndarray, y_ref_on_y: np.ndarray):
    diff = y - y_ref_on_y
    l2 = np.sqrt(np.mean(diff ** 2)) / max(np.sqrt(np.mean(y_ref_on_y ** 2)), 1e-30)
    linf = np.max(np.abs(diff)) / max(np.max(np.abs(y_ref_on_y)), 1e-30)
    return float(l2), float(linf)


def build_table2(matrix_runs: Sequence[RunRecord], out_dir: Path, *, target_over_ter: float, reference_nx: int = 1200, reference_cfl: float = 0.2) -> List[Dict]:
    ref_run = next((r for r in matrix_runs if r.nx == reference_nx and abs(r.cfl - reference_cfl) < 1e-12), None)
    if ref_run is None:
        raise FileNotFoundError("Reference run Nx=1200, CFL=0.2 not found.")
    ref_snap, ref_status = _nearest_snapshot(ref_run.snapshots, target_over_ter, require_reached=True)
    if ref_snap is None or ref_status == "filled_missing_target":
        raise RuntimeError(f"Reference run did not truly reach t/t_er={target_over_ter}.")
    x_ref = np.asarray(ref_snap["x"], dtype=float)
    p_ref = np.asarray(ref_snap["p"], dtype=float)
    r_ref = np.asarray(ref_snap["R"], dtype=float)
    phi_ref = np.asarray(ref_snap["phi"], dtype=float)
    rows = []
    for run in sorted(matrix_runs, key=lambda r: (r.nx, r.cfl)):
        snap, status = _nearest_snapshot(run.snapshots, target_over_ter, require_reached=True)
        row = {"Nx": run.nx, "CFL": run.cfl, "target_over_ter": target_over_ter, "status": status, "p_L2": float("nan"), "p_Linf": float("nan"), "R_L2": float("nan"), "R_Linf": float("nan"), "phi_L2": float("nan"), "phi_Linf": float("nan")}
        if snap is not None and status != "filled_missing_target":
            x = np.asarray(snap["x"], dtype=float)
            p = np.asarray(snap["p"], dtype=float)
            r = np.asarray(snap["R"], dtype=float)
            phi = np.asarray(snap["phi"], dtype=float)
            row["p_L2"], row["p_Linf"] = _relative_errors(p, np.interp(x, x_ref, p_ref))
            row["R_L2"], row["R_Linf"] = _relative_errors(r, np.interp(x, x_ref, r_ref))
            row["phi_L2"], row["phi_Linf"] = _relative_errors(phi, np.interp(x, x_ref, phi_ref))
        rows.append(row)
    suffix = str(target_over_ter).replace(".", "p")
    _write_csv(out_dir / f"table2_{suffix}_relative_errors.csv", rows)
    _write_md(out_dir / f"table2_{suffix}_relative_errors.md", rows)
    return rows


def plot_figure3_profiles(long_run: RunRecord, out_dir: Path) -> Dict[str, str]:
    p = long_run.config
    snaps = long_run.snapshots
    init_snap = _initial_snapshot(snaps)
    tau0 = max(float(np.mean(np.abs(np.asarray(init_snap["tau_b"], dtype=float)))), 1e-30)
    fig, axes = plt.subplots(4, 3, figsize=(13.5, 10), sharex=True)
    labels = ["R/R0", "p/Pin", "phi/phi_soil", "|tau_b|/|tau_b0|"]
    for j, target in enumerate([0.0, 2.2, 4.5]):
        snap, status = (init_snap, "initial") if target == 0.0 else _nearest_snapshot(snaps, target, require_reached=False)
        if snap is None:
            raise RuntimeError(f"Missing snapshot for t/t_er={target}.")
        x = np.asarray(snap["x"], dtype=float) / float(p["L"])
        ys = [
            np.asarray(snap["R"], dtype=float) / float(p["R0"]),
            np.asarray(snap["p"], dtype=float) / float(p["Pin"]),
            np.asarray(snap["phi"], dtype=float) / float(p["phi_soil"]),
            np.abs(np.asarray(snap["tau_b"], dtype=float)) / tau0,
        ]
        actual = _float_or_nan(snap.get("t_actual", snap.get("t", np.nan))) / max(_float_or_nan(snap.get("ter", np.nan)), 1e-30)
        axes[0, j].set_title(f"t/t_er={target:.1f}\nactual={actual:.3f}, {status}")
        for i, y in enumerate(ys):
            axes[i, j].plot(x, y, linewidth=1.8)
            axes[i, j].grid(True, alpha=0.3)
            if j == 0:
                axes[i, j].set_ylabel(labels[i])
    for j in range(3):
        axes[-1, j].set_xlabel("x/L")
    fig.suptitle(f"Figure 3 style — dimensionless profiles\nscenario={long_run.scenario}, Nx={long_run.nx}, CFL={long_run.cfl}", fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    return _save_figure_all_formats(fig, out_dir / "figure3_dimensionless_profiles")


def plot_figure4_mesh_convergence(matrix_runs: Sequence[RunRecord], out_dir: Path, *, preferred_cfl: float = 0.2) -> Dict[str, str]:
    chosen_cfl = preferred_cfl if preferred_cfl in {r.cfl for r in matrix_runs} else min(r.cfl for r in matrix_runs)
    runs = sorted([r for r in matrix_runs if abs(r.cfl - chosen_cfl) < 1e-12], key=lambda r: r.nx)
    if len(runs) < 2:
        raise RuntimeError("Not enough common-CFL runs to plot Figure 4.")
    fig, axes = plt.subplots(3, 2, figsize=(12.5, 9.5), sharex='col')
    row_labels = ["p/Pin", "R/R0", "phi/phi_soil"]
    for j, target in enumerate([2.2, 4.5]):
        for run in runs:
            snap, status = _nearest_snapshot(run.snapshots, target, require_reached=False)
            if snap is None:
                continue
            cfg = run.config
            x = np.asarray(snap["x"], dtype=float) / float(cfg["L"])
            ys = [
                np.asarray(snap["p"], dtype=float) / float(cfg["Pin"]),
                np.asarray(snap["R"], dtype=float) / float(cfg["R0"]),
                np.asarray(snap["phi"], dtype=float) / float(cfg["phi_soil"]),
            ]
            label = f"Nx={run.nx}" + (" (filled)" if status == "filled_missing_target" else "")
            for i, y in enumerate(ys):
                axes[i, j].plot(x, y, linewidth=1.6, label=label)
                axes[i, j].grid(True, alpha=0.3)
                if j == 0:
                    axes[i, j].set_ylabel(row_labels[i])
        axes[0, j].set_title(f"t/t_er={target:.1f}")
        axes[-1, j].set_xlabel("x/L")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    if handles:
        fig.legend(handles, labels, loc="upper center", ncol=max(2, len(labels)))
    fig.suptitle(f"Figure 4 style — mesh sensitivity at CFL={chosen_cfl}", fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    return _save_figure_all_formats(fig, out_dir / "figure4_mesh_sensitivity")


def postprocess_campaign(*, long_run: RunRecord, matrix_runs: Sequence[RunRecord], output_dir: Path) -> Dict[str, object]:
    output_dir.mkdir(parents=True, exist_ok=True)
    figures_dir = output_dir / "figures"
    tables_dir = output_dir / "tables"
    figures_dir.mkdir(parents=True, exist_ok=True)
    tables_dir.mkdir(parents=True, exist_ok=True)
    fig3 = plot_figure3_profiles(long_run, figures_dir)
    fig4 = plot_figure4_mesh_convergence(matrix_runs, figures_dir)
    t1 = build_table1(matrix_runs, tables_dir)
    t2a = build_table2(matrix_runs, tables_dir, target_over_ter=2.2)
    t2b = build_table2(matrix_runs, tables_dir, target_over_ter=4.5)
    manifest = {
        "scenario": long_run.scenario,
        "version_core": long_run.version_core,
        "long_run_dir": str(long_run.run_dir),
        "matrix_run_dirs": [str(r.run_dir) for r in matrix_runs],
        "figures": {"figure3_dimensionless_profiles": fig3, "figure4_mesh_sensitivity": fig4},
        "tables": {
            "table1": str(tables_dir / "table1_cfl_diagnostics.csv"),
            "table2a": str(tables_dir / "table2_2p2_relative_errors.csv"),
            "table2b": str(tables_dir / "table2_4p5_relative_errors.csv"),
        },
        "row_counts": {"table1": len(t1), "table2a": len(t2a), "table2b": len(t2b)},
    }
    with open(output_dir / "postprocess_manifest.json", "w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=2, ensure_ascii=False)
    return manifest


def postprocess_from_output_dir(outdir: Path, output_dir: Path, *, scenario: Optional[str] = None) -> Dict[str, object]:
    records = discover_run_records(outdir)
    long_run, matrix_runs = choose_latest_campaign(records, scenario=scenario)
    return postprocess_campaign(long_run=long_run, matrix_runs=matrix_runs, output_dir=output_dir)
