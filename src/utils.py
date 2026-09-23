import csv
import json
import logging
from datetime import datetime
from pathlib import Path
from typing import Any

import numpy as np

from .config import params_from_dict


def load_scenario(name: str):
    data_path = Path("data") / "scenarios.json"
    with open(data_path, "r", encoding="utf-8") as f:
        all_sc = json.load(f)
    if name not in all_sc:
        raise KeyError(f"Scenario '{name}' not found in {data_path}")
    return params_from_dict(all_sc[name])


def make_run_dir(base_dir: str = "output", version: str = "v0", scenario: str = "scenario") -> Path:
    """
    Creates output/<YYYYMMDD_HHMMSS>_<version>_<scenario>/
    Even if 'output' folder does not exist in the zip, it will be created here.
    """
    base = Path(base_dir)
    base.mkdir(parents=True, exist_ok=True)

    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    run = base / f"{stamp}_{version}_{scenario}"
    run.mkdir(parents=True, exist_ok=True)
    return run


def setup_logger(log_path: Path):
    logger = logging.getLogger("Piping1D")
    logger.setLevel(logging.INFO)
    logger.handlers.clear()

    fmt = logging.Formatter("%(asctime)s | %(levelname)s | %(message)s")

    fh = logging.FileHandler(log_path, encoding="utf-8")
    fh.setFormatter(fmt)
    fh.setLevel(logging.INFO)

    sh = logging.StreamHandler()
    sh.setFormatter(fmt)
    sh.setLevel(logging.INFO)

    logger.addHandler(fh)
    logger.addHandler(sh)
    return logger


def save_config(path: Path, p) -> None:
    with open(path, "w", encoding="utf-8") as f:
        json.dump(p.to_dict(), f, indent=2)


def save_json(path: Path, data: Any) -> None:
    with open(path, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, ensure_ascii=False)


def save_timeseries_csv(path: Path, rows) -> None:
    if len(rows) == 0:
        return
    keys = list(rows[0].keys())
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        for r in rows:
            w.writerow(r)


def _to_npz_value(v: Any) -> np.ndarray:
    if isinstance(v, np.ndarray):
        return v
    if isinstance(v, (np.floating, float)):
        return np.array([float(v)], dtype=float)
    if isinstance(v, (np.integer, int)):
        return np.array([int(v)], dtype=int)
    if isinstance(v, (np.bool_, bool)):
        return np.array([bool(v)], dtype=bool)
    if v is None:
        return np.array([np.nan], dtype=float)
    if isinstance(v, str):
        return np.array([v])
    return np.asarray(v)


def save_snapshots_npz(path: Path, snaps) -> None:
    """
    Store snapshots as a single NPZ with arrays grouped by snapshot index.
    Supports numeric arrays plus scalar strings/bools used by audit metadata.
    """
    out = {}
    for i, s in enumerate(snaps):
        prefix = f"s{i:02d}_"
        for k, v in s.items():
            out[prefix + k] = _to_npz_value(v)
    np.savez_compressed(path, **out)
