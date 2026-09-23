import math
from typing import Any, Dict, Optional, Tuple

import numpy as np


def _find_target_snapshot(snaps, target: float) -> Tuple[Optional[Dict[str, Any]], str]:
    """
    Return the snapshot recorded for a requested target time.
    Status values:
      - "reached": target snapshot exists and was truly reached
      - "filled_missing_target": target snapshot exists but was synthetically filled
      - "missing": no snapshot tagged with this target exists
    """
    for s in snaps:
        tgt = s.get("t_target_over_ter")
        if tgt is None:
            continue
        if abs(float(tgt) - float(target)) < 1e-9:
            if bool(s.get("reached_target", True)):
                return s, "reached"
            return s, "filled_missing_target"
    return None, "missing"


def validate_signatures(snaps, p, logger):
    """
    Implements A/B/C checks (qualitative thresholds).
    IMPORTANT: quantitative metrics only use truly reached target snapshots.
    Synthetic "filled_missing_target" snapshots are never accepted for validation.
    """
    s22, status22 = _find_target_snapshot(snaps, 2.2)
    s157, status157 = _find_target_snapshot(snaps, 15.7)

    ratio_A = math.nan
    Rout_R0 = math.nan
    A_pass = False
    phi_max_22 = math.nan
    x_at_max = math.nan
    B1_pass = False
    ratio_C = math.nan
    C_pass = False
    phi_end_max = math.nan
    B2_pass = False

    if status22 == "reached" and s22 is not None:
        x = s22["x"] / p.L
        R = s22["R"]
        phi = s22["phi"]
        px = s22["px"]

        # --- A: trompette metric ---
        R0 = p.R0
        up = np.mean(R[x > 0.9] - R0)
        mid = np.mean(R[x < 0.5] - R0) + 1e-12
        ratio_A = float(up / mid)
        Rout_R0 = float(R[-1] / R0)
        A_pass = (ratio_A >= 5.0) and (Rout_R0 >= 2.0)

        # --- B1: mud wave at t/t_er=2.2 ---
        phi22 = phi / p.phi_soil
        phi_max_22 = float(np.max(phi22))
        x_at_max = float(x[np.argmax(phi22)])
        B1_pass = (phi_max_22 >= 0.9) and (x_at_max >= 0.8)

        # --- C: pressure kink proxy at t/t_er=2.2 ---
        Nx = len(px)
        rC = (np.max(np.abs(px[int(0.9 * Nx):])) + 1e-12) / (
            np.max(np.abs(px[:int(0.5 * Nx)])) + 1e-12
        )
        ratio_C = float(rC)
        C_pass = (ratio_C >= 1.2)

    if status157 == "reached" and s157 is not None:
        phi157 = s157["phi"] / p.phi_soil
        phi_end_max = float(np.max(phi157))
        B2_pass = (phi_end_max <= 0.1)

    ok = bool(A_pass and B1_pass and B2_pass and C_pass)

    logger.info("=== VALIDATION (signatures type Fig.2) ===")

    if status22 == "reached":
        logger.info(
            f"A trompette: ratio_A={ratio_A:.2f} (>= 5.0) ; "
            f"Rout/R0={Rout_R0:.2f} (>= 2.0) -> {'PASS' if A_pass else 'FAIL'}"
        )
        logger.info(
            f"B1 saturation aval: max(phi/phi_soil)={phi_max_22:.2f} (>= 0.9), "
            f"x@max={x_at_max:.2f} (>= 0.8) -> {'PASS' if B1_pass else 'FAIL'}"
        )
        logger.info(
            f"C cassure pression: ratio_C={ratio_C:.2f} (>= 1.2) -> {'PASS' if C_pass else 'FAIL'}"
        )
    else:
        logger.warning(
            "Target t/t_er=2.2 not truly reached -> A, B1, C marked as NOT EVALUABLE "
            f"(status={status22})."
        )

    if status157 == "reached":
        logger.info(
            f"B2 dilution fin: max(phi/phi_soil)@15.7={phi_end_max:.3f} (<= 0.1) -> "
            f"{'PASS' if B2_pass else 'FAIL'}"
        )
    else:
        logger.warning(
            "Target t/t_er=15.7 not truly reached -> B2 marked as NOT EVALUABLE "
            f"(status={status157})."
        )

    logger.info("========================================")

    return {
        "ok": ok,
        "target_2p2_status": status22,
        "target_15p7_status": status157,
        "ratio_A": ratio_A,
        "Rout_R0": Rout_R0,
        "A_pass": A_pass,
        "phi_max_2p2": phi_max_22,
        "x_at_phi_max_2p2": x_at_max,
        "B1_pass": B1_pass,
        "phi_end_max_15p7": phi_end_max,
        "B2_pass": B2_pass,
        "ratio_C": ratio_C,
        "C_pass": C_pass,
    }
