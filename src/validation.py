"""
Qualitative signature checks retained for regression testing.

These checks are not the scientific verification/validation of the thesis.
They only detect whether a selected reference scenario preserves expected
qualitative signatures. Quantitative verification and validation belong to
Chapter 4.
"""
import math
from typing import Any, Dict, Optional, Tuple

import numpy as np


def _find_target_snapshot(snaps, target: float) -> Tuple[Optional[Dict[str, Any]], str]:
    for s in snaps:
        tgt = s.get("t_target_over_ter")
        if tgt is None:
            continue
        if abs(float(tgt) - float(target)) < 1e-9:
            return (s, "reached") if bool(s.get("reached_target", True)) else (s, "filled_missing_target")
    return None, "missing"


def validate_signatures(snaps, p, logger):
    """
    Qualitative regression checks only.

    The thresholds are retained from the development campaign and must not
    be presented as formal validation criteria in the thesis.
    """
    s22, status22 = _find_target_snapshot(snaps, 2.2)
    s157, status157 = _find_target_snapshot(snaps, 15.7)

    ratio_A = Rout_R0 = phi_max_22 = x_at_max = ratio_C = phi_end_max = math.nan
    A_pass = B1_pass = B2_pass = C_pass = False

    if status22 == "reached" and s22 is not None:
        x = s22["x"] / p.L
        R = s22["R"]
        phi = s22["phi"]
        px = s22["px"]

        up = np.mean(R[x > 0.9] - p.R0)
        mid = np.mean(R[x < 0.5] - p.R0) + 1e-12
        ratio_A = float(up / mid)
        Rout_R0 = float(R[-1] / p.R0)
        A_pass = ratio_A >= 5.0 and Rout_R0 >= 2.0

        phi22 = phi / p.phi_soil
        phi_max_22 = float(np.max(phi22))
        x_at_max = float(x[np.argmax(phi22)])
        B1_pass = phi_max_22 >= 0.9 and x_at_max >= 0.8

        Nx = len(px)
        ratio_C = float(
            (np.max(np.abs(px[int(0.9 * Nx):])) + 1e-12) /
            (np.max(np.abs(px[:int(0.5 * Nx)])) + 1e-12)
        )
        C_pass = ratio_C >= 1.2

    if status157 == "reached" and s157 is not None:
        phi157 = s157["phi"] / p.phi_soil
        phi_end_max = float(np.max(phi157))
        B2_pass = phi_end_max <= 0.1

    ok = bool(A_pass and B1_pass and B2_pass and C_pass)
    logger.info("=== QUALITATIVE REGRESSION SIGNATURES ===")
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
