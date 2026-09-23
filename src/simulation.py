from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, Optional

import numpy as np

from .numerics import advection_conservative
from .physics import (
    barenblatt_fw,
    beta_barenblatt,
    conservative_solid_source,
    fm_julien,
    fm_julien_raw,
    mdot_erosion,
    radius_growth_rate,
    rho_mix,
    rho_soil_sat,
    shear_tau_b,
)


class PressureSolveError(RuntimeError):
    """Raised when the imposed-pressure discharge cannot be solved."""

    def __init__(self, message: str, info: Optional[Dict[str, Any]] = None):
        super().__init__(message)
        self.info = info or {}


@dataclass
class HydraulicState:
    """All hydraulic and erosion fields associated with one state."""

    Q: float
    residual: float
    tau_b: np.ndarray
    mdot: np.ndarray
    fw: np.ndarray
    beta: np.ndarray
    u: np.ndarray
    fm: np.ndarray
    fm_raw: np.ndarray
    solver_info: Dict[str, Any]



def characteristic_time(p) -> float:
    """Return the physical erosion time scale used for normalisation."""
    pressure_drop = p.Pin - p.Pout
    if pressure_drop <= 0.0:
        raise ValueError("Pin-Pout must be strictly positive for t_er")
    soil_density = rho_soil_sat(p.rho_w, p.rho_p, p.phi_soil)
    # t_er = 2 L rho_soil_sat / (k_er Delta_p)
    if p.k_er <= 0.0:
        raise ValueError("k_er must be strictly positive for t_er")
    return 2.0 * p.L * soil_density / (p.k_er * pressure_drop)


def _state_from_evaluation(
    Q: float,
    evaluation,
    solver_info: Dict[str, Any],
) -> HydraulicState:
    residual, tau_b, mdot, fw, beta, u, fm, fm_raw = evaluation
    return HydraulicState(
        Q=float(Q),
        residual=float(residual),
        tau_b=np.asarray(tau_b, dtype=float),
        mdot=np.asarray(mdot, dtype=float),
        fw=np.asarray(fw, dtype=float),
        beta=np.asarray(beta, dtype=float),
        u=np.asarray(u, dtype=float),
        fm=np.asarray(fm, dtype=float),
        fm_raw=np.asarray(fm_raw, dtype=float),
        solver_info=dict(solver_info),
    )


def solve_Q_pressure_imposed(
    R: np.ndarray,
    phi: np.ndarray,
    p,
    dx: float,
    Q_guess: float,
) -> HydraulicState:
    """Solve the quasi-steady discharge for the imposed pressure drop.

    For a trial discharge, the pressure loss is evaluated from the local
    Darcy-Weisbach wall stress and the outlet minor loss.  A bracketed
    bisection is used because the reference model requires a robust scalar
    solve.  A failed bracket is a hard numerical failure; silently continuing
    with an unsatisfied pressure residual would invalidate the physical run.
    """
    R = np.asarray(R, dtype=float)
    phi = np.asarray(phi, dtype=float)
    if R.ndim != 1 or phi.ndim != 1 or R.size != phi.size:
        raise ValueError("R and phi must be one-dimensional arrays of equal size")
    if dx <= 0.0:
        raise ValueError("dx must be strictly positive")

    rho = rho_mix(phi, p.rho_w, p.rho_p)
    fm_raw = fm_julien_raw(phi, p)
    fm = fm_julien(phi, p)

    def evaluate(Q: float):
        Q = float(Q)
        u = Q / (np.pi * R**2 + 1e-30)
        Re = 2.0 * rho * np.abs(u) * R / p.mu_w
        Re = np.maximum(Re, p.Re_min)
        fw_loc = barenblatt_fw(Re)
        beta_loc = beta_barenblatt(Re)
        tau_b = shear_tau_b(rho, fw_loc, fm, u)
        mdot = mdot_erosion(tau_b, p.tau_c, p.k_er)

        px = (2.0 / np.maximum(R, p.R_min)) * tau_b
        pL = (
            p.Pin
            + np.sum(px) * dx
            - p.K_out * 0.5 * rho[-1] * (u[-1] ** 2)
        )
        residual = float(pL - p.Pout)
        result = (residual, tau_b, mdot, fw_loc, beta_loc, u, fm, fm_raw)
        return result

    Q0 = max(float(Q_guess), 1e-12)
    Q_low = 0.0
    eval_low = evaluate(1e-16)
    r_low = float(eval_low[0])
    if not np.isfinite(r_low):
        info = {
            "bracket_found": False,
            "used_fallback": False,
            "n_expand": 0,
            "n_bisect": 0,
            "residual_abs": float("inf"),
            "Q_low": Q_low,
            "Q_high": Q0,
        }
        raise PressureSolveError("Non-finite pressure residual at Q=0", info)

    if abs(r_low) <= p.Q_tol_abs:
        info = {
            "bracket_found": True,
            "used_fallback": False,
            "n_expand": 0,
            "n_bisect": 0,
            "residual_abs": abs(r_low),
            "Q_low": Q_low,
            "Q_high": Q0,
        }
        return _state_from_evaluation(0.0, eval_low, info)

    Q_high = Q0
    eval_high = evaluate(Q_high)
    r_high = float(eval_high[0])
    expansions = 0
    bracket_found = np.sign(r_low) != np.sign(r_high)

    while not bracket_found and expansions < int(p.Q_bracket_max_expand):
        Q_high *= float(p.Q_bracket_growth)
        expansions += 1
        eval_high = evaluate(Q_high)
        r_high = float(eval_high[0])
        if not np.isfinite(r_high):
            break
        bracket_found = np.sign(r_low) != np.sign(r_high)

    if not bracket_found:
        info = {
            "bracket_found": False,
            "used_fallback": False,
            "n_expand": expansions,
            "n_bisect": 0,
            "residual_abs": float(abs(r_high)) if np.isfinite(r_high) else float("inf"),
            "Q_low": float(Q_low),
            "Q_high": float(Q_high),
        }
        raise PressureSolveError(
            "Could not bracket the imposed-pressure discharge", info
        )

    last_eval = eval_high
    for n_bisect in range(1, int(p.Q_bisect_max_iter) + 1):
        Q_mid = 0.5 * (Q_low + Q_high)
        eval_mid = evaluate(Q_mid)
        r_mid = float(eval_mid[0])
        if not np.isfinite(r_mid):
            info = {
                "bracket_found": True,
                "used_fallback": False,
                "n_expand": expansions,
                "n_bisect": n_bisect,
                "residual_abs": float("inf"),
                "Q_low": float(Q_low),
                "Q_high": float(Q_high),
            }
            raise PressureSolveError("Non-finite pressure residual during bisection", info)

        last_eval = eval_mid
        if abs(r_mid) <= p.Q_tol_abs:
            info = {
                "bracket_found": True,
                "used_fallback": False,
                "n_expand": expansions,
                "n_bisect": n_bisect,
                "residual_abs": abs(r_mid),
                "Q_low": float(Q_low),
                "Q_high": float(Q_high),
            }
            return _state_from_evaluation(Q_mid, eval_mid, info)

        if np.sign(r_mid) == np.sign(r_low):
            Q_low = Q_mid
            r_low = r_mid
        else:
            Q_high = Q_mid
            r_high = r_mid

    info = {
        "bracket_found": True,
        "used_fallback": False,
        "n_expand": expansions,
        "n_bisect": int(p.Q_bisect_max_iter),
        "residual_abs": abs(float(last_eval[0])),
        "Q_low": float(Q_low),
        "Q_high": float(Q_high),
    }
    raise PressureSolveError(
        "Pressure bisection reached its iteration limit", info
    )


def run_simulation(p, logger):
    """Run one coupled physical-time simulation."""
    p.validate()
    Nx = int(p.Nx)
    dx = p.L / Nx
    x = (np.arange(Nx) + 0.5) * dx

    ter = characteristic_time(p)
    t_end = p.t_end_factor * ter
    soil_density = rho_soil_sat(p.rho_w, p.rho_p, p.phi_soil)

    R = np.full(Nx, p.R0, dtype=float)
    A = np.pi * R**2
    # S is suspended solid volume per unit axial length: S=A*phi.
    S = np.zeros(Nx, dtype=float)

    hydro = solve_Q_pressure_imposed(R, np.zeros(Nx), p, dx, Q_guess=1e-4)
    timeseries = []
    snapshots = []

    snap_targets = [float(v) for v in p.snap_t_over_ter]
    snap_idx = 0
    target_times_reached = []

    effective_models = {
        "transport_variable": "S=A*phi",
        "phi_scheme_effective": str(p.phi_scheme),
        "phi_form_effective": "conservative_area",
        "phi_limiter_effective": str(p.phi_limiter),
        "clear_water_friction_model": "barenblatt_fw",
        "beta_model": "beta_barenblatt",
        "mixture_multiplier_model": "julien_tanh_capped",
        "pressure_solver": "bracketing_bisection",
        "erosion_rate_definition": "R_t=mdot/rho_soil_sat",
    }

    def current_phi():
        return np.clip(S / np.maximum(A, np.pi * p.R_min**2), 0.0, p.phi_soil)

    def take_snapshot(
        t: float,
        *,
        t_target_over_ter: Optional[float] = None,
        note: str = "",
        reached_target: bool = True,
    ):
        phi = current_phi()
        px = (2.0 / np.maximum(R, p.R_min)) * hydro.tau_b
        pprof = p.Pin + np.cumsum(px) * dx - (px * dx / 2.0)
        split = max(1, int(0.5 * Nx))
        ratio_C = (
            (np.max(np.abs(px[int(0.9 * Nx) :])) + 1e-12)
            / (np.max(np.abs(px[:split])) + 1e-12)
        )
        t_actual = float(t)
        t_plot = (
            t_actual
            if t_target_over_ter is None
            else float(t_target_over_ter) * float(ter)
        )
        return {
            "t": t_plot,
            "t_actual": t_actual,
            "t_target_over_ter": (
                None
                if t_target_over_ter is None
                else float(t_target_over_ter)
            ),
            "reached_target": bool(reached_target),
            "note": str(note),
            "ter": float(ter),
            "x": x.copy(),
            "R": R.copy(),
            "A": A.copy(),
            "S": S.copy(),
            "phi": phi.copy(),
            "Q": float(hydro.Q),
            "u": hydro.u.copy(),
            "tau_b": hydro.tau_b.copy(),
            "mdot": hydro.mdot.copy(),
            "erosion_mass_flux": hydro.mdot.copy(),
            "fw": hydro.fw.copy(),
            "beta": hydro.beta.copy(),
            "fm": hydro.fm.copy(),
            "fm_raw": hydro.fm_raw.copy(),
            "p": pprof.copy(),
            "px": px.copy(),
            "res_pL": float(hydro.residual),
            "ratio_C": float(ratio_C),
        }

    logger.info(f"t_er = {ter:.6f} s")
    logger.info(
        f"u0 = {np.mean(hydro.u):.6f} m/s ; Q0 = {hydro.Q:.6e} m3/s ; "
        f"tau_b0={np.mean(np.abs(hydro.tau_b)):.3f} Pa ; "
        f"tau_c={p.tau_c:.3f} Pa"
    )
    logger.info(
        f"Nx={Nx} dx={dx:.6e} L={p.L:.3f} R0={p.R0:.6e} "
        f"Pin={p.Pin:.3e} Pout={p.Pout:.3e}"
    )
    logger.info(
        f"transport={effective_models['transport_variable']} "
        f"scheme={effective_models['phi_scheme_effective']} "
        f"limiter={effective_models['phi_limiter_effective']}"
    )
    logger.info(
        f"friction={effective_models['mixture_multiplier_model']} "
        f"fm_max={p.fm_max:.6g} pressure_solver={effective_models['pressure_solver']}"
    )

    # Consume target t/ter=0 at the actual initial state.  This avoids the
    # former duplicate zero-time snapshot.
    while snap_idx < len(snap_targets) and snap_targets[snap_idx] <= 1e-14:
        target = snap_targets[snap_idx]
        snapshots.append(
            take_snapshot(
                0.0,
                t_target_over_ter=target,
                note="initial",
                reached_target=True,
            )
        )
        target_times_reached.append(target)
        snap_idx += 1

    t = 0.0
    dt = p.dt_max if p.dt_mode == "fixed" else min(1e-4, p.dt_max)
    stop_reason = "completed"
    stop_message = "Reached t_end."
    step = 0
    clip_dtmax_count = 0
    clip_dtmin_count = 0
    clip_dtend_count = 0
    adaptive_step_count = 0
    state_bound_violation_count = 0
    target_times_missing = []
    t_star_1D = None
    save_ts_every = int(p.save_ts_every)

    while t < t_end - 1e-14:
        if np.max(R) > p.R_break:
            stop_reason = "R_break"
            stop_message = (
                f"R exceeded R_break (maxR={np.max(R):.3f} m > "
                f"{p.R_break:.3f} m)."
            )
            logger.warning(stop_message)
            break

        dt_cfl = np.nan
        active_cfl = np.nan
        dt_clipped_to_max = False
        dt_clipped_to_min = False
        dt_clipped_to_end = False

        if p.dt_mode == "adaptive":
            adaptive_step_count += 1
            a_for_dt = hydro.beta * hydro.u
            amax = float(np.max(np.abs(a_for_dt)) + 1e-12)
            dt_adv = p.CFL * dx / amax
            k_arr = 2.0 * hydro.mdot / (
                np.maximum(R, p.R_min) * soil_density + 1e-30
            )
            k_max = float(np.max(k_arr))
            dt_src = np.inf if k_max <= 0.0 else 0.5 / k_max
            dt_cfl = float(min(dt_adv, dt_src))

            if dt_cfl < p.dt_min:
                stop_reason = "dt_min_violation"
                stop_message = (
                    f"Required dt={dt_cfl:.3e} s is below dt_min="
                    f"{p.dt_min:.3e} s."
                )
                logger.error(stop_message)
                break

            dt = float(min(dt_cfl, p.dt_max))
            dt_clipped_to_max = dt_cfl > p.dt_max
            clip_dtmax_count += int(dt_clipped_to_max)
            active_cfl = float(amax * dt / dx)
        else:
            dt = float(p.dt_max)

        remaining = t_end - t
        if dt > remaining:
            dt = remaining
            dt_clipped_to_end = True
            clip_dtend_count += 1

        if dt <= 0.0 or not np.isfinite(dt):
            stop_reason = "invalid_dt"
            stop_message = f"Invalid time step dt={dt}."
            logger.error(stop_message)
            break

        R_old = R.copy()
        A_old = A.copy()
        S_old = S.copy()

        # The validated scientific law is radial growth directly from the
        # erosion mass flux; no sqrt(1+R_x**2) geometric multiplier is used.
        Rt = radius_growth_rate(hydro.mdot, soil_density)
        R_trial = np.maximum(R_old + dt * Rt, p.R_min)
        A_trial = np.pi * R_trial**2

        # Conservative transport of S=A*phi.  The inlet solid-volume value is
        # based on the current inlet area and the prescribed phi_in.
        a = hydro.beta * hydro.u
        S_adv = advection_conservative(
            S_old,
            a,
            dx,
            dt,
            q_in=float(A_old[0] * p.phi_in),
            scheme=p.phi_scheme,
            limiter=p.phi_limiter,
        )

        # Eroded saturated soil opens an area increment carrying phi_soil of
        # solid volume.  This is the conservative counterpart of the former
        # concentration relaxation source.
        S_trial = S_adv + conservative_solid_source(
            A_old,
            A_trial,
            p.phi_soil,
        )
        upper_S = p.phi_soil * A_trial
        scale_S = max(float(np.max(upper_S)), 1.0)
        bound_tol = 1e-11 * scale_S
        if (
            not np.all(np.isfinite(S_trial))
            or float(np.min(S_trial)) < -bound_tol
            or float(np.max(S_trial - upper_S)) > bound_tol
        ):
            state_bound_violation_count += 1
            stop_reason = "state_bounds"
            stop_message = "Conservative solid-volume update violated physical bounds."
            logger.error(stop_message)
            break
        S_trial = np.clip(S_trial, 0.0, upper_S)
        phi_trial = S_trial / np.maximum(A_trial, np.pi * p.R_min**2)

        try:
            hydro_trial = solve_Q_pressure_imposed(
                R_trial,
                phi_trial,
                p,
                dx,
                Q_guess=hydro.Q,
            )
        except PressureSolveError as exc:
            stop_reason = "pressure_solver_failure"
            stop_message = str(exc)
            logger.error(
                f"{stop_message}; info={getattr(exc, 'info', {})}"
            )
            break

        R = R_trial
        A = A_trial
        S = S_trial
        hydro = hydro_trial
        t += dt
        step += 1

        if step % 50 == 0 or step == 1:
            prog_time = (t / t_end) * 100.0
            prog_R = (np.max(R) / p.R_break) * 100.0
            msg = (
                f"\rProg: {max(prog_time, prog_R):5.1f}% | "
                f"t/ter: {t/ter:.2f} | "
                f"R_max: {np.max(R)*1000:5.2f}mm | "
                f"Q: {hydro.Q:.2e} | dt: {dt:.1e}s | "
                f"Solv: {hydro.solver_info['n_bisect']:2d}"
            )
            print(msg, end="", flush=True)

        if step % p.log_every_steps == 0:
            phi = current_phi()
            logger.info(
                f"t={t:.3f} (t/ter={t / ter:.3f}) dt={dt:.2e} Q={hydro.Q:.3e} "
                f"Rout={1e3 * R[-1]:.2f}mm phi_out={phi[-1] / p.phi_soil:.3f} "
                f"max_phi={np.max(phi) / p.phi_soil:.3f} "
                f"max_fm={np.max(hydro.fm):.6g} "
                f"max_fm_raw={np.max(hydro.fm_raw):.3e} "
                f"max_mdot={np.max(hydro.mdot):.3e} "
                f"pL_res={hydro.residual:.3e} active_CFL={active_cfl}"
            )

        if step % save_ts_every == 0:
            phi = current_phi()
            slenderness = float(R[-1] / p.L)
            max_dRdx = 0.0 if Nx < 2 else float(
                np.max(np.abs(np.gradient(R, dx)))
            )
            if t_star_1D is None and slenderness > 0.1:
                t_star_1D = float(t)
                logger.warning(
                    f"t* (R_out/L > 0.1) reached: t*={t_star_1D:.4f} s "
                    f"| t*/ter={t_star_1D / ter:.4f}"
                )

            u_safe = np.maximum(np.abs(hydro.u), 1e-30)
            R_safe = np.maximum(R, p.R_min)
            epsilon_arr = 2.0 * hydro.mdot / (
                soil_density * R_safe * u_safe
            )
            timeseries.append(
                {
                    "t": float(t),
                    "t_over_ter": float(t / ter),
                    "dt": float(dt),
                    "dt_cfl": float(dt_cfl) if np.isfinite(dt_cfl) else np.nan,
                    "active_cfl": float(active_cfl) if np.isfinite(active_cfl) else np.nan,
                    "dt_clipped_to_max": bool(dt_clipped_to_max),
                    "dt_clipped_to_min": bool(dt_clipped_to_min),
                    "dt_clipped_to_end": bool(dt_clipped_to_end),
                    "Q": float(hydro.Q),
                    "R_out": float(R[-1]),
                    "phi_out": float(phi[-1]),
                    "res_pL": float(hydro.residual),
                    "max_R": float(np.max(R)),
                    "max_phi": float(np.max(phi)),
                    "max_mdot": float(np.max(hydro.mdot)),
                    "max_fw": float(np.max(hydro.fw)),
                    "max_beta": float(np.max(hydro.beta)),
                    "max_fm": float(np.max(hydro.fm)),
                    "mean_fm": float(np.mean(hydro.fm)),
                    "fm_out": float(hydro.fm[-1]),
                    "max_fm_raw": float(np.max(hydro.fm_raw)),
                    "mean_fm_raw": float(np.mean(hydro.fm_raw)),
                    "fm_raw_out": float(hydro.fm_raw[-1]),
                    "slenderness_ratio": slenderness,
                    "max_dR_dx": max_dRdx,
                    "epsilon_incomp_max": float(np.max(epsilon_arr)),
                    "solver_bracket_found": bool(hydro.solver_info["bracket_found"]),
                    "solver_used_fallback": bool(hydro.solver_info["used_fallback"]),
                    "solver_n_expand": int(hydro.solver_info["n_expand"]),
                    "solver_n_bisect": int(hydro.solver_info["n_bisect"]),
                    "solver_residual_abs": float(hydro.solver_info["residual_abs"]),
                }
            )

        while snap_idx < len(snap_targets) and (t / ter) >= snap_targets[snap_idx] - 1e-10:
            target = snap_targets[snap_idx]
            snapshots.append(
                take_snapshot(
                    t,
                    t_target_over_ter=target,
                    note="target",
                    reached_target=True,
                )
            )
            target_times_reached.append(target)
            snap_idx += 1

        if not (
            np.all(np.isfinite(R))
            and np.all(np.isfinite(A))
            and np.all(np.isfinite(S))
            and np.isfinite(hydro.Q)
        ):
            stop_reason = "nan_inf"
            stop_message = "NaN/Inf detected in the coupled state."
            logger.error(stop_message)
            break

        if np.max(R) > p.R_break:
            stop_reason = "R_break"
            stop_message = (
                f"R exceeded R_break after update (maxR={np.max(R):.3f} m)."
            )
            logger.warning(stop_message)
            break

    print()

    target_times_missing = [float(v) for v in snap_targets[snap_idx:]]
    for target in target_times_missing:
        snapshots.append(
            take_snapshot(
                t,
                t_target_over_ter=target,
                note="filled_missing_target",
                reached_target=False,
            )
        )
        logger.warning(
            f"SNAP (filled) t/ter={target:.3f} using last state at "
            f"t/ter={t / ter:.3f} (reason: simulation ended before target)."
        )

    active_cfl_values = [
        row["active_cfl"]
        for row in timeseries
        if np.isfinite(row.get("active_cfl", np.nan))
    ]
    run_manifest = {
        **effective_models,
        "stop_reason": stop_reason,
        "stop_message": stop_message,
        "n_steps": int(step),
        "t_final": float(t),
        "t_final_over_ter": float(t / ter) if ter > 0.0 else np.nan,
        "target_times_requested": [float(v) for v in snap_targets],
        "target_times_reached": [float(v) for v in target_times_reached],
        "target_times_not_reached": target_times_missing,
        "adaptive_step_count": int(adaptive_step_count),
        "clip_dtmax_fraction": (
            float(clip_dtmax_count / adaptive_step_count)
            if adaptive_step_count > 0
            else np.nan
        ),
        "clip_dtmin_fraction": (
            float(clip_dtmin_count / adaptive_step_count)
            if adaptive_step_count > 0
            else 0.0
        ),
        "clip_dtend_count": int(clip_dtend_count),
        "cfl_mean": float(np.mean(active_cfl_values)) if active_cfl_values else np.nan,
        "cfl_max": float(np.max(active_cfl_values)) if active_cfl_values else np.nan,
        "solver_fallback_count": 0,
        "last_pressure_residual_abs": float(abs(hydro.residual)),
        "state_bound_violation_count": int(state_bound_violation_count),
        "t_star_1D": t_star_1D if t_star_1D is not None else None,
        "t_star_1D_over_ter": (
            float(t_star_1D / ter) if t_star_1D is not None else None
        ),
    }

    logger.info(
        "RUN MANIFEST SUMMARY | "
        f"stop_reason={run_manifest['stop_reason']} "
        f"t_final_over_ter={run_manifest['t_final_over_ter']:.3f} "
        f"cfl_mean={run_manifest['cfl_mean']:.3f} "
        f"cfl_max={run_manifest['cfl_max']:.3f}"
    )

    return {
        "params": p,
        "x": x,
        "dx": dx,
        "t_er": ter,
        "t_end": t_end,
        "timeseries": timeseries,
        "snapshots": snapshots,
        "run_manifest": run_manifest,
    }
