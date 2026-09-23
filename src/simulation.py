import numpy as np

from .numerics import advection_phi, ddx_centered
from .physics import (
    barenblatt_fw,
    beta_barenblatt,
    fm_julien,
    mdot_erosion,
    rho_mix,
    rho_soil_sat,
    shear_tau_b,
)


def characteristic_time(p):
    rho_s = rho_soil_sat(p.rho_w, p.rho_p, p.phi_soil)
    # t_er = 2 L rho_soil / (k_er p_in)
    return 2.0 * p.L * rho_s / (p.k_er * p.Pin)


def solve_Q_pressure_imposed(R, phi, p, dx, Q_guess):
    """
    Find Q such that p(L)=Pout given p(0)=Pin and p_x formula.
    Robust: bracketing + bisection.
    Also returns solver diagnostics so each run can state exactly what happened.
    """
    rho = rho_mix(phi, p.rho_w, p.rho_p)

    # Precompute fm depends on phi only
    fm = fm_julien(phi, p)

    def residual(Q):
        # u cell-centered
        u = Q / (np.pi * R**2 + 1e-30)

        Re = 2.0 * rho * np.abs(u) * R / p.mu_w
        Re = np.maximum(Re, p.Re_min)

        fw_loc = barenblatt_fw(Re)
        beta_loc = beta_barenblatt(Re)

        tau_b = shear_tau_b(rho, fw_loc, fm, u)
        mdot = mdot_erosion(tau_b, p.tau_c, p.k_er)

        # Dominant friction term for p_x (strict mode)
        px = (2.0 / np.maximum(R, p.R_min)) * tau_b
        pL = p.Pin + np.sum(px) * dx - p.K_out * 0.5 * rho[-1] * (u[-1] ** 2)
        return pL - p.Pout, tau_b, mdot, fw_loc, beta_loc, u

    # initial guess from previous step or a mild value
    Q0 = max(Q_guess, 1e-12)

    # build bracket
    Q_low = 0.0
    r_low, *_ = residual(Q_low + 1e-16)  # avoid exactly zero
    Q_high = Q0
    r_high, *_ = residual(Q_high)

    max_expand = int(getattr(p, "Q_bracket_max_expand", 60))
    growth = float(getattr(p, "Q_bracket_growth", 2.0))
    expansions = 0
    bracket_found = np.sign(r_low) != np.sign(r_high)

    # expand high until sign change or max expansions
    while (not bracket_found) and expansions < max_expand:
        Q_high *= growth
        expansions += 1
        r_high, *_ = residual(Q_high)
        if not np.isfinite(r_high):
            break
        bracket_found = np.sign(r_low) != np.sign(r_high)

    # fallback: if no bracket, return Q0 with diagnostics
    if not bracket_found:
        r, tau_b, mdot, fw_loc, beta_loc, u = residual(Q0)
        info = {
            "bracket_found": False,
            "used_fallback": True,
            "n_expand": expansions,
            "n_bisect": 0,
            "residual_abs": float(abs(r)),
            "Q_low": float(Q_low),
            "Q_high": float(Q_high),
        }
        return Q0, r, tau_b, mdot, fw_loc, beta_loc, u, info

    tol_abs = float(getattr(p, "Q_tol_abs", 1e-3))
    max_bisect = int(getattr(p, "Q_bisect_max_iter", 80))

    # bisection
    for n_bisect in range(1, max_bisect + 1):
        Qm = 0.5 * (Q_low + Q_high)
        rm, tau_b, mdot, fw_loc, beta_loc, u = residual(Qm)
        if abs(rm) < tol_abs:
            info = {
                "bracket_found": True,
                "used_fallback": False,
                "n_expand": expansions,
                "n_bisect": n_bisect,
                "residual_abs": float(abs(rm)),
                "Q_low": float(Q_low),
                "Q_high": float(Q_high),
            }
            return Qm, rm, tau_b, mdot, fw_loc, beta_loc, u, info
        if np.sign(rm) == np.sign(r_low):
            Q_low, r_low = Qm, rm
        else:
            Q_high, r_high = Qm, rm

    info = {
        "bracket_found": True,
        "used_fallback": False,
        "n_expand": expansions,
        "n_bisect": max_bisect,
        "residual_abs": float(abs(rm)),
        "Q_low": float(Q_low),
        "Q_high": float(Q_high),
    }
    return Qm, rm, tau_b, mdot, fw_loc, beta_loc, u, info


def run_simulation(p, logger):
    Nx = p.Nx
    dx = p.L / Nx
    x = (np.arange(Nx) + 0.5) * dx

    ter = characteristic_time(p)
    t_end = p.t_end_factor * ter

    # R0_init allows decoupling pipe geometry from Julien lm calibration.
    # If p.R0_init is set, use it for the geometric initial radius;
    # p.R0 is then used only for lm = cl * p.R0 inside fm_julien.
    _R0_geom = getattr(p, "R0_init", p.R0)
    R = np.full(Nx, _R0_geom, dtype=float)
    phi = np.zeros(Nx, dtype=float)

    Q = 1e-12
    timeseries = []
    snapshots = []

    snap_targets = list(p.snap_t_over_ter)
    snap_idx = 0

    effective_models = {
        "phi_scheme_effective": str(p.phi_scheme),
        "phi_form_effective": str(p.phi_form),
        "clear_water_friction_model": "barenblatt_fw",
        "beta_model": "beta_barenblatt",
        "mixture_multiplier_model": "fm_julien",
        "pressure_solver": "bracketing_bisection",
    }

    def take_snapshot(
        t,
        Q,
        res_pL,
        tau_b,
        mdot,
        fw,
        beta,
        u,
        *,
        t_target_over_ter=None,
        note="",
        reached_target=True,
    ):
        fm = fm_julien(phi, p)

        # CORRECTION : Profil de pression évalué proprement au centre des mailles
        px = (2.0 / np.maximum(R, p.R_min)) * tau_b
        pprof = p.Pin + np.cumsum(px) * dx - (px * dx / 2.0)

        ratio_C = (
            np.max(np.abs(px[int(0.9 * Nx):])) + 1e-12
        ) / (
            np.max(np.abs(px[:int(0.5 * Nx)])) + 1e-12
        )

        t_actual = float(t)
        t_plot = t_actual if t_target_over_ter is None else float(t_target_over_ter) * float(ter)

        return {
            "t": t_plot,
            "t_actual": t_actual,
            "t_target_over_ter": (None if t_target_over_ter is None else float(t_target_over_ter)),
            "reached_target": bool(reached_target),
            "note": str(note),
            "ter": float(ter),
            "x": x.copy(),
            "R": R.copy(),
            "phi": phi.copy(),
            "Q": float(Q),
            "u": u.copy(),
            "tau_b": tau_b.copy(),
            "mdot": mdot.copy(),
            "fw": fw.copy(),
            "beta": beta.copy(),
            "fm": fm.copy(),
            "p": pprof.copy(),
            "px": px.copy(),
            "res_pL": float(res_pL),
            "ratio_C": float(ratio_C),
        }

    Q, res_pL, tau_b, mdot, fw, beta, u, solver_info = solve_Q_pressure_imposed(
        R, phi, p, dx, Q_guess=1e-4
    )

    tau_b0 = float(np.mean(np.abs(tau_b)))
    logger.info(f"t_er = {ter:.6f} s (cible ~11.11 s)")
    logger.info(
        f"u0 = {np.mean(u):.6f} m/s ; Q0 = {Q:.6e} m3/s ; "
        f"tau_b0=P0 = {tau_b0:.3f} Pa ; tau_c={p.tau_c:.3f} Pa"
    )
    logger.info(
        f"Nx={Nx} dx={dx:.6e} L={p.L:.3f} R0={p.R0:.6e} Pin={p.Pin:.3e} Pout={p.Pout:.3e}"
    )
    logger.info(
        f"phi_scheme_effective={effective_models['phi_scheme_effective']} "
        f"phi_form_effective={effective_models['phi_form_effective']}"
    )
    logger.info(
        f"clear_water_friction_model={effective_models['clear_water_friction_model']} "
        f"beta_model={effective_models['beta_model']} "
        f"mixture_multiplier_model={effective_models['mixture_multiplier_model']}"
    )
    logger.info(
        f"dt_mode={p.dt_mode} dt_max={p.dt_max:.2e} CFL={p.CFL:.2f} K_out={p.K_out:.3e} "
        f"pressure_solver={effective_models['pressure_solver']}"
    )
    if solver_info["used_fallback"]:
        logger.warning(
            "Initial pressure solve used fallback without valid bracket. "
            f"|residual|={solver_info['residual_abs']:.3e} Pa"
        )

    t = 0.0
    dt = p.dt_max if p.dt_mode == "fixed" else min(1e-4, p.dt_max)
    stop_reason = "completed"
    stop_message = "Reached t_end."

    # snapshot at t=0
    snapshots.append(take_snapshot(t, Q, res_pL, tau_b, mdot, fw, beta, u, note="initial", reached_target=True))

    step = 0
    clip_dtmax_count = 0
    clip_dtmin_count = 0
    adaptive_step_count = 0
    solver_fallback_count = 1 if solver_info["used_fallback"] else 0
    target_times_reached = []
    # ---- 1D validity tracking ----
    t_star_1D = None  # premier t où R_out/L > 0.1

    save_ts_every = int(getattr(p, "save_ts_every", 1))  # default 1 (every step)

    while t < t_end:
        if np.max(R) > p.R_break:
            stop_reason = "R_break"
            stop_message = f"R exceeded R_break (maxR={np.max(R):.3f} m > {p.R_break:.3f} m)."
            logger.warning(f"Stop: {stop_message}")
            break

        dt_cfl = np.nan
        active_cfl = np.nan
        dt_clipped_to_max = False
        dt_clipped_to_min = False

        if p.dt_mode == "adaptive":
            adaptive_step_count += 1

            # --- Advective CFL ---
            a_for_dt = beta * u
            amax = float(np.max(np.abs(a_for_dt)) + 1e-12)
            dt_adv = p.CFL * dx / amax

            # --- STABILIZATION (Step 2) : Source-term stiffness limiter ---
            # The erosion relaxation rate k = 2*mdot/(R*rho_s) defines how fast
            # phi converges to phi_soil. Even with the exact exponential update,
            # the *coupling* between phi, fm, tau_b and mdot within a single
            # explicit step can overshoot when k*dt >> 1.
            # We bound dt <= 0.5 / k_max (half the shortest relaxation time).
            rho_s_step = rho_soil_sat(p.rho_w, p.rho_p, p.phi_soil)
            k_arr = 2.0 * mdot / (np.maximum(R, p.R_min) * rho_s_step + 1e-30)
            k_max = float(np.max(k_arr) + 1e-30)
            dt_src = 0.5 / k_max

            dt_cfl = min(dt_adv, dt_src)   # combined constraint
            dt = float(np.clip(dt_cfl, p.dt_min, p.dt_max))
            dt_clipped_to_max = dt_cfl > p.dt_max
            dt_clipped_to_min = dt_cfl < p.dt_min
            active_cfl = float(amax * dt / dx)
            clip_dtmax_count += int(dt_clipped_to_max)
            clip_dtmin_count += int(dt_clipped_to_min)

        # --- update R (erosion) ---
        Rx = ddx_centered(R, dx)
        rho_s = rho_soil_sat(p.rho_w, p.rho_p, p.phi_soil)
        Rt = (mdot / rho_s) * np.sqrt(1.0 + Rx**2)
        R = np.maximum(R + dt * Rt, p.R_min)

        # --- advect phi ---
        a = beta * u
        _limiter = getattr(p, "phi_limiter", "minmod")
        phi_adv = advection_phi(
            phi,
            a,
            dx,
            dt,
            scheme=p.phi_scheme,
            form=p.phi_form,
            limiter=_limiter,
        )

        # --- exact relaxation source for phi ---
        k = 2.0 * mdot / (np.maximum(R, p.R_min) * rho_s)
        phi = p.phi_soil - (p.phi_soil - phi_adv) * np.exp(-k * dt)
        phi = np.clip(phi, 0.0, p.phi_soil * (1.0 - 1e-6))

        # --- solve hydraulics for next step ---
        Q, res_pL, tau_b, mdot, fw, beta, u, solver_info = solve_Q_pressure_imposed(
            R, phi, p, dx, Q_guess=Q
        )
        solver_fallback_count += int(solver_info["used_fallback"])

        t += dt
        step += 1

        # ---- Dashboard Console (Temps réel) ----
        if step % 50 == 0 or step == 1:
            prog_time = (t / t_end) * 100
            prog_R = (np.max(R) / p.R_break) * 100
            prog_max = max(prog_time, prog_R)
            
            dt_status = "WARN" if dt < p.dt_min * 10 else "OK"
            solv_status = "WARN" if solver_info["n_bisect"] > 50 else "OK"
            
            msg = (
                f"Prog: {prog_max:5.1f}% | "
                f"t/ter: {t/ter:.2f} | "
                f"R_max: {np.max(R)*1000:5.2f}mm (Stop: {p.R_break*1000:.0f}mm) | "
                f"Q: {Q:.2e} | "
                f"dt: {dt:.1e}s {dt_status} | "
                f"Solv: {solver_info['n_bisect']:2d} {solv_status}"
            )
            print(f"\r{msg}", end="", flush=True)

        # ---- Standard Logging (Kept to write in run.log) ----
        if step % p.log_every_steps == 0:
            fm = fm_julien(phi, p)
            logger.info(
                f"t={t:.3f} (t/ter={t / ter:.3f}) dt={dt:.2e} Q={Q:.3e} "
                f"Rout={1e3 * R[-1]:.2f}mm phi_out={phi[-1] / p.phi_soil:.3f} "
                f"max_phi={np.max(phi) / p.phi_soil:.3f} max_fm={np.max(fm):.2e} "
                f"max_mdot={np.max(mdot):.2e} pL_res={res_pL:.2e} active_CFL={active_cfl:.3f}"
            )

        # ---- Time series saving (subsampled) ----
        if step % save_ts_every == 0:
            # ---- 1D validity diagnostics ----
            slenderness = float(R[-1] / p.L)       # R_out / L
            max_dRdx    = float(np.max(np.abs(ddx_centered(R, dx))))  # max |dR/dx|

            # ---- t* detection ----
            if t_star_1D is None and slenderness > 0.1:
                t_star_1D = float(t)
                logger.warning(
                    f"t* (R_out/L > 0.1 — limite validité 1D) atteint : "
                    f"t*={t_star_1D:.4f} s | t*/ter={t_star_1D / ter:.4f} | "
                    f"R_out/L={slenderness:.4f}"
                )

            # ---- volumetric incompressibility ratio ----
            u_safe   = np.maximum(np.abs(u), 1e-30)
            R_safe   = np.maximum(R, p.R_min)
            rho_s_val = rho_soil_sat(p.rho_w, p.rho_p, p.phi_soil)
            epsilon_arr = 2.0 * mdot / (rho_s_val * R_safe * u_safe)
            epsilon_max = float(np.max(epsilon_arr))

            # ---- fm recording ----
            fm_arr  = fm_julien(phi, p)

            timeseries.append(
                {
                    "t": float(t),
                    "t_over_ter": float(t / ter),
                    "dt": float(dt),
                    "dt_cfl": float(dt_cfl) if np.isfinite(dt_cfl) else np.nan,
                    "active_cfl": float(active_cfl) if np.isfinite(active_cfl) else np.nan,
                    "dt_clipped_to_max": bool(dt_clipped_to_max),
                    "dt_clipped_to_min": bool(dt_clipped_to_min),
                    "Q": float(Q),
                    "R_out": float(R[-1]),
                    "phi_out": float(phi[-1]),
                    "res_pL": float(res_pL),
                    "max_R": float(np.max(R)),
                    "max_phi": float(np.max(phi)),
                    "max_mdot": float(np.max(mdot)),
                    "max_fw": float(np.max(fw)),
                    "max_beta": float(np.max(beta)),
                    # fm recording
                    "max_fm": float(np.max(fm_arr)),
                    "mean_fm": float(np.mean(fm_arr)),
                    "fm_out": float(fm_arr[-1]),
                    # 1D validity diagnostics
                    "slenderness_ratio": slenderness,
                    "max_dR_dx": max_dRdx,
                    # incompressibility
                    "epsilon_incomp_max": epsilon_max,
                    "solver_bracket_found": bool(solver_info["bracket_found"]),
                    "solver_used_fallback": bool(solver_info["used_fallback"]),
                    "solver_n_expand": int(solver_info["n_expand"]),
                    "solver_n_bisect": int(solver_info["n_bisect"]),
                    "solver_residual_abs": float(solver_info["residual_abs"]),
                }
            )

        # ---- Snapshots at target times ----
        while snap_idx < len(snap_targets) and (t / ter) >= snap_targets[snap_idx] - 1e-6:
            tgt = float(snap_targets[snap_idx])
            snapshots.append(
                take_snapshot(
                    t,
                    Q,
                    res_pL,
                    tau_b,
                    mdot,
                    fw,
                    beta,
                    u,
                    t_target_over_ter=tgt,
                    note="target",
                    reached_target=True,
                )
            )
            target_times_reached.append(tgt)
            snap_idx += 1

        # ---- Safety stop on NaN/Inf ----
        if not (np.all(np.isfinite(R)) and np.all(np.isfinite(phi)) and np.isfinite(Q)):
            stop_reason = "nan_inf"
            stop_message = "NaN/Inf detected in state variables or discharge."
            logger.error("NaN/Inf detected -> abort.")
            break

    print()  # Skips a line in the console at the end of the time loop to avoid overwriting the display

    target_times_requested = [float(v) for v in snap_targets]
    target_times_missing = [float(v) for v in snap_targets[snap_idx:]]

    if target_times_missing:
        for tgt in target_times_missing:
            snapshots.append(
                take_snapshot(
                    t,
                    Q,
                    res_pL,
                    tau_b,
                    mdot,
                    fw,
                    beta,
                    u,
                    t_target_over_ter=tgt,
                    note="filled_missing_target",
                    reached_target=False,
                )
            )
            logger.warning(
                f"SNAP (filled) t/ter={tgt:.3f} using last state at t/ter={t / ter:.3f} "
                f"(reason: simulation ended before reaching target)."
            )

    active_cfl_values = [
        row["active_cfl"] for row in timeseries if np.isfinite(row.get("active_cfl", np.nan))
    ]

    run_manifest = {
        **effective_models,
        "stop_reason": stop_reason,
        "stop_message": stop_message,
        "n_steps": int(step),
        "t_final": float(t),
        "t_final_over_ter": float(t / ter) if ter > 0 else np.nan,
        "target_times_requested": target_times_requested,
        "target_times_reached": [float(v) for v in target_times_reached],
        "target_times_not_reached": target_times_missing,
        "adaptive_step_count": int(adaptive_step_count),
        "clip_dtmax_fraction": (
            float(clip_dtmax_count / adaptive_step_count) if adaptive_step_count > 0 else np.nan
        ),
        "clip_dtmin_fraction": (
            float(clip_dtmin_count / adaptive_step_count) if adaptive_step_count > 0 else np.nan
        ),
        "cfl_mean": float(np.mean(active_cfl_values)) if active_cfl_values else np.nan,
        "cfl_max": float(np.max(active_cfl_values)) if active_cfl_values else np.nan,
        "solver_fallback_count": int(solver_fallback_count),
        "last_pressure_residual_abs": float(abs(res_pL)),
        # 1D validity
        "t_star_1D": t_star_1D if t_star_1D is not None else None,
        "t_star_1D_over_ter": (float(t_star_1D / ter) if t_star_1D is not None else None),
    }

    logger.info(
        "RUN MANIFEST SUMMARY | "
        f"stop_reason={run_manifest['stop_reason']} "
        f"t_final_over_ter={run_manifest['t_final_over_ter']:.3f} "
        f"cfl_mean={run_manifest['cfl_mean']:.3f} "
        f"cfl_max={run_manifest['cfl_max']:.3f} "
        f"clip_dtmax_fraction={run_manifest['clip_dtmax_fraction']:.3f}"
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
