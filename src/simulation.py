import numpy as np

from .numerics import advection_phi, ddx_centered
from .slfv import slfv_advection_phi
from .physics import (
    barenblatt_fw,
    fm_mixture,
    mdot_erosion,
    rho_mix,
    rho_soil_sat,
    shear_tau_b,
)


def characteristic_time(p):
    rho_s = rho_soil_sat(p.rho_w, p.rho_p, p.phi_soil)
    return 2.0 * p.L * rho_s / (p.k_er * p.Pin)


def _validate_parameters(p):
    if p.L <= 0 or p.Nx < 2:
        raise ValueError("L must be positive and Nx >= 2.")
    if p.R0 <= p.R_min:
        raise ValueError("R0 must be greater than R_min.")
    if p.Pin <= p.Pout:
        raise ValueError("Reference solver requires Pin > Pout.")
    if not (0.0 <= p.phi_init < p.phi_soil):
        raise ValueError("phi_init must satisfy 0 <= phi_init < phi_soil.")
    if not (0.0 <= p.phi_in < p.phi_soil):
        raise ValueError("phi_in must satisfy 0 <= phi_in < phi_soil.")
    if p.k_er <= 0 or p.tau_c < 0:
        raise ValueError("k_er must be positive and tau_c non-negative.")
    if p.CFL <= 0 or p.CFL > 1:
        raise ValueError("CFL must satisfy 0 < CFL <= 1.")
    if p.dt_min <= 0 or p.dt_max <= 0 or p.dt_min > p.dt_max:
        raise ValueError("Require 0 < dt_min <= dt_max.")


def solve_Q_pressure_imposed(R, phi, p, dx, Q_guess):
    """
    Solve F(Q)=p_calc(L;Q)-Pout=0 by bracketing and bisection.

    The thesis reference solver never accepts an unbracketed fallback.
    Failure to bracket or failure to reach the pressure tolerance raises
    RuntimeError so that an invalid hydraulic state cannot silently propagate.
    """
    R = np.asarray(R, dtype=float)
    phi = np.asarray(phi, dtype=float)
    rho = rho_mix(phi, p.rho_w, p.rho_p)
    fm = fm_mixture(phi, p)

    def residual(Q):
        u = Q / (np.pi * np.maximum(R, p.R_min) ** 2)
        Re = 2.0 * rho * np.abs(u) * R / p.mu_w
        Re = np.maximum(Re, p.Re_min)
        fw = barenblatt_fw(Re)
        tau_b = shear_tau_b(rho, fw, fm, u)
        mdot = mdot_erosion(tau_b, p.tau_c, p.k_er)

        px = (2.0 / np.maximum(R, p.R_min)) * tau_b
        pL = (
            p.Pin
            + np.sum(px) * dx
            - p.K_out * 0.5 * rho[-1] * u[-1] ** 2
        )
        return pL - p.Pout, tau_b, mdot, fw, u

    Q0 = max(float(Q_guess), 1e-12)
    Q_low = 0.0
    r_low, *_ = residual(1e-16)
    Q_high = Q0
    r_high, *_ = residual(Q_high)

    growth = float(p.Q_bracket_growth)
    expansions = 0
    bracket_found = np.sign(r_low) != np.sign(r_high)

    while not bracket_found and expansions < int(p.Q_bracket_max_expand):
        Q_high *= growth
        expansions += 1
        r_high, *_ = residual(Q_high)
        if not np.isfinite(r_high):
            break
        bracket_found = np.sign(r_low) != np.sign(r_high)

    if not bracket_found:
        raise RuntimeError(
            "Pressure solver failed to bracket Q: "
            f"Q_low={Q_low:.6e}, Q_high={Q_high:.6e}, "
            f"res_low={r_low:.6e}, res_high={r_high:.6e}."
        )

    tol = max(float(p.Q_tol_abs), float(p.Q_tol_rel) * max(abs(p.Pin - p.Pout), 1.0))

    for n_bisect in range(1, int(p.Q_bisect_max_iter) + 1):
        Qm = 0.5 * (Q_low + Q_high)
        rm, tau_b, mdot, fw, u = residual(Qm)

        if not np.isfinite(rm):
            raise RuntimeError("Non-finite pressure residual during bisection.")

        if abs(rm) <= tol:
            return Qm, rm, tau_b, mdot, fw, u, {
                "bracket_found": True,
                "n_expand": expansions,
                "n_bisect": n_bisect,
                "residual_abs": float(abs(rm)),
                "Q_low": float(Q_low),
                "Q_high": float(Q_high),
            }

        if np.sign(rm) == np.sign(r_low):
            Q_low, r_low = Qm, rm
        else:
            Q_high, r_high = Qm, rm

    raise RuntimeError(
        "Pressure solver reached the maximum number of bisection iterations "
        f"without meeting tolerance: residual={abs(rm):.6e} Pa, tol={tol:.6e} Pa."
    )


def run_simulation(p, logger):
    _validate_parameters(p)

    Nx = p.Nx
    dx = p.L / Nx
    x = (np.arange(Nx) + 0.5) * dx

    ter = characteristic_time(p)
    t_end = p.t_end_factor * ter

    R = np.full(Nx, p.R0, dtype=float)
    phi = np.full(Nx, p.phi_init, dtype=float)

    Q = 1e-12
    timeseries = []
    snapshots = []

    snap_targets = list(p.snap_t_over_ter)
    snap_idx = 0

    effective_models = {
        "phi_scheme_effective": str(p.phi_scheme),
        "phi_form_effective": "conservative",
        "clear_water_friction_model": "barenblatt_fw",
        "mixture_multiplier_model": "fm_mixture",
        "pressure_solver": "bracketing_bisection",
        "transport_velocity": "u=Q/A",
        "transport_method": str(p.transport_method),
        "time_step_policy": (
            "min(dt_src,dt_morph,dt_max)" if str(p.transport_method).lower() == "slfv"
            else "min(dt_adv,dt_src,dt_max)"
        ),
        "morph_rel_change": float(p.morph_rel_change),
        "soil_density_definition": "rho_s=rho_w+phi_soil*(rho_p-rho_w)",
    }

    def take_snapshot(
        t, Q, res_pL, tau_b, mdot, fw, u,
        *, t_target_over_ter=None, note="", reached_target=True
    ):
        fm = fm_mixture(phi, p)
        px = (2.0 / np.maximum(R, p.R_min)) * tau_b
        pprof = p.Pin + np.cumsum(px) * dx - 0.5 * px * dx

        ratio_C = (
            np.max(np.abs(px[int(0.9 * Nx):])) + 1e-12
        ) / (
            np.max(np.abs(px[:int(0.5 * Nx)])) + 1e-12
        )

        return {
            "t": float(t if t_target_over_ter is None else t_target_over_ter * ter),
            "t_actual": float(t),
            "t_target_over_ter": None if t_target_over_ter is None else float(t_target_over_ter),
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
            "fm": fm.copy(),
            "p": pprof.copy(),
            "px": px.copy(),
            "res_pL": float(res_pL),
            "ratio_C": float(ratio_C),
        }

    Q, res_pL, tau_b, mdot, fw, u, solver_info = solve_Q_pressure_imposed(
        R, phi, p, dx, Q_guess=1e-4
    )

    logger.info(f"t_er = {ter:.6f} s")
    logger.info(
        f"u0={np.mean(u):.6f} m/s; Q0={Q:.6e} m3/s; "
        f"|tau_b|_mean={np.mean(np.abs(tau_b)):.3f} Pa; tau_c={p.tau_c:.3f} Pa"
    )
    logger.info(
        f"Nx={Nx}; dx={dx:.6e} m; L={p.L:.6e} m; R0={p.R0:.6e} m; "
        f"Pin={p.Pin:.6e} Pa; Pout={p.Pout:.6e} Pa"
    )
    logger.info(
        f"transport={p.transport_method}:{p.phi_scheme}/{p.phi_limiter}/conservative; "
        f"friction=Barenblatt; pressure_solver=bracketing+bisection"
    )

    t = 0.0
    dt = p.dt_max if p.dt_mode == "fixed" else min(1e-4, p.dt_max)
    stop_reason = "completed"
    stop_message = "Reached t_end."

    snapshots.append(
        take_snapshot(t, Q, res_pL, tau_b, mdot, fw, u, note="initial")
    )

    step = 0
    adaptive_step_count = 0
    clip_dtmax_count = 0
    target_times_reached = []
    t_star_1D = None
    save_ts_every = int(p.save_ts_every)

    rho_s = rho_soil_sat(p.rho_w, p.rho_p, p.phi_soil)

    while t < t_end:
        if np.max(R) > p.R_break:
            stop_reason = "R_break"
            stop_message = (
                f"max(R)={np.max(R):.6e} m exceeded R_break={p.R_break:.6e} m."
            )
            logger.warning(stop_message)
            break

        dt_adv = np.inf
        dt_src = np.inf
        dt_candidate = np.inf
        active_cfl = np.nan
        dt_clipped_to_max = False
        dt_morph = np.inf

        if p.dt_mode == "adaptive":
            adaptive_step_count += 1

            amax = float(np.max(np.abs(u)))
            dt_adv = np.inf if amax <= 1e-30 else p.CFL * dx / amax

            k_arr = 2.0 * mdot / (
                np.maximum(R, p.R_min) * rho_s + 1e-30
            )
            k_max = float(np.max(k_arr))

            # The relaxation source is integrated analytically below.
            # Therefore k*dt <= 0.5 is an accuracy criterion, not a
            # stability restriction. It is retained for controlled
            # operator splitting and comparison with the v1 reference.
            dt_src = np.inf if k_max <= 1e-30 else 0.5 / k_max

            # Explicit Euler geometry update: constrain the fractional
            # radius change. This is the physical evolution that remains
            # explicitly time-integrated after the CFL restriction is lifted.
            Rt_est = (mdot / rho_s)
            morph_rate = np.max(Rt_est / np.maximum(R, p.R_min))
            dt_morph = (
                np.inf if morph_rate <= 1e-30
                else p.morph_rel_change / morph_rate
            )

            is_slfv = str(p.transport_method).lower() == "slfv"
            dt_candidate = min(
                dt_src, dt_morph, p.dt_max, t_end - t
            ) if is_slfv else min(
                dt_adv, dt_src, p.dt_max, t_end - t
            )

            limiting_values = (
                [dt_src, dt_morph]
                if is_slfv else [dt_adv, dt_src]
            )
            dt_clipped_to_max = (
                dt_candidate == p.dt_max
                and p.dt_max < min(limiting_values)
            )

            if dt_candidate < p.dt_min:
                raise RuntimeError(
                    "Required stable time step is below dt_min: "
                    f"dt_candidate={dt_candidate:.6e} s < dt_min={p.dt_min:.6e} s."
                )
            dt = float(dt_candidate)
            active_cfl = float(amax * dt / dx) if amax > 0 else 0.0
            clip_dtmax_count += int(dt_clipped_to_max)
        else:
            dt = min(float(dt), p.dt_max, t_end - t)

        # Explicit splitting: all rates below are evaluated from state n.
        Rx = ddx_centered(R, dx)
        Rt = (mdot / rho_s) * np.sqrt(1.0 + Rx ** 2)
        R_new = np.maximum(R + dt * Rt, p.R_min)

        if str(p.transport_method).lower() == "slfv":
            phi_adv = slfv_advection_phi(
                phi,
                u,
                dx,
                dt,
                phi_in=p.phi_in,
            )
        else:
            phi_adv = advection_phi(
                phi,
                u,
                dx,
                dt,
                scheme=p.phi_scheme,
                form="conservative",
                phi_in=p.phi_in,
                limiter=p.phi_limiter,
            )

        k = 2.0 * mdot / (
            np.maximum(R, p.R_min) * rho_s + 1e-30
        )
        phi_new = p.phi_soil - (
            p.phi_soil - phi_adv
        ) * np.exp(-k * dt)
        phi_new = np.clip(phi_new, 0.0, p.phi_soil * (1.0 - 1e-6))

        R = R_new
        phi = phi_new

        Q, res_pL, tau_b, mdot, fw, u, solver_info = solve_Q_pressure_imposed(
            R, phi, p, dx, Q_guess=Q
        )

        t += dt
        step += 1

        if step % p.log_every_steps == 0:
            fm = fm_mixture(phi, p)
            logger.info(
                f"t={t:.6e} s; t/ter={t/ter:.4f}; dt={dt:.3e} s; Q={Q:.3e}; "
                f"Rout={1e3*R[-1]:.4f} mm; phi_out/phi_soil={phi[-1]/p.phi_soil:.4f}; "
                f"max_fm={np.max(fm):.3e}; max_mdot={np.max(mdot):.3e}; "
                f"|pL residual|={abs(res_pL):.3e} Pa; active_CFL={active_cfl:.4f}"
            )

        if step % save_ts_every == 0:
            slenderness = float(R[-1] / p.L)
            max_dRdx = float(np.max(np.abs(ddx_centered(R, dx))))

            if t_star_1D is None and slenderness > 0.05:
                t_star_1D = float(t)
                logger.warning(
                    f"Nominal 1D validity threshold crossed: R_out/L={slenderness:.5f} > 0.05; "
                    f"t*={t_star_1D:.6e} s; t*/ter={t_star_1D/ter:.5f}."
                )

            u_safe = np.maximum(np.abs(u), 1e-30)
            R_safe = np.maximum(R, p.R_min)
            epsilon_arr = 2.0 * mdot / (rho_s * R_safe * u_safe)
            fm_arr = fm_mixture(phi, p)

            timeseries.append({
                "t": float(t),
                "t_over_ter": float(t / ter),
                "dt": float(dt),
                "dt_adv": float(dt_adv) if np.isfinite(dt_adv) else np.nan,
                "dt_src": float(dt_src) if np.isfinite(dt_src) else np.nan,
                "dt_candidate": float(dt_candidate) if np.isfinite(dt_candidate) else np.nan,
                "dt_morph": float(dt_morph) if np.isfinite(dt_morph) else np.nan,
                "active_cfl": float(active_cfl) if np.isfinite(active_cfl) else np.nan,
                "dt_clipped_to_max": bool(dt_clipped_to_max),
                "Q": float(Q),
                "R_out": float(R[-1]),
                "phi_out": float(phi[-1]),
                "res_pL": float(res_pL),
                "max_R": float(np.max(R)),
                "max_phi": float(np.max(phi)),
                "max_mdot": float(np.max(mdot)),
                "max_fw": float(np.max(fw)),
                "max_fm": float(np.max(fm_arr)),
                "mean_fm": float(np.mean(fm_arr)),
                "fm_out": float(fm_arr[-1]),
                "slenderness_ratio": slenderness,
                "max_dR_dx": max_dRdx,
                "epsilon_incomp_max": float(np.max(epsilon_arr)),
                "solver_n_expand": int(solver_info["n_expand"]),
                "solver_n_bisect": int(solver_info["n_bisect"]),
                "solver_residual_abs": float(solver_info["residual_abs"]),
            })

        while snap_idx < len(snap_targets) and (t / ter) >= snap_targets[snap_idx] - 1e-9:
            tgt = float(snap_targets[snap_idx])
            snapshots.append(
                take_snapshot(
                    t, Q, res_pL, tau_b, mdot, fw, u,
                    t_target_over_ter=tgt,
                    note="target",
                    reached_target=True,
                )
            )
            target_times_reached.append(tgt)
            snap_idx += 1

        if not (np.all(np.isfinite(R)) and np.all(np.isfinite(phi)) and np.isfinite(Q)):
            stop_reason = "nan_inf"
            stop_message = "NaN/Inf detected in state variables or discharge."
            logger.error(stop_message)
            break

    target_times_requested = [float(v) for v in snap_targets]
    target_times_missing = [float(v) for v in snap_targets[snap_idx:]]

    if target_times_missing:
        for tgt in target_times_missing:
            snapshots.append(
                take_snapshot(
                    t, Q, res_pL, tau_b, mdot, fw, u,
                    t_target_over_ter=tgt,
                    note="filled_missing_target",
                    reached_target=False,
                )
            )

    active_cfl_values = [
        row["active_cfl"] for row in timeseries
        if np.isfinite(row.get("active_cfl", np.nan))
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
            float(clip_dtmax_count / adaptive_step_count)
            if adaptive_step_count > 0 else np.nan
        ),
        "cfl_mean": float(np.mean(active_cfl_values)) if active_cfl_values else np.nan,
        "cfl_max": float(np.max(active_cfl_values)) if active_cfl_values else np.nan,
        "last_pressure_residual_abs": float(abs(res_pL)),
        "t_star_1D": t_star_1D if t_star_1D is not None else None,
        "t_star_1D_over_ter": (
            float(t_star_1D / ter) if t_star_1D is not None else None
        ),
    }

    logger.info(
        "RUN MANIFEST | "
        f"stop_reason={run_manifest['stop_reason']} | "
        f"t_final/ter={run_manifest['t_final_over_ter']:.5f} | "
        f"cfl_mean={run_manifest['cfl_mean']:.5f} | "
        f"cfl_max={run_manifest['cfl_max']:.5f}"
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
