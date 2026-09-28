"""
run_revision.py
===============
Revision runner: regenerate Tables 3, 4 and 5 of the manuscript with the
retained parameter set.

- tau_b coefficient: 1/8 (already in shear_tau_b)
- dp = 75 um, fm_max = 5.0 (defaults of config.py), K_in = 0 (no inlet loss)
- R/L validity check (already in the simulation loop)

ALIGNEMENT AVEC LE CHAPITRE 4 (a lire avant de modifier quoi que ce soit)
-------------------------------------------------------------------------
Ce script utilisait k_er = 0.01 (Sol A) et 0.001 (Sol B) et phi_soil = 0.62,
alors que le jeu de parametres retenu (chapitre 4, tableau 4.5, et article) est

    Sol A : k_er = 1.0e-4 s/m, tau_c = 1.0 Pa, phi_soil = 0.60
    Sol B : k_er = 1.0e-5 s/m, tau_c = 10.0 Pa, phi_soil = 0.60

avec 2 L rho_s / (k_er * Delta p) donnant t_er = 949.4 s (Sol A) et 9493.6 s
(Sol B), soit les valeurs de l'article.

Les valeurs du script n'etaient pas fausses pour autant : elles correspondent au
jeu retenu multiplie par un facteur d'acceleration temporelle

    k_er_solver = ALPHA * k_er_retenu,   ALPHA = 100 par defaut

Le probleme etait qu'aucun document ne le disait, et que les temps affiches
n'etaient donc pas les temps dimensionnels de l'article.

POURQUOI L'ACCELERATION EST VALIDE
----------------------------------
Les equations sont auto-similaires en temps sous k_er -> ALPHA * k_er : la
geometrie R(x,t) verifie R(x, t ; ALPHA k_er) = R(x, ALPHA t ; k_er), car (i)
l'hydraulique est quasi-stationnaire (aucune derivee temporelle), (ii) la mise a
jour du rayon est lineaire en k_er, et (iii) la source de concentration est une
relaxation exacte, phi = phi_s - (phi_s - phi_adv) exp(-k dt), dont l'argument
k*dt est proportionnel a k_er. La seule quantite qui ne se rescale pas est le
temps de transit advectif L/(beta u) = 0,026 s, negligeable devant t_er meme
dans le calcul accelere (9,5 s). Le champ de concentration est donc en regime
quasi-stationnaire et le retour rheologique reste inactif (phi/phi_s < 4e-4),
ce que le chapitre etablit par ailleurs pour la campagne A.

Verification numerique (R_out/R0 a t/t_er = 0,05, Nx = 50, reproduite par
--verifier-echelle) :

    ALPHA         1        10        100       1000
    R_out/R0  1.050459  1.050476  1.050616  1.051407
    ecart         --    +0,002 %  +0,015 %  +0,090 %
    duree      385,7 s   41,8 s     4,1 s     0,4 s

ALPHA = 100 est retenu. Valeur de controle reproduite (Nx = 100, 1,5 t_er) :
t_2 = 665,3 s, soit 0,701 t_er, contre 668,6 s (0,704 t_er) a Nx = 600 dans le
chapitre, et 666 s dans l'article : ecart 0,5 %, dans la dispersion de maillage
documentee (<= 1,1 % entre Nx = 100 et Nx = 600, section 4.10 du chapitre). Tous les temps restitues (t_er, t2, t_end) sont
multiplies par ALPHA : les valeurs imprimées et ecrites dans summary.json sont
donc directement comparables a l'article. Utiliser --verifier-echelle pour
reproduire la table ci-dessus.

COUT
----
Le cout est domine par Nx et par t_end_factor (dt ~ CFL dx / u, donc ~2x plus de
pas a Nx = 100 qu'a Nx = 50). Valeurs mesurees :
    Nx = 50,  0,2 t_er  ->  13 s par cas                       (--rapide)
    Nx = 100, 1,5 t_er  ->  ~5 min par cas Sol A
    Nx = 100, 5,0 t_er  ->  ~15 min par cas Sol A, ~2,5 h pour Sol B
    Nx = 50,  0,05 t_er ->  de 386 s (ALPHA = 1) a 4,1 s (ALPHA = 100)

Usage
-----
    python run_revision.py                     # jeu retenu, ALPHA = 100
    python run_revision.py --case SolA_Kout0   # un seul cas (les noms sont ceux
                                               # de SCENARIOS, ci-dessous)
    python run_revision.py --rapide            # verification de bon fonctionnement
    python run_revision.py --nx 200            # maillage plus fin (plus lent)
    python run_revision.py --verifier-echelle  # revalide l'acceleration temporelle
    python run_revision.py --alpha 1           # calcul non accelere (reference)
"""

import argparse
import json
import os
import sys
import time as time_module

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from src.config import Params                          # noqa: E402
from src.simulation import run_simulation, characteristic_time   # noqa: E402
from src.utils import (make_run_dir, setup_logger, save_config,   # noqa: E402
                       save_timeseries_csv, save_snapshots_npz, save_json)


# =====================================================================
# 1. JEU DE PARAMETRES RETENU (chapitre 4, tableau 4.5 ; article)
# =====================================================================
K_ER_RETENU = {"SolA": 1.0e-4, "SolB": 1.0e-5}   # s/m
TAU_C = {"SolA": 1.0, "SolB": 10.0}               # Pa
PHI_SOIL_RETENU = 0.60                            # etait 0.62 dans ce script

# Facteur d'acceleration temporelle (voir l'en-tete). 1.0 = calcul non accelere.
ALPHA_DEFAUT = 100.0


def rho_soil_sat(rho_w, rho_p, phi_soil):
    return rho_w * (1.0 - phi_soil) + rho_p * phi_soil


def t_er_retenu(p, k_er_retenu):
    """t_er = 2 L rho_s / (k_er Delta p), avec le k_er RETENU (pas l'accelere)."""
    return 2.0 * p.L * rho_soil_sat(p.rho_w, p.rho_p, p.phi_soil) / (k_er_retenu * p.Pin)


# =====================================================================
# 2. COMMUN AUX CAS  (identique au jeu retenu, k_er mis a part)
# =====================================================================
COMMON = dict(
    L=0.117, R0=0.003, Nx=100,
    rho_w=1000.0, rho_p=2650.0,
    mu_w=1e-3, phi_soil=PHI_SOIL_RETENU, phi_in=0.0,
    dp=75e-6, fm_max=5.0, cB=0.2, cl=0.07,
    julien_lambda_power=2.0, lm_mode="clR0", fm_cap_mode="hard",
    K_in=0.0, inlet_kinetic=False,
    Re_min=1.0, R_min=1e-6, R_break=0.5,
    phi_scheme="muscl", phi_form="conservative", phi_limiter="vanleer",
    dt_mode="adaptive", dt_min=1e-12, dt_max=10.0, CFL=0.35,
    t_end_factor=5.0,                      # run to 5 * t_er (accelere)
    Q_tol_abs=1e-3, Q_bisect_max_iter=80,
    Q_bracket_max_expand=60, Q_bracket_growth=2.0,
    snap_t_over_ter=[0.0, 2.0, 3.0, 4.0, 5.0],
    log_every_steps=100,
    save_ts_every=10,                      # 1 saturerait la memoire sur les runs longs
)

SCENARIOS = {
    # Table 3: Soil A (erodible), K_out=0
    "SolA_Kout0": dict(
        Pin=4905.0, Pout=0.0, K_out=0.0, sol="SolA",
        description="Soil A (erodible silt-clay), K_out=0, Delta_p=0.5m"
    ),
    # Table 3: Soil B (resistant), K_out=0
    "SolB_Kout0": dict(
        Pin=4905.0, Pout=0.0, K_out=0.0, sol="SolB",
        description="Soil B (resistant silt-clay), K_out=0, Delta_p=0.5m"
    ),
    # Table 4: Soil A, K_out=10
    "SolA_Kout10": dict(
        Pin=4905.0, Pout=0.0, K_out=10.0, sol="SolA",
        description="Soil A (erodible), K_out=10, Delta_p=0.5m"
    ),
    # Table 5: Soil A, K_out=0, Delta_p=0.25m
    "SolA_Kout0_dp0.25": dict(
        Pin=2452.5, Pout=0.0, K_out=0.0, sol="SolA",
        description="Soil A, K_out=0, Delta_p=0.25m"
    ),
    # Table 5: Soil A, K_out=0, Delta_p=1.0m
    "SolA_Kout0_dp1.00": dict(
        Pin=9810.0, Pout=0.0, K_out=0.0, sol="SolA",
        description="Soil A, K_out=0, Delta_p=1.0m"
    ),
}


# =====================================================================
# 3. VERIFICATION DE L'ACCELERATION TEMPORELLE
# =====================================================================
def verifier_echelle(alphas=(1.0, 10.0, 100.0, 1000.0), t_cible=0.05, Nx=50):
    """Revalide l'auto-similarite : R_out/R0 a t/t_er doit etre independant de ALPHA."""
    import io
    import contextlib
    print("\n" + "=" * 78)
    print("VERIFICATION DE L'ACCELERATION TEMPORELLE  (Sol A, t/t_er = %.2f, Nx = %d)"
          % (t_cible, Nx))
    print("=" * 78)
    print("  Si k_er -> ALPHA k_er, alors R_out/R0 a t/t_er doit etre identique.")
    print("  Le cas ALPHA = 1 est la reference non acceleree.\n")
    print(f"  {'ALPHA':>8s} {'k_er_solver':>13s} {'t_er_acc (s)':>13s} "
          f"{'R_out/R0':>11s} {'ecart':>9s} {'duree':>8s}")
    ref = None
    for a in alphas:
        cfg = dict(COMMON, Nx=Nx, k_er=K_ER_RETENU["SolA"] * a, tau_c=TAU_C["SolA"],
                   Pin=4905.0, Pout=0.0, K_out=0.0, t_end_factor=t_cible,
                   snap_t_over_ter=[0.0], save_ts_every=1)
        p = Params(**{k: v for k, v in cfg.items() if k in Params.__dataclass_fields__})
        for k, v in cfg.items():
            if not hasattr(p, k):
                setattr(p, k, v)
        logger = setup_logger("output_revision/verif_echelle_alpha%g.log" % a)
        t0 = time_module.time()
        with contextlib.redirect_stdout(io.StringIO()):
            res = run_simulation(p, logger)
        el = time_module.time() - t0
        val = res["timeseries"][-1]["R_out"] / p.R0
        if ref is None:
            ref = val
        print(f"  {a:8.0f} {cfg['k_er']:13.1e} {characteristic_time(p):13.3f} "
              f"{val:11.6f} {100 * (val / ref - 1):+8.3f}% {el:7.1f}s")
    print("\n  -> ALPHA = %.0f est retenu : ecart mesure +0,015 %% a t/t_er = 0,05."
          % ALPHA_DEFAUT)


# =====================================================================
# 4. EXECUTION
# =====================================================================
def main():
    ap = argparse.ArgumentParser(
        description="Regenere les tableaux 3, 4 et 5 du manuscrit avec le jeu de "
                    "parametres retenu (chapitre 4).")
    ap.add_argument("--nx", type=int, default=COMMON["Nx"],
                    help=f"nombre de cellules (defaut {COMMON['Nx']}). Le chapitre montre que "
                         "les grandeurs integrales (t2, R_out/R0) sont insensibles au maillage "
                         "a mieux de 1,1 %% entre Nx = 100 et Nx = 600.")
    ap.add_argument("--alpha", type=float, default=ALPHA_DEFAUT,
                    help=f"facteur d'acceleration temporelle (defaut {ALPHA_DEFAUT:.0f}). "
                         "Les temps affiches sont multiplies par ce facteur.")
    ap.add_argument("--facteur-fin", type=float, default=COMMON["t_end_factor"],
                    help=f"horizon en multiples de t_er (defaut {COMMON['t_end_factor']}).")
    ap.add_argument("--case", type=str, default=None,
                    help="n'executer qu'un seul cas (nom de SCENARIOS) : "
                         + ", ".join(sorted(SCENARIOS)))
    ap.add_argument("--rapide", action="store_true",
                    help="test de bon fonctionnement : Nx = 50, 0,2 t_er.")
    ap.add_argument("--verifier-echelle", action="store_true",
                    help="revalide l'acceleration temporelle puis quitte.")
    args = ap.parse_args()

    if args.verifier_echelle:
        os.makedirs("output_revision", exist_ok=True)
        verifier_echelle()
        return

    nx = 50 if args.rapide else args.nx
    facteur_fin = 0.2 if args.rapide else args.facteur_fin
    alpha = args.alpha

    print("=" * 80)
    print("RUN_REVISION -- jeu de parametres retenu (chapitre 4 / article)")
    print("=" * 80)
    print(f"  Sol A : k_er = {K_ER_RETENU['SolA']:.1e} s/m, tau_c = {TAU_C['SolA']} Pa")
    print(f"  Sol B : k_er = {K_ER_RETENU['SolB']:.1e} s/m, tau_c = {TAU_C['SolB']} Pa")
    print(f"  phi_soil = {PHI_SOIL_RETENU} ; L = {COMMON['L']} m ; R0 = {COMMON['R0']} m")
    print(f"  acceleration temporelle ALPHA = {alpha:.0f} "
          f"(temps affiches = temps calcules x {alpha:.0f})")
    print(f"  Nx = {nx} ; horizon = {facteur_fin} t_er ; sauvegarde tous les "
          f"{COMMON['save_ts_every']} pas")

    noms = sorted(SCENARIOS)
    if args.case:
        if args.case not in SCENARIOS:
            raise SystemExit(f"cas inconnu : {args.case} (choix : {', '.join(noms)})")
        noms = [args.case]

    results = {}
    for name in noms:
        overrides = dict(SCENARIOS[name])
        print(f"\n{'=' * 60}")
        print(f"Running: {name}")
        print(f"  {overrides['description']}")
        print(f"{'=' * 60}")

        sol = overrides.pop("sol")
        cfg = dict(COMMON, **overrides)
        cfg["Nx"] = nx
        cfg["t_end_factor"] = facteur_fin
        cfg["phi_soil"] = PHI_SOIL_RETENU
        cfg["k_er"] = K_ER_RETENU[sol] * alpha      # <- seul le solveur voit l'acceleration
        cfg["tau_c"] = TAU_C[sol]

        p = Params(**{k: v for k, v in cfg.items() if k in Params.__dataclass_fields__})
        for k, v in cfg.items():
            if not hasattr(p, k):
                setattr(p, k, v)

        # t_er et t2 dimensionnels : calcules avec le k_er RETENU, par mise a l'echelle
        ter_vrai = t_er_retenu(p, K_ER_RETENU[sol])
        ter_acc = characteristic_time(p)

        run_dir = make_run_dir(base_dir="output_revision",
                               version="v2_retenu", scenario=name)
        logger = setup_logger(run_dir / "run.log")
        save_config(run_dir / "config.json", p)

        t0 = time_module.time()
        try:
            res = run_simulation(p, logger)
            elapsed = time_module.time() - t0

            info = res["run_manifest"]
            ts = res["timeseries"]
            snaps = res["snapshots"]
            t_end_acc = res["t_end"]

            R_out_final = float(ts[-1]["R_out"]) if len(ts) > 0 else p.R0
            R_out_R0 = R_out_final / p.R0

            # temps de doublement : premier passage de R_out a 2 R0, puis rescale
            t2_acc = None
            for row in ts:
                if row.get("R_out", 0) >= 2.0 * p.R0:
                    t2_acc = row["t"]
                    break

            # --- temps DIMENSIONNELS (comparables a l'article) -----------
            t2 = t2_acc * alpha if t2_acc is not None else None
            t_end_vrai = t_end_acc * alpha

            results[name] = {
                "sol": sol,
                "k_er_retenu": K_ER_RETENU[sol],
                "k_er_solver": cfg["k_er"],
                "alpha": alpha,
                "ter_analytic_s": ter_vrai,
                "ter_accelerated_s": ter_acc,
                "t_end_s": t_end_vrai,
                "t_end_over_ter": t_end_vrai / ter_vrai,
                "t2_s": t2,
                "t2_over_ter": (t2 / ter_vrai) if t2 else None,
                "R_out_final": R_out_final,
                "R_out_R0": R_out_R0,
                "n_steps": info.get("n_steps", 0),
                "stop_reason": info.get("stop_reason", "unknown"),
                "R_out_L": R_out_final / p.L,
                "t_star_1D_s": (info.get("t_star_1D", None) * alpha
                                if info.get("t_star_1D", None) is not None else None),
                "elapsed_s": elapsed,
            }

            save_timeseries_csv(run_dir / "timeseries.csv", ts)
            save_snapshots_npz(run_dir / "snapshots.npz", snaps)
            save_json(run_dir / "run_manifest.json", {
                "name": name,
                "description": overrides["description"],
                "results": results[name],
                "params": cfg,
                "note_temps": ("Tous les temps du bloc 'results' sont dimensionnels : ils "
                               "integrent le facteur d'acceleration ALPHA. Le solveur a "
                               "utilise k_er_solver = ALPHA * k_er_retenu."),
            })

            print(f"  Done in {elapsed:.1f}s")
            print(f"  t_er (retenu) = {ter_vrai:.1f} s  |  solveur accelere : "
                  f"t_er = {ter_acc:.2f} s")
            print(f"  t_end = {t_end_vrai:.1f} s ({t_end_vrai / ter_vrai:.3f} t_er)")
            print(f"  R_out_final = {R_out_final * 1000:.3f} mm, R_out/R0 = {R_out_R0:.2f}")
            print(f"  t2 = {t2:.1f} s ({t2 / ter_vrai:.3f} t_er)" if t2 else "  t2: NOT REACHED")
            print(f"  R_out/L = {R_out_final / p.L:.3f}, t*_1D = "
                  f"{results[name]['t_star_1D_s']}")

        except Exception as e:
            elapsed = time_module.time() - t0
            print(f"  FAILED after {elapsed:.1f}s: {e}")
            import traceback
            traceback.print_exc()
            results[name] = {"error": str(e)}

    # ---------------- resume -------------------------------------------
    print("\n" + "=" * 92)
    print("RESUME -- jeu de parametres retenu, temps dimensionnels (secondes)")
    print("=" * 92)
    print(f"{'Case':<22} {'t_er(s)':>10} {'t2(s)':>10} {'t2/t_er':>9} "
          f"{'R_out/R0':>10} {'R_out/L':>9} {'Stop':>12}")
    print("-" * 92)
    for name, r in results.items():
        if "error" in r:
            print(f"{name:<22} {'ERROR: ' + r['error'][:50]}")
        else:
            t2s = f"{r['t2_s']:.1f}" if r.get("t2_s") else "N/A"
            t2r = f"{r['t2_over_ter']:.3f}" if r.get("t2_over_ter") else "N/A"
            print(f"{name:<22} {r['ter_analytic_s']:>10.1f} {t2s:>10} {t2r:>9} "
                  f"{r['R_out_R0']:>10.2f} {r['R_out_L']:>9.3f} {r['stop_reason']:>12}")

    os.makedirs("output_revision", exist_ok=True)
    with open("output_revision/summary.json", "w", encoding="utf-8") as f:
        json.dump(results, f, indent=2, ensure_ascii=False)
    print("\nSummary saved to output_revision/summary.json")
    print("Rappel : t_er = 2 L rho_s / (k_er_retenu Delta p) ; la colonne t_er est donc")
    print("analytique, les autres proviennent du calcul accelere puis rescale.")


if __name__ == "__main__":
    main()
