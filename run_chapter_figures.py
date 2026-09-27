"""
run_chapter_figures.py
======================
Genere les QUATRE figures attendues par le chapitre 4
(`docs/chapitre4_cadre_modelisation_numerique.tex`), avec les noms de
fichiers exactement tels que le chapitre les appelle :

    docs/figure1.png   -> \\includegraphics{...}{figure1.png}  (fig:conduit-1D)
    docs/figure2.png   -> \\includegraphics{...}{figure2.png}  (fig:boucle)
    docs/figure3.png   -> \\includegraphics{...}{figure3.png}  (fig:profils-spatiaux)
    docs/figure4.png   -> \\includegraphics{...}{figure4.png}  (fig:convergence-maillage)

Le chapitre compile depuis `docs/`, les fichiers sont donc ecrits a cote du
`.tex`. Les `.png` sont volontairement ignores par le depot (`.gitignore`) :
ce sont des artefacts de build, a regenérer avant `pdflatex`.

Correspondance avec les legendes du chapitre
--------------------------------------------
figure1  schema du conduit effectif 1D : pressions imposees a l'entree et a la
         sortie, vitesse moyenne u, concentration phi, contrainte parietal
         tau_b et elargissement induit par l'erosion.
figure2  structure de retroaction du modele couple : hydraulique -> frottement
         du melange -> erosion -> evolution geometrique -> transport, plus les
         deux retours (rheologique et geometrique).
figure3  4 panneaux empiles (R/R0, p/Pin, phi/phi_soil, f_m) le long du
         conduit, pour la configuration a perte aval de la campagne B, a
         t/ter = 0 / 0.5 / 1.0 (front raide, monotone et convergé en maillage ;
         cf. "reserve numerique" ci-dessous). La ligne tiretee du premier
         panneau materialise la limite de validite geometrique R/L = 0.05.
figure4  convergence spatiale : champs couples pour Nx = 150 / 300 / 600 au
         meme instant t/ter = 1.0, plus un panneau quantitatif donnant l'ecart
         de concentration au maillage le plus fin (echelle logarithmique) ;
         l'ecart est localise au front advectif et decroit avec le
         raffinement, tandis que pression et rayon sont deja superposes.

Usage
-----
    cd <racine du depot>
    python run_chapter_figures.py                 # les 4 figures
    python run_chapter_figures.py --only 3 4      # simulations seules
    python run_chapter_figures.py --Nx 300        # figure 3 moins couteuse

Duree indicative (machine de bureau, sans parallelisation) :
    figures 1 et 2 :   < 1 s
    figure 3       :   ~1-2 min   (Nx = 600, jusqu'a t = 5.3 t_er)
    figure 4       :   ~1 min     (Nx = 150 / 300 / 600, jusqu'a t = 2.4 t_er)
"""

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

from src.config import Params                      # noqa: E402
from src.simulation import run_simulation          # noqa: E402
from src.utils import setup_logger                 # noqa: E402
import logging                                     # noqa: E402

OUT_DIR = SCRIPT_DIR / "docs"                      # le chapitre compile depuis docs/
LOG_DIR = SCRIPT_DIR / "Figures_HD"
LOG_DIR.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# Charte graphique (aspect "these")
# ---------------------------------------------------------------------------
plt.rcParams.update({
    "font.size": 11,
    "axes.labelsize": 12,
    "axes.titlesize": 11,
    "legend.fontsize": 9.5,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "lines.linewidth": 1.8,
    "axes.grid": True,
    "grid.alpha": 0.25,
    "grid.linestyle": "--",
    "mathtext.fontset": "cm",
    "figure.facecolor": "white",
    "savefig.facecolor": "white",
})

COLORS = ["#1f4e79", "#c0504d", "#4f8a3d", "#7030a0", "#bf8f00"]


# ---------------------------------------------------------------------------
# Configuration physique des campagnes (identique a data/scenarios.json)
# ---------------------------------------------------------------------------
# Campagne B (terrain) : utilisee par les figures 3 et 4, car c'est la seule
# ou le front de concentration, la saturation de f_m et la localisation aval
# du gradient de pression sont effectivement visibles (campagne A :
# phi/phi_soil ~ 1e-4 et f_m ~ 1.000, profils quasi cylindriques).
CAMPAGNE_B = dict(
    L=100.0, R0=0.003, Pin=3330000.0, Pout=0.0,
    rho_w=1000.0, rho_p=2700.0, phi_soil=0.50, phi_in=0.0, mu_w=1.0e-3,
    tau_c=10.0, k_er=1.0e-2,
    dp=1.0e-3, cl=0.10, cB=0.01, julien_lambda_power=2.0,
    fm_max=2000.0, lm_mode="clR0", Re_min=3000.0, R_min=1.0e-6, R_break=0.5,
    phi_scheme="muscl", phi_form="conservative", phi_limiter="vanleer",
    dt_mode="adaptive", dt_min=1.0e-6, dt_max=10.0, CFL=0.35,
    Q_tol_abs=1.0e-3, Q_bisect_max_iter=60, Q_bracket_max_expand=40,
    K_in=0.0, inlet_kinetic=False, log_every_steps=10 ** 9, save_ts_every=1,
)

# Configuration a perte aval : celle dont le chapitre discute l'effet
# (section 4.10.4, "une resistance aval ne peut que retarder l'instabilite").
# K_out = 10 (perte moderee) ou 500 (perte forte, retenue par defaut).
# RESERVE NUMERIQUE (a lire avant de changer ces instants)
#
# Le plafond fm <= fm_max est applique dans le solveur, conformement a
# l'equation publiee ; c'est ce qui permet de retrouver les temps d'echec
# publies a mieux de 0.5 %. En contrepartie, une fois phi sature a phi_soil
# sur une portion significative du conduit, le champ fm sature a son tour et
# le front de concentration developpe des oscillations dont l'amplitude CROIT
# avec le maillage. Mesure de la variation totale de phi/phi_soil, qui vaut 1
# pour un saut raide monotone :
#
#   t/ter   = 1.6 :  Nx=300 -> 0.927   Nx=600 -> 0.952   (converge, ecart 2.7 %)
#   t/ter   = 1.7 :  Nx=300 -> 1.173   Nx=600 -> 1.424   (divergent, ecart 21 %)
#   t/ter   = 2.2 :  Nx=150 -> 2.59    Nx=300 -> 4.38    Nx=600 -> 5.17
#            (sans plafond, le meme indicateur vaut ~1.2 quel que soit le maillage)
#
# L'artefact est donc TARDIF (il apparaît entre t/ter = 1.6 et 1.7, quand phi
# sature) et LOCAL (il ne compromet pas les grandeurs integrales : le temps
# d'echec reste stable a mieux de 1.1 % pour Nx = 100 / 300 / 600, et les
# rapports de retard par K_out = 500 sont reproduits a mieux de 0.2 %).
#
# Les instants retenus ci-dessous sont donc places AVANT le seuil, la ou le
# front est raide, monotone et converge en maillage. Voir la section
# "Contrepartie numerique du plafonnement" du chapitre 4 (tableau 4.12).
FIG3_KOUT = 500.0
FIG3_SNAPS = [0.0, 0.5, 1.0]
FIG3_NX = 600

FIG4_KOUT = 500.0
FIG4_T_OVER_TER = 1.0
FIG4_NX = [150, 300, 600]

# Critere de validite geometrique du chapitre : R/L < 0.05
VALIDITE_RL = 0.05


# ---------------------------------------------------------------------------
# Utilitaires
# ---------------------------------------------------------------------------
def make_params(**kwargs):
    """Construit un Params en tolerant les cles non declarees du dataclass."""
    known = set(Params.__dataclass_fields__)
    p = Params(**{k: v for k, v in kwargs.items() if k in known})
    for k, v in kwargs.items():
        if k not in known:
            setattr(p, k, v)
    return p


def run_case(tag, t_end_factor, snap_t_over_ter, **overrides):
    """Lance une simulation et renvoie (params, resultats)."""
    cfg = dict(CAMPAGNE_B, **overrides)
    cfg["snap_t_over_ter"] = list(snap_t_over_ter)
    cfg["t_end_factor"] = float(t_end_factor)
    p = make_params(**cfg)
    logger = setup_logger(LOG_DIR / f"run_chapter_fig_{tag}.log")
    # les simulations tracent leur progression sur stdout ; on ne garde que le
    # resume final du script, le detail restant disponible dans le .log
    logger.propagate = False
    for _h in list(logger.handlers):
        if isinstance(_h, logging.StreamHandler) and not isinstance(_h, logging.FileHandler):
            logger.removeHandler(_h)
    t0 = time.time()
    res = run_simulation(p, logger)
    elapsed = time.time() - t0
    ts = res["timeseries"]
    n_steps = res["run_manifest"]["n_steps"]
    stop = res["run_manifest"]["stop_reason"]
    if ts:
        t_fin, t_over = ts[-1]["t"], ts[-1]["t_over_ter"]
    else:                       # serie non echantillonnee : lire les snapshots
        s_last = res["snapshots"][-1]
        t_fin = float(s_last.get("t_actual", s_last["t"]))
        t_over = t_fin / res["t_er"]
    print(f"    [{tag}] t_fin = {t_fin:.2f} s ({t_over:.2f} t_er), "
          f"{n_steps} pas, arret = {stop}, {elapsed:.1f} s")
    return p, res


def pick_snapshot(snaps, target_t):
    """Instant le plus proche d'une cible (les cibles sont t/t_er)."""
    best, best_d = snaps[0], 1e30
    for s in snaps:
        d = abs(float(s.get("t_actual", s["t"])) - target_t)
        if d < best_d:
            best, best_d = s, d
    return best


def save(fig, name):
    fig.savefig(OUT_DIR / name, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  [OK] docs/{name}")


# ===========================================================================
# Figure 1 -- schema du conduit effectif unidimensionnel
# ===========================================================================
def figure1_schema_conduit(out="figure1.png"):
    fig, ax = plt.subplots(figsize=(9.0, 3.4))
    ax.set_xlim(-0.12, 1.14)
    ax.set_ylim(-0.66, 0.80)
    ax.axis("off")

    # --- profil de rayon : cylindrique puis "trompette" en aval -----------
    x = np.linspace(0.0, 1.0, 400)
    R = 0.13 + 0.20 * (1.0 - np.exp(-3.2 * x))
    ax.plot(x, R, color="black", lw=1.6, zorder=3)
    ax.plot(x, -R, color="black", lw=1.6, zorder=3)
    ax.fill_between(x, -R, R, color="#cfe2f3", alpha=0.85, zorder=1)

    # axe du conduit
    ax.plot([0.0, 1.0], [0.0, 0.0], ls=(0, (6, 5)), color="0.35", lw=1.0, zorder=2)

    # --- conditions aux limites : pressions imposees ----------------------
    ax.add_patch(FancyArrowPatch((-0.095, 0.0), (-0.015, 0.0), arrowstyle="-|>",
                                 mutation_scale=16, lw=2.0, color=COLORS[1], zorder=4))
    ax.text(-0.105, 0.15, r"$P_{\mathrm{in}}$", color=COLORS[1],
            ha="center", va="bottom", fontsize=13)
    ax.add_patch(FancyArrowPatch((1.095, 0.0), (1.015, 0.0), arrowstyle="-|>",
                                 mutation_scale=16, lw=2.0, color=COLORS[1], zorder=4))
    ax.text(1.105, 0.15, r"$P_{\mathrm{out}}$", color=COLORS[1],
            ha="center", va="bottom", fontsize=13)

    # --- ecoulement : u et phi -------------------------------------------
    for xa in (0.16, 0.42, 0.68):
        ax.add_patch(FancyArrowPatch((xa, 0.0), (xa + 0.16, 0.0), arrowstyle="-|>",
                                     mutation_scale=14, lw=1.6, color=COLORS[0], zorder=4))
    ax.text(0.46, 0.045, r"$u(x)$,  $\phi(x)$", color=COLORS[0],
            ha="center", va="bottom", fontsize=12)

    # --- contrainte parietale (sur la paroi superieure) -------------------
    for xa in (0.20, 0.50, 0.80):
        yl = 0.13 + 0.20 * (1.0 - np.exp(-3.2 * xa))
        ax.add_patch(FancyArrowPatch((xa + 0.07, yl + 0.02), (xa - 0.07, yl + 0.02),
                                     arrowstyle="-|>", mutation_scale=13,
                                     lw=1.8, color=COLORS[2], zorder=4))
    ax.text(0.50, 0.50, r"$\tau_b$ (contrainte pariétale)", color=COLORS[2],
            ha="center", va="bottom", fontsize=11.5)

    # --- cotation des rayons ---------------------------------------------
    # R0 : cote placee sur la frontiere amont ; le libelle est decale vers
    # l'interieur du conduit pour ne pas chevaucher la fleche de pression.
    ax.annotate("", xy=(0.0, 0.13), xytext=(0.0, -0.13),
                arrowprops=dict(arrowstyle="<|-|>", lw=1.2, color="black"))
    ax.text(0.012, 0.135, r"$R_0$", ha="left", va="bottom", fontsize=12)

    ax.annotate("", xy=(0.90, R[-1]), xytext=(0.90, -R[-1]),
                arrowprops=dict(arrowstyle="<|-|>", lw=1.2, color="black"))
    ax.text(0.878, -R[-1] - 0.03, r"$R(x,t)$", ha="right", va="top", fontsize=12)

    # --- elargissement induit --------------------------------------------
    ax.annotate("élargissement induit\npar l'érosion",
                xy=(0.97, 0.315), xytext=(0.70, 0.72),
                ha="center", va="center", fontsize=10.5, color=COLORS[3],
                arrowprops=dict(arrowstyle="-|>", lw=1.3, color=COLORS[3],
                                connectionstyle="arc3,rad=-0.25"))

    # --- axe x ------------------------------------------------------------
    ax.annotate("", xy=(1.0, -0.56), xytext=(0.0, -0.56),
                arrowprops=dict(arrowstyle="-|>", lw=1.2, color="black"))
    ax.text(0.0, -0.605, r"$x=0$", ha="center", va="top", fontsize=11)
    ax.text(1.0, -0.605, r"$x=L$", ha="center", va="top", fontsize=11)

    fig.tight_layout()
    save(fig, out)


# ===========================================================================
# Figure 2 -- structure de retroaction du modele couple
# ===========================================================================
def figure2_boucle(out="figure2.png"):
    fig, ax = plt.subplots(figsize=(8.4, 6.8))
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(-0.02, 1.10)
    ax.axis("off")

    nodes = [
        (0.20, 0.78, "Hydraulique\n" r"$Q,\ u(x)$", COLORS[0]),
        (0.63, 0.88, "Frottement du mélange\n" r"$\tau_b=\frac{1}{8}f_D\,f_m\rho u^2$", COLORS[2]),
        (0.88, 0.52, "Érosion\n" r"$\dot m=k_{er}\left(|\tau_b|-\tau_c\right)^{+}$", COLORS[1]),
        (0.63, 0.14, "Géométrie du conduit\n" r"$\partial_t R=\dot m/\rho_s$", COLORS[3]),
        (0.18, 0.24, "Transport en suspension\n" r"$\partial_t\phi+\partial_x(\beta u\phi)=S$", COLORS[4]),
    ]
    W, H = 0.30, 0.14

    boxes = []
    for (cx, cy, label, color) in nodes:
        box = FancyBboxPatch((cx - W / 2, cy - H / 2), W, H,
                             boxstyle="round,pad=0.012,rounding_size=0.02",
                             linewidth=1.5, edgecolor=color, facecolor=color,
                             alpha=0.13, zorder=2)
        ax.add_patch(box)
        ax.text(cx, cy, label, ha="center", va="center", fontsize=9.6,
                color=color, zorder=3)
        boxes.append(box)

    def connect(i, j, label, *, color="0.25", ls="-", lw=1.7, rad=0.0,
                label_off=(0.0, 0.0)):
        a = FancyArrowPatch(nodes[i][:2], nodes[j][:2], arrowstyle="-|>",
                            mutation_scale=18, lw=lw, color=color, linestyle=ls,
                            connectionstyle=f"arc3,rad={rad}", zorder=1)
        a.set_patchA(boxes[i])
        a.set_patchB(boxes[j])
        ax.add_patch(a)
        mx = (nodes[i][0] + nodes[j][0]) / 2 + label_off[0]
        my = (nodes[i][1] + nodes[j][1]) / 2 + label_off[1]
        ax.text(mx, my, label, ha="center", va="center", fontsize=9.5,
                color="0.15" if color == "0.25" else color,
                bbox=dict(boxstyle="round,pad=0.18", fc="white", ec="none",
                          alpha=0.9), zorder=4)

    # cycle principal : hydraulique -> frottement -> erosion -> geometrie
    #                   -> transport
    connect(0, 1, r"$u(x),\ \rho$")
    connect(1, 2, r"$\tau_b$")
    connect(2, 3, r"$\dot m$")
    connect(3, 4, r"$\dot m,\ R\ \rightarrow\ S$", label_off=(-0.02, 0.03))

    # retour rheologique : transport -> hydraulique (le couplage fort)
    a = FancyArrowPatch(nodes[4][:2], nodes[0][:2], arrowstyle="-|>",
                        mutation_scale=20, lw=2.4, color=COLORS[1],
                        linestyle="--", connectionstyle="arc3,rad=-0.30", zorder=1)
    a.set_patchA(boxes[4])
    a.set_patchB(boxes[0])
    ax.add_patch(a)
    ax.text(0.015, 0.52, "rétroaction\nrhéologique\n" r"$\rho(\phi),\ f_m(\phi)$",
            ha="left", va="center", fontsize=10.0, color=COLORS[1],
            bbox=dict(boxstyle="round,pad=0.22", fc="white", ec=COLORS[1],
                      lw=0.9, alpha=0.95), zorder=4)

    # retour geometrique : geometrie -> hydraulique
    connect(3, 0, r"$R(x,t)$", color="0.35", ls=(0, (5, 3)), lw=1.5,
            rad=0.26, label_off=(0.0, 0.02))

    ax.text(0.5, 1.06, "Couplage hydraulique–érosion–transport",
            ha="center", va="bottom", fontsize=12.5, color="0.15")

    fig.tight_layout()
    save(fig, out)


# ===========================================================================
# Figure 3 -- profils spatiaux et fronts multiphysiques
# ===========================================================================
def figure3_profils(out="figure3.png", Nx=FIG3_NX, K_out=FIG3_KOUT):
    print(f"  simulation : campagne B, K_out = {K_out:.0f}, Nx = {Nx}, "
          f"t/ter = {FIG3_SNAPS}")
    p, res = run_case("fig3", t_end_factor=FIG3_SNAPS[-1] + 0.3,
                      snap_t_over_ter=FIG3_SNAPS, Nx=Nx, K_out=K_out)

    ter = res["t_er"]
    x = res["x"] / p.L
    chosen = [pick_snapshot(res["snapshots"], t * ter) for t in FIG3_SNAPS]

    fig, axes = plt.subplots(4, 1, figsize=(7.4, 10.4), sharex=True)
    lim = VALIDITE_RL * p.L / p.R0          # R/L = 0.05  ->  R/R0

    for k, s in enumerate(chosen):
        t_over = float(s.get("t_actual", s["t"])) / ter
        c = COLORS[k % len(COLORS)]
        lab = rf"$t/t_{{er}} = {t_over:.1f}$"
        axes[0].plot(x, s["R"] / p.R0, color=c, label=lab)
        axes[1].plot(x, s["p"] / p.Pin, color=c, label=lab)
        axes[2].plot(x, s["phi"] / p.phi_soil, color=c, label=lab)
        axes[3].plot(x, s["fm"], color=c, label=lab)

    # Limite de validite geometrique R/L = 0.05 : elle correspond a
    # R/R0 = 1667 pour la campagne B (L = 100 m, R0 = 3 mm), donc bien au-dela
    # de la plage utile de la figure : la tracer serait illisible. L'echelle
    # logarithmique est conservee pour montrer que l'elargissement reste
    # faible (R/R0 < 3) sur l'intervalle represente.
    axes[0].set_yscale("log")
    axes[0].set_ylim(1.0, 6.0)

    axes[0].set_ylabel(r"$R/R_0$")
    axes[1].set_ylabel(r"$p/P_{\mathrm{in}}$")
    axes[2].set_ylabel(r"$\phi/\phi_{\mathrm{soil}}$")
    axes[3].set_ylabel(r"$f_m$")
    axes[3].set_xlabel(r"$x/L$")

    axes[2].set_ylim(-0.02, 1.05)
    for ax in axes:
        ax.set_xlim(0.0, 1.0)
        ax.legend(loc="upper left", framealpha=0.92)
    axes[0].set_title(
        "Campagne B — configuration à perte aval "
        rf"($\Delta p = {p.Pin/1e6:.2f}$ MPa, $K_{{\mathrm{{out}}}} = {K_out:.0f}$)",
        fontsize=11)

    fig.tight_layout(h_pad=0.35)
    save(fig, out)


# ===========================================================================
# Figure 4 -- convergence spatiale (sensibilite au maillage)
# ===========================================================================
def figure4_maillage(out="figure4.png", Nx_list=None, K_out=FIG4_KOUT,
                     t_over_ter=FIG4_T_OVER_TER):
    """
    Convergence spatiale des champs couples.

    Trois panneaux de comparaison (pression, rayon, concentration) et un
    quatrieme panneau quantitatif : l'ecart de concentration au maillage le
    plus fin, en echelle logarithmique. C'est ce dernier panneau qui montre
    que l'ecart est localise au front advectif et qu'il decroit avec le
    raffinement, les champs lisses (pression, rayon) etant deja superposes.
    """
    Nx_list = sorted(Nx_list or FIG4_NX)
    print(f"  simulations : Nx = {Nx_list}, t/ter = {t_over_ter}, "
          f"K_out = {K_out:.0f}")

    fig, axes = plt.subplots(4, 1, figsize=(7.4, 10.6), sharex=True)
    styles = ["-", "--", ":"]
    Nx_ref = Nx_list[-1]

    profs = {}
    for k, Nx in enumerate(Nx_list):
        p, res = run_case(f"fig4_Nx{Nx}", t_end_factor=t_over_ter + 0.2,
                          snap_t_over_ter=[t_over_ter], Nx=Nx, K_out=K_out)
        # le dernier instantane est celui de fin de calcul : selectionner la cible
        s = pick_snapshot(res["snapshots"], t_over_ter * res["t_er"])
        t_plot = float(s.get("t_actual", s["t"])) / res["t_er"]
        x = res["x"] / p.L
        profs[Nx] = (x, s["p"] / p.Pin, s["R"] / p.R0, s["phi"] / p.phi_soil)
        st, col = styles[k % len(styles)], COLORS[k % len(COLORS)]
        lab = rf"$N_x = {Nx}$" + ("  (référence)" if Nx == Nx_ref else "")
        axes[0].plot(x, profs[Nx][1], st, color=col, label=lab)
        axes[1].plot(x, profs[Nx][2], st, color=col, label=lab)
        axes[2].plot(x, profs[Nx][3], st, color=col, label=lab)

    # --- panneau 4 : ecart de concentration au maillage de reference -------
    x_ref, _, _, phi_ref = profs[Nx_ref]
    print(f"    ecart L1 / L_inf sur phi/phi_soil, reference Nx = {Nx_ref} :")
    for k, Nx in enumerate(Nx_list[:-1]):
        x, _, _, phi = profs[Nx]
        phi_i = np.interp(x_ref, x, phi)          # interpolation sur la grille fine
        err = np.abs(phi_i - phi_ref)
        L1 = float(np.trapezoid(err, x_ref)) if hasattr(np, "trapezoid") else float(np.trapz(err, x_ref))
        Linf = float(err.max())
        print(f"      Nx = {Nx:4d} : L1 = {L1:.3e}   L_inf = {Linf:.3e}")
        axes[3].plot(x_ref, np.maximum(err, 1e-8), styles[k % len(styles)],
                     color=COLORS[k % len(COLORS)],
                     label=rf"$N_x = {Nx}$  ($L_1 = {L1:.1e}$)")

    axes[0].set_ylabel(r"$p/P_{\mathrm{in}}$")
    axes[1].set_ylabel(r"$R/R_0$")
    axes[2].set_ylabel(r"$\phi/\phi_{\mathrm{soil}}$")
    axes[3].set_ylabel(r"$|\Delta\phi|/\phi_{\mathrm{soil}}$")
    axes[3].set_xlabel(r"$x/L$")
    axes[2].set_ylim(-0.02, 1.05)
    axes[3].set_yscale("log")
    axes[3].set_ylim(1e-6, 1e-1)
    for ax in axes:
        ax.set_xlim(0.0, 1.0)
        ax.legend(loc="upper left", framealpha=0.92)
    axes[0].set_title(
        "Convergence spatiale des champs couplés — campagne B, "
        rf"$t/t_{{er}} = {t_plot:.2f}$, $K_{{\mathrm{{out}}}} = {K_out:.0f}$",
        fontsize=11)

    fig.tight_layout(h_pad=0.35)
    save(fig, out)


# ===========================================================================
def main():
    ap = argparse.ArgumentParser(
        description="Genere docs/figure1.png ... docs/figure4.png du chapitre 4.")
    ap.add_argument("--only", nargs="+", type=int, default=[1, 2, 3, 4],
                    help="numeros des figures a generer (defaut : 1 2 3 4)")
    ap.add_argument("--Nx", type=int, default=FIG3_NX,
                    help=f"maillage de la figure 3 (defaut {FIG3_NX})")
    ap.add_argument("--Kout", type=float, default=FIG3_KOUT,
                    help=f"K_out des figures 3 et 4 (defaut {FIG3_KOUT:.0f})")
    args = ap.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    print("=" * 74)
    print("FIGURES DU CHAPITRE 4  ->  docs/figure1.png ... docs/figure4.png")
    print("=" * 74)

    t0 = time.time()
    if 1 in args.only:
        print("\n[1/4] Schéma du conduit effectif 1D")
        figure1_schema_conduit()
    if 2 in args.only:
        print("\n[2/4] Structure de rétroaction du modèle couplé")
        figure2_boucle()
    if 3 in args.only:
        print("\n[3/4] Profils spatiaux et fronts multiphysiques")
        figure3_profils(Nx=args.Nx, K_out=args.Kout)
    if 4 in args.only:
        print("\n[4/4] Sensibilité au maillage")
        figure4_maillage(K_out=args.Kout)

    print("\n" + "=" * 74)
    print(f"Termine en {time.time() - t0:.1f} s")
    for i in (1, 2, 3, 4):
        f = OUT_DIR / f"figure{i}.png"
        if f.exists():
            print(f"  {f.relative_to(SCRIPT_DIR)}  ({f.stat().st_size / 1024:.0f} ko)")
    print("\nCompiler ensuite le chapitre depuis docs/ (pdflatex).")
    print("Rappel : les .png sont gitignores par le depot (artefacts de build).")


if __name__ == "__main__":
    main()
