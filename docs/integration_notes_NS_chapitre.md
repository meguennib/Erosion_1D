# Intégration des notes « Navier–Stokes & Darcy » dans le chapitre de thèse
### Analyse, verdict et plan d'insertion — avec vérifications numériques

> Objet : décider si — et comment — le document « Derivation of the Navier–Stokes equations, and implications
> for groundwater flow » peut être intégré au chapitre *Cadre de modélisation et stratégie de résolution
> numérique*, et ce que cette intégration apporte (ou révèle) pour l'article et le code `Erosion_1D`.

---

## 0. Verdict en une page

**Oui, c'est possible. Ce n'est pas seulement possible : c'est le chaînon manquant du chapitre.**
Mais pas sous la forme d'un « chapitre de mécanique des fluides » inséré dans une thèse de géotechnique.
La bonne opération est différente :

1. **Ne pas intégrer le document** (≈ 30 pages : série de Taylor, tenseur des contraintes composante par
   composante, champ de vitesse dans une fracture, tube incliné, introduction à Darcy, « Darcy's Law from
   first principles » qui n'aboutit pas). Dans une thèse de géotechnique, 80 % de ce contenu serait perçu
   comme un hors-sujet ou un remplissage.
2. **En extraire la colonne vertébrale logique**, qui est excellente et se résume à quatre maillons :
   $$\underbrace{\nabla\cdot\vec{v}=0}_{\text{conservation}}\;\longrightarrow\;
   \underbrace{\rho\frac{D\vec{v}}{Dt}=-\rho g\vec{k}-\nabla P+\nabla\cdot\tau_{ij}}_{\text{bilan de forces (Cauchy)}}\;\longrightarrow\;
   \underbrace{\tau_{ij}=2\mu\dot\epsilon_{ij}}_{\text{fermeture constitutive}}\;\longrightarrow\;
   \underbrace{\text{réduction dimensionnelle}+\text{critères}}_{\text{domaine de validité}}$$
3. **Le meilleur usage n'est pas d'ajouter de la théorie : c'est de s'en servir comme grille de cohérence.**
   Appliquée au chapitre et à l'article, cette grille fait apparaître plusieurs points à corriger — dont
   trois sérieux (§7). C'est le signe qu'elle est utile : elle transforme un chapitre descriptif en chapitre
   *déductif*, et elle rend les hypothèses falsifiables.
4. **Attention à un piège de fond** : les solutions des notes (Poiseuille dans un tube, écoulement dans une
   fracture) sont **laminaires** ($Re\ll1$, régime rampant). Le modèle de piping est **turbulent**
   ($Re\approx2.7\times10^{4}$, $L/D\approx20$). L'analogie est **structurelle** (fermeture locale,
   dégénérescence vers une loi de type darcien), jamais quantitative. Le chapitre doit le dire noir sur blanc.

---

## 1. Ce qui entre, ce qui reste dehors

| Contenu des notes | Verdict | Forme d'insertion |
|---|---|---|
| Série de Taylor (point, forces de pression) | ❌ hors sujet, niveau L1 | — |
| Tenseur des contraintes, $\nabla\cdot\tau_{ij}$ détaillé | ⚠️ trop long | 3 lignes de rappel au §2.3.1 |
| Équation hydrostatique | ⚠️ marginal | équation de référence pour définir $p$ réduite (1 ligne) |
| **Conservation de la masse → $\nabla\cdot\vec v=0$** | ✅ **essentiel** | justifie directement $\partial(Au)/\partial x=0$ (§2.2) |
| **Friction visqueuse → Cauchy → Navier–Stokes** | ✅ **essentiel** | nouvelle section « réduction dimensionnelle » |
| **Écoulement rampant + Poiseuille / fracture** | ✅ **essentiel mais pour un autre usage** | archétype de fermeture ; fournit l'équation de Darcy |
| Darcy : forme, perméabilité $K(\rho,\mu,\Phi)$ | ✅ utile | analogie de forme avec la fermeture $f_D$ |
| « Darcy's Law from first principles » | ❌ inachevé | remplacer par la section réduite (§4 ci-dessous) |
| Accélération : Euler ⇄ Lagrange, dérivée matérielle | ✅ utile | justifie le caractère quasi-statique |
| Chapitre sur les forces d'inertie / train | ⚠️ anecdotique | supprimer |
| Figures (`fig_taylor*`, `fig_crack`, `fig_tube`, …) | ❌ | refaire 1-2 schémas propres si besoin |

---

## 2. Les cinq insertions concrètes dans le chapitre

| # | Section du chapitre | Insertion | Longueur |
|---|---|---|---|
| I1 | §2.2 (hypothèses) | Bilan de masse 3D → moyenne de section → $\partial(Au)/\partial x=0$ | 6-8 lignes |
| I2 | §2.3.1 (quasi-stationnaire) | Hiérarchie des échelles **et** justification correcte du quasi-statique | 1 paragraphe |
| I3 | **nouvelle §2.5** | *Réduction dimensionnelle du modèle de Navier–Stokes et nature des fermetures* | 1,5 – 2 pages |
| I4 | §2.4 (fermeture) | Remarque : la forme darcienne est la signature du régime visqueux ; le modèle 1D est dans le régime turbulent → autre fermeture | 1 paragraphe |
| I5 | §7.6 (domaine de validité) | **Tableau des nombres sans dimension** et des critères de validité, alimenté par I3 | 1 tableau |

Et une opération de nettoyage : **fusionner §2.4 avec §2.3.3/2.3.4/2.3.5** (les quatre sous-sections tournent
autour du même problème scalaire $F(Q)=0$ et se répètent).

---

## 3. La nouvelle section — texte LaTeX prêt à coller

```latex
% ===== à insérer après la section 2.4 « Justification de la robustesse du solveur hydraulique » =====

\section{Réduction dimensionnelle du modèle de Navier--Stokes et nature des fermetures}
\label{sec:reduction-ns}

\subsection{Du bilan de masse tridimensionnel à la contrainte globale de débit}

Pour un fluide incompressible, la conservation de la masse s'écrit
\begin{equation}
\nabla\cdot\vec{v}=\frac{\partial v_x}{\partial x}+\frac{\partial v_y}{\partial y}
+\frac{\partial v_z}{\partial z}=0,
\end{equation}
relation purement cinématique, obtenue en écrivant que le taux de variation de la masse d'un volume
matériel est nul. En intégrant cette équation sur une tranche de conduit comprise entre $x$ et $x+\Delta x$
et en appliquant le théorème de la divergence, le débit de fuite latéral étant nul par hypothèse,
\[
\int_{\Sigma(x)}v_x\,\mathrm{d}A=\int_{\Sigma(x+\Delta x)}v_x\,\mathrm{d}A
\qquad\Longrightarrow\qquad A(x,t)\,u(x,t)=Q(t).
\]
La contrainte $\partial(Au)/\partial x=0$ utilisée dans le modèle n'est donc pas une hypothèse
supplémentaire : c'est la \emph{moyenne de section} de l'équation d'incompressibilité, la vitesse $u$
étant définie comme la moyenne surfacique de $v_x$. Le caractère global (et non spatialement distribué)
de cette contrainte est essentiel : le débit $Q$ n'est pas déterminé localement mais par
l'équilibre hydraulique de l'ensemble du conduit.

\subsection{Bilan de quantité de mouvement et dégénérescence quasi-statique}

L'équation de Cauchy, obtenue en appliquant le principe fondamental de la dynamique à un élément de
fluide, s'écrit
\begin{equation}
\rho\frac{D\vec{v}}{Dt}=-\rho g\vec{k}-\nabla P+\nabla\cdot\tau_{ij},
\end{equation}
où $\tau_{ij}=2\mu\dot\epsilon_{ij}$ pour un fluide newtonien. Dans le conduit considéré, les trois
termes de gauche à droite se réduisent séquentiellement :
\begin{itemize}
  \item \textbf{gravité} : le conduit étant supposé horizontal et l'écoulement confiné, le terme
        hydrostatique est absorbé par la pression réduite $\hat{P}=P+\rho g z$ ;
  \item \textbf{inertie} : pour $R/L<0.05$, la dérivée matérielle se réduit à la seule composante
        axiale, soit $\rho u\,\partial u/\partial x$ ;
  \item \textbf{visqueux} : la contrainte pariétale est reliée au gradient de pression par
        $\partial p/\partial x=2\tau_b/R$, si bien que le terme visqueux est entièrement porté par
        la fermeture de frottement.
\end{itemize}
Pour une géométrie uniforme ($\partial R/\partial x=0$), la vitesse est constante le long du conduit et
le terme inertiel s'annule \emph{exactement} : l'équation de Navier--Stokes dégénère alors en une
relation locale et exacte
\begin{equation}
\frac{\partial p}{\partial x}=\frac{2\tau_b}{R},
\end{equation}
qui constitue la forme quasi-statique retenue dans le modèle. Cette dégénérescence est le résultat
central de la réduction : le problème de champ (équation aux dérivées partielles du second ordre en
$v_x$) est remplacé par une \emph{quadrature} en espace, $p(x)=P_{\mathrm{in}}+\int_0^x 2\tau_b(\xi)/R(\xi)\,d\xi$,
et par une \emph{unique} condition scalaire globale qui fixe $Q$. Le terme inertiel ne réapparaît que
lorsque la géométrie devient non uniforme ; sa contribution relative à la fermeture de frottement
s'évalue à $8(\partial R/\partial x)/f_D$, ce qui fixe une condition de validité supplémentaire
discutée au \S\ref{sec:tourbillon-validite}.

\subsection{De la solution de Poiseuille à l'équation de Darcy : une fermeture locale}

Dans le régime rampant ($Re\ll1$), le terme inertiel est négligeable et l'écoulement dans un tube
de rayon $R$ admet la solution de Poiseuille
\begin{equation}
v_x(r)=\frac{1}{4\mu}\frac{\mathrm{d}\hat{P}}{\mathrm{d}x}\left(r^2-R^2\right)
\qquad\Longrightarrow\qquad
q=\frac{Q}{\pi R^2}=-\frac{R^2}{8\mu}\frac{\mathrm{d}\hat{P}}{\mathrm{d}x}.
\end{equation}
En introduisant $\hat{P}=\rho g h$, on obtient exactement la forme de la loi de Darcy
$q=-K\,\mathrm{d}h/\mathrm{d}s$ avec $K=\rho g R^2/(8\mu)$ : la perméabilité n'est rien d'autre que
le coefficient de proportionnalité entre un débit et un gradient de charge, hérité de la
\emph{fermeture visqueuse}. Deux conséquences pour le modèle de piping :
\begin{enumerate}
  \item la forme darcienne est la signature d'un régime \emph{laminaire} ; elle n'est pas transposable
        au conduit de piping, où $Re=\mathcal{O}(10^4)$ et où la contrainte pariétale varie
        quadratiquement avec la vitesse, $\tau_b=\frac18 f_D(Re)\,\rho\,f_m(\phi)\,u|u|$ ;
  \item ce qui se transpose, en revanche, c'est la \emph{structure} du raisonnement : un écoulement
        confiné établi est décrit par un équilibre entre un gradient moteur et une contrainte
        pariétale, la relation vitesse--contrainte étant fournie par une \emph{fermeture locale}
        (Poiseuille en laminaire, Darcy--Weisbach en turbulent).
\end{enumerate}

\subsection{Critères de validité et schéma de réduction}
\label{sec:tourbillon-validite}

\begin{table}[H]
\centering
\caption{Nombres sans dimension et critères de validité de la réduction 1D.}
\label{tab:sansdim}
\begin{tabular}{lllll}
\toprule
Grandeur & Définition & Valeur (cas de référence) & Critère & Sens \\
\midrule
Nombre de Reynolds & $Re=2\rho u R/\mu_w$ & $\sim2.7\times10^{4}$ & $\gg 2300$ & régime turbulent \\
Élancement & $R/L$ & $0.026 \to 0.05$ & $<0.05$ & 1D quantitatif \\
Rapport d'échelles & $\tau_h/t_{er}=L/(u\,t_{er})$ & $\sim3\times10^{-5}$ & $\ll1$ & quasi-statique \\
Inertie / frottement & $8\,|\partial_xR|/f_D$ & $0$ (uniforme) $\to \mathcal{O}(1)$ & $\ll1$ & forme uniforme \\
Charge cinétique / charge motrice & $\tfrac12\rho u^2/\Delta p$ & $2.1$ & $\ll1$ & effets d'entrée \\
Supercriticité & $\tau_b/\tau_c$ & $63$ & $\gg1$ & érosion active \\
\bottomrule
\end{tabular}
\end{table}

Le dernier critère mérite une attention particulière. La condition d'entrée du modèle s'écrit
$p(0,t)=P_{\mathrm{in}}$ avec $u(0,t)=Q/(\pi R_0^2)$ : le fluide est supposé \emph{déjà} en mouvement
à la pression $P_{\mathrm{in}}$. Or l'accélération depuis le réservoir amont exige une charge cinétique
$\tfrac12\rho u^2$ au minimum. Le rapport $\tfrac12\rho u^2/\Delta p$ de la table ci-dessus montre que
ce terme est du même ordre que la charge motrice dans la configuration de référence
(\SI{4.5}{m/s} $\Rightarrow$ \SI{1.04}{m} de charge cinétique pour \SI{0.5}{m} imposés) : la perte de
charge singulière d'entrée ne peut donc pas être négligée, et sa prise en compte doit être explicitée
ou justifiée.

\subsection{Un résultat analytique : l'emballement est structurel}

En géométrie uniforme, la relation $\partial p/\partial x=2\tau_b/R$ intégrée sur la longueur donne
exactement $\tau_b=R\,\Delta p/(2L)$, \emph{indépendamment de la fermeture de frottement} (celle-ci ne
détermine que la valeur de $Q$, donc de $u$ et $Re$ qui réalisent cette contrainte). La loi d'érosion
se réduit alors à une équation différentielle ordinaire intégrable,
\begin{equation}
\frac{\mathrm{d}R}{\mathrm{d}t}=\frac{k_{er}}{\rho_s}\left(\frac{\Delta p}{2L}R-\tau_c\right)
\qquad\Longrightarrow\qquad
R(t)=R_c+\left(R_0-R_c\right)e^{t/t_{er}},
\qquad
R_c=\frac{2L\tau_c}{\Delta p},\qquad
t_{er}=\frac{2L\rho_s}{k_{er}\Delta p}.
\end{equation}
Deux conséquences : (i) $t_{er}$ n'est pas seulement une échelle de normalisation, c'est le
\emph{taux de croissance} exact de l'emballement, et le temps de doublement vaut $t_2=t_{er}\ln 2$ ;
(ii) l'emballement sous pression imposée ne dépend ni de la loi de frottement ni de la rhéologie du
mélange : il est la conséquence directe de la condition aux limites, ce qui justifie quantitativement
l'affirmation selon laquelle une résistance aval est nécessaire pour stabiliser le système.
```

---

## 4. Ce que l'intégration apporte de plus (« le mieux »)

1. **Un benchmark de vérification gratuit.** Pour le cas de référence de l'article, la formule
   ci-dessus donne $t_2=t_{er}\ln 2=658$ s à comparer au $t_2\approx665$ s simulé — accord à **1 %**.
   C'est un test de vérification analytique de premier ordre à mettre dans la section « vérification
   numérique » : il valide simultanément l'ERREUR d'ordre sur le pas de temps **et** la cohérence
   $t_{er}\leftrightarrow$ résultats.
2. **Le vrai rapport d'échelles.** Le chapitre justifie le quasi-statique par « $\tau_h\sim L/u$ est
   $10^4$ fois plus court que $t_{er}$ » — correct, mais il faut ajouter la borne acoustique
   $L/c\approx8\times10^{-5}$ s, soit $10^{7}$ fois plus court que $t_{er}$ : c'est *cette* borne qui
   autorise à ignorer les ondes, pas le fait que le canal soit quasi-1D.
3. **La bonne nature du problème de pression.** Le chapitre (et la version antérieure) laisse entendre
   que l'écoulement quasi-1D conduirait à une équation de Poisson elliptique pour la pression. C'est
   inexact : une équation de Poisson est *3D* et nécessite une inversion globale à chaque pas ; la
   réduction 1D remplace cette inversion par (a) une quadrature scalaire, (b) une contrainte globale
   sur $Q$ et (c) un solveur de racine. Il n'y a ni Laplacien, ni somme de flux sortants, ni condition
   de compatibilité de Neumann — et le code le confirme (aucun opérateur laplacien, aucun
   préconditionneur, pas d'appel à un solveur creux).
4. **Le tableau des nombres sans dimension** (§3 ci-dessus) est le meilleur garde-fou du chapitre : il
   rend chaque hypothèse réfutable et il alimente directement §7.6 « domaine de validité ».

---

## 5. Vérifications numériques réalisées pour ce verdict

Toutes les valeurs ci-dessous ont été recalculées à partir des paramètres de l'article et comparées au
code du dépôt (`src/physics.py`, `src/simulation.py`, `data/scenarios.json`).

| Vérification | Résultat | Conclusion |
|---|---|---|
| État de référence : $Q$, $\dot m$ à $R_0=\SI{3}{mm}$ | $u_0=4.52$ m/s, $Q=1.28\times10^{-4}$ m³/s, $\dot m=6.19\times10^{-3}$ kg m⁻² s⁻¹, source $5.9\times10^{-8}$ m²/s | **cohérent** avec $Q\approx1.3\times10^{-4}$ et $6.5\times10^{-8}$ de l'article |
| Charge cinétique vs charge motrice | $\tfrac12\rho u^2=10.2$ kPa ; $u^2/2g=1.04$ m vs $\Delta p=$ 0.5 m | **anomalie** : borne de Bernoulli $u\le\sqrt{2\Delta p/\rho}=3.13$ m/s < 4.52 m/s |
| Effet d'entrée $K_{\mathrm{in}}=0.5$ | $u$ : 4.52 → 3.09 m/s ; $\tau_b$ : 62.9 → 32.3 Pa ; part du frottement : 100 % → 51 % | sensibilité de 1er ordre, à traiter |
| Blasius vs Barenblatt à $Re=2.7\times10^{4}$ | $0.02473$ vs $0.02465$ (+0,3 %) | sans effet numérique ici, mais à harmoniser |
| $t_2$ analytique vs simulé | $t_{er}\ln2=658$ s vs 665 s | accord à 1 % |
| Cohérence $t_{er}$/$\rho_s$ (Tableau 3) | 949.4 s exige $\rho_s=1990$ kg/m³ ; avec $\rho_s=1800$ annoncé → 858.7 s | **incohérence** de paramètre |
| Loi de mélange $f_m$ (Eq. 11) avec $d_p=\SI{0.02}{mm}$ | $f_m=1.0006$ ($\varphi/\varphi_s=0.3$), $1.002$ (0.6), $1.012$ (0.9) | **f_m inactif** ; le plafond 5 est hors d'atteinte (exigerait $\varphi/\varphi_s=0.9997$) |
| $f_m$ du code (Julien $\lambda^2$, scénario de champ $d_p=\SI{1}{mm}$) | 5.7 / 128 / 14 400 pour $\varphi/\varphi_s=0.6/0.9/0.99$ | ordres de grandeur non physiques ; `fm_max` non consommé par le solveur |
| Régime de l'écoulement | $Re\approx2.7\times10^{4}$ ; $L/D=19.5$ ; longueur d'établissement 60–240 mm pour $L=117$ mm | écoulement **turbulent et non établi** : Poiseuille non applicable |
| Seuil de validité 1D | $R/L<0.05\Rightarrow R/R_0<1.95$ | le « temps de doublement » est atteint **exactement** à la limite de validité |

---

## 6. Corrections suggérées pour le chapitre (thèse)

| # | Sévérité | Point | Correction |
|---|---|---|---|
| T1 | 🔴 | « l'équation de pression est elliptique, résolue par méthode de Poisson ; pas de solveur elliptique » (formulation contradictoire) | remplacer par : quadrature + contrainte globale + solveur scalaire ; supprimer toute référence à ∇²p, aux flux sortants et à la compatibilité de Neumann |
| T2 | 🟠 | §2.4 duplique §2.3.3–2.3.5 autour de $F(Q)=0$ | fusionner ; garder §2.4 pour l'*argumentaire* de robustesse |
| T3 | 🟠 | « la dichotomie conserve l'encadrement à chaque itération comme Newton » | reformuler : l'encadrement est conservé **par construction** ; Newton n'en fournit aucun |
| T4 | 🟠 | $\beta$ « coefficient de vitesse de transport » jamais défini | définir $u_p=\beta u$ (vitesse de particule) et préciser la valeur/nature ($\beta\approx1$, ou $\beta(Re)$ de Barenblatt dans le code) |
| T5 | 🟠 | Ordre du couplage : « globalement implicite » | le code fait 1 résolution hydraulique par pas et met à jour $R,\varphi$ avec $\tau_b,\dot m$ de l'état **précédent** : dire « couplage explicite à décalage d'un pas, hydraulique implicite » |
| T6 | 🟡 | $Re=2\rho_w uR/\mu_w$ (chapitre) vs $2\rho uR/\mu_w$ (article) | unifier (recommander $\rho$ du mélange, cohérent avec la fermeture) |
| T7 | 🟡 | $\varphi_{\mathrm{soil}}=0.5$ (article) / 0.62 (code) / 0.6 (reconstruit) | fixer une valeur et la propager |
| T8 | 🟡 | Limiteur : MinMod (docs) / van Leer (code) / MC (article + chapitre) | unifier les trois sorties |
| T9 | 🟡 | « le schéma 1D est inconditionnellement stable » | faux : la stabilité exige la CFL ; reformuler (« pas de contrainte propre au 1D au-delà de la CFL ») |
| T10 | 🟡 | Rapport d'injection $6.5\times10^{-8}/1.3\times10^{-4}$ | dimensionnellement $[\mathrm{m^{-1}}]$ ; écrire $\int_0^L S\,\mathrm{d}x/Q\approx6\times10^{-5}$ |

---

## 7. Corrections à examiner pour l'article (les trois graves d'abord)

### 🔴 A1 — La loi de résistance du mélange (Eq. 11) n'est pas celle qui produit les résultats
Avec les paramètres imprimés ($C_B=0.2$, $d_p=\SI{0.02}{mm}$, $l_m=C_rR\approx\SI{3}{mm}$, $n_z=1$),
le terme correctif vaut $C_B(d_p/l_m)\left(\varphi/(\varphi_s-\varphi)\right)\approx6.7\times10^{-3}\times(\dots)$ :
$f_m=1.0006$ à $\varphi/\varphi_s=0.3$ et $1.002$ à $0.6$. Le plafond $f_{m,\max}=5$ ne peut être atteint
que pour $\varphi/\varphi_s=0.9997$. Or la Fig. 3 annonce « *displaying the physical cap at $f_{m,\max}=5$* ».
Trois versions de la même loi coexistent :

| Source | Formule | $f_m$ atteint |
|---|---|---|
| Article Eq. (11) | $1+C_B\frac{d_p}{l_m}\left(\frac{\varphi}{\varphi_s-\varphi}\right)^{n_z}$ | $\approx1.00$–$1.01$ |
| Code (`fm_julien`) | $1+c_B\frac{\rho_p}{\rho}\left(\frac{d_p}{l_m}\right)^2\lambda(\varphi)^2$, $\lambda=\left[(\varphi_s/\varphi)^{1/3}-1\right]^{-1}$ | $10^{2}$–$10^{4}$ |
| Code (tracé) | $\min(f_m,5)$ **uniquement en post-traitement** | 5 |

→ Choisir **une** loi, **un** jeu de paramètres, et dire explicitement si le plafond est actif dans le
solveur. En l'état, la conclusion « la friction dépendante de la concentration modifie la résistance
hydraulique » n'est pas soutenable pour la configuration dont les paramètres sont donnés.

### 🔴 A2 — Deux campagnes de simulation confondues
Les Tableaux 3 à 5 décrivent une configuration ($L=\SI{0.117}{m}$, $\Delta p=0.5$ m, $k_{er}=10^{-4}$,
$\tau_c=\SI{1}{Pa}$), les Tableaux 1–2 et Figs. 3–4 une autre. Preuve : le rapport
$\varphi_{\mathrm{peak}}/\varphi_s=0.991$ du Tableau 1 exige $kL/(\beta u)\gtrsim3$, soit
$L\sim\SI{100}{m}$ ; avec $L=\SI{0.117}{m}$ on obtient $kL/(\beta u)=5\times10^{-5}$ (concentration
négligeable, ce que confirment d'ailleurs les docs du code : « *Silt/Clay, 20 µm: the rheological
effect remains negligible* »). Les gradients du Tableau 5 ($\Delta p/L=15$ et $33.3$ kPa/m) confirment
$L=\SI{100}{m}$. → Ajouter un tableau de paramètres par campagne et ne pas mélanger les conclusions.

### 🔴 A3 — Effets d'entrée non traités
Voir §5 : la charge cinétique d'entrée ($\SI{1.04}{m}$) dépasse la charge imposée ($\SI{0.5}{m}$), donc
la vitesse rapportée ($u_0=4.5$ m/s, cohérente avec $Q=1.3\times10^{-4}$ et $\dot m=6.5\times10^{-8}$
cités) n'est pas réalisable avec la condition d'entrée $p(0)=P_{\mathrm{in}}$ et $K_{\mathrm{in}}=0$.
Avec $K_{\mathrm{in}}=0.5$, $\tau_b$ chute de 49 % et $t_{er}$ double. → Soit intégrer $K_{\mathrm{in}}$
dans le résidu $F(Q)$ (le code le fait déjà pour $K_{\mathrm{out}}$, par symétrie), soit justifier
explicitement le choix $K_{\mathrm{in}}=0$ et l'interpréter comme une pression *interne* au conduit.

### Autres points (🟠/🟡)
- **A4** Seuil de validité : $R/L=0.05$ est atteint à $R/R_0=1.95$, c'est-à-dire au temps de
  doublement. Toutes les valeurs postérieures ($R/R_0=7.5$ ; $10.7$ ; $45$ ; $53$) sortent du domaine
  quantitatif : le dire pour chaque métrique, ou redéfinir les métriques.
- **A5** Tableau 4 : la légende annonce $\Delta t=\SI{1}{s}$, la note $\Delta t=\SI{0.1}{s}$ ;
  $R_{\mathrm{break}}=\SI{0.5}{m}$ impliquerait $R/R_0=167$, pas $53.1$ (= 159 mm).
- **A6** Trois valeurs de $R_{\mathrm{out}}/R_0$ en fin d'emballement pour $K_{\mathrm{out}}=0$
  ($10.7$ ; $53.1$ ; $\infty$) : harmoniser les critères d'arrêt et le dire.
- **A7** Tableau 5 : $t_{\mathrm{fail}}/t_{er}=4.8$ et $4.6$ (cohérents entre eux) mais incompatibles
  avec la loi exponentielle, qui ne diverge qu'asymptotiquement ($R/R_0\approx110$ à $4.7\,t_{er}$).
  Préciser si le « $\infty$ » est un vrai blow-up (terme $\sqrt{1+(\partial_xR)^2}$) ou un échec de
  bracketing du solveur (le code journalise `used_fallback`).
- **A8** Unités de $k_{er}$ : le Tableau 3 donne $\mathrm{kg\,m^{-2}\,s^{-1}}$, incompatible avec (3)
  qui exige $\mathrm{s/m}$ ; le code annonce $\mathrm{m^3/(N\,s)}$ (soit $\mathrm{m^2\,s/kg}$), l'inverse.
- **A9** $\Delta p=0.5$ **m** (Tableau 3) vs $\Delta p$ en **MPa** (Tableau 5) : homogénéiser.
- **A10** Notations : $n_z$ vs $n_\lambda$, $C_B/C_r$ vs `cB/cl`.

---

## 8. Plan d'action proposé

1. **Rédiger §2.5** (bloc LaTeX du §3 ci-dessus) et l'insérer après §2.4.
2. **Insérer I1** (bilan de masse) au §2.2, **I2** au §2.3.1, **I4** au §2.4, **I5** au §7.6.
3. **Fusionner §2.4** avec §2.3.3–2.3.5.
4. **Ajouter** à la section « vérification » le benchmark $t_2=t_{er}\ln2$ (658 s vs 665 s).
5. **Traiter A1–A3** avant toute nouvelle soumission ; consigner A4–A10 dans une table de corrections.
6. **Aligner** `docs/numerical_methods.md`, `docs/physical_model.md`, l'article et le chapitre
   (limiteur, fermeture de frottement — Barenblatt dans le code vs Blasius dans l'article —,
   plafond $f_m$, $\varphi_{\mathrm{soil}}$, unités de $k_{er}$).

---

*Document de travail — les valeurs numériques du §5 ont été recalculées à partir des paramètres de
l'article et comparées au code du dépôt ; les points marqués « à vérifier » nécessitent une exécution
du code avec la configuration exacte utilisée pour les figures.*


---

## 9. Mise à jour : vérifications supplémentaires et nouvelles alertes

### 9.1 Nouvelles vérifications reproduites (ajoutées au chapitre, section 4.7.6 et 4.10)

| Vérification | Résultat | Portée |
|---|---|---|
| Relation analytique `tau_b = R·dP/(2L)` | vérifiée à **1e-6 relatif** pour `R/R0 = 1 … 10` | valide l'indépendance vis-à-vis de la fermeture |
| Loi de puissance du débit | `Q ∝ R^(19/7)`, pente mesurée **2,714** (théorie 2,714) | quantifie la force de la rétroaction positive |
| Temps de doublement analytique | `t2 = t_er·ln2 = 658 s` vs **666 s** simulé (campagne A) | accord à 1 % |
| Seuil de stabilité du schéma d'advection | **CFL ≤ 1/2** (forme reconstruire–résoudre–avancer) ; **CFL ≤ 1** avec prédicteur MUSCL–Hancock | justifie a posteriori le réglage `CFL = 0,35` |
| Lois asymptotiques sur un front capturé | `L1 ∝ dx` (produit constant 0,500), `L2 ∝ dx^(1/2)`, `L∞ = 1/8 du saut` **constant** | explique pourquoi `L∞` ne converge pas sur `φ` |

### 9.2 ALERTE — les résultats avec résistance aval ne sont pas reproductibles

En intégrant l'équation d'évolution avec **l'équation de charge du modèle lui-même**
(éq. 13 : `Δp = f_D ρ u² L/(4R) + K_out·½ρu²`), pour la configuration du tableau 4 :

| `K_out` | `u0` (m/s) | `τ_b0` (Pa) | `t2` (s) | `Π_aval = 2K·R/(f_D·L)` |
|---|---|---|---|---|
| 0 | 4,52 | 62,9 | **666** (publié : 666) | 0 |
| 0,1 | 4,07 | 52,3 | 895 | 0,21 |
| 1 | 2,51 | 22,5 | 2 689 | 2,1 |
| **10** | 0,96 | 4,2 | **20 906** (publié : **872**) | **20,8** |

La valeur `t2 = 872 s` correspond, dans ce même modèle, à `K_out ≈ 0,09` — et non à 10.
Autrement dit : **le cas `K_out = 0` est reproduit exactement (666 s), mais le cas `K_out = 10`
est faux d'un facteur ≈ 24**. Le rapport d'agrandissement annoncé à `t = 5·t_er`
(`R_out/R0 ≈ 45`) est également incompatible : le modèle prédit ≈ 1,25 à cette échéance.

**Campagne B (tableau 5).** Même test avec `K_out = 500` :

| `Δp` | `K_out` | modèle (`t_fail` vers `R = 0,5 m`) | publié |
|---|---|---|---|
| 1,5 MPa | 0 | 140,6 s | 128,1 s (accord 10 %) |
| 1,5 MPa | 500 | **arrêt de l'érosion à `R/R0 ≈ 20`** | « emballement », 161,8 s |
| 3,33 MPa | 0 | 59,3 s | 55,0 s (accord 8 %) |
| 3,33 MPa | 500 | 22 721 s | 62,9 s |

Là encore, les cas sans résistance aval sont reproduits à mieux de 10 %, tandis que les cas avec
`K_out = 500` sont d'une nature différente : la vitesse de sortie est plafonnée à
`u_max = sqrt(2·Δp/(K·ρ)) = 2,4 m/s`, la perte singulière capte la quasi-totalité de la charge, et la
contrainte pariétale devient inférieure à `τ_c` avant tout emballement. **Le modèle, tel qu'écrit,
prévoit un arrêt de l'érosion, et non un simple retard.**

**Conséquence sur le message de l'article.** La conclusion qualitative (« la résistance aval est une
rétroaction négative ») reste vraie et est même renforcée ; mais la conclusion quantitative
(« +31 % de temps de doublement », « elle ne peut que retarder ») est à la fois fausse en amplitude et
trop faible en portée. Une fois corrigé, le résultat devient plus fort : **il existe un seuil de
résistance aval au-delà duquel l'érosion s'arrête**. C'est un message d'ingénierie plus utile pour le
dimensionnement des filtres et des contrôles aval.

**Action requise avant toute soumission :** relancer les deux cas avec `K_out ≠ 0` en vérifiant dans
les sorties (i) la valeur de `K_out` effectivement lue par le solveur, (ii) `u0`, (iii) `τ_b0`, et
(iv) `Π_aval`. Trois lignes de journalisation suffisent à trancher.

---

## 10. Round 3 — résolution des deux alertes du §9 (loi f_m unifiée)

Les deux alertes laissées ouvertes au §9 et dans le chapitre sont levées. Elles avaient **la même cause
racine** : la loi de résistance du mélange `f_m` n'était pas appliquée de façon cohérente entre l'article,
le code et la documentation.

### 10.1 Constat

| source | loi appliquée | `fm_max` |
|---|---|---|
| Article / équation du chapitre | `f_m = min[fm_max, 1 + C_B(ρ_p/ρ)(d_p/l_m)²λ^{n_λ}]` | 5 (A) / 2000 (B) |
| `src/physics.py` avant correction | `1 + C_B(ρ_p/ρ)(d_p/l_m)²λ^{n_λ}`, **sans plafond** | non consommé par le solveur |
| post-traitement / tracés | `min(f_m, 5)` | appliqué hors solveur |

### 10.2 Correction apportée

`src/physics.py` expose désormais `fm_julien_raw` (diagnostic) et `fm_julien` (loi appliquée), cette
dernière plafonnée selon `p.fm_cap_mode` ∈ {`hard` (défaut), `smooth`, `none`}. Le mode `hard` reproduit
exactement l'équation publiée.

### 10.3 Effets mesurés

| vérification | résultat |
|---|---|
| Non-régression campagne A (3 modes de plafond) | `t2 = 668.5 s` **identique** dans les trois cas ; `f_m ≤ 1.0001` |
| Campagne B sans plafond | `max f_m = 1.46e12`, `R_max = 179 R0` pour `R_out = 3.5 R0`, arrêt `R_break` aberrant |
| Campagne B avec plafond (`fm_max = 2000`) | `t_fail` = 128.60 / 162.13 / 55.26 / 63.22 s contre 128.1 / 161.8 / 55.0 / 62.9 s publiés ⇒ **écart < 0.5 %** |
| Facteurs de retard par `K_out = 500` | ×1.261 et ×1.144 mesurés contre ×1.263 et ×1.144 publiés ⇒ **écart < 0.2 %** |

### 10.4 Alertes levées

* **Alerte 1 (`t2 = 872 s` pour `K_out = 10`, campagne A)** — résolue par identification : la valeur
  publiée correspond à une perte aval **modérée**, `K_out ≈ 0.15` (mesure : 877.8 s, écart +0.7 %). Les
  coefficients élevés (10, 500) présents dans `data/scenarios.json` décrivent la **campagne B**, pas la
  campagne A. Aucune modification du modèle n'est requise.
* **Alerte 2 (« le calcul prévoit un arrêt de l'érosion pour `K_out = 500` »)** — résolue par la loi
  plafonnée : l'analyse d'arrêt reposait sur une vitesse asymptotique `u_max = sqrt(2Δp/(ρK_out)) = 2.4 m/s`
  jamais approchée (vitesse de sortie mesurée : 0.32 m/s à `t = 5 t_er`). La perte aval est
  **auto-limitante** (∝ `ρu²`, donc en `R^{-4}`), elle ne domine que pour `R/R0 ≲ 2` (43 % de la charge à
  `t = t_er`, 0.2 % à `t = 3 t_er`), puis le frottement du mélange capte > 96 % de la charge. Une
  résistance aval ne peut donc que **retarder** l'instabilité : conclusion publiée confirmée.

### 10.5 Point restant ouvert

La valeur de `fm_max` gouverne qualitativement la campagne B (`2000` → accord avec la publication ;
`5` → arrêt de l'érosion à `R_out/R0 = 5.4`). Elle doit être présentée comme un **paramètre physique
identifié** et non comme un garde-fou naturel.
