# Numerical Methods and Algorithms

This document describes the numerical formulation currently implemented by the
reference solver.

## 1. Conservative finite-volume state

The conduit is discretised on a uniform, cell-centred one-dimensional grid:

\[
\Delta x=L/N_x,\qquad x_i=(i+1/2)\Delta x.
\]

The radius and area are stored at cell centres:

\[
A_i=\pi R_i^2.
\]

The transported conservative variable is the suspended solid volume per unit
axial length:

\[
S=A\phi.
\]

The concentration used by the constitutive laws is reconstructed as

\[
\phi=S/A.
\]

This distinction is essential when the conduit widens. Transporting `phi`
alone does not conserve the total solid volume when `A` varies in space or
time.

## 2. MUSCL finite-volume advection

The conservative transport equation is written as

\[
\frac{\partial S}{\partial t}
+\frac{\partial}{\partial x}(aS)=\mathcal{S}_S,
\qquad a=\beta u.
\]

The production configuration uses a MUSCL reconstruction with the MinMod
limiter and an upwind Godunov flux. A first-order upwind option and the van
Leer/MC limiters remain available for controlled numerical comparisons, but
MinMod is the reference limiter.

The inlet condition is imposed on the conservative variable:

\[
S_{in}=A_{in}\phi_{in},
\]

and the outlet uses a zero-gradient outflow condition.

The reconstruction is spatially TVD under the usual CFL conditions. The
complete coupled method is not claimed to be globally second order in time:
the morphology, conservative transport and source coupling are advanced by
an explicit operator-split update.

## 3. Hydraulic solve

For a trial discharge `Q`, the cell velocity and Reynolds number are

\[
u_i=\frac{Q}{\pi R_i^2},
\qquad
Re_i=\frac{2\rho_i|u_i|R_i}{\mu_w}.
\]

The Barenblatt law provides the full Darcy factor `f_D` and the transport
coefficient `beta`. The wall stress is

\[
\tau_b=-\frac18 f_D\rho f_m u^2.
\]

The pressure gradient contains the wall-friction contribution

\[
p_x=\frac{2\tau_b}{R}.
\]

The outlet minor loss is included in the pressure residual. The discharge is
found with a bracket search followed by bisection. A failed bracket or a
bisection that does not reach the requested pressure tolerance is treated as
a failed simulation step; the solver no longer continues with an arbitrary
fallback discharge.

The current reference implementation is quasi-steady. It does not include a
separate transient momentum equation or an axial Bernoulli/inertial term.

## 4. Julien friction multiplier

The raw Julien multiplier is

\[
f_{m,raw}=1+c_B\frac{\rho_p}{\rho}
\left(\frac{d_p}{\ell_m}\right)^2
\lambda(\phi)^n,
\]

with

\[
\lambda(\phi)=
\frac{1}{(\phi_{soil}/\phi)^{1/3}-1},
\qquad \ell_m=c_lR_0.
\]

The multiplier used by the hydraulic solver is the smooth bounded form

\[
f_m=1+(f_{m,max}-1)
\tanh\left(\frac{f_{m,raw}-1}{f_{m,max}-1}\right).
\]

The hyperbolic tangent provides a smooth transition to the physical upper
bound `fm_max`. Endpoint clipping of `phi` is retained only as a numerical
protection against the singular Julien expression and is documented
separately from the rheological smoothing.

## 5. Erosion and source coupling

The erosion law returns a mass flux per unit wall area:

\[
\dot m=k_{er}\max(|\tau_b|-\tau_c,0)
\qquad [\mathrm{kg\,m^{-2}\,s^{-1}}].
\]

The radial growth law is strictly

\[
R_t=\frac{\dot m}{\rho_{soil,sat}},
\]

with no artificial morphological acceleration factor and no additional
\(\sqrt{1+R_x^2}\) multiplier.

When the wall opens an area increment \(\Delta A\), the conservative source
adds the solid volume carried by the eroded saturated soil:

\[
\Delta S_{source}=\phi_{soil}\Delta A.
\]

This is the conservative counterpart of the former concentration-relaxation
term and ensures that the geometry change does not create or remove solid
volume numerically.

## 6. Physical time stepping

The adaptive step is constrained by the advective speed and by the local
area-growth/source rate:

\[
\Delta t_{adv}=CFL\frac{\Delta x}{\max|\beta u|},
\]

\[
k_A=\frac{2\dot m}{R\rho_{soil,sat}},
\qquad
\Delta t_{src}=\frac{0.5}{\max k_A}.
\]

The reference step is

\[
\Delta t=\min(\Delta t_{adv},\Delta t_{src}),
\]

subject to the configured upper bound. If the required step is below
`dt_min`, the run stops explicitly rather than silently violating the
stability constraint. The last step is truncated so that the requested final
time is not overshot.

The characteristic time is based on the imposed pressure drop:

\[
t_{er}=\frac{2L\rho_{soil,sat}}
{k_{er}(P_{in}-P_{out})}.
\]

It is a physical normalisation and stopping-time scale, not a morphological
acceleration factor.
