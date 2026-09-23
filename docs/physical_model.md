# Physical and Morphological Dynamics Model

The model represents water and suspended eroded soil flowing through an
initially cylindrical, erodible conduit. The conduit radius changes because
of wall shear, while the suspended solid volume changes through advection and
wall-erosion injection.

## 1. Mixture and geometry

The local mixture density is

\[
\rho(\phi)=\rho_w+\phi(\rho_p-\rho_w).
\]

The conduit area and the conservative suspended-solid variable are

\[
A=\pi R^2,
\qquad
S=A\phi.
\]

Here `S` is the solid volume per unit axial length. The concentration used in
the rheology is reconstructed from `phi=S/A`. This area-weighted formulation
is required once the radius varies longitudinally or evolves in time.

## 2. Hydrodynamics and erosion

For the imposed pressure drop, the solver determines the quasi-steady
constant discharge `Q`. The local velocity is

\[
u=Q/A.
\]

The Darcy-Weisbach wall stress is

\[
\tau_b=-\frac18 f_D\rho f_m u^2.
\]

The erosion mass flux is

\[
\dot m=k_{er}\max(|\tau_b|-\tau_c,0),
\]

and the radius evolves according to the strict radial law

\[
R_t=\frac{\dot m}{\rho_{soil,sat}}.
\]

There is no artificial morphological acceleration factor.

## 3. Conservative solid transport

The transported variable satisfies the area-conservative equation

\[
\frac{\partial S}{\partial t}
+\frac{\partial}{\partial x}(\beta u S)
=\mathcal S_S.
\]

The area increment caused by erosion is

\[
\Delta A=A^{n+1}-A^n.
\]

The newly opened saturated soil carries its intact solid fraction, so the
source update is

\[
\Delta S_{source}=\phi_{soil}\Delta A.
\]

This formulation is equivalent to a local concentration relaxation only in
the special case where the area evolution and source are represented through
\(\phi=S/A\). It additionally guarantees that the total suspended solid
volume is tracked consistently as the conduit widens.

## 4. Julien rheological feedback

The raw Julien multiplier scales as

\[
f_{m,raw}-1
\propto
\left(\frac{d_p}{\ell_m}\right)^2
\lambda(\phi)^n.
\]

It is converted into the effective multiplier used in the stress law through

\[
f_m=1+(f_{m,max}-1)
\tanh\left(\frac{f_{m,raw}-1}{f_{m,max}-1}\right).
\]

Thus the suspension increases the hydraulic resistance downstream while the
multiplier remains bounded and smoothly varying with the raw constitutive
value.

## 5. Emergence of the trumpet profile

The downstream widening is not imposed as an initial condition or an
explicit downstream source. It emerges through the following feedback:

1. water enters with the prescribed inlet concentration;
2. wall erosion injects solid volume into the conduit;
3. the conservative variable `S=A*phi` is advected downstream;
4. the concentration and the Julien multiplier increase where the solid load
   accumulates;
5. the larger multiplier increases local wall stress and erosion;
6. the radius therefore grows faster in the downstream region;
7. the larger area lowers the local velocity, providing a stabilising
   hydraulic feedback.

The particle diameter acts quadratically in the raw Julien coefficient. Larger
particles therefore amplify the rheological feedback, while small particles
can leave the multiplier close to unity. The domain length controls both the
integrated pressure loss and the residence distance available for suspended
solid accumulation.

The one-dimensional validity indicator remains important: once `R_out/L`
becomes large, the result should be interpreted as an extrapolation of the
1D model rather than as a fully resolved three-dimensional conduit geometry.
