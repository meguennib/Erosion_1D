# Scientific Freeze v1 — Thesis Reference Solver

Branch: `thesis-v1-freeze`

This branch defines the reference numerical solver to be described in Chapter 3.

## Reference equations

The evolved state is (R(x,t)) and (phi(x,t)). The reference transport equation is

[
\partial_t\phi+\partial_x(u\phi)
=
\frac{2\dot m}{R\rho_s}(\phi_{soil}-\phi),
qquad
u=Q/A,quad A=\pi R^2.
]

No additional (eta) multiplier is used in the reference transport velocity.

The wall shear stress is

[
\tau_b=-\frac18 f_D(Re)f_m(\phi)\rho_m u|u|.
]

The erosion law is

[
\dot m=k_{er}\max(|\tau_b|-\tau_c,0),
]

with ([k_{er}]=\mathrm{s/m}).

The saturated bulk density is

[
\rho_s=\rho_w+\phi_{soil}(\rho_p-\rho_w).
]

## Numerical method

The concentration equation is discretized with cell-centred finite volumes, MUSCL reconstruction, Van Leer limiting, and an upwind/Godunov flux. The inlet concentration is Dirichlet and the outlet is zero-gradient.

The relaxation source is integrated analytically:

[
\phi^{n+1}
=
\phi_{soil}-
(\phi_{soil}-\phi^*)e^{-k^n\Delta t}.
]

The adaptive time step uses

[
\Delta t_{adv}=CFL\frac{\Delta x}{\max |u|}
]

and

[
\Delta t_{src}=\frac{0.5}{\max k},
qquad
k=\frac{2\dot m}{R\rho_s}.
]

The selected step is the minimum of the active constraints, (dt_{max}), and the remaining simulation time. A required step below (dt_{min}) is treated as a numerical failure rather than being forced.

## Pressure-imposed solve

At each time level,

[
F(Q)=p_{calc}(L;Q)-P_{out}=0
]

is solved by bracket expansion followed by bisection. There is no silent fallback to the previous discharge.

## Mixture resistance

The reference closure is the model-specific phenomenological relation

[
f_m=
1+C_B\frac{\rho_p}{\rho_m}
\left(\frac{d_p}{c_lR_0}\right)^2\lambda^2,
]

with

[
\lambda=
\left[
\left(\frac{\phi_{soil}}{\phi}\right)^{1/3}-1
\right]^{-1}.
]

No cap is applied to (f_m) during the solve. `fm_max` is retained only for post-processing compatibility.

## 1-D validity diagnostic

The nominal geometric criterion is

[
D/L<0.1
\Longleftrightarrow
R/L<0.05.
]

The code records the first time (t^*) for which (R_{out}/L>0.05). The simulation is not automatically stopped at (t^*); the diagnostic marks the limit of the nominal reduced-model domain.

## Deliberately excluded from the reference solver

- (eta)-corrected transport velocity;
- hard/smooth (f_m) cap during the solve;
- silent pressure-solver fallback;
- the old (R/L>0.1) validity threshold;
- MinMod as the reference limiter;
- stale experimental parameter branches.

Quantitative verification, convergence studies and physical validation are reserved for Chapter 4.
