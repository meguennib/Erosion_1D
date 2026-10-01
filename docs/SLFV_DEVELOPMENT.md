# Conservative semi-Lagrangian FV prototype (thesis-v2)

This branch is experimental. The frozen thesis reference solver is kept on
`thesis-v1-freeze`.

## Purpose

The v1 diagnostic showed that the explicit transport CFL is the dominant
time-step restriction. This prototype tests a conservative finite-volume
semi-Lagrangian transport operator without changing the continuous model

[
\partial_t\phi+\partial_x(u\phi)=0
]

and keeps the relaxation source in exact form

[
\partial_t\phi=k(\phi_{soil}-\phi).
]

## Algorithm

For each target cell, both faces are traced backward along the frozen velocity
field. The departure interval is then integrated against the MUSCL/Van-Leer
reconstruction of the old cell averages.

For positive velocity, the characteristic map is represented through the
travel-time coordinate

[
T(x)=\int_0^x \frac{d\xi}{u(\xi)}.
]

The departure face satisfies

[
T(x_d)=T(x_f)-\Delta t.
]

The current implementation assumes positive velocity, which is consistent with
the reference pressure ordering (P_{in}>P_{out}). Variable velocity is
represented as piecewise linear between face velocities, so travel-time
inversion can cross many cells without CFL subcycling.

## Source coupling

For the first coupled prototype, use Strang splitting:

[
S_{\Delta t/2}\;T_{\Delta t}\;S_{\Delta t/2}.
]

The source operator is integrated analytically for frozen (k), preserving the
same constitutive law as v1.

## What has been demonstrated

The local verification suite checks:

1. pure periodic advection at CFL up to 50 with non-integer total displacement;
2. open-boundary transport with Dirichlet inflow and exact mass balance;
3. variable positive velocity with conservation at CFL about 50;
4. advection plus relaxation against a closed-form solution.

The prototype is **not yet the production erosion solver**. Before integration into
`src/simulation.py`, it needs a coupled A/B comparison against the frozen v1
solver and a timestep/accuracy study using the actual erosion model.
