# Experimental SLFV solver

This branch keeps the thesis v1 physical model unchanged and replaces only the concentration transport operator by a conservative semi-Lagrangian finite-volume (SLFV) remapping.

## Transport

For the transport substep,

\( \partial_t\phi+\partial_x(u\phi)=0 \),

the velocity is frozen over one step and characteristics are traced backward:

\( dX/dt=u(X) \).

The cell average is reconstructed from the swept departure interval:

\( \phi_i^*=(1/\Delta x)\int_{X^-_i}^{X^+_i}\phi^n(x)\,dx \).

MUSCL reconstruction with Van-Leer limiting is used. The present implementation targets the reference positive-flow configuration \(u>0\).

## Why the CFL restriction is lifted

The Eulerian v1 transport scheme requires \(\Delta t_{adv}=CFL\,\Delta x/\max|u|\). The SLFV operator does not use this stability restriction. Its characteristics may cross several cells during one time step.

Therefore, on the SLFV branch, the transport CFL is monitored but is not an active time-step limiter.

## Remaining temporal constraints

The erosion source is integrated analytically:

\( \phi^{n+1}=\phi_{soil}-(\phi_{soil}-\phi^*)e^{-k^n\Delta t} \).

Consequently \(\Delta t_{src}=0.5/k_{max}\) is interpreted as an accuracy criterion for operator splitting, not as a source stability condition.

The radius equation remains explicitly integrated:

\( R^{n+1}=R^n+\Delta t\,\dot R^n \).

The SLFV branch therefore uses the morphological accuracy constraint

\( \Delta t_{morph}=\eta_R\min_i(R_i/|\dot R_i|) \),

implemented with the actual slope factor in \(\dot R\), and controlled by the parameter `morph_rel_change` (default 0.01).

The active SLFV step is

\( \Delta t=\min(\Delta t_{src},\Delta t_{morph},\Delta t_{max},t_{end}-t) \),

while the v1 Eulerian branch retains \(\Delta t=\min(\Delta t_{adv},\Delta t_{src},\Delta t_{max},t_{end}-t) \).

This separation makes it possible to quantify which physical process controls computational cost after the CFL restriction has been removed.

## Required verification sequence

1. Compare v1 and SLFV at identical time steps.
2. Verify mass conservation and boundedness of \(\phi\).
3. Run SLFV with the adaptive policy above.
4. Record \(\Delta t_{adv}\), \(\Delta t_{src}\), \(\Delta t_{morph}\) and the active time-step constraint.
5. Repeat with `morph_rel_change` = 0.02, 0.01 and 0.005.
6. Compare \(Q\), \(R_{out}\), \(\phi\), \(\tau_b\), \(\dot m\) and the 1-D validity time against the reference solver.
7. Only after this convergence study should the SLFV strategy be considered for the thesis reference implementation.
