# Numerical Methods and Algorithms

This document details the numerical algorithms utilized to simulate the 1D internal erosion processes.

## 1. 1D Finite Volume Solver (Hydrodynamics & Transport)
The core of the simulation relies on a 1D transient Finite Volume Method (FVM) to solve the coupled fluid and sediment transport equations.
*   **Grid:** A structured, cell-centered 1D grid is employed (e.g., `Nx = 400` uniform cells representing the longitudinal domain).
*   **Variable Storage:** Conservative variables (cross-sectional area $A$, sediment volume fraction $\phi$) are stored at the cell centers.
*   **Flux Computation:** The Monotonic Upstream-centered Scheme for Conservation Laws (MUSCL) is utilized to ensure 2nd-order spatial accuracy.
*   **Limiter:** A `MinMod` flux limiter prevents unphysical oscillations near sharp gradients, thereby guaranteeing a Total Variation Diminishing (TVD) scheme.

## 2. Dynamic Time-Stepping (CFL Condition)
To maintain strict numerical stability, the time step ($\Delta t$) is dynamically calculated at every iteration.
*   **Advective CFL:** $\Delta t_{adv} = \text{CFL} \times \frac{\Delta x}{U_{max}}$, where $U_{max}$ is the maximum fluid velocity.
*   **Erosive CFL:** Limits the radial change per time step: $\Delta t_{er} = 0.01 \times \frac{R}{\dot{R}}$.
*   **Global Time Step:** $\Delta t = \min(\Delta t_{adv}, \Delta t_{er})$. The CFL number is strictly maintained below the theoretical stability limit of $1.0$ (typically around $0.45$).

## 3. Hydraulic Algorithm
At each time step, the flow rate $Q(t)$ is determined such that the total pressure drop across the conduit matches the prescribed boundary conditions ($P_{in} - P_{out}$).
*   **Algorithm:** Brent's method (`scipy.optimize.brentq`) is used to find the root of the function $\Delta P_{calc}(Q) - \Delta P_{target} = 0$.
*   **Pressure Integration:** For a given $Q$, the pressure is numerically integrated backward from the outlet to the inlet using a 4th-order Runge-Kutta (RK4) or Euler integration scheme. This integration accounts for fluid acceleration (Bernoulli) and continuous friction losses (Darcy-Weisbach).

## 4. Rheological Smoothing
To avoid numerical shocks that could disrupt the spatial TVD scheme, the abrupt capping of the friction multiplier $f_m$ has been replaced with a $C^1$-continuous hyperbolic tangent smoothing function:
$$f_{m, smooth} = 1.0 + (f_{m, max} - 1.0) \times \tanh\left( \frac{f_{m, raw} - 1.0}{f_{m, max} - 1.0} \right)$$
This ensures a smooth evolution of the rheological parameters over space and time.
