# Performance profiling and acceleration

## Profiling method

The reproducible profiler entry point is:

```bash
python tools/profile_simulation.py --nx 100 --t-end-factor 0.001
```

A binary profile can be saved with:

```bash
python tools/profile_simulation.py --stats profile/erosion.prof
```

The production campaign should be profiled with a short representative case
first. Profiling the complete paper campaign directly makes it difficult to
distinguish startup, I/O and numerical costs.

## Baseline bottleneck

A `cProfile` run on a representative `Nx=100` case with 109 time steps showed
that the hydraulic pressure solve dominated the execution:

- `solve_Q_pressure_imposed`: approximately 89% of cumulative runtime;
- its scalar `evaluate` function: approximately 82%;
- the repeated Barenblatt friction and beta evaluations inside the bisection;
- the advection and morphology update were comparatively small contributors.

The root solver was evaluating complete temporary arrays at every bisection
point even though only the scalar pressure residual was needed until the
accepted discharge had been found.

## Implemented acceleration

`src/kernels.py` contains two isolated hot kernels:

1. `pressure_residual_kernel`: computes only the scalar pressure residual for
   a trial discharge;
2. `hydraulic_fields_kernel`: reconstructs `tau_b`, `mdot`, `f_D`, `beta` and
   `u` once for the accepted discharge.

When Numba is installed these kernels are compiled with `njit(cache=True)`.
The code retains a Python fallback so that the scientific solver remains
usable in environments without Numba.

The outer bisection, state transitions, logging and I/O remain in Python. This
keeps the acceleration isolated from the physical model and avoids compiling
objects such as `Params`, logger instances or dictionaries.

## Measured result

On the profiling environment used during development, the representative
case took approximately:

| Configuration | Runtime |
|---|---:|
| Original vectorised residual with field construction at every trial | 0.225 s |
| Accelerated kernel, first process including JIT compilation | 0.369 s |
| Accelerated kernel, warmed process | 0.070 s |

The first Numba call includes compilation and is therefore not a fair measure
of a long campaign. For the large reference campaign, the compilation cost is
amortised over thousands of hydraulic evaluations. The exact speed-up remains
hardware- and problem-size-dependent.

The optimized kernels were checked against the original NumPy constitutive
expressions to machine precision for the residual and all reconstructed
hydraulic fields.

## Scope and constraints

No parallel execution, morphological acceleration, change to the physical
model, or change to the pressure root criterion was introduced. The intended
next performance checks are:

- benchmark at the actual paper resolutions (`Nx=300`, `600`, `1200`);
- measure cold-start and warmed-process times separately;
- verify that the pressure residual and final fields remain unchanged within
  numerical round-off;
- profile I/O only after the hydraulic kernel cost has been measured.
