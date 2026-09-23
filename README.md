# 1D Internal Erosion Model

This repository contains the source code for a one-dimensional (1D) internal
erosion simulation. It includes the numerical algorithms, the physical model,
and the scripts required to generate the results and figures associated with
the study.

The reference implementation transports the conservative suspended-solid
volume variable

\[
S=A\phi,
\qquad A=\pi R^2,
\]

and uses the MinMod MUSCL limiter, the active C1-smoothed Julien multiplier,
and physical-time erosion without a morphological acceleration factor.

For the scientific context and the project objectives, see
[PROJECT_FRAMEWORK.md](PROJECT_FRAMEWORK.md). The current physical and
numerical formulations are described in:

- [docs/physical_model.md](docs/physical_model.md)
- [docs/numerical_methods.md](docs/numerical_methods.md)

## Project Structure

- `src/`: physical laws, numerical operators, simulation orchestration,
  validation and post-processing utilities;
- `data/`: simulation configurations;
- `docs/`: scientific and numerical documentation;
- `main.py`: official entry point;
- `run_revision.py`: reference campaign runner;
- `run_figures_data.py`: figure-data campaign runner.

## Prerequisites

Python 3.8+ and the scientific Python stack are required:

- NumPy;
- Matplotlib and Pandas for the legacy figure generator;
- Pytest for the test suite.

Install the runtime dependencies from `requirements.txt` and the development
dependencies from `requirements-dev.txt` when available.

## Execution

Display the available reference cases:

```bash
python main.py --help
```

Run one reference case:

```bash
python main.py --case SolA_Kout10
```

Run the complete reference campaign:

```bash
python main.py
```

The runner writes each case into `output_revision/` and records the effective
transport variable, limiter, friction model, pressure residuals and validation
metrics in the run manifest.

## Performance profiling

The hydraulic pressure residual is accelerated with optional Numba kernels.
Profile a representative case with:

```bash
python tools/profile_simulation.py --nx 100 --t-end-factor 0.001
```

See [docs/performance.md](docs/performance.md) for the profiling results and
the warm-up interpretation of the Numba timings.
