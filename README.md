# 1D Internal Erosion Model

This repository contains the source code for a one-dimensional (1D) internal erosion simulation.
It includes the numerical algorithms, the physical model, and the scripts required to generate the results and figures associated with the study.

For details regarding the scientific context, publication details, and achieved objectives, please see the [Project Framework](PROJECT_FRAMEWORK.md) document.

## Project Structure

- `src/`: Contains the main source code (physical models, numerical methods, simulation utilities).
- `data/`: Contains the simulation configurations (e.g., `scenarios.json`).
- `docs/`: Contains technical documentation on the numerical algorithms and the physical model.
- `main.py`, `run_revision.py`, `run_figures_data.py`: Main execution scripts for the simulation.

## Prerequisites
- Python 3.8+
- Standard scientific libraries (NumPy, SciPy, Matplotlib)

## Execution
To launch the main simulation:
```bash
python main.py
```
