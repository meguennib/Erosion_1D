# Project Framework and Objectives

This document outlines the scientific framework, the initial goals, and the major milestones achieved during the development of this 1D internal erosion model.

## Publication Details
This codebase was developed in the context of the following scientific publication:

*   **Title:** [A coupled one-dimensional model for piping erosion under constant pressure drop]
*   **Authors:** [meguenni & all]
*   **Journal:** International Journal of Geotechnical Engineering (IJGE)
*   **Status:** [DOI - 10.1080/19386362.2026.2737123]

*If you use this code in your research, please cite the article above.*

## 1. Initial Objectives
The primary goal of this research project was to develop a robust, physics-based one-dimensional (1D) model capable of simulating internal erosion within a soil conduit. The specific initial objectives were:
*   To strictly couple the hydrodynamics (fluid flow, pressure gradients) with sediment transport and morphological evolution (erosion of the conduit walls).
*   To accurately model the rheological feedback: how the concentration of eroded particles suspended in the fluid increases the flow resistance.
*   To investigate the physical mechanisms leading to non-uniform geometric evolution, specifically the downstream widening commonly referred to as the "trumpet" effect.

## 2. Challenges Encountered
During the development phase, the explicit coupling of the highly non-linear erosion source terms with the advective transport equations introduced severe numerical stiffness. Initial iterations suffered from high-frequency spatial oscillations ("checkerboarding") and required artificial morphological acceleration factors (`morfac`) to observe significant erosion at laboratory scales, which masked the true physical time scales.

## 3. Achieved Objectives and Milestones
The current version of the code successfully overcomes these challenges, achieving all initial objectives through rigorous physical and numerical stabilization:

*   **Numerical Stabilization:** The model was successfully stabilized using a 2nd-order Finite Volume Method with a MUSCL scheme and a MinMod flux limiter. This guarantees a Total Variation Diminishing (TVD) property, eliminating unphysical oscillations.
*   **Strict Physical Time-Stepping:** We implemented a dynamic, dual-constrained CFL condition (advective and erosive). The simulation now runs entirely in real physical time, entirely eliminating the need for artificial morphological acceleration.
*   **Continuous Rheological Coupling:** The implementation of a $C^1$-continuous smoothing function for the Julien friction multiplier prevents numerical shocks while accurately representing the increased dissipation caused by the sediment load.
*   **Recovery of the "Trumpet" Geometry:** We successfully demonstrated that the downstream widening of the pipe is a natural consequence of the rheological feedback. The model proves that this geometric asymmetry is directly triggered by the physical properties of the sediment (specifically the particle diameter) and the spatial scale of the domain (laboratory vs. field scale).
