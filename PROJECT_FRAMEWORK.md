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
The current reference implementation addresses these challenges through the following physical and numerical choices:

*   **Conservative transport:** The transported state is the suspended-solid volume per unit axial length, $S=A\phi$, rather than concentration alone. This preserves the solid-volume balance as the conduit area evolves.
*   **Numerical stabilization:** The reference solver uses a MUSCL finite-volume reconstruction with a MinMod flux limiter for the conservative transport equation. The complete coupled method is explicitly documented as operator-split in time.
*   **Strict physical time-stepping:** The radius is advanced directly with $R_t=\dot m/\rho_{soil,sat}$ using a dynamic advective/source constraint. No artificial morphological acceleration factor is used.
*   **Continuous rheological coupling:** The Julien multiplier is regularised with the documented hyperbolic-tangent law and actively bounded by $f_{m,max}$ during the hydraulic solve.
*   **Trumpet geometry:** Downstream widening emerges from conservative solid transport, concentration-dependent rheology, and the dependence of wall stress on the evolving radius. Particle diameter and spatial scale control the strength of this feedback.
