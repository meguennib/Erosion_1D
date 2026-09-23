# Physical and Morphological Dynamics Model

This document outlines the physical principles implemented in the 1D internal erosion simulation code, and details the dynamic processes leading to asymmetric geometric profiles (the "trumpet" effect).

## 1. Physical Context of Internal Erosion
The model simulates water flow through an erodible conduit. The geometry of the conduit evolves over time due to the shear stress exerted by the fluid on the walls.
Two main dynamics are strongly coupled:
1. **Hydrodynamics:** Calculation of the fluid velocity and pressure along the conduit for a given pressure gradient.
2. **Sediment Transport:** The erosion of the walls increases the suspended sediment concentration ($\phi$) within the flow.

## 2. Rheological Feedback
A key aspect of the model is the feedback of the sediment concentration on the fluid rheology.
As water flows and erodes the walls, the suspended sediment concentration increases from upstream to downstream. This denser suspension modifies the effective fluid friction against the wall.

A friction multiplier, $f_m(x)$, is utilized to model this additional resistance.
*   **Upstream:** Clear water enters with zero sediment concentration ($f_m \approx 1.0$).
*   **Downstream:** The accumulated sediment concentration generates collisional dissipation within the boundary layer, increasing the friction ($f_m > 1.0$).

## 3. Mechanism of the "Trumpet" Effect (Downstream Widening)
The rheological feedback is the primary driver of the observed geometric asymmetry (the "trumpet" profile). The mechanism is detailed as follows:

1. **Spatial Accumulation:** The concentration $\phi$ increases from upstream to downstream.
2. **Friction Gradient:** The rheological multiplier increases downstream, raising the flow resistance.
3. **Differential Stress:** The heightened friction generates a significantly higher local shear stress on the downstream walls compared to the upstream section.
4. **Differential Erosion:** The conduit erodes more rapidly downstream, resulting in the flared geometric profile.
5. **Regulation:** Once the downstream conduit is sufficiently widened, the fluid velocity drops locally, diluting the concentration and progressively reducing the erosion rate, thereby stabilizing the overall shape.

## 4. Importance of the Spatial Scale
The emergence of this effect is highly dependent on the fluid's **residence time** within the conduit.
*   **Laboratory Scale (e.g., 10 cm):** Water traverses the conduit too rapidly to accumulate a significant sediment concentration. The fluid behaves essentially as clear water over the entire length, leading to perfectly cylindrical erosion.
*   **Field Scale (e.g., 10 m):** The water has sufficient time to accumulate a substantial sediment load along its path, naturally triggering the rheological feedback loop without relying on unphysical morphological acceleration factors.

## 5. Influence of Particle Size
The rheological parameter (based on Julien's model) is proportional to the square of the particle diameter ($d_p$).
*   **Very Fine Particles (e.g., Silt/Clay, 20 µm):** The rheological effect remains negligible. The conduit erodes uniformly into a cylinder.
*   **Coarse Particles (e.g., Fine Sand, 200 µm):** Larger grains generate significant dissipation. The feedback is strongly activated, creating the trumpet geometry.

The code is architected to handle these varying spatial and granulometric scales in a completely physical and stable manner, enforced by strict time-step control (CFL).
