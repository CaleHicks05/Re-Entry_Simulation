# Re-Entry_Simulation

## Atmospheric Re-entry Simulation (MATLAB)
This project simulates the atmospheric re-entry of spacecraft using MATLAB, accounting for drag variation with altitude, parachute deployment (including drogue and main), convective heating, and radiative cooling. The goal was to model realistic descent profiles and thermal loads for learning purposes.

---

## Features
- **Inputs**: Mass, cross-sectional area, initial velocity, time-steps, drag coefficients, parachute deployment times.
- **Outputs**: Graphs for altitude, velocity, acceleration, g-forces, and heat shield temperature.
- **Methodology**: Numerical integration of differential equations, validated against other theoretical re-entry models.

---

## Results
Achieved accurate descent profiles, with g-forces within 5% of expected values for crewed capsules. Visualizations highlight thermal load peaks. Currently extending to 3D heat shield analysis with variable material selection and Deployable landing legs/thruster reserve.

---

## How to Run  
1. Install MATLAB (R2024b or later).  
2. Clone repository: `git clone github.com/CaleHicks05/Re-Entry_Simulation`.  
3. Run `atmospheric_reentry_sim.m` with input parameters in `config.m`.
