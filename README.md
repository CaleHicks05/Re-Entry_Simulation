# Re-Entry_Simulation

## Atmospheric Re-entry Simulation (MATLAB)
This project simulates atmospheric re-entry of a spacecraft using MATLAB, accounting for drag variation with altitude, parachute deployment (drogue and main), convective heating, and radiative cooling. The goal was to model realistic descent profiles and thermal loads for learning purposes.

---

## Features
- Variable air density based on altitude
- Two-stage parachute deployment with adjustable areas and deployment times
- Heat shield modeling with drag and ambient temperature as inputs
- G-force tracking to analyze loads during descent
- Graph outputs for velocity, altitude, heat shield temperature, and acceleration profiles

---

## How it works
- **Variable inputs**: mass, cross-sectional area, velocity, drag coefficients, parachute deployment times, altitude, etc.
- **Simulation Process**: code runs through initial conditions, arrays, and time setup. Then runs through main loop for desired time and given conditions. Finally, produces graphs based on arrays produced.
- **Outputs**: Graphs for altitude, velocity, G-force, and temperature. 

---

## Future Improvements
- Heat shield 3D-temperature modeling
- Mach effects on drag
- Deployable landing legs / final touchdown velocity modeling

---

## How to Run
**For MATLAB:**
```matlab
% Open the .m file in MATLAB and run:
reentry_simulation
