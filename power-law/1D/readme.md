One-dimensional SW simulations.
1. Modified saint-venant solver (ignore topography term and gravity variable G).
2. Modified HLLC Riemann solver for the new eigenstructure.
3. Simulation file.

TODO:
- [x] Kurganov's free-of-Riemann solver approach (central-upwind).
- [x] Central-upwind scheme with numerical dissipation reduction (Kurganov & Lin 2007).
  - [This link](https://github.com/MGYBY/Gerris_practice/blob/main/power-law/2D/damBreak_onSlope/explicit_scheme/riemann-power-law_sharp.h)
- [ ] Passive scalar transport equation.
