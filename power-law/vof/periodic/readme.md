Periodic power-law roll waves.

The most updated one with consistent gravity is stored here: *[link](https://github.com/MGYBY/Basilisk_practice/tree/main/power-law/herschel-bulkley/periodic)*.

* Validated 2D D2 calculation.
* Controllable top location.
* Serial or parallel run.
* Gravity body-force & surface tension.
* Droplet & bubble removal (for robustness).
* Tuned adaptivity criteria.
* Customized multiple (>1) filters for VOF tracer, see [code here](https://github.com/MGYBY/Basilisk_practice/blob/main/power-law/vof/periodic/two-phasePL_multipleFilters.h).
* More reasonable initialization. Solid placement. Avoid refine the solid part.
* Permanent refinement for the boundary layer. (Quite necessary for guarantee the problem is physical, and the convergence of the solver.)
* Inclusion (in a reasonable way) of surface tension: follow [code here](https://github.com/MGYBY/Basilisk_practice/tree/main/power-law/vof/normalFlow/fixedSurfaceTension).
