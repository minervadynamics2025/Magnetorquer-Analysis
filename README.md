# Magnetorquer-Analysis
Here we share primarily MATLAB codes developed for simulating magnetorquer magnetic fields and processing test data

## Revised analysis (manuscript revision R1)

Run `mainDipoleMoment.m`. It

1. solves the axisymmetric magnetostatic problem **with the ferromagnetic core included**
   (`solveAxisymmetricField.m`, finite-volume method for psi = r*A_theta) and computes the
   transfer factor kappa = <B_z>_rod / B_z(sensor) and the demagnetizing factor N_d
   (`getFieldMetrics.m`);
2. converts the teslameter readings to the volume-averaged internal flux density;
3. fits the Langevin anhysteretic model by simulated annealing + Levenberg-Marquardt
   (`optimizeSimulatedAnnealingLangevinModel.m`, `fitLMLangevin.m`, `langevinRodModel.m`):
   once with M_s free (effective parameters) and once with M_s fixed to the HyperCo 50A
   datasheet value, plus a constant-mu_r reference fit (`constantPermeabilityFit.m`);
4. quantifies the sensitivity of kappa and of the fitted mu_r to the sensor position (+-0.5 mm)
   and builds the uncertainty budget (fit scatter, probe calibration, probe length, winding position);
5. saves the four manuscript figures (`rev_FEMCenterline`, `rev_MagneticField`, `rev_Sensitivity`,
   `rev_RelativePermeability`; styled by `newFig.m`, `styleAxes.m`, `styleMarkers.m`, `saveFig.m`)
   and `results_macros.tex` into the manuscript folder.

**Edit the geometry block at the top of `mainDipoleMoment.m` (total rod length, coil radial
build, sensor position) before using the results.**

Why the change: about 99% of the field at the sensor comes from the core magnetization, whose
spatial distribution differs from that of the coil current. The former air-core conversion factor
(r_s0*r_volume = 25.54, from `getAmpFactor.m`, `Bz_volume_avg_full.m`, `Bz_analytical.m`)
is therefore biased and is no longer used; those files are kept only for reference.
`Bz_on_axis.m` now implements the correct finite-solenoid on-axis field.
`permeabilityModelFit.m` (exponential mu_r model) is no longer used in the manuscript.

A Python version of the same analysis (`revised_analysis.py`) is kept with the manuscript and
gives identical numbers (all macros agree).
