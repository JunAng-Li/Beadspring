# Beadspring
This folder include all the data and codes necessary to do the analysis in the paper and to reproduce all the figures in the paper 'Quantifying dissipation using fluctuating currents'.
It has three separate folders, code, data and figures.
The code folder includes codes to do the bead-spring simulation, Kernel estimation, calculate the TUR Lower bound and MC optimization.
-multiBeads_data.m
This code is used to generate the simulation trajectories for different number of beads at arbitrary temperature baths from the discrete-time scheme.
-twoBeads_Gillespie.nb
This code is used to generate the simulation trajectories for two beads at arbitrary temperature baths from the discrete-space scheme.
-multiThermodynamicFroce_CurFluc.m
This code is used to estimate the thermodynamic force along the trajectory using Kernel estimation and calculate the accumulated current which will be used for both the temporal estimator and TUR lower bound later.
-twoBeads_KDE_Spatial.m
This code is using kernel estimator to do the spatial integral.
-EstthermoForceCurFluc.nb
This code takes the files(accumulated current) generated from multiThermodynamicFroce_CurFluc.m to calculate the TUR lower bound.
-TemperatureOpt.nb
This code is for generating the optimized lower bound for two beads at different temperature ratios using the Monte Carlo random sampling from the thermodynamic force.
-OptimizationThermodunamicsFroce.nb
Doing Monte Carlo random sampling from the thermodynamic force and plot out the optimized configuration of the d vector.
-OptimizationRandom.nb
Doing Monte Carlo random sampling from a random d vector and plot out both the initial random d vector and the optimized configuration of the d vector.

The data folder includes the raw data to do analysis and all the data to generate the figures.

The figure folder includes all the codes to make the figures. Along with the data file, it can reproduce the figures directly.
