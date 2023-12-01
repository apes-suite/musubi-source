title: Channel 2D VelNonEqExpol PressNonEqExpol
@warning WORK IN PROGRESS @endwarning

## Channel 2D - Non-Equilibrium Extrapolation BCs ## {#ex_C2D_BC_VelNonEqExpol_PressNonEqExpol}

### Test Case Description ###
Simple channel flow driven by pressure drop across the channel.
Using velocity nonequilibrium extrapolation boundary as inlet and pressure non-
equilibrium extrapolation as outlet.
In order to run the test case with the same settings as the other C2D_BC cases,
double the number of elements in height  (nHeight = 40) in seeder.lua to use 40
instead of 20 elements and thus a finer resolution.

### Plots ###
For the nightly regression check, the tracking is mostly deactivated. Please
make sure to remove the comments to activate it for usual runs.

### Templates for Advanced Studies ###
The templates directory contains a shepherd config file and all necessary
template files to run convergence studies for mesh refinement and omega.
For numerical stability the Reynolds number has been decreased to 10 (instead of
100) to ensure that the simulations stay below limits of the LBM method.
