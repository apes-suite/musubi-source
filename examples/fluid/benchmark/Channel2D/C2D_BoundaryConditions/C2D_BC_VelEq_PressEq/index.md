title: Channel 2D VelEq PressEq
@warning WORK IN PROGRESS @endwarning

Test case description:
Simple channel flow using velocity equilibrium inlet boundary and pressure equilibrium outlet boundary.

using symmetry at y-Axis. So computing only the lower part of the channel.
The upper part of the channel is symmetric to the lower part.

The templates directory contains a shepherd config file and all necessary template files to run convergence studies for mesh refinement and omega.
For numeric stability the Reynolds number has been decreased to 10  to ensure the simulations stay below limits of LBM method.
