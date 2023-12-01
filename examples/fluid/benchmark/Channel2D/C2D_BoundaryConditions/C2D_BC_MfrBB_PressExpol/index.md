title: Channel 2D MfrBB PressExpol
@warning WORK IN PROGRESS @endwarning

Test case description:
Simple channel flow using mass flow rate bounce back as inlet and pressure extrapolation as outlet boundary.
Be aware: Mass flow rate is given as scalar value, so the inflow condition does not match the analytical solution.

using symmetry at y-Axis. So computing only the lower part of the channel.
The upper part of the channel is symmetric to the lower part.

The templates directory contains a shepherd config file and all necessary template files to run convergence studies for mesh refinement and omega.
For numeric stability the Reynolds number has been decreased to 10  to ensure the simulations stay below limits of LBM method.
