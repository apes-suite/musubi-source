This directory contains simulations on the trajectories of particles in the
steady flow over a backward-facing step. The setup is identical to that in
earlier numerical research by:

[1] I. E. Barton, "Computation of particle tracks over a backward-facing step,"
Journal of Aerosol Science, vol. 26, no. 6, pp. 887–901, Sep. 1995, doi:
10.1016/0021-8502(95)00018-8.

and the experimental work by

[2] B. F. Armaly, F. Durst, J. C. F. Pereira, and B. Sch, "Experimental and
theoretical investigation of backward-facing step flow," p. 25.

Four cases were simulated with a difference in Stokes number and fluid-particle
coupling (one or two-way) in our publication:

[3] Vlogman, T. G., & Jain, K. (2024). Efficient coupled lattice Boltzmann and
Discrete Element Method simulations of small particles in complex geometries.
Computers & Mathematics with Applications, 175, 313-329.

This testcase simulates one of the four cases.
The Stokes number is set to 0.01 and two way coupling has been employed. 
Particle trajectories are written from the simulations as .dat files whereas
the flow field can be obtained by tracking .vtu files.

The exact setup for simulations is given in the .lua scripts.
1. It is recommended to first execute musubi_init.lua that simulates the flow
field and writes a restart file.
Particle computations are expensive and they can be initiated after a pure flow
simulation.
2. Then musubi_particles.lua can be executed that incorporates particles in the
already simulated flow field by reading the restart.
