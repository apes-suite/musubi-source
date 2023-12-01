Definition of Initial Conditions {#tut_6initial}
========

Navigate: \ref tut_5restart "&larr; Restarting Simulations"
| \ref mus_tutorials_overview "Overview"
| \ref tut_7convergence "Setting up a Convergence Criterium &rarr;"

In this tutorial, we cover the definition of initial conditions. 
They can be defined either in lattice units or in physical units.
Lattice units can be confusing, especially to people just starting in this field.
The relationships between lattice and physical units can be found here.

> [Here](\ref "wiki.palabos.org/_media/howtos:lbunits.pdf"), you can have a look at a paper about unit conversion by Jonas Latt.

In this tutorial, a two-dimensional plain channel is set up. 
Not only are the boundaries specified to obtain a defined pressure drop over the channel length, 
but also are the initial conditions set in a consistent manner.

For the viscous, laminar two-dimensional plain channel flow, an analytical solution
of the incompressible Navier-Stokes equation can be derived. From the analytical solution
the pressure drop, the velocity profile and the shear stress distribution can be computed.

Before starting, we need to define the flow regime and physical reference values.

\snippet testcases/channel/musubi.lua reference physical values

For the lattice Boltzmann simulation, basic simulation parameters such as a characteristic velocity, 
a resolution and a relaxation parameter omega need to be specified.

\snippet testcases/channel/musubi.lua reference LB values
Depending on the model used, the reference pressure differs. For the incompressible model, 
the reference pressure is 0, while for the compressible model the reference pressure is rho0*cs2
\snippet testcases/channel/musubi.lua reference pressure

From the solution of the Navier-Stokes equation, the following relations for the
velocity distribution across the channel height can be obtained
\snippet testcases/channel/musubi.lua velocity function
Similarly for the pressure drop along the channel length
\snippet testcases/channel/musubi.lua pressure function
and the shear stress across the channel height
\snippet testcases/channel/musubi.lua shear stress function

Now the physics table establishes the connection between the lattice reference values and the 
physical values and gives Musubi means of transferring between these two unit systems
\snippet testcases/channel/musubi.lua physics table

For the Lattice Boltzmann algorithm, a reference density and the relaxation rate omega need to be defined
\snippet testcases/channel/musubi.lua fluid table

Now the initial conditions for each element in the simulation domain is defined by setting
each physical quantity and connecting it to a lua function, which we defined above.
\snippet testcases/channel/musubi.lua initial conditions

There are still the in- and outlet, which also have to be given consistent conditions.
Here, we make use of pressure boundaries, as we know the analytical solution for the pressure drop along
the channel length.
The inlet is impingned with the pressure at the inlet position, which is in our definition minus one half of 
the total channel length -length*0.5. Similiarly the pressure at the outlet is defined by the pressure function
evaluated at the outlet position at +length*0.5
\snippet testcases/channel/musubi.lua boundary conditions

Next chapter: \ref tut_7convergence "Set up a Convergence Criterium &rarr;"
