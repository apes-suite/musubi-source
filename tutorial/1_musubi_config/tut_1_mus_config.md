Set up your first simulation with Musubi {#tut_1mus_config}
========

Navigate: \ref tut_0prerequisites "&larr; Download, build and run Musubi"
| \ref mus_tutorials_overview "Overview"
| \ref tut_2toolchain "Using the full toolchain &rarr;"

## The musubi.lua file ##

As a user, your best way to communicate with Musubi is a *Lua*-Script that
tells Musubi what you want to simulate and which options you want to choose.
The Lua-syntax is fairly easy, and for this tutorial, we will use only a
small subset of it, so don't be afraid if you hear of Lua for the first time.

By default, Musubi looks for a file called `musubi.lua` inside your current
working directory. If you want, you can use different names (for example if
you want to keep more than one script) and pass the name of the script you
want the application to use as a command-line argument.

We are now going to build a working, almost-minimal `musubi.lua` file
step-by-step.
Please navigate to your `musubi` directory, and create a new file called
`musubi.lua` with your favorite text editor.

\note Most configuration options have some sane defaults, which allows you to
leave them alone, if they are not relevant to you.

\note Before you run the following simulation setup, make sure that you have
created the directories `restart` and `tracking` in the same path where 
`musubi.lua` is located.

## General simulation name ##

Moreover, every simulation should have a name,
which is specified by a variable called `simulation_name`.

In this example, we will simulate a
[gaussian pulse](http://en.wikipedia.org/wiki/Gaussian_function)
in pressure, hence the name *Gausspulse* seems adequate.
You can use just anything as a simulation name, just make it clear and
descriptive, as it will help you to identify produced results from the
simulation.

\note The simulation name will appear in all output files, so if you use good
and unique names, you always know where a particular output came from.

\snippet testcases/gaussian_pulse/musubi.lua general settings

## Physics table ##

Before we go to the next step, lets learn about the unit system used in 
Musubi to define simulation parameters: `kinematic_viscosity` in `fluid` table 
and macroscopic variables: `pressure`, `velocity` in `initial_condition` and 
`boundary_condition` tables. These parameters/variables can be defined either in 
physical SI units or lattice (non-dimensional) units but not mixed i.e. some in 
physical and other in lattice. It is important that they are consistent. 
Base SI units used to define the physical variabes are: kg (mass), m (length),
s (time), mol (amount of substance), K (thermodynamic temperature) and A (electric 
current). In Musubi, the following reference values 
are used used to convert parameters/variables in physical units to lattice
units and vice-versa: 

| Measure      | Unit     | Reference          | Notation. | Name in physics table | Default value |
|:------------:|:--------:|:------------------:|:---------:|:---------------------:|:-------------:|
| Mass         | $kg$     | fluid mass         | $m_{ref}$    | `mass0`            | -             |
| Density      | $kg/m^3$ | mean fluid density | $\rho_{ref}$ | `rho0`             | $mass0/dx^3$   |
| Length       | $m$      | coarsest element   | $dx_{ref}$   | Not required | is taken from mesh  |
| Time         | $s$      | coarsest time step size | $dt_{ref}$   | `dt`             | -                     |
| Amount of substance | $mol$ | Inverse of Avagadro's number  | $N_A$ | `mole0`       | $1/6.02214129e^{-23}$ |
| Mole density | $mol/m^3$ | mean mole density | $c_{ref}$    | `moleDens0`           | $1/(N_A dx^3)$        |
| Molecular weight | $kg/mol$ | Largest molecular weight of species| $M_{ref}$ | `molWeight0` | $mass0/mole0$ |
| Temperature  | $K$      | fluid temperature  | $T_{ref}$    | `temp0`               | 1.0                   |
| Electric charge | $C$    | fundamental electric charge | $q_{ref}$ | `coulomb0`     | $1.60217657e^{-19}$   |

As shown in the table some of them have default values and coarsest element size 
is taken from the mesh. So, the only two reference units that are mandatory to 
define in the `physics` table are reference for $kg$ and $s$. Most common approach
to define $kg$ is using reference density $\rho_{ref}$: $m_{ref}=\rho_{ref} m^3$.
Here is an example `physics` table used in a channel flow testcase.
\snippet testcases/channel/musubi.lua physics table

The non-dimensional lattice units for the base SI units are computed
with physical values = reference values so they are usually 1.0.
such that
- Lattice density, $\rho_{l} = \frac{\rho_{p}}{\rho_{ref}} = 1.0$
- Lattice element size, $dx_{l} = \frac{dx_p}{dx_{ref}} = 1.0$
- Lattice time step, $dt_{l} = \frac{dt_p}{dt_{ref}} = 1.0$ and so on.

All other derive variable in lattice units are converted to physical SI units
using reference variables in base SI units. 
For simplicity, underscore 'ref' are dropped in the equations below for the reference values.
- Lattice pressure $p_l$ to physical pressure, $p_p = p_l * \frac{\rho dx^2}{dt^2}$ $[kg/m/s^2]$
- Lattice velocity $u_l$ to physical velocity, $u_p = u_l * \frac{dx}{dt}$ $[m/s]$
- Lattice kinematic viscosity to physical kinematic viscosity, $\nu_p = \nu_l * \frac{dx^2}{dt}$ $[m^2/s]$
- Lattice force $f_l$ to physical force, $f_p = f_l * \frac{\rho dx^4}{dt^2}$ $[kg m/s^2]$
- Lattice current $I_l$ to physical current, $I_p = I_l * \frac{q}{dt}$ $[C/s]$
- Lattice potential $\phi_l$ to physical potential, $\phi_p = \phi_l * \frac{\rho dx^5}{q dt^2}$ $[kg m^2 /(C s^2)=V]$

For more unit conversion for other derive variables refer to source/mus_physics_table.f90

Now you know how to define variables in `fluid`/`species`/`mixture` tables and 
macroscopic variables in `initial_conditon` and `boundary_condition`.
Remember that the variables must be defined in consistent units.
Here is an example which defines all variables in physical units
\snippet testcases/channel/musubi.lua reference physical values

Finally, the physical time step size required in the `physics` table is computed
either by acoustic or diffusive scaling.
\snippet testcases/channel/musubi.lua time step determination

When physical units are not of interest and you would like to just run a simulation
by defining lattice units then just remove `physics` table from the musubi
configuration file. In this case, the variables are then defined in lattice units
\snippet testcases/channel/musubi.lua reference LB values 

## The fluid table ##

Next, we define the fluid properties, required for the LBM simulations.
As mentioned in previous session, they can be defined either in physical or
in lattice units. When defining in physical units, DO NOT FORGET to define
`physics` table. The properties to be defined in the fluid table depends
on the type of equations being solved. 

Here is the example for `fluid_incompressible` kind with `bgk` collision
\snippet testcases/channel/musubi.lua fluid table 
which only requires kinematic viscosity of the fluid. For a weakly-compressible 
`fluid` kind, `bulk_viscosity` must also be defined.
Both of them are inside the table `fluid` (Lua has a single neat concept
[tables](http://www.lua.org/pil/2.5.html) for representing multiple datas in
a structure. It allows you to treat variables a bit like trees: they can
contain more variables, other tables or even references to functions).
Tables are denoted by curly braces in Lua.

Here is a simple configuration file with `fluid` table in lattice units.
\snippet testcases/gaussian_pulse/musubi.lua general settings
With a single string variable describing the simulation name, and a table
containing two reals of names `kinematic_viscosity` and 
`bulk_viscosity` describing the fluid.
`viscLB` can be computed from the relaxation parameter of the collision operator
`omega` as
\snippet testcases/gaussian_pulse/musubi.lua general settings

## Time settings ##

Now we should specify the number of time steps for our simulation.
We do this by opening a `time` section with at least the variable `max`
given to specify the maximal number of time steps you want to simulate.
Additionally you might provide the `interval`-setting, that controls the
debug output: after `[interval]` time steps, the total density of the system
is calculated and written to the console, this defaults to 1.
Since we have conservation of mass, the total density should not change
much, so if it does, it is a hint for the user that the results of the
simulation are probably wrong. Setting the interval to something around
`[max/10]` is reasonable.
\snippet testcases/gaussian_pulse/musubi.lua time_control settings

## identification ##

It should be followed by a part that identifies the nature of the simulation.
To this part belong the following aspects: the `label`, the `layout`, the `kind` and the
`relaxation`. For now, it is not necessary to know these settings in detail.
You can choose between some values that you can look up in the 
[mus_scheme_header_module](@ref mus_scheme_header_module).
If these settings are not defined, musubi will use  default values.
As a first example, you can set up the identify table as:
\snippet testcases/gaussian_pulse/musubi.lua identify

## Geometry ##

The geometry is usually a more complicated thing to define. For anything
but the simplest case, we need *Seeder*, a tool that will be dealt with
in the next tutorial chapter.
For now, we will use a pre-defined geometry, which is just a simple cube
with periodic boundaries and no obstacles. We can specify the `length`
and position (`origin`) of that cube. The position of the mesh is only
relevant if you specify positions of other objects, too, such as
initial conditions or *Trackers* (we will explain what that is in the next
section).

Last but not least, the `refinementLevel` tells Musubi how fine the space
is discretized. A refinement level of `4` means that the space is cut into
\f$ 2^4 \f$ portions in each dimension, which leaves us with
\f$ 2^{3\cdot4} \f$ cells in `3` dimensions.
\snippet testcases/gaussian_pulse/musubi.lua geometry

\attention Lua is case-sensitve, and some Musubi
options (like `refinementLevel`) have to be written in camelCase.

## Tracking ##

Of course you do not only want to simulate, you also want to see results.
The `tracking` section is your way to control some subset output from Musubi.
The idea behind it is that for any real-life scale problem, you just
cannot save all the information of every point at every time step,
since the shear amount of data would be too much to handle. Thus, we
tell Musubi what we want to see (like pressure, velocity...), at what
time, and in what position.

Each output is handled by a Tracker. Let us define a Tracker now, step
by step.
First, a Tracker has a name (specified in `label`) that will appear in
every output file created by this tracker. Anything that is precise and
meaningful will do.
\snippet testcases/gaussian_pulse/musubi.lua tracking basics

Next, we should tell the Tracker what we are interested in. The `variable`
section will do just that for us. In this example, we will ask the tracker
to store density and velocity information for us.

The output will be stored on your hard drive, in the location you specify
in `folder`.
\attention You have to create the folder yourself
*before* you start the simulation. If it doesn't exist, Musubi will crash
with a Fortran runtime error.

\snippet testcases/gaussian_pulse/musubi.lua tracking vars

The most complicated part is the [shape](@ref tem_shape_module) variable.
It defines *where* inside
the simulation domain (which is a cube in our example) you want to observe the
variables. You can take samples at a point, along a line, everywhere on a plane,
or even everywhere inside a given box. These options will be discussed in more
depth in the next chapter.
For now, we are happy with a point at position `(1,1,1)`.
\snippet testcases/gaussian_pulse/musubi.lua tracking shape

The `format` option lets you choose the file format in which the data
should be stored. The options you have are `ascii` and `harvester`.
Ascii is only really usable for point trackers. For anything more
sophisticated (starting in the next chapter), we will use the tool
`Harvester` for post-processing the data.
\snippet testcases/gaussian_pulse/musubi.lua tracking format

Finally, we have to specify at which time steps we want to save data.
Inside the `time` section, you can specify the first (`min`) and last
(`max`) time step that is of interest to you. If `[max]` is negative,
every time step till the end will be considered. In `interval`, you can
choose if you want to save every time step (set it to `1` then) or only
every `[n]`-th timestep (set `interval` to `n` in this case).
\snippet testcases/gaussian_pulse/musubi.lua tracking time_control
Don't forget to close the `tracking` section with another curly brace.

## Restart ##

In the `restart` part of your musubi.lua file you can create so-called restart files
in order to save the results of your simulation. With these restart files you can create vtk files 
with Harvester to view them in Paraview. In addition to that, you can also read the restart files
to go on with your simulation beginning at the saved point written in the restart files.
In the "restart" section you are also able to use `time_control` to define the `min`, `max` and the
`interval` time, when Musubi writes a restart file to the `restart/` folder. 
Be sure that you set up useful definitions for that.
For example, you can set up the restart settings like this:
\note You have to create the `restart` folder within your simulation folder at first.

\snippet testcases/gaussian_pulse/musubi.lua restart settings
So far, these are only basic information about `restart`. If you want to know more about
the possibilities in that case \ref restart_usage "click here."

## Initial Conditions ##

Initial conditions can be specified in the `initial_condition` section.
First, let us define that the velocity at `t=0` should be `0` everywhere.

But if every point is equal, not much is going to happen. So for the
initial pressure, we will do something more fancy: We will define a function
`gausspulse` that will set the initial pressure such that it has a peak
in the middle of the domain, and decreases quickly towards the sides.
With this, we will create two waves running from the center towards both
sides. The `initial_condition`-section will point to this function for the
density:

\snippet testcases/gaussian_pulse/musubi.lua initial conditions
The function can be defined like this in the Lua script:
\snippet testcases/gaussian_pulse/musubi.lua ic function

\attention Define the function before the `initial_condition` table or you will receive an error whilst
running Musubi.

Moreover,define the p0 as global variable in the function.
So that it will return the function like `p0+amplitude*math.exp(-0.5/halfwidth^2)*
(x-originX)^2`:
\snippet testcases/gaussian_pulse/musubi.lua ref pressure 

\note A last short remark on variables: If you do something more complicated
in your `musubi.lua` script, you are of course free to use functions,
variables, loops, branches, whatever you want.
The only thing that counts is what is defined after running your script.
Be careful, however, when using your own variables, as you might set Musubi
options that you are not aware of.
It might be wise to use variables only inside your own functions,
or at least give them some prefix so you don't confuse your own variables
with Musubi options.

## Result ##

After doing all this, your file should look a lot like
`musubi/tutorial/testcases/gaussian_pulse/musubi.lua`. You should be able to run
Musubi now with your own configuration.
After that, you will find an output file called `Gausspulse__track_pressure_p00000.res` in your
`tracking`-folder. 
If you open it, you will find densities and velocities for each and every
time step.

\note The Ascii output is used only in this first tutorial. You can visualize
this data with any program you like, for example with
[Gnuplot](http://www.gnuplot.info), if you
have it installed. If not, wait for the next chapter, as we are going to
introduce different tools there, anyway.
Simply run Gnuplot with
\code
gnuplot -p -e "plot \"tracking/Gausspulse__track_pressure_p00000.res\" using 1:2"
\endcode
to get the following plot of the pressure: 
\image html images/pulse_point.png "Several pressure pulses hitting the point tracker"
The picture shows the pulse running through your point tracker several
times. This happens because of the periodic boundaries. The intensity of
the pulse decreases due to dissipation.

In the next chapter, we will learn to define more complex geometries, and
to create more sophisticated outputs.

Next chapter: \ref tut_2toolchain "Using the full toolchain &rarr;"
