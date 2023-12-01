Steady Flow Around Cylinder 2D {#tut_steady2D}
=======

Navigate: \ref tut_FlowAroundCyl2D "&larr; Test case Flow Around Cylinder 2D"
| \ref tut_unsteady2D "Unsteady Flow Around Cylinder 2D &rarr;"

# Simulation of the steady flow

We would like to use the whole apes suite. This means that we are going through
a few steps:

* We create a mesh in Seeder, containing the channel and the obstacle (sphere
  or cylinder, which depends on the testcase).

* We visualize the mesh with Harvester and have a look at it in Paraview.

* We create a simulation file for Musubi and find out when the simulation gets 
  stable. Plot the results with Gleaner.
  \note We use Gleaner inside our group which is based upon the free available
  matplotlib. You can use

* If we have reached a steady status, we track the needed variable values in 
  order to get results.

* The results can be plotted with Gleaner and shall be sum up at the end.

> This test case is available in the repository in the 
> `tutorial/testcases/FlowAroundCylinder/`. 

### Conditions provided by the test case

In the 2D test case for a steady flow around a cylinder we use a Reynolds number
of \f$ Re=20 \f$ . We have a kinetic viscosity of 
\f$ \nu = 0.001\frac{m^{2}}{s} \f$ and a fluid density of 
\f$ \rho = 1\frac{kg}{m^3} \f$ .

Test begins here:
\latexonly \f$ \nu = 0.001\frac{m^{2}}{s} \f$ \endlatexonly and a fluid density of 
\latexonly \f$ \rho = 1\frac{kg}{m^3} \f$ \endlatexonly.
Test ends here.

The mean velocity is given with \f$ U_{m}=0.3\frac{m}{s} \f$ . 
The geometry information is given in the above image that displays a sphere 
inside of a channel.
Although the length is given by 2.2 m, we will run the simulation with different
ratios between length and height of the channel in order to have a look at the 
influence on the tracked variable values. We will call it `l_h` in the code.


We calculate the mean velocity \f$ \bar{U} = Re\frac{\nu}{D} \f$ and call it 
`u_mean_phy` and the incoming velocity 
\f$ U_{m} = \frac{3}{2}\bar{U} \f$ and call it `u_in_phy`.

We also define a variable to check the Reynolds number. In this testcase we 
can make use of either acoustic scaling where dt is proportional to \f$dx\f$ or
diffusive scaling where dt is proportional to \f$ dx^{2} \f$ . 
We define a boolian variable for this. The calculation of lattice units depends
on the use of different scaling methods. 
After that, we calculate the physical pressure which depends on the physical 
density in the solver.
The function for the incoming velocity of the flow is given in the paper 
"Benchmark Computations of Laminar Flow Around a Cylinder" by M. Schaefer and 
S. Turek for the steady test case with:

\f$ U(0,y,t)=\frac{4*U_{m}*y*(H-y)}{H^{2}}, V=0 \f$

Therefore we might start, creating a file containing all the information about 
the physical conditions of the fluid and the convertion between physical and 
lattice units and basics that can be used in both files, seeder.lua and 
musubi.lua. This is done in the following.

# common.lua, units.lua, printSimParam.lua

In order to separate things and to get a better overview of the simulation files
it is recommend to create a common.lua file containing all the information 
dealing with physical and lattice boltzmann units. Here, we can define variables
and functions that are needed for the whole package of seeder, harvester and 
musubi.

In this test case we have to set up variables concerning refinement, geometry, 
physics, fluid and scaling method.

We write down everything we know from the provided conditions.

\snippet testcases/FlowAroundCylinder/common.lua geometry

In the above code, nHeight and l_h defines number of elements along height and 
length to height ratio respectively.

In the steady flow simulation (Re = 20), l_h = 5 (fixed) and 
nHeight = 50,100,150,200,300,400. A higher `nHeight` refines the mesh so that 
the simulation will be more exact. 

Then we write down the flow parameters which are known and which have to be 
calculated.

\snippet testcases/FlowAroundCylinder/common.lua flow parameters

There is a difference between acoustic scaling and diffusive scaling. 
This is why we have to define two different cases.

\snippet testcases/FlowAroundCylinder/common.lua scaling

For acoustic scaling we fix the incoming lattice velocity and compute omega.

\snippet testcases/FlowAroundCylinder/common.lua acoustic scaling

And for diffusive scaling it is the other way round, we fix omega and compute 
the lattice incoming velocity.

\snippet testcases/FlowAroundCylinder/common.lua diffusive scaling

> To get better results, we could also change the Mach number which is a 
> division by the lattice velocity over the square root of the lattice stream of 
> sound.
> Therefore we have to change the omega value. `omega` leads to `nu_L`, `nu_L` 
> leads to `dt` and `dt` leads to `u_in_L`. So changing omega means changing 
> lattice velocity. 
> This is for diffusive scaling.

\snippet testcases/FlowAroundCylinder/common.lua mach

The initial pressure is dependent on the density and the computed values for dx
and dt.

\snippet testcases/FlowAroundCylinder/common.lua pressure

In the paper mentioned above, the velocity of the flow is given by this space 
time function that is the last line of common.lua:

\snippet testcases/FlowAroundCylinder/common.lua velocity input function

There is another file that we call `units.lua`. There we write down the lattice units:

\include testcases/FlowAroundCylinder/units.lua 

We also could write a lua script that prints all the necessary 
variable values for a documentation of the results later on.

We could call it `printSimParams.lua`
It just prints all the values to the screen. It looks like this:

\include testcases/FlowAroundCylinder/printSimParam.lua 

After finishing this additional file we can go on with the main files for seeder
and musubi. 

# Seeder
   
We use Seeder in order to generate the mesh for our simulation of the flow 
around a cylinder. 
The goal is to create a channel with a cylinder inside. The cylinder is regarded
as an obstacle and is placed vertically to the sidewalls of the channel. Because 
of the 2d case, we can use a sphere.

We have to build a seeder.lua file containing the following aspects:

* folder
* debug table
* bounding_cube table
* spatial_object table

These are the basic things we need in this file. Because of the fact, that we 
created a file common.lua with physical and technical information about the 
simulation we need to make sure that seeder.lua gets its information.

We do this with: 

\snippet testcases/FlowAroundCylinder/seeder.lua input

Normally, we make a directory called "mesh" for the seeder output. This is the 
location for the mesh which we define with:

\snippet testcases/FlowAroundCylinder/seeder.lua folder

For debugging cases we are able to define a debug table. We have to activate the
debugMode and the debugFiles.
Before we run seeder.lua later on, we make sure that we create a directory 
called 'debug' in the same folder where seeder.lua is located.

\snippet testcases/FlowAroundCylinder/seeder.lua debug

Now, we can define the bounding_cube table. The bounding_cube is a single 
element of the mesh. We need to define the origin and the length of a single 
edge. After that we have done this, we set up a variable named `minlevel` that 
we need later on in musubi.lua.

\snippet testcases/FlowAroundCylinder/seeder.lua bounding cube

Then we go on with the spatial_object table where we define the seed, the 
outside boundaries and the obstacle that is also a boundary.
A spatial_object is defined by its attribute table with kind and label and by 
its geometry table containing the kind and the object.

To have a look at the seed, it is important to place it inside of a fluid 
element. It cannot be placed onto a boundary. The seed is a simple point in the
fluid mesh where the flood begins to spread in every direction.
At first, we open the spatial_object table and go on with every attribute. For 
example, the seed is an attribute.

So we can define the seed like that:

\snippet testcases/FlowAroundCylinder/seeder.lua seed

After that, we define the four boundaries of the 2D channel. To have a good 
overview, we can give them labels like 'north', 'south', 'west' and 'east'. 
These boundaries are planes which are defined by the origin and two vectors.
Planes are also of kind 'canoND'. Each plane is one another attribute in the 
spatial_object table.

\snippet testcases/FlowAroundCylinder/seeder.lua boundaries

Because of the fact, that there are the same boundary conditions along the 
z-axis, we make use of periodic boundaries ('usePeriodic')which is defined in 
common.lua. This is the reason why we define an attribute which is periodic. The
object is defined by two planes.

\snippet testcases/FlowAroundCylinder/seeder.lua periodic

The last attribute in this case is the definition of the obstacle which is a 
cylinder in 3D and a sphere in 2D. Therefore, we need to define the object by 
origin and radius. The geometry kind is the special kind 'sphere'.

\snippet testcases/FlowAroundCylinder/seeder.lua cylinder

Run Seeder with 

\code
~/apes/seeder/build/seeder seeder.lua
\endcode

to get a mesh in the folder you have chosen. You can check your mesh using 
Paraview. 
Therefore write a harvester.lua script. 

## Visualisation of the mesh

To have a look onto the mesh in Paraview, we create a harvester.lua file that 
converts the mesh information into VTU format. Therefore, we need to define 
certain information. These are 

* simulation_name
* input
* output

We use 'FlowAroundCylinder2D' as simulation_name for example. To define the 
input, we have to set `mesh` to the location of the mesh files. To define the 
output we need the location where harvester shall dump the VTU and PVTU files. 
In addition to that we define the format as VTU and give a label like 'mesh' for
example.

\snippet testcases/FlowAroundCylinder/harvester.lua mesh

Now run Harvester with
\code
~/apes/harvester/build/harvester harvester.lua
\endcode
to get the VTU file and open it with Paraview.
The mesh looks like this:
\image html images/FlowAroundCylinder2D_mesh.png 
"The generated mesh visualised with Paraview"

# Musubi

You have a generated mesh now that can be used in musubi.lua. While you write 
the script for Musubi, you might have a look at the tutorials for Musubi if you
feel unsure about the context.

## general information

We want to simulate the pressure drop, the lift and drag coefficient and the 
recirculation length as mentioned above. 
Therefore we need to track the pressure over time at two points, one at the 
beginning and one at the ending of the obstacle. 
The difference between these pressure values will be plotted over time. 

## Which tracking variable values are needed for each quantity pressure 
## difference, recirculation zone, drag coefficient, lift coefficient?

The following refers to Gleaner which is software we use in our group.
Instead you can use Gnuplot, Matlab or matplotlib of python.

* pressure difference
  
  We need pressure_phys over time at two points. It will be plotted with 
  Gleaner.

* recirculation zone

  We need velocityX over length. Using the resulting ascii spatial file 
  in Gleaner we can compute it there.
  
* drag coefficient and lift coefficient

  We need to track the forces at the obstacle (sphere). It is directly possible
  to compute the coefficients in the tracking table. The results are dumped in 
  ascii file format.

A stable simulation is reached, if the pressure over time will reach a straight
line. Therefore plot it with Gleaner first.
In addition to that or instead, you can use the steady state criteria in 
musubi.lua that is described later on.

## preparations for the simulation

So we start with creating directories for tracking, output and restart. The 
names are free to choose, but in this tutorial we call them 'tracking', 
'output' and 'restart'. Make sure that these directories are located in the same
level like 'musubi.lua'. We have to do this every time we start a new simulation
so that it might be helpful to create a template which automatically creates 
directories that are needed to run the simulation progress. Therefore we use 
Shepherd.

\code
mkdir debug tracking restart output mesh
\endcode

## general musubi settings


As we make use of different inlet and outlet kinds, we define the currently used
kind in the beginning of the file as a variable which is called in the tracking 
table.

In addition to that it is necessary to distinguish between the used test cases 
with different file names. We also define the tracking, the mesh and the restart
folder in the beginning. 
In the following code, we can see the basic information for musubi.

\snippet testcases/FlowAroundCylinder/musubi.lua basic info

## sim_control

After that, we define the time step settings which is done inside the 
sim_control table. We have to define the maximum time, the interval between 
different time steps and an abort_criteria in order to finish the simulation as
a steady_state has been reached.

\snippet testcases/FlowAroundCylinder/musubi.lua time

## restart

This is followed by the *restart* table where we have the possibility to save 
our simulated results and resume a simulation from a certain time step. 
For the first run, the restart files are not existing so far. 
We have to write them first. This is why in the following code, the read statement is deactivated.

\snippet testcases/FlowAroundCylinder/musubi.lua restart

You can find more information about the restart table in the [Documentation tutorials](@ref tut_5restart).

## physics

In the next part, we have to give information about the *physics* and the 
*fluid*. In the physics table we define the physical reference values in SI unit
system. You can find all the possible variables in the Documentation for Musubi.
They belong to the mus_physics_module. It is the mus_physics_type. 
Have a look [here](@ref mus_physics_module::mus_physics_type).

In this test case we have information about the reference time and the reference
density. We can call it from common.lua.

\snippet testcases/FlowAroundCylinder/musubi.lua physics

## interpolation_method

As *interpolation_method* we choose linear interpolation. 
[Here](@ref intp_methods) we find information about the interpolation methods that are possible so far.

\snippet testcases/FlowAroundCylinder/musubi.lua interpolation

## initial_condition

Going on with the initial_condition table, here we define the initial pressure and the initial velocity. The velocity is a 
three-dimensional vector. That is why we have to give velocity values for each component.
There is a [tutorial](@ref tut_6initial) in the Documentation about initial_condition table that gives more information with links to papers etc.

\snippet testcases/FlowAroundCylinder/musubi.lua initial

## identify

The next step is to define the nature of the simulation. Therefore, define a 
table called `identify` with information about label, layout (stencil) and kind. 
 
\snippet testcases/FlowAroundCylinder/musubi.lua identify

There are more information about the identify table in the [tutorial section](@ref tut_1mus_config) of the Documentation. 

## boundary_condition

As shown in the tutorial for [Definition of Boundary Conditions](@ref tut_4boundaries) we define for each boundary which was set up
in seeder.lua kind and information about variable values if this is necessary. 
Therefore we open the bnd.lua file from the mesh folder to see the labels of 
the existing boundaries.

The *boundary conditions* used at the inlet and outlet are

* inflow BC: inlet_ubb
* outflow BC: outlet_expol, outlet_eq, outlet_dnt, outlet_pab
* outlet_expol - extrapolation boundary
* outlet_dnt - do-nothing
    * both outlet_expol and outlet_dnt from "Asymptotic Analysis of Lattice 
      Boltzmann Outflow Treatments, Junk, Michael Yang, Zhaoxia,Communications 
      in Computational Physics, 2011.
* outlet_eq - equilibrium boundary from "Analysis of open boundary effects in 
  unsteady lattice Boltzmann simulations, Izquierdo, Salvador, Martínez-Lera, 
  Paula Fueyo, Norberto, Computers & Mathematics with Applications, 2009.
* outlet_pab - pressure anti-bounce back
* inlet_ubb - velocity bounce back
    * both outlet_pab and inlet_ubb from "Characteristic nonreflecting boundary
      conditions for open boundaries in lattice Boltzmann methods, Izquierdo, 
      Salvador Fueyo, Norberto,Physical Review E 2008.

\note For group members that have access to the Collaboration Server:
      Four different boundary conditions are used for this testcase and their 
      comparison results are shown in the corresponding tickets steady flow #915 and 
      unsteady flow #916.

The first boundary we define is the inlet which is located in the west. As 
mentioned in the beginning, we can use the variable for the chosen inlet kind 
here. The velocity profile at the inlet looks like a parabol with its maximum 
at the middle of the channel. We use the incoming velocity that we have defined 
in common.lua here.

\snippet testcases/FlowAroundCylinder/musubi.lua inlet

The opposite boundary states the outlet of the channel where we can use several 
outlet kinds. At this boundary the pressure is known and has to be defined.

\snippet testcases/FlowAroundCylinder/musubi.lua outlet

The other boundaries represent walls. They are not part of the fluid domain. 
These are:

\snippet testcases/FlowAroundCylinder/musubi.lua wall

## variable

In order to compute lift and drag coefficient with Musubi, we have to define a 
variable table that is used by the tracking table in the next step. So we 
define the formulas for the coefficients and multiply the formula with the 
tracked forces later on. We start with the formulas now:

\snippet testcases/FlowAroundCylinder/musubi.lua coeff

Next, we define the variable table:

\snippet testcases/FlowAroundCylinder/musubi.lua variables

In the variable table, we define the coefficient vector itself and then we 
define the result as a multiplication by force and coefficient vector.
The next step would be to define the corresponding tracking table. 
We will do this in the next chapter.

## tracking

Now we have to define the tracking table. 

For each variable value we would like to track, we must define a table 
containing label, variable, shape, time_control, format and folder.
In the first simulation run, we only track the global shape and the pressure 
over time in order to have a look if a steady state is reached.

The first thing is to define the label of the tracked variable value. The label
will be part of the filename so that it is necessary to give a meaningful and 
recognizable name. We can use `..` to combine the values of different variables 
to a long character. Next, we define the variable we want to track. Here it is 
pressure. 

To track the pressure difference over time, we have to place two points in the 
channel, one at the beginning and another one at the ending of the sphere. 
This is right the middle of the channel height. Therefore we have to give its 
coordinates that are defined in the shape table. It looks quite similar to the 
geometry definition in seeder.lua.

> There is a tutorial for the tracking table in the [Documentation](@ref tut_3tracking)

### pressure over time at two points

After that, we have to make clear at which time steps we want musubi to track the pressure. In the time_control
table we define the min, max and the interval. The last things are to define the tracking format and the 
folder. Here is the code for pressure difference over time in ascii format.

\snippet testcases/FlowAroundCylinder/musubi.lua pressure

As result, we get the pressure in physical SI units at two points over time. 
So for all time steps, there are two pressure values.

### velocity over time 

To track the velocity over time at the middle of the channel, we use ascii 
format. The results are values for each direction, velocityX, velocityY, 
velocityZ at each time step. Here is the code to do that:

\snippet testcases/FlowAroundCylinder/musubi.lua velocity over time

### lift and drag coefficient at the obstacle

The lift and drag coefficient are very important characteristics of a flow 
around an obstacle. While tracking the forces at the obstacle we can directly 
compute them using the variable table which is explained above.
In the tracking table we use ascii format and the obstacle (sphere) as a shape. 

\snippet testcases/FlowAroundCylinder/musubi.lua Re20

### pressure and velocity over time along the global shape

To have a look at the flow around the cylinder in Paraview, we have to track the pressure and the velocity over
the global shape. We make use of harvester format to create a VTU file. The output is a binary file that is needed for Harvester.

\snippet testcases/FlowAroundCylinder/musubi.lua global flow

### velocity along a vertical line at last time step

We can track the velocity at the entry of the channel over the y-axis. The format is asciiSpatial. Our results will be 
values for each directory velocityX, velocityY and velocityZ at different locations which are defined by `segments`.

\snippet testcases/FlowAroundCylinder/musubi.lua velocity over height

### mean velocity along a vertical line at last time step

We can do the same once again for the mean velocity. But instead of getting values for each 
coordinate, we reduce the results to one single vector which is the average along the y-axis.
Therefore, we can use `reduction = average`:

\snippet testcases/FlowAroundCylinder/musubi.lua mean velocity

### pressure and velocity along the channel length at last time step

To have a look at the velocity and pressure expansion along the channel,
we might track them at the middle of the channel along the x-axis at the
last time step.

\snippet testcases/FlowAroundCylinder/musubi.lua velocity over length

### offending forces at the obstacle at the last time step

We can also compute the forces at the obstacle.
In this case, we need harvester format and the definition of the shape in order to visualize the forces. 
Forces can only be tracked using obstacles and q-values. So we make sure to activate q-values in
seeder.lua and in musubi.lua as described above.

\snippet testcases/FlowAroundCylinder/musubi.lua forces at the obstacle 

# Computing the quantities for the steady flow around the cylinder

After running the simulation by 

\code
~/apes/musubi/build/musubi musubi.lua
\endcode

we are able to visualize the flow around the cylinder.

Therefore we can use the tracking output in harvester format.
Every tracking object with output format == 'vtk' produces 
files in the output folder that can be used with Paraview.

Here is a screenshot from Paraview.

\image html images/global_flow.png "Flow Around Cylinder 2D, steady state after 30 s"

At the picture, we can see the velocity in x direction on the global shape. 
This image was taken from the last time step,
30 seconds in this case.

## Is a steady_state reached? 

\note In the following we create plots with Gleaner, a tool that is part of the APES suite
that is used in our group. You can use instead the free tools matplotlib (Python) or GNUplot.

We have to prove that a stable state is reached. Therefore, we could use the 
`probe_velocity` which is velocity at a certain point over time in ascii format.
We use Gleaner to plot the velocity over time. If a straight line is computed at
the end of the plot, a steady state has been reached.

So let us have a look at the Gleaner input file. We use a python script to do 
that. We could create a file `params_plot.py` containing all the plots. 
For now, we start only with this one and continue with the other 
plots later on, after a steady state has been reached.

At the very beginning of the file we write the following:

\snippet testcases/FlowAroundCylinder/params_plot.py header 

After that we place a table with our plotting information.

\snippet testcases/FlowAroundCylinder/params_plot.py velocityX over time

Now we run Gleaner and have a look at the plot. 

We run Gleaner with:

\code
~/apes/harvester/gleaner/gleaner.py params_plot.py
\endcode

If it has reached a steady state it looks similar to this image:

\image html images/steady_state.png "Steady state"

## Pressure difference

To have a look at the results, we have to run musubi and look into the tracking
folder. We recognize the label we defined before and see one *.lua file with 
information about the tracked pressure over time and one *.res file that 
contains three columns. The first column is filled with all the time steps. 
The second one contains the corresponding pressure values at coordinates 
(0.15, 0.2) and the last one at the coordinates (0.25, 0.2).

To have a look at the pressure difference over time we use Gleaner again.
Now, we add the pressure difference plot to `params_plot.py`.

\snippet testcases/FlowAroundCylinder/params_plot.py pressure difference

You can see that there are two additional lines. They are guidelines from the 
benchmark table. The results should be near these guidelines.

We run Gleaner with:

\code
~/apes/harvester/gleaner/gleaner.py params_plot.py
\endcode

The following image shows the results from Gleaner.

\image html images/deltaP.png "Pressure difference"

For the computation we made use of the formula given in the paper mentioned 
above: 

\f$ \Delta P = \Delta P(t)=P(x_{a},y_{a},t)-P(x_{e},y_{e},t)\f$

If we have a look at our simulation directory we will find another directory 
called `filedata`. This is the place where output from Gleaner is dumped. 
There we find the pressure difference results in a ascii file. 
Here you can see the content of this resulting file:

\code
#filaname 	 average 	 min 	 max 	 last 
tracking/FlowAroundCyl_2D_probe_pressure_p00000.res	0.111775987288	0.0	0.417952364225	0.112177193548
\endcode

## Recirculation length

In order to compute the recirculation length, we can use Gleaner once again. 
There is a special Gleaner `kind = rLength` for that.
We will get only the computed value but no plot in this case.
We can enter our predefined *params_plot.py* and add a few lines there:

\snippet testcases/FlowAroundCylinder/params_plot.py recirculation length

There are some values needed. `re_start` is the beginning point of the 
recirculation zone. This value is given by the paper with

\f$x_{e}=0.25\f$. 

It is the place right behind the obstacle. `dx` is computed by our file 
`common.lua`. We need to print the variable first, before we get its value.
> `printSimParams.lua`

The results of Gleaner are dumped in the shell. It will give us two values, 
one with and one without using `dx`. 
It could look like this:

\code
kind = rLength
recirculation length: 0.06979999999999992
recirculation length using dx: 0.06569999999999991
\endcode

Once again we will find the results dumped under `filedata/` in a ascii file.

## Drag and lift coefficient

As described in the tracking chapter, we compute the lift and drag coefficient 
during the simulation run. Therefore we just have a look at the .res file for 
cyl_coeff in the tracking folder. The results could look like the following.

\code
# Rank of the process:       0
#                   time             coeff_red_01             coeff_red_02             coeff_red_03
 0.3000024666665035E+002  0.4937947108860087E+001  0.1379254135365926E-001  0.0000000000000000E+000
\endcode

`coeff_red_01` represents the drag coefficient and `coeff_red_02` is equivalent
to the lift coefficient. The values are taken at the last time step.

The results in this tutorial are created with a refinement level of 6 and a 
small number of unkwons of 12407.

Let us go on with the next test case, \ref tut_unsteady2D "Unsteady Flow Around Cylinder 2D &rarr;".
