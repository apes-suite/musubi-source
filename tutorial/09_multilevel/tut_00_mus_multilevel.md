Multilevel Simulations {#tut_00multilevel}
========

Navigate: \ref tut_7convergence "&larr; Setting up a Convergence Criterium"
| \ref mus_tutorials_overview "Overview"

In some simulations there are regions of special interest. In these regions 
strong gradients might occur, so that a higher resolution is required. 
You then have basically two possibilities to resolve this problem.
1. increase the global resolution of your simulation. This is generally 
   not a good idea because you very quickly might run into memory problems.
   A high number of fluid cells also means long computation times. 
2. Locally refine the mesh in regions of high interest
   You can adaptively refine the mesh in regions where it is desired.
   This focuses the computation cost and memory expenses on the regions
   where it actually is required. Why waste computation time on flow regions
   where very little happens? A general problem are the interfaces between 
   grid levels. You can easily get reflections or non-discontinuities there
   which is why you usually place interfaces in regions where the flow does 
   not fluctuate considerably.

> In mus_interpolate_module you can get some background information on the implemented 
> interpolation methods and the general workflow


## Test Case Description ##

In this tutorial we cover a channel test case with a cylinder inside. 
The area around and behind the cylinder is refined by a refinement 
patch with a higher resolution.


## Mesh Generation ##

Before starting your simulation you have to set up a mesh again.
Again we involve Seeder to generate the required data structures for the
solver.
It is very advisable to define a variable for the most important variables.
Later on you will be able to easily change the resolution of your simulation 
or some aspect ratio.
One very important parameter for setting the resolution of the simulation
is a reference \ref tem_topology_module::tem_levelof "tree level". 
In our case this will be the tree level of the 
channel region. Let's name it `level` and set it to 6
\snippet 99_multilevel/seeder.lua reference level
> Note: For a Description of levels and the layout of the tree have 
> a look at \ref tem_distributed_octree .
Another very important definition is the size of your bounding 
cube. Again, let's define a variable and set all other length 
definitions to depend on this variable. 
Here we set it to the length of the bounding cube and define the length
and the height of what later will become the channel
\snippet 99_multilevel/seeder.lua reference length
The bounding box, which is the universe in which all elements of the tree will live in, 
is defined as
\snippet 99_multilevel/seeder.lua bounding box

The region around the channel is set to the desired tree level of 6. 
\snippet 99_multilevel/seeder.lua initial box
The rest of the bounding cubic domain is not of interest and hence the discretization
level is not of interest. 

### Defining Geometry ###

We need to specify the walls, the inlet, outlet and the cylinder. 
Let's start with the walls at the north, south position.
The walls at east and west direction will later on become the in- and outlet.
\snippet 99_multilevel/seeder.lua geometry walls
 
The walls at the top and bottom (in z-direction) are added 
by extending the `geom` table 
\snippet 99_multilevel/seeder.lua front back planes

For defining the cylinder, we use an STL file. This file was before generated
by Blender, but you can basically use any software you like for generating
such geometry data. Let's include that STL file with 
\snippet 99_multilevel/seeder.lua geometry cylinder

One very important action is the placement of the seed.
The seed determines the contiguous flow domain. It basically 
defines for the shapes, what is inside and what is outside.
We palce it right in the middle of the channel, which is simply 
\snippet 99_multilevel/seeder.lua seed definition

Ok. Now we have defined all geometric constraints.
Let's continue with refining some parts of the simulation domain.

### Defining Refined Regions ###

The region around and behind the cylinder is defined by a refinement box.
Let's define first, by how many levels we want to refine the elements
inside this refined area. For simplicity, let's say that all elements should
be on one level higher than the rest of the channel. 
\snippet 99_multilevel/seeder.lua refinement level
It is also good to have the information about the maximum level in your
simulation available as a parameter. You will later on see why.
\snippet 99_multilevel/seeder.lua maxlevel

You first need to specify the origin of this box and its extents.
Again, we use some variables depending on the total bounding cube.
\snippet 99_multilevel/seeder.lua patch size
The patch will then lie around and behind the cylinder. 

 \snippet 99_multilevel/seeder.lua refinement box

Once you followed through the above explanations, you can visualize the 
generated mesh.

\image html tut_multilevel_cyl_mesh.jpg


## Setting up the Musubi Configuration ##

After generating the mesh above, we need to tell Musubi that it
should use the above generated mesh. The mesh was generated in the
folder `mesh`. Let's define that along with a name for the simulation
\snippet 99_multilevel/musubi.lua general settings

### Initial and Boundary Conditions ###

Initial conditions are set to a medium at rest with constant pressure
\snippet 99_multilevel/musubi.lua initial conditions
The boundary conditions are a little bit more complex, as we have 
solid walls, the cylinder, and also an inlet and an outlet.
Let's start with the wall boundaries. These include among the walls, which we named 
according to the directions and the cylinder, which we simply called `obs`.
They all get the property of a simple wall.
\snippet 99_multilevel/musubi.lua wall boundary conditions
For the outlet, we would like to have a 
\ref mus_bc_fluid_module::outlet_pab "simple pressure" boundary condition.
All you need to specify is a reference pressure. 

\snippet 99_multilevel/musubi.lua outlet condition
The inlet is slightly more complicated. Here we would like to set a 
velocity inlet. We would like to have a spatial distribution of the
x-velocity component according to a 
\ref tem_spatial_module::tem_spatial_parabol3d_for_treeids "three-dimensional parabola".
Also, the velocity should be slowly 
\ref tem_transient_module::load_transient_predef_fun "ramped" 
up to get a smooth transition from the initial state to the steady state.
Afterall, the definition for the inlet boundary condition looks like
\snippet 99_multilevel/musubi.lua inlet condition

### Other simulation parameters ###

\snippet 99_multilevel/seeder.lua iterations
\snippet 99_multilevel/musubi.lua time settings


### Tracking ###

Pressure / Density probe
\image html tut_multilevel_probedensity.png

### Visualize the Results ###

For visualisation of ascii or ascii-spatial format files you can use for example
Gnuplot or matplotlib module of python.

Navigate: \ref tut_7convergence "&larr; Setting up a Convergence Criterium"
| \ref mus_tutorials_overview "Overview"
