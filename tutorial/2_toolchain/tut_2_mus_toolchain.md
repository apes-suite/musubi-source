Using the full toolchain {#tut_2toolchain}
========

Navigate: \ref tut_1mus_config "&larr; Set up your first simulation with Musubi"
| \ref mus_tutorials_overview "Overview"
| \ref tut_3tracking "Tracking quantities during the simualtion run &rarr;"

In this tutorial, you will learn how to use the complete \ref toolchain "tool-chain" of the APES suite to deploy
simulations. This includes the generation of the mesh with Seeder, 
Running the actual simulation with Musubi and finally post-process the results with Harvester
which generates files to be visualized with standard tools like Gnuplot and Paraview.
Therefore we use a channel as an example.

## Preparing the Simulation Folder ##

Generate the folders
* mesh
* tracking
* restart
* harvest


## Generate the Mesh ##

Make sure that Seeder is compiled and you know how to call it.

At first, we define the common header of the file containing general information like the name and
the refinement level of the mesh. 

\snippet testcases/channel/seeder.lua common header

In the next part we have some variables with boolian values. We need them to activate or deactivate 
some parts of the script. With this script we can have a look at the effect of qvalues, periodic boundaries
and the use of different STL files (3D model).  

\snippet testcases/channel/seeder.lua geometry definition

After that, we can specify some basic relations of the channel's geometry.

\snippet testcases/channel/seeder.lua geometry relations

You might remember the `fixHeight` value from above. If it is set to `true` the `dx` is computed by
the specified height of the channel. Else `dx` is computed defining the length of the channel first.

\snippet testcases/channel/seeder.lua dx computed by height

\snippet testcases/channel/seeder.lua dx computed by length

Here, we have some reference values that have effect on the geometry of the channel. 
You might change them and compare the effects later on in Paraview.

\snippet testcases/channel/seeder.lua reference values

Next, we define the bounding cube. 
The bounding cube, i.e. the box that includes the complete simulation domain.
It must be large enough, otherwise it just truncates all geometries that go beyond this cube.

\snippet testcases/channel/seeder.lua bounding cube

This part is followed by the `spatial_object` table. Here we define the seed, the channel itself
and obstacles. As you see, we have only the seed defined here. The missing parts depend on the use of 
periodic boundaries, obstacles and more that are activated and deactivated above.

\snippet testcases/channel/seeder.lua spatial objects

In this script we make use of an area with special refined mesh. Therefore we defined the levels above.
You can enter a lua table from any place after its definition. So we are able to deactivate
the following refinement boxes if we think that they are not necessary. To do this, you just make use of `--`
in front of each line.

\snippet testcases/channel/seeder.lua refinement box

Here we can define the inlet and the outlet of the channel. They are boundaries and belong to 
the `spatial_object` table. 

\snippet testcases/channel/seeder.lua inlet and outlet boundary conditions

It is possible for the rectangular channel to make use of periodic boundaries. Therefore you need to define 
two planes close to each other that are repeated along the channel length (`usePeriodic=true`) or the height
(`periodicX=true`). There exist these scripts:

\snippet testcases/channel/seeder.lua use periodicX boundaries

\snippet testcases/channel/seeder.lua use periodic boundaries

Of course you do not have to use periodic boundaries. Instead you could use a script that defines specific
boundaries for `top` and `bottom` of the channel. 

\snippet testcases/channel/seeder.lua non-periodic

In order to simulate a flow around an obstacle you might activate `useObstacle` (`=true`). Here you can choose
between two stl files for a sphere (2D) or a cylinder (3D). You have to make sure that the files are located in the 
defined directory. Now it makes sense to activate qvalues. This allows to refine the mesh at i.e. circular
boundaries. This is done with `calc_dist = true` in the attribute table of a spatial object.

\snippet testcases/channel/seeder.lua use of obstacles

If you do not want to use qvalues you need this script:

\snippet testcases/channel/seeder.lua cylinderObj

Now you can start Seeder. If it successfully finishes, the simulation mesh was generated.

## Control the generated Mesh ##

Harvester can process the mesh and generate a VTU file, which can then be visualized by Paraview.
Let's prepare the harvester.lua for visualization.
\snippet testcases/channel/harvester.lua input mesh visualisation
How to output the mesh and what
\snippet testcases/channel/harvester.lua output
Visualize the mesh with paraview. If you set up `usePeriodic = false` the VTU file should look like this:
\image html images/use_periodic_false.png "3D channel"
Otherwise you can build a two-dimensional channel like this:
\image html images/use_periodic_true.png "2D channel"

## Run Musubi ##

Let's have a look at how to configure Musubi. 
\note Many things were mentioned in the last tutorial

\note Describe here the most important parameters for running the simulation.


In order to repeat the most important things, you can make use of a checklist. 
For the simulation you need to define the 
  1. `simulation_name`
  2. `mesh` The location of the directory is needed.
  3. `identify = { label = "..", layout = "..", kind = "..", relaxation = ".."}`
     There exist default values, which are explained in the [last tutorial](@ref tut_1mus_config)
  4. time settings that describe the time interval, the time maximum and an abort criteria
  5. `fluid` table with information about physics and LBM units, i.e. the fluid density
  6. `initial_condition` table with the initial pressure and velocity
  7. `boundary_condition` This table aims to tell Musubi how to handle the boundaries that were set up
      in seeder.lua. It reads the boundaries' labels and converts them into a wall, an inlet or an outlet.
  8. `tracking` table for information about the target output. It exist a description module
      for tracking that you find [here](@ref tracking_mod)
      


## Post-process the Results ##

Run Harvester on the tracked quantities:

\snippet testcases/channel/harvester.lua input mesh false

\snippet testcases/channel/harvester.lua output

Here, you can see the output in Paraview:

\image html images/cylinder_pressure_screenshot_2D.png "2D channel with tracked pressure"

Now you are free to go on with the next chapter that deals with Tracking.

Next chapter: \ref tut_3tracking "Tracking quantities during the simulation run &rarr;"
