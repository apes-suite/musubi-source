Tracking quantities during the simulation run {#tut_3tracking}
========

Navigate: \ref tut_2toolchain "&larr; Using the full toolchain"
| \ref mus_tutorials_overview "Overview"
| \ref tut_4boundaries "Definition of Boundary Conditions &rarr;"

In this tutorial we deal with tracking, that is tracking flow quantities 
during the simulation run and visualizing them.
> An extensive documentation of the tracking features can be found in 
> \ref tem_tracking_module 

In the `musubi.lua` file, there has to be a section called `tracking`.
In order to get a simulation output for harvester you set up some
information for each tracked quantity. 
These are in general:
* label 

This defines the name of the tracked quantity.

* variable 

What variable values are tracked? (i.e. pressure, velocity)

* shape

How does the geometry look like in the mesh? Shall it be a point
or a line? Where is it placed? These are the questions you are
supposed to answer here.

* time_control

This behaves similar to the general time settings in the
musubi file. You have to define a starting and an end point
for the tracked quantity. In addition to that, you must 
give musubi information about the time interval, when
i.e. the pressure gets its values during the simulation.

* format

   You can choose between these three formats:

  * vtk

    This format creates results in the treelm format. It is used to process 
    input files for Paraview. When using this format, make sure there is a
    directory called `output` in the simulation path.
    This has been done in the last tutorial as well. Here you can see another example:

    \image html screenshot_cylinder_shearStress.png "Shear stress in a channel"

  * ascii

    This output format creates a .res file that contains the tracked quantities
    like pressure and velocity over certain time steps. It is normally used for 
    1D shapes (a single point). Nevertheless, it is possible to track the variable values
    over a 2D shape. In this case, a reduction is used. There are picked up
    certain points of the shape and they are reduced to a single result.
    The table in the resulting file is build up as following.
    In the first column, you can find the time (t_0, t_1, t_2, ..., t_n). It is
    followed by the variable values. Every row replaces one timestep with values for 
    each variable. 
    Here you can see an example that is plotted by GNUplot:

    \image html images/pressureOverTime.png "Pressure over time"

  * asciiSpatial      
    
    This kind of format creates not only one res. file. Moreover it creates 
    one .res file for each time step. This means, that you get a table with 
    the x, y and z coordinates of the shape in the first three columns. These
    are followed by the different variable values like pressure and sheer stress.
    In summary, you get all these information for each time step in another
    file. You can plot these values as well, so that you get something like this:

    \image html images/pressureOverLength.png "Pressure over length"

* folder

  Here, you have to define the location where your resulting files shall be dumped.

Here is a demonstration example of how a tracking section in the musubi.lua file 
can look like:

\snippet testcases/channel/musubi.lua Tracking example

Of course, you have to define your variable values first before you can track them.

Here ends the tutorial for tracking quantities during the simulation. You can go on for
more information about boundary conditions in the next tutorial.

Next chapter: \ref tut_4boundaries "Definition of Boundary Conditions &rarr;"
