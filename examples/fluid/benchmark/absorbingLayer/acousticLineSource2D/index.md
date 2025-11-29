title: Acoustic Line Source 2D

@warning WORK IN PROGRESS @endwarning

Navigate: [&larr; Test case absorbing layer](../index.md)
| [Test case acoustic pulse 2D &rarr;](../acousticPulse2D/index.md)
| [Test case acoustic pulse 3D &rarr;](../acousticPulse3D/index.md)
| [Test case acoustic cylinder 2D &rarr;](../acousticCylinder2D/index.md)

## Acoustic 2D line source for plane absorbing layer ## {ex_fluid_abslayer_lineSrc}

This test case provides an example for using absorbing layer at outlet boundary.
The ambient pressure is distrubed by a 2D acousic line source. The total
pressure field (ambient plus fluctuations) is defined at the west boundary and
ambient pressure is defined at the east boundary. The planar absorbing layer
is defined at the east boundary with a thickness of 0.8.
The objectives of this example is to introduce how to:

* Create a mesh with boundaries using Seeder.
* Post-process the mesh using Seeder-harvester and visualize it Paraview.
* Simulate a 2D acoustic line source in a channel using Musubi.
* Set up an absorbing layer.
* Validate the numerical results.
* Visualize the flow in Paraview.
* Create a 2D plot using Gleaner tool. Gleaner is a Python tool which
extracts data from Musubi ascii output and uses matplotlib in python library
to create a plot.

## Generating mesh ##
### Define geometry information ###
### Define spatial objects ###

## Running simulation ##
### Define flow parametes ###
### Define collision parameter ###
### Define boundary condition ###
## Post-processing ##

Here are the results from the simulation.

Pressure across the length of the channel:
![Pressure_Profile](PressureAlongLine.png)

To create these plots, run <tt>python plot_track.py<tt> to create the plots.
Before running the plot script, open 'plot_track.py' and update path to
Gleaner script in 'glrPath'.
Download Gleaner script using
<tt>hg clone https://geb.sts.nt.uni-siegen.de/hg/gleaner</tt>
