title: Taylor-Green Vortex with WALE Smagorinsky (PDF) model at Re = 1600
@warning WORK IN PROGRESS @endwarning

# Taylor-Green Vortex LES-Smagorinsky (using PDF) # {#eg_TGV_SmagPDF}

In this example, we will investigate the creation and evolution of vortices in a
simple pre-defined cube with periodic boundaries. For this case, the Reynolds
number (*Re*) is set to 1600.

A detailed description of test case can be found in the parent directory
[Description of test case TGV](../../index.html).

A description about the LES cases can be found here:
[Description of test case TGV_LES](../index.html).

The objectives of this example is to introduce how to:
* Use a pre-defined geometry instead of creating a mesh with Seeder.
* Use a turbulence model.
* Simulate the Taylor-Green Vortex in the cube using Musubi.
* Create 2D plots using the Gleaner tool. Gleaner is a Python tool that extracts
  data from Musubi ASCII output and uses the plotting library Matplotlib in
  Python to create a plot.
* Post-process the results by calculating the volume average of the tracked
  kinetic energy.
* Calculate the dissipation rate by means of the tracked quantities.
* Compare the tracked and calculated dissipation rate with a reference solution.

## Running simulation ##
For the Smagorinksy turbulence model, the strain rate can be computed from
the particle distribution function (PDF) as in this example here or from
velocity gradient. For the setup of latter one see:

* [TGV-LES Smagorinksy from Velocity Gradient ](../TGV_SmagGradU/index.html)

### Define flow parametes ###

### Define collision parameter ###

## Post-processing ##

Here are the results from the simulation.

Kinetic energy over time compared:
![Kinetic Energy](media/KineticEnergy.png)

Dissipation rate over time compared to reference solution from Brachet:
![Dissipation Rate](media/DissipationRate.png)

To create these plots, run <tt>python plot_track.py<tt> to create the plots.
Before running the plot script, open 'plot_track.py' and update path to
Gleaner script in 'glrPath'. Download Gleaner script using
<tt>hg clone https://geb.sts.nt.uni-siegen.de/hg/gleaner</tt>
