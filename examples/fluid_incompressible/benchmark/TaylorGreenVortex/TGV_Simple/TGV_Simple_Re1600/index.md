title: Taylor-Green Vortex at Re = 1600
@warning WORK IN PROGRESS @endwarning

# Taylor-Green Vortex # {#eg_TGV_Re1600}

In this example, we will investigate the creation and evolution of vortices in a
simple pre-defined cube with periodic boundaries. For this simple geometry we
even do not need seeder.
The case is based on the well-known Taylor-Green vortex benchmark case, which is
subject of the following papers:
[1] Brachet, M. E., Meiron, D. I., Nickel, B. G., Morf, R. H., Frisch, U., &
Orszag, S. A. (1983). Small-scale structure of the taylor-green vortex. Journal
of Fluid Mechanics, 130, 411–452.
[2] Brachet, M. E. (1991). Direct simulation of three-dimensional turbulence in
the Taylor-Green vortex. Fluid Dynamics Research, 8(1–4), 1–8.
[3] DeBonis, J. R. (2013). Solutions of the Taylor-Green vortex problem using
high-resolution explicit finite difference methods. 51st AIAA Aerospace Sciences
Meeting Including the New Horizons Forum and Aerospace Exposition 2013,
(February 2013).

The objectives of this example is to introduce how to:

* Create the mesh internally with periodic boundaries using Musubi.
* Simulate the evolution if the vortices in the cube using Musubi.
* Visualize the simulation results in Paraview.
* Create a 2D plot using the Gleaner tool. Gleaner is a Python tool which
  extracts data from Musubi ascii output and uses matplotlib in python library
  to create a plot.
