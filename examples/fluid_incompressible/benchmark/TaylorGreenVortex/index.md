title: Taylor Green Vortex
@warning WORK IN PROGRESS @endwarning

Progress of this test case collection for TGV is tracked in the following ticket:
Action #2246: Test case - Benchmark: Taylor-Green Vortex

# Taylor-Green Vortex # {#eg_TGV}

In these examples, we will investigate the creation and evolution of vortices in
simple pre-defined cubes with periodic boundaries. For such a simple geometry we
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

The objective of the TGV examples are to introduce the following features in
*Musubi*:

* Define functions (TGV) to use them for boundary conditions.
* Calculate the kinetic energy after runtime and compare it to reference data.
* Use of different turbulence models.

Here is the list of test cases to learn these features in Musubi:

* [Simple TGV](TGV_Simple/index.md)
* [TGV with different LES models](TGV_LES/index.md)

In these test cases, in addition to the aforementioned features, you will also
learn how to:

* Post-processing nd visualization in Paraview.
* Create a 2D plot using Gleaner tool. Gleaner is a Python tool which
extracts data from Musubi ascii output and uses matplotlib in python library
to create a plot.
