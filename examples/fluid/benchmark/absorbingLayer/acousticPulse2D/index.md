title: Acoustic pulse 2D

## Acoustic 2D pulse test cases with different spatial absorbing layers ##
{ex_fluid_abslayer_pulse2D}

The objective of the acoustic 2D pulse examples are to introduce the following
features in Musubi:

* Define a function (pulse) to use it for boundary conditions.
* Define absorbing layer as source term.
* Setup of different sponges as listed below.

Here are the list of test cases to learn these features in Musubi:

* [Absorbing layer plane](absLayerPlane/index.md)
* [Absorbing layer box](absLayerBox/index.md)
* [Absorbing layer radial](absLayerRadial/index.md)

In these test cases, in addition to the aforementioned features, you will also
learn how to:

* Create a mesh with boundaries using Seeder.
* Post-processing nd visualization in Paraview.
* Create a 2D plot using Gleaner tool. Gleaner is a Python tool which
extracts data from Musubi ascii output and uses matplotlib in python library
to create a plot.
