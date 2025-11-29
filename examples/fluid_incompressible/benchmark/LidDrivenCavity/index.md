title: Lid driven cavity
@warning WORK IN PROGRESS @endwarning

Navigate: [&larr; Benchmark Test cases](../index.html)
| [Test case lid driven cavity simple &rarr;](LDC_Simple/index.html)

# Lid driven cavity in 2D #

Lid-driven cavity is a simple 2D problem with simple boundary conditions.
In a square domain, dirichlet boundary conditions are set on all sides.
Fixed wall (\( \vec{u} = 0\)) boundary conditions on 3 sides (west, east and
south) and moving wall (\( u_x = u_c, u_y = 0\)) boundary condition on north
side.

Here are the list of test cases to learn these features in Musubi:

* [Simple](LDC_Simple/index.html)
* [Multi-level](LDC_MultiLevel/index.html)

In these test cases, in addition to the aforementioned features, 
you will also learn how to:

* Create a mesh with boundaries using Seeder.
* Post-processing nd visualization in Paraview.
* Create a 2D plot using Gleaner tool. Gleaner is a Python tool which
extracts data from Musubi ascii output and uses matplotlib in python library
to create a plot.
