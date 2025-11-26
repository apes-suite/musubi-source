title: Channel 2D
@warning WORK IN PROGRESS @endwarning

## Channel 2D ## {#ex_channel2D}

The objective of the channel 2D examples are to introduce the following 
features in Musubi:

* Boundary conditions to define pressure and velocity at the boundaries.
* Q-Values at boundaries to improve the accuracy of boundary conditions.  
* Multilevel simulations
* Non-Newtonian models.
* External forces. 

Here is the list of test cases to learn these features in Musbi:

* [Simple](C2D_Simple/index.md)
* [Boundary conditions](C2D_BoundaryConditions/index.md)
* Flow around the cylinder
    * [Single-level](C2D_Cylinder_SingleLevel/index.md)
    * [Multi-level](C2D_Cylinder_MultiLevel/index.md)
* [Non-Newtonian](C2D_NonNewtonian/index.md)
* [External force](C2D_Force/index.md)    

In these test cases, in addition to the aforementioned features, 
you will also learn how to:

* Create a mesh with boundaries using Seeder.
* Post-processing nd visualization in Paraview.
* Create a 2D plot using Gleaner tool. Gleaner is a Python tool which
extracts data from Musubi ascii output and uses matplotlib in python library
to create a plot.
