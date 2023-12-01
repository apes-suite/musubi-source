Flow around cylinder {#tut_FlowAroundCyl2D}
========

# Introduction

In the paper "Benchmark Computations of Laminar Flow Around a Cylinder" by M. Schaefer and S. Turek a few test
cases concerning the flow around a cylinder in 2D and 3D are demonstrated. This test case here has been created
with those information from the paper and its results has been compared with other solver's results.

\image html images/channel2D_geometry.png "Geometry of the 2D channel"

In this image we can see the geometry definitions of the channel containing the cylinder in 2D. Exactly this 
channel has to be created then with Seeder. 

With Musubi we aim to get results about the pressure drop from the front to the end of the cylinder, 
recirculation length, the drag coefficient and the lift coefficient.

There are three 2D test cases. The first one is a steady flow using the Reynolds number Re=20 and an incoming velocity of 
0.3 m/s. The second one describes the unsteady flow with a higher Reynolds number Re=100 and a higher incoming velocity of
1.5 m/s. 

We start with the \ref tut_steady2D "Steady Flow Around Cylinder 2D &rarr;".
