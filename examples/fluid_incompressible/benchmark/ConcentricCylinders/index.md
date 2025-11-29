title: Concentric Cylinders
@warning WORK IN PROGRESS @endwarning

Navigate: [&larr; Benchmark Test Cases](../index.html)
| [Test case concentric cylinder Couette flow &rarr;](COC_CoeutteFlow/index.html)

# Concentric Cylinders-Couette Flow #

Couette flow between two circular cylinders is used to investigate the behaviour
of velocity boundary condition for complex geometries with q-values.
In this flow, the inner cylinder with radius \(r_1\) rotates with a constant
angular velocity \(\omega\) and the outer cylinder with radius \(r_2\) is kept
stationary. This Couette flow has the following analytical solution
\begin{equation}
 u_\theta(r) = \frac{u_0 \beta}{1-\beta^2}
 \left( \frac{r_2}{r} - \frac{r}{r_2}\right)
\end{equation}
where \(u_0=\omega r_1\) and \(\beta=r_1/r_2\).

Here are the list of test cases:

* [Single-level](COC_CoeutteFlow/index.html)
* [Multi-level](COC_CoeutteFlow_MultiLevel/index.html)

The geometry of the this cylindrical channel is created with the help of the
STL in `cylinder_Dia1m.stl`.

In these test cases, in addition to the aforementioned features,
you will also learn how to:

* Create a mesh with boundaries with q-values using Seeder.
* Post-processing and visualization in Paraview.
* Create a 2D plot using Gleaner tool. Gleaner is a Python tool which
extracts data from Musubi ascii output and uses matplotlib in python
to create a plot.
