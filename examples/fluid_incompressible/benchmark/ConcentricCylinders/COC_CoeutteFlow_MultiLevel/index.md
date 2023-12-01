title: Concentric Cylinders 
@warning WORK IN PROGRESS @endwarning

Navigate: [&larr; Test Case Concentric Cylinder Couette Flow](../index.html)
| [Test Case Concentric Cylinder Couette Flow Single-Level &rarr;](COC_CoeutteFlow/index.html)

# Concentric Cylinders-Couette Flow Multi-level#

Couette flow between two circular cylinders is used to investigate the behaviour
of velocity boundary condition for complex geometries with q-values.
In this flow, the inner cylinder with radius \(r_1=0.25/2 m\) rotates with a
constant angular velocity \(\omega=0.213 rad/s\) and the outer cylinder with
radius \(r_2=0.5 m\) is kept stationary.

This Couette flow has the following analytical solution
\begin{equation}
 u_\theta(r) = \frac{u_0 \beta}{1-\beta^2}
 \left( \frac{r_2}{r} - \frac{r}{r_2}\right)
\end{equation}
where \(u_0=\omega r_1\) and \(\beta=r_1/r_2\).

This analytical solution is provided as variable `vel_an` in the configuration
of the testcase.

The mesh created here is multiple resolutions.
Distance refinement is used to refine towards both inner and outer cylinder
with level=minlevel+2. The element size near the cylinder is similar to 64
elements in inner cylinder of a single level mesh.

To decide the convergence, the x axis in the right side of the cylinder is
used. This line is defined as the `track_line` table and used in the abort
criteria, a tracking over time and the spatial tracking for the final output.
Once the average velocity magnitude over this line doesn't change more than
the provided threshold in the abort section anymore or the `tmax_phy` time is
reached, the simulation will stop.

The `run.sh` runs the complete simulation, from mesh generation to creating an
image of the velocity profile over the `track_line`.
You need to set the executable paths there accordingly and have gleaner for it
to work.
