Unsteady Flow Around Cylinder 2D {#tut_unsteady2D}
=======

Navigate: \ref tut_steady2D "&larr; Steady Flow Around Cylinder 2D" 
| \ref tut_FlowAroundCyl2D "Test case Flow Around Cylinder 2D"

# Simulation of the unsteady flow

In this test case, the inflow conditions are not the same as for the steady 
flow. This is why we have to change parts in our created files that are already
existing from the last test case.

> Please make sure that you have gone through the Steady Flow Around Cylinder 
> test case. We need those files here.

## inflow conditions

* velocity
  \f$
  U(0,y,t)=4U_{m}y\frac{(H-y)}{H^{2}}, V = 0
  \f$
* Reynolds number
  \f$
  Re = 100
  \f$
* physical incoming velocity
  \f$
  U_{m} = 1.5 \frac{m}{s}
  \f$

The other conditions like mean velocity and geometry remain the same.

## What shall be computed?

* drag coefficient \f$c_{D}\f$ and its maximum
* lift coefficient \f$c_{L}\f$ and its maximum
* pressure difference \f$\Delta P(t_{1},t_{2})\f$ at \f$[t_{0},t_{0}+\frac{1}{f}]\f$
* pressure difference \f$\Delta P(t)\f$ at \f$t=t_{0}+\frac{1}{2f}\f$
* Strouhal number  
  \f$
  St=\frac{Df}{\bar{U}}=\frac{D}{T\bar{U}}
  \f$

## common.lua

Now we go on with the common.lua file. Here, we just have to change the Reynolds
number: 

instead of
\code
Re = 20
\endcode
we write
\code
Re = 100
\endcode

## musubi.lua

In the musubi file we have to change the tracking table. The other ones remain 
the same.
So we have to think about which variables have to be tracked in order to compute
the above mentioned quantities.

* pressure difference
 
 In this case, the pressure difference is needed for a period of time at one 
 point in the channel.

### tracking table

What we need to track:

* velocity over time
* global flow 
* pressure over time
* lift and drag coefficient over time
* velocity over length

#### velocity over time

At first we have a look if our simulation reaches a steady state. Therefore, 
we need velocity over time.
We deactivate our old tracking table from the last test case and build a new one
with our first entry using ascii format for Gleaner.

\snippet testcases/FlowAroundCylinder/musubi.lua velocity over time

#### global flow

The next one will be the velocity and pressure over time regarding the global 
shape using harvester format for Harvester.

\snippet testcases/FlowAroundCylinder/musubi.lua global flow

#### pressure over time

In order to compute the pressure difference later on, we add pressure over time
at two points to tracking table:

\snippet testcases/FlowAroundCylinder/musubi.lua pressure

#### lift and drag coefficient over time

To go on with the lift and drag coefficients we add the coefficient over time to
the tracking table:

\snippet testcases/FlowAroundCylinder/musubi.lua Re100

#### velocity over length

To have a look at the development of velocity vector over length we need to 
track velocity in ascii spatial format. It is tracked at the last time step.

\snippet testcases/FlowAroundCylinder/musubi.lua velocity over length

What we need to visualize:

* global flow as animation

In order to visualize the global flow as an animation we have to create a file 
with pvd format.
Musubi generates the file for you and stores it in the output folder with all 
vtu files. Have a look at your tracking or output folder for those vtu-, vtk-
and pvd-files.

Now open the PVD file with Paraview.

## params_plot.py

What we need to plot:

* drag coefficient over time
* lift coefficient over time
* pressure difference over time
* frequency
* velocityX over time
* velocityX over length
* velocityY over time
* velocityY over length

This is the input file for Gleaner that will plot the tracked quantities of 
Musubi. It belongs to the post-processing part.

In the beginning we define some basic settings for Gleaner like if the plot 
shall be shown after running Gleaner, or if the output shall be written to 
files, or just the font type.

\snippet testcases/FlowAroundCylinder/params_plot.py header

After that, we can list each plot step by step.
So we do it step by step. The first one are velocity over time which are 
separated to x and y components.
We use a simple `xy` plot. We can see when a steady state is reached by the 
periodic signal.

> How we use Gleaner is shown in the Gleaner tuturial 

### velocityX over time 

\snippet testcases/FlowAroundCylinder/params_plot.py header

### velocityY over time

\snippet testcases/FlowAroundCylinder/params_plot.py velocityY over time

Here we can see both plots:

\image html images/velX_over_time.png "velocity X over time"
\image html images/velY_over_time.png "velocity Y over time"

The next thing is to plot the lift and drag coefficient over the time. Therefore
we can use again an `xy` plot.

### drag coefficient over time

\snippet testcases/FlowAroundCylinder/params_plot.py drag coefficient over time

### lift coefficient over time

\snippet testcases/FlowAroundCylinder/params_plot.py lift coefficient over time

Here we can see the results:

\image html images/cD_over_time.png "drag coefficient over time"
\image html images/cL_over_time.png "lift coefficient over time"

We go on with the pressure difference over time. It is computed by using two 
points which were defined in the tracking table. It has a special `kind = deltaP`.

### pressure difference over time

\snippet testcases/FlowAroundCylinder/params_plot.py pressure difference over time

We get two results. One image with the plot and one text file under `filedata/`
containing the values.

\code
#filaname 	 average 	 min 	 max 	 last 
tracking/Re100/FlowAroundCyl_2D_probe_pressure_p00000.res	2.4689231451	0.0	4.72413262507	2.57346840196
\endcode

\image html images/deltaP_over_time.png "pressure difference over time"

Then, we plot the x and y components of velocity over the length. Therefore, 
we use once again a simple `kind = xy` plot.

### velocityY over length

\snippet testcases/FlowAroundCylinder/params_plot.py velocityY over length

### velocityX over length

\snippet testcases/FlowAroundCylinder/params_plot.py pressure difference over time

Here are the resulting plots:

\image html images/velX_over_length.png "velocity X over length of the channel"
\image html images/velY_over_length.png "velocity Y over length of the channel"

Last but not least we plot the frequency which is needed to compute the Strouhal
number. We use a special `kind = FFT` which means 'Fast Fourier Transform' to do 
so. 

> It is important to refer to velocityY instead of the pressure difference. If 
> we use pressure difference we will get the right solution multiplied by 2.

### frequency

\snippet testcases/FlowAroundCylinder/params_plot.py pressure difference over time

Again we get two results, a textfile and a plot. The plot shows the frequency on 
the x-axis and velocityY on the y-axis.
In the text file, we get the exact value for the frequency which we use to 
compute the Strouhal number.

\image html images/frequency.png "velocity Y over frequency"
Here is the resulting text file:
\code
#filaname 	 Frequency 
# tracking/Re100/FlowAroundCyl_2D_probe_velocity_p00000.res	3.0989387494
\endcode

What we have to compute after that:

* Strouhal number with results of frequency
* maximum values for lift and drag coefficient with results from tracking folder

### computing Strouhal number and lift and drag coefficients

The frequency is used then for computing the Strouhal number in Python for 
example. It is defined as

\code
St = Dia*f/u_mean_phy
\endcode
So we get our result for Strouhal number like:
\code
>>> 3.0989387494*0.1
0.30989387494
\endcode

The maximum of lift and drag coefficient can be computed by a simple python 
script. We call it `get_coeff_max.py` for example.

\include testcases/FlowAroundCylinder/get_coeff_max.py

The results are:

\code
t_0 =  5.64226175439
cD_max =  3.19177160837
cL_max =  0.802261936509
\endcode

Here ends the tutorial for the Flow Around Cylinder test case.

Navigate: \ref tut_steady2D "&larr; Steady Flow Around Cylinder 2D" 
| \ref tut_FlowAroundCyl2D "Test case Flow Around Cylinder 2D"
