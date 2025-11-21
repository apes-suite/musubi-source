title: Advection-diffusion of a Gaussian Hill

# Advection-diffusion of a Gaussian Hill # {#eg_GPP}

In this example, we will simulate the time evolution of a Gaussian hill
in an infinite space with periodic boundaries.
The initial form of the Gaussian pulse gives the form:
  $$ C(\mathbf{x}, t=0) = C_0 \exp\left(-\frac{(\mathbf{x}
  - \mathbf{x}_0)^2}{2\sigma_0^2}\right) + C_{base} $$
The analytical solution with time is:
    $$ C(\mathbf{x}, t) = \frac{\sigma_0^2}{\sigma_0^2 + \sigma_D^2} C_0 \exp \left(
    -\frac{(\mathbf{x} - \mathbf{x}_0)^2}{2(\sigma_0^2 + \sigma_D^2)} \right) + C_{base} $$

In our simulation, we define $nelem = 500$ which generates a 1000\*1000 2D mesh. The mesh is recognized  
 as an infinite area for $\sigma_0=40$. The D2Q9 layout is used. Boundaries from all directions are
 periodic. As the area is much larger than the $ \sigma_0 $, the effect from the boundaries
 is neglected. Results of t=200 are used to calculate the diffusion factor.

![The comparison of concentration between simulation and analytical solution](media/compare_sim_ana_C.png)

In <tt>./musubi.lua<tt>, define 'shape' in 'tracking' block as

```
shape = {
  -- kind = 'all'
  kind  = 'canoND',
  object  = {
    origin = {-nelem-1, 0, 0},
    vec = { {2.0*(nelem+1), 0., 0.0} }
  }
},
```

To generate a midline profile comparison between the analytical and computed solution, run <tt>profile.sh<tt> 
to create the plot.

![The error of diffusion factor when u\_bg=0](media/err_diffusion_factor.png)

In <tt>./musubi.lua<tt>, define 'shape' in 'tracking' block as

```
shape = {
  kind = 'all'
  -- kind  = 'canoND',
  -- object  = {
  --   origin = {-nelem-1, 0, 0},
  --   vec = { {2.0*(nelem+1), 0., 0.0} }
  -- }
},
```

To generate the above plot, simply run <tt>test_stability.sh<tt>. Notice that when $ \tau = 1 $,
the error of computed diffusion factor is exactly 0.

![Absolute error of D with u\_bg and order of equilibriums](media/d_abs_err_tau0.501.jpg)

In <tt>test_order.sh<tt>, define $ \tau=0.501 $, and run <tt>test_order.sh. The chart shows
the error of diffusion factor changing with the background velocity. The following 
charts can be obtained with different relaxation time $ \tau $.

![Absolute error of D with u\_bg and order of equilibriums](media/d_abs_err_tau0.503.jpg)
![Absolute error of D with u\_bg and order of equilibriums](media/d_abs_err_tau0.8.jpg)
![Absolute error of D with u\_bg and order of equilibriums](media/d_abs_err_tau2.jpg)
![Absolute error of D with u\_bg and order of equilibriums](media/d_abs_err_tau5.jpg)

The objectives of this example is to introduce how to:
* Simulate time evolution of the advection-diffusion process of a 2D
Gaussian Hill.
* Compare the profile between first and second order bgk equilibria
* Compute the diffusion parameter and compare it to the theoritical value for different
background velocities

