title: Advection-diffusion of a Gaussian Hill
@warning WORK IN PROGRESS @endwarning

# Advection-diffusion of a Gaussian Hill # {#eg_GPP}

In this example, we will simulate the time evolution of a Gaussian hill 
in an infinite space with periodic boundaries.
The initial form of the Gaussian pulse gives the form: 
  $$ C(\mathbf{x}, t=0) = C_0 \exp\left(-\frac{(\mathbf{x} 
  - \mathbf{x}_0)^2}{2\sigma_0^2}\right) + C_{base} $$
The analytical solution with time is:
    $$ C(\mathbf{x}, t) = \frac{\sigma_0^2}{\sigma_0^2 + \sigma_D^2} C_0 \exp \left( 
    -\frac{(\mathbf{x} - \mathbf{x}_0)^2}{2(\sigma_0^2 + \sigma_D^2)} \right) + C_{base} $$

In our simulation, $\sigma_0=40$ is used as a default value. A 1000*1000 2D mesh is recognized
 here as an infinite area, and the D2Q9 layout is used. Boundaries from all directions are 
 periodic. As the area is much larger than the $ \sigma_0 $, the effect from the boundaries 
 is neglected. Results of t=200 are used to calculate the diffusion factor.

![The comparison of concentration between simulation and analytical solution](media/compare_sim_ana_C.jpg)
![The error of diffusion factor when u\_bg=0](media/err_numerical_diffusion.jpg)
![Absolute error of D with u\_bg and order of equilibriums](media/d_abs_err_tau0.501.jpg)
![Absolute error of D with u\_bg and order of equilibriums](media/d_abs_err_tau0.503.jpg)
![Absolute error of D with u\_bg and order of equilibriums](media/d_abs_err_tau0.8.jpg)
![Absolute error of D with u\_bg and order of equilibriums](media/d_abs_err_tau2.jpg)
![Absolute error of D with u\_bg and order of equilibriums](media/d_abs_err_tau5.jpg)
![The model error for first equilibrium](media/d_1st_model_error.jpg)

The objectives of this example is to introduce how to:
* Simulate time evolution of the advection-diffution process of a 2D 
Gaussian Hill.
* Compare the profile between first and second order bgk equilibria
* Compare the profile between bgk and trt schemes
* Compute the diffusion parameter and compare it to the theoritical value for different 
background velocities

