title: Brinkman Channel 2D Example with TRT operator

# Brinkman Channel 2D Example with TRT operator ## {#ex_brinkman_channel2d_trt}

Brinkman flow in a channel is a classical benchmark problem for validating 
the Navier--Stokes--Brinkman (NSB) solver, as it combines viscous diffusion 
and Darcy drag effects characteristic of porous media. 
The problem considers a steady laminar flow through a channel filled with a porous matrix, 
where the additional Darcy drag term significantly alters the velocity distribution. 
The governing one-dimensional Brinkman equation reads

$$\frac{\partial^2 u_x}{\partial y^2} = F_0 u_x, \label{nsb_eq_1}$$

where $u_x$ is the velocity in the streamwise ($x$) direction, and 
$F_0 = \tfrac{\mu_0}{k\mu}$ is the normalized Darcy coefficient, with 
$\mu$ the fluid viscosity, $\mu_0$ the effective viscosity, and $k$ the permeability. 

The boundary conditions are imposed as a no-slip condition at the lower wall 
and a symmetric velocity condition at the upper wall,

  $$u(y=0) = 0, \qquad u(y=h) = u_m. \label{nsb_bc}$$

The analytical solution of this boundary-value problem is

  $$u(y) = k_1 e^{\sqrt{F_0} y} + k_2 e^{-\sqrt{F_0} y}, \label{nsb_solution}$$

with $k_1 = \tfrac{u_m}{e^{\sqrt{F_0} h} - e^{-\sqrt{F_0} h}}$ and $k_2 = -k_1$. 
The analytical solution is written in `func.lua`.

The values of $F_0$ is set in `args.lua`. Readers are encouraged to modify it
to see how the Brinkman flow profile changes with different $F_0$ values.

Readers can run the workflow by executing the script `multi_f0.sh` to simulate 
multiple $F_0$ values in one go. 
The post-process of results `plot_contour.py` is to generate contour plots of 
the velocity profiles for different $F_0$ values.
To run the python script, make sure you have `matplotlib` and `numpy`
installed in your Python environment.

The example uses the following settings in `params.lua`:

$$
  resolution = 2 ^ 3
  tmax       = 20000 * resolution^2
  interval   = 2000 * resolution^2
$$

The reference profiles for different $F_0$ values are provided in
`media/`.