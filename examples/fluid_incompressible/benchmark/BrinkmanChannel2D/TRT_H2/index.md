title: Brinkman Channel 2D Example with TRT operator and H2 Brinkman force

# Brinkman Channel 2D Example with TRT operator and H2 Brinkman force ## {#ex_brinkman_channel2d_trt_h2}

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

  $$u(y) = k_1 e^{\sqrt{F_0} y} + k_2 e^{-\sqrt{F_0} y}$$

with $k_1 = \tfrac{u_m}{e^{\sqrt{F_0} h} - e^{-\sqrt{F_0} h}}$ and $k_2 = -k_1$. 

The analytical solution is written in `func.lua`. The values of $F_0$ is and 
other parameters are set in `params.lua`. 

In this example, TRT collision model is adopted. Meanwhile, the source term is expanded 
up to the second Hermite order. Separating it according to lattice parity gives

  $$ S_i=S_i^{-}+S_i^{+} $$
where the first-order, antisymmetric contribution is

$$
  S_i^{-}=w_i\omega^{-}
  \left(
  \frac{1}{\omega^{-}}-\frac{\Delta t}{2}
  \right)
  \frac{\mathbf{c}_i\cdot\mathbf{F}_{\mathrm{tot}}}{c_s^2},
$$
and the second-order, symmetric contribution is
$$
  S_i^{+}=w_i\omega^{+}
  \left(
  \frac{1}{\omega^{+}}-\frac{\Delta t}{2}
  \right)
  \left[
  \frac{
  (\mathbf{c}_i\cdot\mathbf{u})
  (\mathbf{c}_i\cdot\mathbf{F}_{\mathrm{tot}})
  }{c_s^4}-\frac{\mathbf{u}\cdot\mathbf{F}_{\mathrm{tot}}}{c_s^2}
  \right].
$$