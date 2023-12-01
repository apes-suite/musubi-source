Non-newtonian fluid simulations {#nonNewtonian}
========

# Theory # {#nonNewtonianTheory}

For this type of fluid, its viscosity is not constant anymore, but dependent on
the strain rate.  So the relaxation parameter in LBM needs to be adjusted at
every time step, bacause viscosity is related to relaxation.

In the following two common models are described, namely the Casson and the
Carreau-Yasuda model.

A detailed description can be found in \cite Boyd2007 .

## Common Basis ## {#comBasis}

First the strain rate tensor \f$ S_{\alpha\beta} \f$ needs to be calculated:

\f[
  S_{\alpha\beta} =
    -\frac{3\omega}{2\rho}\sum{i}{}f^{(1)}_{i}\vec{c}_{i\alpha}\vec{c}_{i\beta}
\f]

where \f$ f^{(1)}_{i} \f$ is approxiamated by \f$ f^{neq}_{i} \f$.

Then the second invariant of the strain rate tensor is defined as

\f[
  D_{II} = \sum_{\alpha,\beta=1}^{l}S_{\alpha\beta}
\f]
where \f$ l=3 \f$ in the case of 3D and the shear rate \f$ \dot{\gamma} \f$ (\f$ [\dot{\gamma}]=\textnormal{s}^{-1} \f$) as

\f[
  \dot{\gamma} = 2\sqrt{D}_{II}
\f]

## The Casson Model ## {#nonNewt_Casson}

The Casson model is one of those used to describe the shear thinning behavior
of blood. Its viscosity updating rule is:

\f[
    \nu(\dot{\gamma}) = 
      \frac{1}{\dot{\gamma}} [k_0(c) + k_1(c)\sqrt{\dot{\gamma}}]^2
\f]

where \f$ k_0(c) \f$ and \f$ k_1(c) \f$ are functions of the hematocrit.
One set of these two parameters can be:

- \f$ k_0(c) = 0.1937 \textnormal{ Pa}^{1/2} \f$
- \f$ k_1(c) = 0.055 \textnormal{ (Pa s)}^{1/2} \f$

## The Carreau-Yasuda Model ## {#nonNewt_CarYas}

After getting shear rate, a non-Newtonian model can be applied to update local
dynamic viscosity.  The most common one used in hemodynamics is called
Carreau-Yasuda model, which models for the shear thinning behavior of blood.
The equation is following, where \f$ a, n \f$ and \f$\lambda\f$  are empirically determined
constant parameters.  \f$a\f$ and \f$n\f$ are dimensionless; \f$\lambda\f$ has units of s.

\f[
    \mu(\dot{\gamma}) = 
      \mu{\infty} + (\mu_0-\mu_{\infty})[1+(\lambda\dot{\gamma})^a]^{(n-1)/a}
\f]

one set of parameters is:

- \f$ \nu_{0} = 0.16 \textnormal{ Pa s} \f$
- \f$ \nu_{\infty} = 0.0035 \textnormal{ Pa s} \f$
- \f$ \lambda = 8.2 \textnormal{ s} \f$
- \f$ a = 0.64 \f$
- \f$ n = 0.2128 \f$

Next convert dynamic viscosity to kinematic viscosity:

\f[
  \nu = \frac{\mu}{\rho}
\f]

As the viscosity is in physical units, it has to be convert into LB units:
\f[
  \nu_{LB} = \nu \cdot \frac{dt}{dx^2}
\f]

finally, the relaxation parameter - \f$\omega\f$ can be obtained from \f$\nu_{LB}\f$:

\f[
  \omega = \frac{2}{6\nu_{LB} + 1}
\f]
