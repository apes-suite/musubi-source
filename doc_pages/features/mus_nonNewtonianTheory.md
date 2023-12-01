title: Non-Newtonian Fluid Simulations


# Theory # {#nonNewtonianTheory}

For this type of fluid, its viscosity is not constant anymore, but dependent on
the strain rate. So the relaxation parameter in LBM needs to be adjusted at
every timestep, bacause viscosity is related to relaxation.

In the following two common models are described, namely the Casson and the
Carreau-Yasuda model.

A detailed description can be found in \cite Boyd2007 .

## Common Basis ## {#comBasis}

First, the strain rate tensor \( S_{\alpha\beta} \) needs to be calculated:

\[
  S_{\alpha\beta} =
    -\frac{3\omega}{2\rho}\sum{i}{}f^{(1)}_{i}\vec{c}_{i\alpha}\vec{c}_{i\beta}
\]

where \( f^{(1)}_{i} \) is approxiamated by \( f^{neq}_{i} \).

Then, the second invariant of the strain rate tensor is defined as

\[
  D_{II} = \sum_{\alpha,\beta=1}^{l}S_{\alpha\beta}
\]
where \( l=3 \) in the case of 3D and the shear rate \( \dot{\gamma} \) (\( [\dot{\gamma}]= {s}^{-1} \)) as

\[
  \dot{\gamma} = 2\sqrt{D}_{II}
\]

## The Casson Model ## {#nonNewt_Casson}

The Casson model is one of those used to describe the shear thinning behavior
of blood. Its viscosity updating rule is:

\[
    \nu(\dot{\gamma}) = 
      \frac{1}{\dot{\gamma}} [k_0(c) + k_1(c)\sqrt{\dot{\gamma}}]^2
\]

where \( k_0(c) \) and \( k_1(c) \) are functions of the hematocrit.
One set of these two parameters can be:

- \( k_0(c) = 0.1937  { Pa}^{1/2} \)
- \( k_1(c) = 0.055  { (Pa s)}^{1/2} \)

## The Carreau-Yasuda Model ## {#nonNewt_CarYas}

After getting shear rate, a non-Newtonian model can be applied to update local
dynamic viscosity. The most common one used in hemodynamics is called
Carreau-Yasuda model, which models for the shear thinning behavior of blood.
The equation is following, where \( a, n \) and \(\lambda\)  are empirically
determined constant parameters.  \(a\) and \(n\) are dimensionless;
\(\lambda\) has units of s.

\[
    \mu(\dot{\gamma}) = 
      \mu{\infty} + (\mu_0-\mu_{\infty})[1+(\lambda\dot{\gamma})^a]^{(n-1)/a}
\]

one set of parameters is:

- \( \nu_{0} = 0.16  { Pa s} \)
- \( \nu_{\infty} = 0.0035  { Pa s} \)
- \( \lambda = 8.2  { s} \)
- \( a = 0.64 \)
- \( n = 0.2128 \)

Next, convert dynamic viscosity to kinematic viscosity:

\[
  \nu = \frac{\mu}{\rho}
\]

As the viscosity is in physical units, it has to be convert into LB units:
\[
  \nu_{LB} = \nu \cdot \frac{dt}{dx^2}
\]

Finally, the relaxation parameter - \(\omega\) can be obtained from \(\nu_{LB}\):

\[
  \omega = \frac{2}{6\nu_{LB} + 1}
\]
