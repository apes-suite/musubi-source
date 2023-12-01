title: Multispecies

Navigate: [&larr; Features](index.html) | [Multispecies test case tutorial](../../page/examples/tutorials/tut_09_mus_multilevel.html)

Multispecies approach implemented in the code is based on the paper
"Multi-species Lattice Boltzmann Model and Practical Examples. Short Course
material Pietro Asinari PhD."

Equlibrium distribution function for multispecies is given in the paper as
 \[  f^{\sigma(eq)}_\alpha(\rho,\mathbf{u^*_{\sigma}}) =
 \omega_\alpha\cdot\Big[s^{\sigma}_{\alpha}+\frac{1}{c^2_s}
 (\mathbf{e}_\alpha \cdot\mathbf{u^*_{\sigma}})+\frac{1}{2c^4_s}(\mathbf{e}_\alpha
 \cdot\mathbf{u^*_{\sigma}})^2-\frac{1}{c^2_s}(\mathbf{u^*_{\sigma}}\cdot
 \mathbf{u^*_{\sigma}})\Big] \]

 where,
 \(s^{\sigma}_0 = (9-5 \phi^\sigma)/4 \)
 \(s^{\sigma}_{\alpha} = \phi^\sigma\) for \(1\leq\alpha\leq8\) and
 \(\phi^\sigma=min_\varsigma(m^\varsigma)/m^\sigma\),
 \( p_\sigma=\rho_\sigma\phi^\sigma/3 \)
 \(m^\sigma\) is molecular weight for species \(\sigma\).

 \(u^*_{\sigma} \) is given as
 \[\mathbf{u}^*_{\sigma} = \mathbf{u}_\sigma + \sum_{\varsigma}
 {\frac{m^2}{m^\sigma m^{\varsigma}}\frac{B_{\sigma \varsigma}}{B_{mm}}
 x_\varsigma (\mathbf{u}_\varsigma-\mathbf{u}_\sigma) }  \]
 \(x_\sigma= \rho_\sigma/\rho\)

 Relaxation time is given as

 \( \lambda_\sigma = \frac{pB_{mm}}{\rho } \)

# Multispecies: Variable Transformation

A semi-implicit Lattice Boltzmann equation is given as,

\[ f^{\sigma,+}(\mathbf{x}+\mathbf{e},t+1) = f^\sigma(\mathbf{x},t) + (1-\frac{1}{2})
\lambda^\sigma[f^{\sigma(eq)}-f^\sigma]
+ \frac{1}{2} \lambda{^\sigma,+}[f^{\sigma(eq),+}-f^{\sigma,+}]  \]

Variable transformation presented in above paper involves three steps

## Step 1. Transforming f -> g

\( g^\sigma = f^\sigma-\frac{1}{2}\lambda^\sigma[f^{\sigma(eq)}-f^\sigma] \)

## Step 2. Stream and Collide i.e g -> g^+

\( g^{\sigma,+} = g^\sigma + \frac{\lambda^\sigma}{1+\frac{1}{2}\lambda^\sigma}
[f^{\sigma(eq)}-g^\sigma] \)

## Step 3. Back Transormation to f i.e g^+ -> f^+

\[ f^{\sigma,+} = \frac{g^{\sigma,+} + \frac{1}{2}\lambda^{\sigma,+}
f^{\sigma(eq),+}}{1+\frac{1}{2}\lambda^{\sigma,+}} \]

In back tranformation, to compute feq we need \(\rho\) and \( \mathbf{u}_\sigma \)
\( \rho \) can be computed directly from g
\(  \rho^+_\sigma = <g^{\sigma,+}> \)

where as the \( \mathbf{u}_\sigma \) computed by solving the linear
system of equation given below

\[  <\mathbf{e}_i g^{\sigma,+}> = \Bigg[ 1 + \frac{1}{2} \lambda^{\sigma,+} \sum_\varsigma
(\chi_{\sigma \varsigma} x^+_\varsigma) \Bigg] q^+_{\sigma,i} - \frac{1}{2}
\lambda^{\sigma,+}x^+_\sigma \sum_\varsigma (\chi_{\sigma \varsigma}q^+_{\varsigma,i})
\]

where, \( \chi_{\sigma \varsigma} \) is,
 \[ \chi_{\sigma \varsigma} = \frac{m^2}{m_\sigma m_\varsigma}
 \frac{B_{\sigma \varsigma}}{B_{\sigma \sigma}} \]

@todo use for optimized routine
```fortran
       !West
       u_n(1) = - uxstar(s)
       !south
       u_n(2) =             - uystar(s)
       !bottom
       u_n(3) =                         - uzstar(s)
       !east
       u_n(4) =   uxstar(s)
       !north
       u_n(5) =               uystar(s)
       !top
       u_n(6) =                           uzstar(s)
       !bottom south
       u_n(7)
       !top south
       u_n(7)
       !bottom north
       u_n(7)
       !top north
       u_n(7)
       !bottom west
       u_n(7)
       !bo
       u_n(7)
       u_n(7)
       u_n(7)
       u_n(7)
```

Navigate: [&larr; Features](index.html) | [Multispecies test case tutorial](../../page/tutorial/taylorgreen/index.html)
