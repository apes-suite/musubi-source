title: Advection-diffusion-reaction benchmark: Semi-infinite 1D domain

# Advection-diffusion-reaction benchmark: Semi-infinite 1D domain ## {#ex_semi_infinite_1d}

The movement of a tracer in a semi-infinite column is a classical benchmark problem in transport phenomena and porous media hydrology.

References:
  Jacob Bear. Dynamics of Fluids in Porous Media. New York: Elsevier,
  1972. isbn: 978-0-486-65675-5.

It serves to validate the accuracy of the LBM framework for solving advection-diffusion-reaction (ADR) equations. 
The setup considers a tracer initially injected at the inlet of a semi-infinite column, 
which then advects, diffuses, and reacts along the flow direction. 

The governing ADR equation is

$$
  \partial_t c + U \,\partial_x c = D \,\partial_{xx} c - \lambda c, \label{adr}
$$

where $c(x,t)$ is the concentration, $U$ the advective velocity, $D$ the molecular diffusion coefficient, 
and $\lambda$ the reaction rate. 

The boundary and initial conditions are

$$
  c(x,0) = 0, \qquad c(0,t) = C_0, \qquad c(\infty,t) = 0, \label{adr_bc}
$$

corresponding to a constant concentration injection at the inlet. 

The analytical solution for this problem is given in \cite{Bear1972a} as

$$
\begin{split}
  C(x,t) = \frac{C_0}{2} \exp\!\left(\frac{Ux}{2D}\right)
  \Bigg[
  & \exp(-\beta x)\, \mathrm{erfc}\!\left(\frac{x - \sqrt{U^2 + 4 \lambda D}\, t}{2\sqrt{Dt}}\right) \\
  +\, & \exp(\beta x)\, \mathrm{erfc}\!\left(\frac{x + \sqrt{U^2 + 4 \lambda D}\, t}{2\sqrt{Dt}}\right)
  \Bigg],
\end{split}
\label{adr_solution}
$$

with $\beta^2 = \tfrac{U^2 + 4\lambda}{4D}$ and $\mathrm{erfc}(x) = 1 - \mathrm{erf}(x)$ 
the complementary error function. 

In this benchmark, the simulation parameters are set in `params.lua`, and the analytical solution is implemented in `func.lua`. 
To view the profile of the concentration along the column at a specific time, run the example and plot with GNUplot:

```gnuplot
gnuplot -p -e "
    plot $(printf "'%s' " tracking/*.res) using 1:4 title 'calculated' with points,
         $(printf "'%s' " tracking/*.res) using 1:5 title 'analytical' with lines
"
```
This will generate a plot comparing the calculated concentration profile from the LBM simulation against the analytical solution at the specified time.

If your parameters are set as in the provided `params.lua`, you should observe the same results as shown in `media/reference.png`.