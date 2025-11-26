title: Examples
@warning WORK IN PROGRESS @endwarning

This is a collection of examples that illustrates the usage of the various
capabilities of Musubi.

Musubi is configured via Lua scripts which need to define several variables.
Most of these variables are tables with multiple components. See any of the
specific examples below for a complete configuration.

At least the following need to be defined:

* `sim_control`, see [[tem_simControl_module]], note that there are additional
                 abort options defined for Musubi,
                 see [[mus_abortCriteria_module]]
* `mesh`, see [[treelmesh_module]]
* `physics`, see [[mus_physics_module]]
* `identify`, see [[mus_scheme_header_module]]
* `initial_condition`, see [[mus_flow_module(module):mus_init_byIC]]
* `boundary_condition`, see [[mus_bc_header_module]]

Some other variables may be set to enable optional features or override
defaults:

* `tracking`, see [[mus_tracking_module]]
* `restart`, see [[mus_restart_module]]
* `logging`, see [[tem_logging_module]]

Treelm also provides various general settings that may be specified in the
configuration, see [[tem_general_module]].

Please note that you can include other Lua scripts with
[require](https://www.lua.org/pil/8.1.html).
And you can access table components with a dot notation like `equation.name`.

The Lua script will be executed by Musubi and in the end the defined variables
will be used as configuration for the simulation.

A configuration that shows the typical parameters for a flow simulation is
provided in [Nozzle Flow inside](fluid/application/Nozzle/NOZ_FlowInside).

## Examples ##  {#mus_examples}

The example setups are structured according to the physics supported by Musubi:

- Fluid flows
    - [Weakly compressible flows](fluid/index.md)
    - [Incompressible flows](fluid_incompressible/index.md)
    - [Isothermal acoustics flows](isothermal_acoustics/index.md)
- Multicomponent flows
    - [Gas mixtures](multispecies_gas/index.md)
    - [Liquid mixtures](multispecies_liquid/index.md)
- Advection-Diffusion equations (Scalar transport equations)
    - [Nernst-Planck equation for ionic species transport](nernst_planck/index.md)
- [Poisson equation](poisson/index.md)
- Fluid Particle coupling 
    - [fully resolved model](particles/index.md)
    - [unresolved model](particles/index.md)
