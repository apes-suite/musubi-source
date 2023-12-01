title: Nonreflecting Boundary Conditions

@warning WORK IN PROGRESS, see also issue #2286 @endwarning

This setup shows the use of non-reflecting boundary conditions.

Here is an example for the definition of a boundary condition of this type:

```lua
  {
    label = 'east',
    kind = 'outlet',
    pressure = 0,
    kappa = 1,
    sigma = 0.1,
    length = 256,
    mach_lat = 0.05 / math.sqrt(1./3.)
  }
```

Look in [[mus_bc_header_module]] for details on configuring boundary conditions
and [[mus_bc_fluid_module]] for their actual implementation.

The configuration file `musubi.lua` provides a complete example for a simulation
setup with this boundary condition:

```lua
{!examples/fluid/application/Nonreflecting_BC/musubi.lua!}
```
