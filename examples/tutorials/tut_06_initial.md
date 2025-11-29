title: Initial Conditions
@warning WORK IN PROGRESS @endwarning

Navigate: [&larr; Restart](tut_05_restart.md)
| [Overview](index.md)
| [Abort Criteria &rarr;](tut_07_convergence.md)

# Initial Conditions
Definition of Initial Conditions

In this tutorial, we cover the definition of initial conditions.
They can be defined either in lattice units or in physical units.
Lattice units can be confusing, especially to people just starting in this field.
The relationships between lattice and physical units can be found in this paper:

> [unit conversion by Jonas Latt](http://wiki.palabos.org/_media/howtos:lbunits.pdf).

In this tutorial, a two-dimensional plain channel is set up.
Not only are the boundaries specified to obtain a defined pressure drop over
the channel length, but also are the initial conditions set in a consistent
manner.

@note To run this test case use the generic channel test case and set in
`seeder.lua`:
```lua
case2d = true
usePeriodic = false
qValues = false
useObstacle = false
fixHeight = true
useRefine = false
```
@endnote

For the viscous, laminar two-dimensional plain channel flow, an analytical
solution of the incompressible Navier-Stokes equation can be derived. From the
analytical solution the pressure drop, the velocity profile and the shear stress
distribution can be computed.

Before starting, we need to define the flow regime and physical reference values.

```lua
-- Density of the fluid [kg/m^3]
rho0_phy = 1.0
-- Reynolds number of the flow
Re = 60
-- Inflow velocity [m/s]
vel_phy = 1.0
-- Kinematic viscosity of the fluid calculated from Re [m^2/s]
nu_phy = vel_phy * height_phy / Re
```

For the Lattice_Boltzmann simulation, basic simulation parameters such as a
lattice velocity, the timestep and the lattice viscosity need to
be specified.

```lua
if (scaling == 'acoustic') then
  -- Lattice Mach number
  Ma_lat = 0.05
  -- Lattice velocity
  vel_lat = Ma_lat * cs_lat
  -- Physical timestep computed from physical and lattice velocity
  dt  = vel_lat / vel_phy*dx
  -- Lattice viscosity
  nu_lat  = nu_phy*dt /dx /dx
  -- Relaxation parameter
  omega   = 1.0/(3.0*nu_lat + 0.5)
else
  -- Diffusive scaling
  -- Relaxation parameter
  omega   = 1.7
  -- Lattice viscosity
  nu_lat  = ( 1.0/omega - 0.5 ) / 3.0
  -- Physical timestep computed from physical and lattice velocity
  dt      = nu_lat/nu_phy*dx*dx
  -- Lattice velocity
  vel_lat     = vel_phy*dt/dx
  -- Lattice Mach number
  Ma_lat = vel_lat * cs_lat
end
--------------------------------------------------------------------------------

--! [Reference LB values]
-- Square of lattice speed of sound
cs2      = 1.0/3.0
-- Lattice density
rho0_lat = 1.0
-- Zero lattice density
rho0_lat0 = 0.0
--! [Reference LB values]
```

Depending on the model used, the reference pressure differs. For the
incompressible model, the reference pressure is 0, while for the compressible
model the reference pressure is `rho0*cs2`

```lua
--! [Reference pressure dependent on fluid kind]
if physicsModel == 'fluid_incompressible' then
  p0 = 0.0
else
  p0 = rho0_lat*cs2*dx*dx/dt/dt
end
--! [Reference pressure dependent on fluid kind]
```

From the solution of the Navier-Stokes equation, the following relations for the
velocity distribution across the channel height can be obtained

```lua
--! [Velocity function]
function velX(x,y,z)
  velX_phy = vel_phy * ( 1.0 - ( 2.0*y/height_phy )^2 )
  return velX_phy
end
--! [Velocity function]
```

Similarly for the pressure drop along the channel length

```lua
--! [Pressure function]
function pressureRef(x,y,z)
  press_drop = vel_phy*8.*nu_phy*rho0_phy/height_phy^2*length
  return p0 + press_drop*0.5 - press_drop/length*x
end
--! [Pressure function]
```

and the shear stress across the channel height

```lua
--! [Shear stress function]
function Sxy(x,y,z)
  tauxy= -nu_phy*rho0_phy*8./height_phy^2*vel_phy*y
  S_xy = tauxy/nu_phy/rho0_phy
  return S_xy
end
--! [Shear stress function]
```

Now the physics table establishes the connection between the lattice reference
values and the physical values and gives *Musubi* means of transferring between
these two unit systems. See [[mus_physics_module]] for more information.

```lua
physics = {
  dt = dt,
  rho0 = rho0_phy
}
```


For the Lattice_Boltzmann algorithm, a reference density and the kinematic
viscosity (for compressible also the bulk viscosity) need to be defined.
See [[mus_fluid_module]] for more information.

```lua
  fluid = {
    kinematic_viscosity = nu_phy,
    bulk_viscosity = bulk_visc
  }
```

Now the initial conditions for each element in the simulation domain is defined
by setting each physical quantity and connecting it to a lua function, which we
defined above.

```lua
--! [Initial conditions]
initial_condition = { pressure  = ic_pressure,
                      velocityX = ic_velX,
                      velocityY = 0.0,
                      velocityZ = 0.0,
                      Sxx = ic_Sxx,
                      Syy = ic_Syy,
                      Sxy = ic_Sxy
                      }
--! [Initial conditions]
```
@note The whole code of musubi.lua is shown in the [chapter_02](tut_02_mus_toolchain.md).

Next chapter: [Abort Criteria &rarr;](tut_07_convergence.md)
