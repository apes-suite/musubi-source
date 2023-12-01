title: Scheme Implementation

The concept of schemes provides the user with a bigger flexibility. It is now
possible to run multiple simulations with different layouts on the same mesh.
This is the basis for LBM passive scalar transport simulations in *Musubi*.


## Scheme Definition

The scheme is part of the individual solvers and contains all information about
initial conditions, boundary conditions, the [[mus_scheme_layout_module]], fluid
and flow properties.

## Usage
To start a flow simulation one can add the table species with the following
quantities:

- the scheme header which includes kind, relaxation, layout:
  take a look at [[mus_scheme_header_module]]
- the initial conditions, take a look at [[mus_field_module:mus_set_ic_states]]
  for a list of variables that are supported for the different schemes
- the boundary conditions, take a look at [[mus_bc_header_module:mus_load_bc]]
  for the different cases of possible boundary conditions in *Musubi*
- the fluid quantities. For all supported input variables check [[mus_load_fluid]]

```lua
--! [Scheme identifier]
identify = {
  layout = stencil,
  relaxation = relaxationModel,
  kind = physicsModel
}
--! [Scheme identifier]


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


--! [Boundary conditions]
-- Label is a boundary identifier and it should be same as in seeder
-- configuration.
boundary_condition = {
  {  label = 'north',
     kind  = 'wall'
  },
  {  label = 'south',
     kind  = 'wall'
  },
  {  label = 'west',
     kind  = 'pressure_expol',
     pressure = pressureIn
  },
  {  label = 'east',
     kind  = 'pressure_expol',
     pressure = pressureOut
  }
}
--! [Boundary conditions]


--! [Fluid]
-- For both, incompressible and compressible kinematic viscosity has to be
-- defined. While for the first one, default values are stored for bulk
-- viscosity,, the user has to explicitly give them for compr. fluid.
if (physicsModel == 'fluid') then
  fluid = {
    kinematic_viscosity = nu_phy,
    bulk_viscosity = bulk_visc
  }
else
  fluid = {
      kinematic_viscosity = nu_phy,
    }
end
--! [Fluid]
```
