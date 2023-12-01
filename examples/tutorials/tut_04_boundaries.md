title: Boundary Conditions
@warning WORK IN PROGRESS @endwarning

Navigate: [&larr; Tracking](tut_03_tracking.html)
| [Overview](index.html)
| [Restart &rarr;](tut_05_restart.html)

# Boundary Conditions

In this tutorial we use the 2D channel case to explain how to define boundary
conditions that are needed for the simulation run.

@note To run this case use the generic channel test case and set in
`seeder.lua`:
```lua
case2d = true
usePeriodic = true
qValues = false
useObstacle = false
fixHeight = true
useRefine = true
```
@endnote

@note Multi-level meshes at boundaries lead to fluctuation in the boundary layer.
Also avoid single refinement layers since they could lead to crashes.
At least four refined elements are recommended.


## Where the Information for the Boundaries are Taken From

In *Seeder*, boundaries are defined by `spatial_object` attribute kind "boundary".
The output of *Seeder* is normally placed in the `./mesh` folder. There you can
find files called `bnd.lua` and `bnd.lsb`.
`bnd.lua` file contains basic information about boundaries like number of
boundaries in a mesh and list of boundary labels in ascii format.
`bnd.lsb` file contains boundary IDs for 26 directions for each fluid element
which has boundary neighbor in binary format.
These two files are important for *Musubi* in case of defining the boundary
conditions.

If you open `bnd.lua` file, you will probably find something like this:
```lua
nSides = 26
nBCtypes = 4
bclabel = {
    'south',
    'north',
    'east',
    'west'
}
```

In the example, there are four different boundaries we have set up with *Seeder*.
In this case they have the labels `north`, `south`, `east` and `west`. This
gives you a feeling of how to view this channel.

These boundaries are defined in the mesh, but it is not clear which function
each boundary has at the moment. So this is the part you have to do in *Musubi*.
As an example, we want to simulate a channel. A channel in 2D needs two walls,
one inlet and one outlet.

> If you would like to know in detail how the boundaries are defined you might
> have a look at [[tem_bc_module]] and [[tem_load_bc_state]].

## How to Define Boundary Conditions in Musubi?

Now, we will have a look on the `boundary_conditions` table in *Musubi* for the
channel test case.

```lua
boundary_condition = {
{  label = 'north',
   kind = 'wall' },
{  label = 'south',
   kind = 'wall' },
{  label = 'west',
   kind = inlet_kind, -- local variable inlet_kind
   pressure= 'pressureIn'}, -- refers to pressure defined in variable table
                            -- (see last tutorial)
{  label = 'east',
   kind = outlet_kind, -- local variable outlet_kind
   pressure= 'pressureOut'}} -- local variable from variable table
```

In this basic example you can see the function of each boundary as outlined
above.

### boundary_condition Table

In `boundary_condition` table, for each boundary label from *Seeder*, a kind
must be defined in *Musubi* which defines what to do with that boundary. The
order of boundary definition in *Musubi* does not depend on the order in
*Seeder* (bnd.lua). It just has to exist in both files.
```lua
boundary_condition = {
...
```

### Label

For every boundary which you have created in *Seeder*, you have to set its status.
Therefore, you call the boundary with the exact name (`label = ...`) that you
can see in the bnd.lua file and give it a certain **kind** that is explained in
the next section.
```lua
boundary_condition = {
  label = 'south',
...
```

### Kind

You can choose between some basic boundary kinds. They define the use of the
boundary in the simulation run. Some of these kinds are described below.

* **wall**

  A wall means that the fluid is not able to penetrate through this boundary.
  It has to regard the wall as an obstacle. Moreover, the wall is seen as a
  **no-slip** boundary. If you would like to observe slip as well you have to
  use the separate kind `slip_wall`.
```lua
boundary_condition = {
  label = 'south',
  kind = 'wall'
}
```

* **slip_wall**

  If slip shall be defined as well, you will have to set `kind` to
  `kind = slip_wall`. Slip means that the normal and the tangential velocity in
  normal direction equal zero. The pressure gradient along the normal direction
  is equal to zero as well. The degree of slip can be defined by the
  multiplication of a slip-factor called `fac` and the velocity. Special cases
  are on the one hand `fac = 1` which means that there is free-slip or full-slip
  and on the other hand, there is `fac = 0` which is used for no-slip. This case
  is the default case for the `kind = wall` which is mentioned above.
```lua
boundary_condition = {
  label = 'north',
  kind = 'slip_wall',
  fac = 0.4
}
```
  > More information can be found in the [[mus_bc_fluid_wall_module]].

* **wall_libb**

There is a possibility to make your simulation more efficient: Instead of
making use of a higher refinement level, you could use linear interpolation.
The obstacle in the fluid will have a higher resolution if the distance
between the barycentric centre position of the fluid element to the obstacle
is calculated which is "q". The total distance between the center positions of
each, the fluid element's and the obstacle's is `q=1`. There is a case
differentiation between `q<0.5` and `q>=0.5`.

  For \(q<\frac{1}{2}\) this formula is used:
    \[
    {\displaystyle f_{i^{\prime}}(\mathbf{r}_{l},t+1)}
    {\displaystyle =2qf_{i}^{c}(\mathbf{r}_{l},t)+(1-2q)f_{i}^{c}(\mathbf{r}_{l}-\mathbf{c}_i{i},t)}
    \]
  For \(q\geq\frac{1}{2}\) this formula is used:
    \[
    {\displaystyle f_{i^{\prime}}(\mathbf{r}_{l},t+1)}
    {\displaystyle =\frac{1}{2q}f_{i}^{c}(\mathbf{r}_{l},t)+\frac{2q-1}{2q}f_{i^{'}}^{c}(\mathbf{r}_{i},t)}
    \]

Therefore, *Seeder* has to be configured. For each spatial_object which has
boundary kind, `calc_dist = true` has to be added to the attributes right
behind kind and label. Here is the `seeder.lua` block before the changes:

```lua
  table.insert(spatial_object,  {
    attribute = {
      kind = 'boundary',
      label='north'
    },
    geometry = {
      kind = 'canoND',
      object = {
        origin = {-length*0.5, height*0.5+dxDash, -length*0.5},
        vec = {{length, 0.0, 0.},
              {0.,0.0, length}}
      }
    }
  })
  table.insert(spatial_object,  {
    attribute = {
      kind = 'boundary',
      label='south'
    },
    geometry = {
      kind = 'canoND',
      object = {
        origin = {-length*0.5,-height*0.5-dxDash, -length*0.5},
        vec = {{length, 0.0, 0.},
              {0.,0.0, length}}
      }
    }
  })
```
  >  In the `spatial_object` tables for the boundaries, you have to add to the
  >  attributes `calc_dist = true`.

  Then there should be written as an example for one boundary:
```lua
table.insert(spatial_object,  {
        attribute = {
          kind = 'boundary',
          level = maxLevel, --local variable
          label = stlLabel, --local variable
          calc_dist = true,
        },
        geometry = {
          kind = 'stl',
          object = { filename = stlfile --local variable
  --          { origin = {-length*0.3,-0.01*height,0.},
  --            radius = radius }
          }
        },
       transformation = {
          deformation =  2.0,
          translation =  {-2., 0., 0. }
        }
        })
```
  After that, you are able to set the kind of the boundary in `musubi.lua` to
  `wall_libb` instead of `wall`. You have to use the same syntax
  as shown above, otherwise it will not work.
```lua
  boundary_condition = {
    label = 'south',
    kind = 'wall_libb',
  }
```

  > The information about linear interpolation are taken from M. Bouzidi,
  > M. Firdaouss, and P. Lallemand, "Momentum transfer of a Boltzmann-lattice
  > fluid with boundaries," Physics of Fluids, vol. 13, no. 11,
  > pp. 3452–3459, Nov. 2001 Equations 5 (linear interpolation).
  > More information about no-slip wall linear interpolation can be found in
  > [[mus_bc_fluid_wall_module]].

### Velocity Boundary Conditions

* **velocity_eq**

`kind = velocity_eq` makes use of the equilibrium function concerning
the Lattice Boltzmann Method (LBM). This function gets the density (rho) from
the fluid. The velocity (u) has to be defined for each coordinate x, y and z.
Therefore you can place the definition of the velocity after the kind of the
boundary. Example:
```lua
boundary_condition = {
{ label = 'inlet',
  kind = 'velocity_eq',
  velocity = 'inlet_vel' }
```

> You can get more information about the equilibrium function currently in
> the documentation for the subroutine [[velocity_eq]].

* **velocity_bounceback**

In addition to that you can set `kind = velocity_bounceback`. You can imagine
bounceback like this: A fluid particle gets near the boundary. If it reaches the
boundary, the particle will bounce back in the same angle as it gets there.
For this kind, you have to set the velocity values, too.

```lua
boundary_condition = {
  { label = 'west',
    kind = 'velocity_bounceback',
    velocity = 'inlet_vel' }
```

> More details can be found in the documentation for subroutine [[velocity_bounceback]].

* **mfr_bounceback**

Like for `velocity_bounceback` there is an "Inlet Velocity Bounce Back" boundary
condition. But in this case, `mass flow rate` (`mfr`) is used as an input
argumet as well. It is used like that:

```lua
boundary_condition = {
  { label = 'inlet',
    kind = 'mfr_bounceback',
    massflowrate = 0.1 }
```

> For more information visit [[mfr_bounceback]].

* **mfr_eq**

In this case, the mass flow rate is used for the equilibrium boundary condition
and the velocity is taken from the configuration file.
```lua
boundary_condition = {
  { label = 'inlet',
    kind = 'mfr_eq',
    massflowrate = 0.1 }
```

> The corresponding Documentation can be found in [[mfr_eq]].

### Pressure Boundary Conditions

* **pressure_expol**

The variable values are extrapolated during the simulation.
```lua
boundary_condition = {
  { label = 'east',
    kind = 'pressure_expol',
    pressure = 2.0 }
```

> Detail information can be found in [[pressure_expol]].

* **pressure_antiBounceBack**

This is the outlet pressure anti-bounce back boundary condition kind. The
velocity is extrapolated by two of its neighbors. The pressure has to be given
as well.
```lua
boundary_condition = {
  { label = 'east',
    kind = 'pressure_antiBounceBack',
    pressure = 2.0 }
```

> More information [[pressure_antiBounceBack]].

* **pressure_eq**

The incoming densities are set to the equilibrium distribution with macroscopic
velocity and pressure.
```lua
boundary_condition = {
  { label = 'east',
    kind = 'pressure_eq',
    pressure = 1.0 }
```

> Detail information in the [[pressure_eq]].

* **symmetry**

For symmetric test cases this boundary condition can be used to mirror the
flow-parameters at this border. This allows to reduce the computational effort.
The missing part can be easily mirrored when post-processing.

@warning Be aware: The symmetric boundary has to be the last boundary created
by *Seeder*!
@end warning

```lua
boundary_condition = {
  { label = 'east',
    kind = 'symmetry'
  }
```


### Physical Conditions

If you define the inlet and the outlet boundary for the shape, you will have to
give *Musubi* further information about the variable values. This depends on the
used test case. They are defined in the [[bc_states_type]]. You can use the
following:

* velocity

* massflowrate

* pressure

* molefrac

* moleflux

* moleDens

* molediff_flux

* pdf

* potential

* surChargeDens


For example, to simulate the channel, you are free to give information about the
variable values with **space time functions**.

### Space Time Functions

> The documentation for these functions can be found in the
> [[tem_spacetime_fun_module]] in subroutine [[tem_load_spacetime_single]].

Variable values like `pressure` and `velocity` are defined as
**space time functions**. Space time functions contain the x-, y- and
z-coordinates and the time as arguments. There are some different ways to get
to a value for i.e. pressure and velocity at different timesteps that are
described shortly in the following.

Before they are used in the `boundary_condition` table they are defined in the
`variable` table.
@note See this [documentation](|temurl|/page/features/variables/index.html)
for more information.

> Inside the `variable = {..}` table space time functions can be defined as a
> constant, a lua function, a predefined fortran function or a combination of
> spatial and transient functions that are themselves defined as a constant, a
> lua function or a predefined fortran function.

Here are a few examples on how to implement space time functions as a boundary
condition.

* **lua function**

```lua
function lua_function(x,y,z,t)
  ...
  return
end

variable = {
  {
    label = 'var_label',
    ncomponents = 1,
    kind = 'st_fun',
    st_fun = lua_function
  },
  ...
}

boundary_condition = {
  {
    label = 'bc_label',
    kind = 'bc_kind'
    bc_variable = 'var_label'
  },
  ...
}
```

* **constant**

```lua
boundary_condition = {
  label = 'bc_label',
  kind = 'bc_kind',
  bc_variable = 1.0
}
```

* **variables with more than one component**

Do not forget to give combination variables to make a vector out of each
velocity, moleflux and molediff_flux component.

Here is an example on how to implement velocity as a boundary variable:

```lua
variable = {
-- example for a predefined space time function
  { name = 'inlet_vel',
    ncomponents = 3,
    vartype = 'st_fun',
    st_fun = {
      predefined = 'combined',
      temporal  = {
        predefined = 'smooth',
        min_factor = 0, max_factor = u_in_L,
        from_time = 0, to_time = tmax/4
      },
      spatial = {
        const = {1.0,0.0,0.0}
      }
    }
  },
-- example for a combination of 3 functions
  { name = 'velX',
    ncomponents = 1,
    vartype = 'st_fun',
    st_fun = lua_function
  },
  { name = 'velY',
    ncomponents = 1,
    vartype = 'st_fun',
    st_fun = 0.0
  },
  { name = 'velZ',
    ncomponents = 1,
    vartype = 'st_fun',
    st_fun = 0.0
  },
  { name = 'inlet_vel2',
    ncomponents = 3,
    vartype = 'operation',
    operation = {
      kind = 'combine',
      input_varname = {'velX','velY','velZ'}
    }
  }
}
```

> You can find more examples in the subroutine [[tem_load_bc_state]] and on
> [variable system page](|temurl|/page/features/variables/index.html).

Next chapter: [Restart &rarr;](tut_05_restart.html)
