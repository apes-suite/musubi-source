Definition of Boundary Conditions {#tut_4boundaries}
========

Navigate: \ref tut_3tracking "&larr; Tracking quantities during the simulation run"
| \ref mus_tutorials_overview "Overview"
| \ref tut_5restart "Restarting Simulations &rarr;"

In this tutorial we explain how to define boundary conditions that are needed for the simulation run.

## Where the information for the boundaries are taken from

In Seeder, boundaries are defined by `spatial_object` attribute kind "boundary". The output of 
Seeder is normally placed in the `./mesh` folder. There you can find files called `bnd.lua` and 'bnd.lsb'.
**bnd.lua** file contains basic information about boundaries like number of boundaries in 
a mesh and list of boundary labels in ascii format. **bnd.lsb** file contains boundary IDs for 26 directions
for each fluid element which has boundary neighbor in binary format.
These two files are important for Musubi in case of defining the boundary conditions. 

If you open bnd.lua file, you will probably find something like this:

\code
 nSides = 26
 nBCtypes = 2
 bclabel = {
    'south',
    'north',
    'west',
    'east'
}
\endcode

In the example, there are four different boundaries we have set up with Seeder. In this case they have
the labels `north`, `south`, `east` and `west`. This gives you a feeling of how to view this channel.

These boundaries are defined in the mesh, but it is not clear which function each boundary has 
at the moment. So this is the part you have to do in Musubi. As an example, we want to simulate 
a channel. A channel in 2D needs two walls, one inlet and one outlet. 

> If you would like to know in detail how the boundaries are defined you might have a look at
> this [page](@ref tem_bc_module).

## How to define boundary conditions in Musubi?

Now, we will have a look on the `boundary_conditions` table in Musubi for the channel testcase.

\snippet testcases/channel/musubi.lua boundary conditions

In this basic example you can see the function of each boundary as outlined above.

### boundary_condition table

In `boundary_condition` table, for each boundary label from Seeder, a kind must be defined in Musubi 
which defines what to do with that boundary. The order of boundary definition in Musubi does not dependent
on the order in seeder (bnd.lua). It just has to be existent in both files.
\code
boundary_condition = {
}
\endcode

### label

For every boundary which you have created in Seeder, you have to set its status.
Therefore, you call the boundary with the exact name (`label = ...`) that you can see in the bnd.lua
file and give it a certain **kind** that is explained in the next section.
\code
boundary_condition = {
  label = 'south',
}
\endcode

### kind

You can choose between some basic boundary kinds. They define the use of the boundary in the simulation run.
Some of these kinds are described below.

* **wall**

  A wall means that the fluid is not able to penetrate through this boundary. It has to regard the wall as 
  an obstacle. Moreover, the wall is seen as a **no-slip** boundary. If you would like to observe slip as
  well you have to use the separate kind `slip_wall`. 
  \code
  boundary_condition = {
    label = 'south',
    kind = 'wall'
  }
  \endcode

* **slip_wall**

  If slip shall be defined as well, you will have to set `kind` to `kind = slip_wall`. Slip means that
  the normal and the tangential velocity in normal direction equal zero. The pressure gradient along the
  normal direction is equal to zero as well. The degree of slip can be defined by the multiplication
  of a slip-factor called `fac` and the velocity. Special cases are on the one hand `fac = 1` which means 
  that there is free-slip or full slip and on the other hand, there is `fac = 0` which is used for no-slip.
  This case is the default case for the `kind = wall` which is mentioned above.
  \code
  boundary_condition = {
    label = 'north',
    kind = 'slip_wall',
    fac = 0.4 
  }
  \endcode
  > More information can be found in the [mus_bc_fluid_wall_module](@ref mus_bc_fluid_wall_module).

* **wall_linearInterpolation**

  There is a possibility to make your simulation more efficient: Instead of making use of a higher 
  refinement level, you could use linear interpolation. The obstacle in the fluid will have a higher
  resolution if the distance between the barycentric centre position of the fluid element to the obstacle
  is calculated which is "q". The total distance between the center positions of each, the fluid element's
  and the obstacle's is q=1. There is a case differentiation between q<0.5 and q>=0.5.
  Actually, it is not possible to display a formula in Latex here. Maybe this can be fixed.
  For q<0.5 this formula is used:  
  \f${ 
  {\displaystyle f_{i^{\prime}}(\mathbf{r}_{l},t+1)} 
  {\displaystyle =2qf_{i}^{c}(\mathbf{r}_{l},t)+(1-2q)f_{i}^{c}(\mathbf{r}_{l}-\mathbf{c}_i{i},t)}
  \f$}
  
  For q>=0.5 this formula is used:

  Therefore, Seeder has to be configured. For each spatial_object which has boundary kind
  `calc_dist = true` has to be added to the attributes right behind kind and label. Here is the seeder.lua
  block before the changes:  
  \snippet testcases/channel/seeder.lua walls
  >  In the `spatial_object` tables for the boundaries, you have to add to the attributes `calc_dist = true`. 
  Then there should be written as an example for one boundary:
  \snippet testcases/channel/seeder.lua use of obstacles
  After that, you are able to set the kind of the boundary in musubi.lua to `wall_linearInterpolation` 
  instead of `wall`. You have to use the same syntax as shown above, otherwise it will not work.  
  \code
  boundary_condition = {
    label = 'south',
    kind = 'wall_linearInterpolation',
  }
  \endcode

  > The information about linear interpolation are taken from M. Bouzidi, M. Firdaouss, and P. Lallemand, 
  > "Momentum transfer of a Boltzmann-lattice fluid with boundaries," Physics of Fluids, vol. 13, no. 11, 
  > pp. 3452–3459, Nov. 2001 Equations 5 (linear interpolation). 
  > More information about no-slip wall linear interpolation can be found [here](@ref mus_bc_fluid_wall_module).

### velocity boundary conditions

Velocity boundary conditions need a given velocity as input. This could be a constant or a space
time function for example. But as shown in the following examples, the velocity is called as a
variable name. This means that you have to define the velocity in the variable table first.
For each component we need to define a variable. After that we combine them to 'inlet_vel'.

\code
variable = {
  { name = 'vel_x',
    ncomponents = 1,
    vartype = 'st_fun',
    st_fun = {
      predefined = 'combined',
      transient = {
        predefined='smooth', min_factor = 0.25, 
        max_factor=1.0, from_time=0, to_time=1000*dt
      },
      patial = { -- example for predefined function
        predefined='parabol', 
        shape = { 
          kind = 'canoND',
          object = {
            origin  = {0.0,0.0,zpos},
            vec = {0.0,height,0.0}
          },
        }, -- shape
        amplitude = 1.0
      },
      spatial = u_inflow -- u_inflow is an example for a lua space time function
    }, -- vel_x
  },
  { name = 'vel_y',
    ncomponents = 1,
    vartype = 'st_fun',
    st_fun = 0.0
  },-- vel_y
  { name = 'vel_z',
    ncomponents = 1,
    vartype = 'st_fun',
    st_fun = 0.0
  },--vel_z
  { name = 'inlet_vel',
    ncomponents = 3,
    vartype = 'operation',
    operation = {
      kind = 'combine',
      input_varname = {'vel_x','vel_y','vel_z'}
    }  
  },--inlet_vel
}
\endcode

* **inlet_eq**

  `kind = inlet_eq` makes use of the equilibrium function concerning the Lattice Boltzmann method (LBM).
  This function gets the density (rho) from the fluid. The velocity (u) has to be defined for 
  each coordinate x, y and z in the variable table. 
  Example:
  \code 
  boundary_condition = {
  { label = 'inlet',
    kind = 'inlet_eq',
    velocity = 'inlet_vel' -- 'inlet_vel' must be set in variable table
  }
  \endcode 
  > You can get more information about the equilibrium function currently in the documentation for the 
  > subroutine [inlet_eq](@ref mus_bc_fluid_module::inlet_eq).

* **inlet_ubb**

  In addition to that you can set `kind = inlet_ubb`. Ubb in this case stands for "velocity bounce back".
  You can imagine this like this: A fluid particle gets near the boundary. If it reaches the boundary, the 
  particle will bounce back in the same angle as it gets there. For this kind, you have to set the velocity
  values, too. 
  \code
  boundary_condition = {
    { label = 'west',
      kind = 'inlet_ubb',
      velocity = 'someVelocity'
  \endcode
  > More details can be found in the documenation for subroutine [inlet_ubb](@ref mus_bc_fluid_module::inlet_ubb_qval).

* **inlet_mfr**

  Like for `inlet_ubb` there is an "Inlet Velocity Bounce Back" boundary condition. But in this case, mass flow
  rate is used as an input argumet as well. The velocity is given by `me%ubb%velocity`. It is used like that.
  \code
  boundary_condition = {
    { label = 'inlet',
      kind = 'inlet_mfr',
      mass_flowrate = 'mfr' } -- 'mfr' must be set in variable table
  \endcode
  > For more information visit the [Documentation](@ref mus_bc_fluid_module::inlet_mfr).

* **inlet_mfr_eq**

  In this case, the mass flow rate is used and the velocity is taken from the configuration file.
  \code
  boundary_condition = {
    { label = 'inlet',
      kind = 'inlet_eq',
      mass_flowrate = 'mfr' } -- 'mfr' must be set in variable table
  \endcode
  > The corresponding Documentation can be found [here](@ref mus_bc_fluid_module::inlet_mfr_eq).

### pressure boundary conditions

Like for velocity boundary conditions we need to define the pressure as variable first.
Here is an example:
\code
variable = {
  { name = 'p0',
    ncomponents = 1,
    vartype = 'st_fun',
    st_fun = 1.0, 
  },-- p0
}
\endcode

* **outlet_expol**

  The variable values are extrapolated during the simulation.
  \code
  boundary_condition = {
    { label = 'east',
      kind = 'outlet_expol',
      pressure = 'p0' } -- must be set in variable table
  \endcode
  > Detail information can be found in the [Documentation](@ref mus_bc_fluid_module::outlet_expol).

* **outlet_pab**
  
  This is the outlet pressure anti-bounce back boundary condition kind. The velocity is extrapolated by
  two of its neighbors. The pressure has to be given as well.
  \code
  boundary_condition = {
    { label = 'east',
      kind = 'outlet_pab',
      pressure = 'p0' } -- must be set in variable table
  \endcode
  > More information [here](@ref mus_bc_fluid_module::outlet_pab).

* **outlet_eq**

  The incoming densities are set to the equilibrium distribution with macroscopic velocity and pressure.
  \code
  boundary_condition = {
    { label = 'east',
      kind = 'outlet_eq',
      pressure = 'p0' } -- must be set in variable table
  \endcode
  > Detail information in the [Documentation](@ref mus_bc_fluid_module::outlet_eq).

### Physical conditions

If you define the inlet and the outlet boundary for the shape, you will have to give Musubi further 
information about the variable values. This depends on the used testcase. They are defined in the 
[bc_states_type](@ref mus_bc_header_module::bc_states_type). The variables that you need to define in
the variable table are named with number of its components in subroutine [mus_store_bcVarPos](@ref mus_variable_module::mus_store_bcVarPos). 
You can use the following:

* velocity (3 components)
  
* mass_flowrate (1 component)

* pressure (1 component)

* mole_fraction (1 component)

* mole_flux (3 components)

* mole_diff_flux (3 components)

For example, to simulate the channel, you are free to give information about the variable values with
**space time functions**.

### Space time functions

> The documentation for these functions can be found in the [tem_spacetime_fun_module](@ref tem_spacetime_fun_module).

Variable values like `pressure` and `velocity` are defined as **space time functions**. 
Space time functions contain the x-, y- and z-coordinates and the time as arguments. 
There are some different ways to get to a value for i.e. pressure and velocity at different 
time steps that are described shortly in the following.

One example is to have a constant value over the time. 
Therefore you are free to give a constant value like:
\code
boundary_condition = {
  label = 'west',
  kind = 'outlet_expol',
  pressure = 'p0'
}

variable = {
  { name = 'p0',
    ncomponents = 1,
    vartype = 'st_fun' -- space time function
    st_fun = 2.0 -- constant
  }
}
\endcode

Or like in this example, the pressure is defined by a spatial function. This means that the pressure is a 
constant over time. In this case, the transient function would be "1".

\code
function pressureRef(x,y,z) -- lua spatial function
  dp = u0Phys*8.*viscPhys*rho0Phys/heightPhys^2*length
  return p0 + dp*0.5 - dp/length*x
  end

boundary_condition = {
  label = 'west',
  kind = 'outlet_expol',
  pressure = 'pressureRef'
}

variable = { -- variable table 
  { name='pressureRef', 
    ncomponents=1, 
    vartype = 'st_fun', 
    st_fun = pressureRef(1.0, 0.0, 0.0) -- uses lua space time function 
  }
}
\endcode

Another one is a combination of a [transient function](@ref tem_transient_module), 
that replaces the time as an argument, and a [spatial function](@ref tem_spatial_module), 
that has the x-, y, and z-coordinate as an argument. 
Each, the transient and the spatial function can be a constant, 
an anywhere predefined function or a function that is given directly in this block.

Here is an example for `velocityX` which is defined by a combination of a 
spatial and a transient function.

\code
variable = {
  { name = 'vel_x',
    ncomponents = 1,
    vartype = 'st_fun',
    st_fun = { 
      style = 'dirchlet',
      transient = {
        predefined ="linear", 
        min_factor = 0.0,
        max_factor = 1.0, 
        from_time = 0, 
        to_time = 1000
      },
-- example for line as shape:

--      spatial = {                  
--        predefined='parabol', 
--        shape = {
--          object = 'line',
--          line = { 
--            origin={-2.0,0.0,0.0},
--            direction={0.0,1.0,0.0}
--          }
--        }
--      },
-- example for plane as shape:
      spatial = {
        predefined='parabol', 
        shape = {
          object = 'plane',
          plane = { 
            origin={-2.0,0.0,0.0},
            vecA={0.0,1.0,0.0},
            vecB={0.0,0.0,1.0} 
          } 
        } 
      }
    }
  }
}
\endcode

\note Do not forget to give functions for each coordinate for variable values like velocity, mole_flux and
mole_diff_flux that are vectors. Please make use of the same syntax like in the examples for velocityX, 
velocityY and velocityZ.

> The examples are taken from the Documentation for the subroutine [tem_load_bc_state](@ref tem_bc_module::tem_load_bc_state).

Next chapter: \ref tut_5restart "Restarting Simulations &rarr;"
