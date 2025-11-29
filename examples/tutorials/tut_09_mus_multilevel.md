title: Multilevel Simulations
@warning WORK IN PROGRESS @endwarning

Navigate: [&larr; Source Terms](tut_08_source.html)
| [Overview](index.html)

# Multi-level Simulations

In some simulations there are regions of special interest. In these regions
strong gradients might occur, so that a higher resolution is required.
You then have basically two possibilities to resolve this problem:

1. Increase the global resolution of your simulation. This is generally
   not a good idea because you very quickly might run into memory problems.
   A high number of fluid cells also means long computation times.
2. Locally refine the mesh in regions of high interest.
   You can refine the mesh in regions where it is desired.
   This focuses the computation cost and memory expenses on the regions
   where it actually is required. Why waste computation time on flow regions
   where very little happens? A general problem are the interfaces between
   grid levels. You can easily get reflections or non-discontinuities there
   which is why you usually place interfaces in regions where the flow does
   not fluctuate considerably.

> In [[mus_interpolate_module]] you can get some background information on the
> implemented interpolation methods and the general workflow.


## Test Case Description ##

In this tutorial we cover a channel test case with a cylinder inside.
The area around and behind the cylinder is refined by a refinement
patch with a higher resolution.

@note To run this case use the generic channel test case and set in
`seeder.lua`:
```lua
case2d = false
usePeriodic = false
qValues = true
useObstacle = true
fixHeight = true
useRefine = true
```
@endnote


## Mesh Generation ##

Before starting your simulation you have to set up a mesh again.
Again we involve *Seeder* to generate the required data structures for the
solver.
It is very advisable to define a variable for the most important variables.
Later on you will be able to easily change the resolution of your simulation
or some aspect ratio.
One very important parameter for setting the resolution of the simulation
is a reference tree level that is done in [[tem_topology_module:tem_levelof]].

In our case this will be the tree level of the channel region. Let's name it
`level` and set it based on the computed `dx`. In our default setting this
results in a level of 8. In addition to that, we define the folder (create it
via `mkdir mesh`), the refinement level and the resulting level for the
refinement box.

The region around and behind the cylinder is defined by a refinement box.
Let's define first, by how many levels we want to refine the elements
inside this refined area. For simplicity, let's say that all elements should
be on one level higher than the rest of the channel (`refinementLevel`).

It is also good to have the information about the maximum level in your
simulation available as a parameter (`maxlevel`). You will later on see why.

```lua
if useRefine then
  refinementLevel  = refinementLevel + 1
  -- Refinement level 2 can only be one larger than 1, otherwise level jump is
  -- too big.
  refinementLevel2 = refinementLevel + 1
end
[...]
maxLevel = level+math.max(refinementLevel, refinementLevel2)
[...]
folder = 'mesh/'
```

> Note: For a description of levels and the layout of the tree have
> a look at the [octree page](|temurl|/page/octree.html).


```lua
height = 1.0         -- Height of the channel [m]
l_h = 8.0            -- Length to height ratio
length = l_h*height  -- Length of the channel [m]
depth = length       -- Depth of the channel [m]
level =  8           -- General level of mesh refinement
```

The bounding box, which is the universe in which all elements of the tree will live in,
is defined as

```lua
bounding_cube = {
  origin = { -0.5*length_bnd, -0.5*length_bnd, -0.5*length_bnd+0.5*dx},
  length = length_bnd
}
```

The region around the channel is set to the desired tree level of 9 and 10.

```lua
spatial_object = {
  [...]
  {
    attribute = {
      kind = 'refinement',
      level = refinementLevel+level,
      label ='box1'
    },
    geometry = {
      kind = 'canoND',
      object = {
        origin = { start_x, start_y, start_z },
        vec = {
          { size_x, 0.0, 0.0 },
          { 0.0, size_y, 0.0 },
          { 0.0, 0.0, size_z }
        }
      }
    }
  },
  {
    attribute = {
      kind = 'refinement',
      level = refinementLevel2+level,
      label ='box2'
    },
    geometry = {
      kind = 'canoND',
      object = {
        origin = { start2_x, start2_y, start2_z },
        vec = {
          { size2_x, 0.0, 0.0 },
          { 0.0,size2_y, 0.0 },
          { 0.0, 0.0, size2_z }
        }
      }
    }
  },
  [....]
}
--! [Refinement box 1 & 2]
```

The rest of the bounding cubic domain is not of interest and hence the
discretization level is not of interest.

### Defining Geometry ###

We need to specify the walls, the inlet, outlet and the cylinder.
Let's start with the walls at the north, south position.
The walls at east and west direction will later on become the in- and outlet.

```lua
[...]
    {
      attribute = {
        kind = 'boundary',
        label ='north'
      },
      geometry = {
        kind = 'canoND',
        object = {
          origin = { -0.5*length-dx_eps, 0.5*height+dx_eps, -0.5*depth-dx_eps },
          vec = {
            { length + 2*dx_eps, 0.0, 0.0 },
            { 0.0, 0.0, depth + 2*dx_eps }
          }
        }
      }
    },
    {
      attribute = {
        kind = 'boundary',
        label ='south'
      },
      geometry = {
        kind = 'canoND',
        object = {
          origin = { -0.5*length-dx_eps, -0.5*height-dx_eps, -0.5*depth-dx_eps },
          vec = {
            { length + 2*dx_eps, 0.0, 0.0 },
            { 0.0, 0.0, depth + 2*dx_eps }
          }
        }
      }
    },
[...]
  {
    attribute = {
      kind = 'boundary',
      label ='east' -- outlet
    },
    geometry = {
      kind = 'canoND',
      object = {
        origin = { 0.5*length+dx_eps, -0.5*height-dx_eps, -0.5*depth-dx_eps },
        vec = {
          { 0.0, height + 2*dx_eps, 0.0 },
          { 0.0, 0.0, depth + 2*dx_eps }
        }
      }
    }
  },
  {
    attribute = {
      kind = 'boundary',
      label ='west' -- inlet
    },
    geometry = {
      kind = 'canoND',
      object = {
        origin = { -0.5*length-dx_eps, -0.5*height-dx_eps, -0.5*depth-dx_eps },
        vec = {
          { 0.0, height + 2*dx_eps, 0.0 },
          { 0.0, 0.0, depth + 2*dx_eps }
        }
      }
    }
  },
```

For defining the cylinder or sphere, we use an STL file. This file was before
generated by Blender, but you can basically use any software you like for
generating such geometry data. The files are located in the test case directory:
`stl/*.stl`.
Let's include the sphere STL file with:

```lua
   stlfile  = 'stl/sphere.stl'
   stlLabel = 'sphere'
```
then add is to `spatiatl_object` tabel.

```lua
  table.insert(spatial_object,
    {
      attribute = {
        kind = 'boundary',
        level = maxLevel,
        label = stlLabel,
        calc_dist = qValues,
      },
      geometry = {
        kind = 'stl',
        object = {
          filename = stlfile,
        }
      },
      transformation = {
        -- Deformation factor to scale the obstacle. Here: 1 --> remains the same.
        deformation =  1,
        -- In this case, we do not have to move the geometry.
        translation =  { 0.0, 0.0, 0.0 }
      }
    }
  )
```

One very important action is the placement of the seed.
The seed determines the contiguous flow domain. It basically
defines for the shapes, what is inside and what is outside.
We palce it right in the middle of the channel, which is simply

```lua
spatial_object = {
  {
    attribute = {
      kind = 'seed'
    },
    geometry = {
      kind = 'canoND',
      object = {
        origin = { 0.0, 0.0, dx_half*0.5 },
      }
    }
  },
[...]
```

Ok. Now we have defined all geometric constraints.
Let's continue with refining some parts of the simulation domain.

### Defining Refined Regions ###

You first need to specify the origin of this box and its extents.
Again, we use some variables depending on the total bounding cube.

```lua
--! [Refinement box 1]
-- Size of refinement box 1 in x-direction
size_x = 0.75*length
-- Start of refinement box 1 (x)  based on length
start_x = -0.425*length
-- Size of refinement box 1 in y-direction
size_y = 0.5*height
-- Start of refinement box 1 (y)  based on height
start_y = -0.25*height
-- Size of refinement box 1 in y-direction
size_z = 0.5*height
-- Start of refinement box 1 (z)  based on height
start_z = -0.25*height
--! [Refinement box 1]


--! [Refinement box 2]
-- Size of refinement box 2 in x-direction and start of it
start2_x = -0.375*length
size2_x  = 0.5*size_x
-- Size of refinement box 2 in y-direction and start of it
size2_y  = 0.5 * size_y
start2_y = -0.5* size2_y
-- Size of refinement box 2 in z-direction and start of it
size2_z  = size2_y - dx_half
start2_z = start2_y + dx_half/2.0
--! [Refinement box 2]
```

The patch will then lie around and behind the cylinder.

```lua
spatial_object = {
[...]
--! [Refinement box 1 & 2]
  {
    attribute = {
      kind = 'refinement',
      level = refinementLevel+level,
      label ='box1'
    },
    geometry = {
      kind = 'canoND',
      object = {
        origin = { start_x, start_y, start_z },
        vec = {
          { size_x, 0.0, 0.0 },
          { 0.0, size_y, 0.0 },
          { 0.0, 0.0, size_z }
        }
      }
    }
  },
  {
    attribute = {
      kind = 'refinement',
      level = refinementLevel2+level,
      label ='box2'
    },
    geometry = {
      kind = 'canoND',
      object = {
        origin = { start2_x, start2_y, start2_z },
        vec = {
          { size2_x, 0.0, 0.0 },
          { 0.0,size2_y, 0.0 },
          { 0.0, 0.0, size2_z }
        }
      }
    }
  },
--! [Refinement box 1 & 2]
[...]
}
```

Once you followed through the above explanations, you can visualize the
generated mesh. Therefore use *Seeder Harvesting*.

First, we prepare the config file, defining the folder of the mesh files,
the output folder for the vtk file (create the folder first) and the output
itself. We name the file `sdr_harvester.lua`.

```lua
mesh          = 'mesh/'
output_folder = 'mesh/'
output = {
   format   = 'vtk',
   dataform = 'binary',
   write_pvd = false
}
```
Then we can run *Seeder Harvesting* with:
`~/apes/seeder/build/sdr_harvesting harvester.lua`

## Setting up the Musubi Configuration ##

After generating the mesh above, we need to tell *Musubi* that it
should use the above generated mesh. The mesh was stored in the
folder `mesh`. Let's define that along with a name for the simulation

```lua
mesh = './mesh/' -- Mesh information
simulation_name = 'channelRefine'
```

### Initial and Boundary Conditions ###

Initial conditions are set to a medium with constant pressure and velocity
`ic_velX`.

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

Next, we define the variable table for the boundary conditions:

```lua
--! [User defined variables]
-- Mainly used for tracking.
-- This variable can be refered to as variable in boundary condition and source
variable = {
  -- Reference pressure dependent on physicsModel
  {
    name='pressureRef',
    ncomponents=1,
    vartype = 'st_fun',
    st_fun = pressureRef
  },
  -- Pressure at inlet
  {
    name='pressureIn',
    ncomponents=1,
    vartype = 'st_fun',
    st_fun = pressureIn
  },
  -- Pressure at outlet
  {
    name='pressureOut',
    ncomponents=1,
    vartype = 'st_fun',
    st_fun = pressureOut
  },
  -- Reference shear stress
  {
    name='stressRef',
    ncomponents=1,
    vartype = 'st_fun',
    st_fun = stressRef },
  -- Difference between numerical pressure and reference pressure
  {
    name='press_diff',
    ncomponents=1,
    vartype = 'operation',
    operation = {
      kind = 'difference',
      input_varname = { 'pressure_phy', 'pressureRef',},
    },
  },
}
--! [User defined variables]
```


The boundary conditions are a little bit more complex, as we have
solid walls, the cylinder, and also an inlet and an outlet.
Let's start with the wall boundaries. These include among the walls, which we named
according to the directions and the cylinder, which we simply called `obs`.
They all get the property of a simple wall.

```lua
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
  [...]
```

For the outlet, we would like to have a
[[mus_bc_fluid_module:pressure_expol]] "simple pressure" boundary
condition, with pressure set to `pressureOut` a variable we defined in
the variable table above.

```lua
[...]
  {  label = 'east',
     kind  = 'pressure_expol',
     pressure = 'pressureOut'
  }
[...]
```

For the inlet we also use `pressure_expol` boundary condition.

```lua
[...]
  {  label = 'west',
     kind  = 'pressure_expol',
     pressure = 'pressureIn'
  },
[...]
```

If we set `usePeriodic=false` we need to define `top` and `bottom` as wall.
If we use periodic *Seeder* will handle the periodic relationship for these
planes.

```lua
-- If we deactivated usePeriodic in seeder.lua, we have to add boundaries for
-- top and bottom.
if usePeriodic ==false then
  table.insert( boundary_condition,
  {  label = 'top',
     kind = 'wall'
  }
  )
  table.insert( boundary_condition,
  {
    label = 'bottom',
    kind = 'wall'
  }
end
```

### Other Simulation Parameters ###

```lua
tmax_phy = 10.0
tmax_iter =  math.ceil(tmax_phy/dt)
interval_phy = tmax_phy/10.0
trac_start = 0.0
rest_start = tmax_phy/4.0
sim_control = {
  time_control = {
    max = { sim = tmax_phy },
    interval = { sim = interval_phy },
    clock = 3600 --s
   },
  [...]
}
```

### Tracking ###

We have several different trackers defined. These can be found inside the
`tracking` table. Tracking has been explained in [[chapter 03]](tut_03_tracking.html)

```lua
tracking = {
  -- Track pressure at the center of the channel over time.
  {
    label = 'probeAtCenter',
    --label = 'probePressure',
    folder = 'tracking/',
    variable = {'pressure_phy'},
    shape = {
      kind = 'canoND',
      object = {
        origin = { 0.0, 0.0, 0.0 }
      }
    },
    time_control = { min = {iter=1}, interval = {iter=1} },
    output = { format = 'ascii' },
  },
  [...]
```

### Visualize the Results ###

For visualisation of ascii or ascii-spatial format files you can use for example
Gnuplot or matplotlib module of python. An example script for postprocessing
(`plot_track.py`) can be found inside the example directory.

Navigate: [Overview](index.html)
