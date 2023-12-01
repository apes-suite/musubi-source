-- In the following boundary conditions (north, east, south and west), the
-- z-direction is treated with length. By this, we make sure that the depth is
-- large enough, to enable all cases (2D/3D w/o periodicity).

-- Use this file as template. Do not modify this file for running some testcases

recheck = true
if recheck then
  stlFolder = '$!stl_path!$'
else
  stlFolder = './'
end

--! [Global variables]
-- Height of the channel [m]
height = 1.0
-- Length to height ratio
l_h = 8.0
-- Length of the channel [m]
length = l_h*height
-- Depth of the channel [m]
depth = length
-- General level of mesh refinement
level =  8
-- Special levels for local refinement. Set to 0 to make refinement optional.
-- If it should be used, activate useRefine in the 'Geometry definition'. The
-- refinmentLevels will be increased then.
refinementLevel  = 0
refinementLevel2 = 0
--outputpreview = true
--! [Global variables]


--! [Geometry definition]
case2d = true           -- Creates a 2D channel if set to true
usePeriodic = true       -- Activates periodic boundaries
qValues = true           -- Activates qvalues
useObstacle = false      -- Activates obstacle in the channel
fixHeight = true         -- Channel's height is fixed and is used to compute dx
useRefine = false        -- Activate refinement, increase levels by 1 and 2.
smoothlevels = false     -- False avoids refinement of single elements.
smoothbounds = false     -- False avoids refinement of single elements.
--! [Geometry definition]


--! [dx computed by height]
if fixHeight then
  -- Number of elements in the height
  nHeight = 2^(level-4)
  -- Element size [m]
  dx = height/nHeight
  -- Number of elements in bounding cube
  -- = number of elements in channel + inlet + outlet
  nLength = nHeight*l_h + 2
  -- Level required to reach computed dx
  level = math.ceil( math.log(nLength)/math.log(2))
  -- Length of the bounding cube
  length_bnd = (2^level)*dx
--! [dx computed by height]


--! [dx computed by length]
else
  -- Length of the bounding cube
  length_bnd = length/(1.0-2.0/2^level)
  -- Element size [m]
  dx  = length_bnd/(2.0^level)
  -- Number of elements in height
  nHeight = 2.0*math.floor(length/(dx*2.0*l_h))
  -- Height of the channel
  height = nHeight*dx
end
--! [dx computed by length]


--! [Reference values]
-- A minimum level, by which all parts in the computational domain should at
-- least be resolved with. Default is 0. Here, set it to value calcualted above.
minlevel=level
-- A maximum level, which is reached on finest level of refinement.
maxLevel = level+math.max(refinementLevel, refinementLevel2)
-- Smallest element size (on finest level) [m]
dxMin  = length_bnd/(2^maxLevel)
-- Half of element size on finest level [m]
dx_half = 0.5*dxMin
-- Smallest possible element size
dx_eps = length_bnd/2^20

-- Increase refinement levels, if useRefine is active.
if useRefine then
  refinementLevel  = refinementLevel + 1
  -- Refinement level 2 can only be one larger than 1, otherwise level jump is
  -- too big.
  refinementLevel2 = refinementLevel + 1
end

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
--! [Reference values]

---------------------------Seeder configurations -------------------------------
-- Directory to write mesh
folder = 'mesh/'

-- File to write measurement time of mesh generation each step
timing_file = 'sdr_timing.res'

-- How detailed should the output be?
-- The higher the level, the more output you'll get from seeder.
logging = { level = 3 }

-- Debug outputs
NOdebug = {
  debugMode = true,
  debugFiles = true,
  debugMesh = 'mesh/' 
}

--! [Bounding cube]
-- Bounding_cube: two entries (origin and length in this order, if no keys used)
bounding_cube = {
  origin = { -0.5*length_bnd, -0.5*length_bnd, -0.5*length_bnd+0.5*dx},
  length = length_bnd
}
--! [Bounding cube]


--! [Spatial objects]
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
--! [Inlet and outlet boundary conditions]
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
--! [Inlet and outlet boundary conditions]

--! [Walls]
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
--! [Walls]
}
--! [Spatial objects]


--! [Use periodic boundaries]
-- If usePeriodic is active, periodic boundary conditions are used in depth,
-- which means in z-direction in this case. This is one way of creating a quasi-
-- 2D domain.
if usePeriodic == true then
  depth = dx
  table.insert(spatial_object,
    {
      attribute = {
        kind = 'periodic',
        label = 'periodic'
      },
      geometry = {
        kind = 'periodic',
        object = {
          plane1 = { -- top plane
            origin = {
              -0.5*length - dx_eps,
              -0.5*height - dx_eps,
              0.5*depth + dx_half/4.0
            },
            vec = {
              { length + 2*dx_eps, 0.0, 0.0 },
              { 0.0, height + 2*dx_eps, 0.0}
            }
          },
          plane2 = { -- bottom plane
            origin = {
              -0.5*length + dx_eps,
              -0.5*height + dx_eps,
              -0.5*depth - dx_half/4.0
            },
            vec = {
              { 0.0, height + 2*dx_eps, 0.0 },
              { length + 2*dx_eps, 0.0, 0.0 }
            }
          },
        }
      }
    }
  )
--! [Use periodic boundaries]

--! [Non-periodic]
-- If usePeriodic is deactived, but case2d is selected, we use another wayt to
-- create a quasi-2D channel: a channel with 1 elements in depth. If not, we
-- create a 3D channel.
else
  if case2d then
    depth = dx
  else
    depth = height
  end
  table.insert(spatial_object,
    {
      attribute = {
        kind = 'boundary',
        label = 'top'
      },
      geometry = {
        kind = 'canoND',
        object = {
          origin = {
            -0.5*length - dx_eps,
            -0.5*height - dx_eps,
             0.5*depth + dx_eps
          },
          vec = {
            { length + 2*dx_eps, 0.0, 0.0 },
            { 0.0, height + 2*dx_eps, 0.0}
          }
        } -- object
      } -- geometry
    } -- spatial object
  )
  table.insert(spatial_object,
    {
      attribute = {
        kind = 'boundary',
        label = 'bottom'
      },
      geometry = {
        kind = 'canoND',
        object = {
          origin = {
            -0.5*length - dx_eps,
            -0.5*height - dx_eps,
            -0.5*depth - dx_eps
          },
          vec = {
            { length + 2*dx_eps, 0.0, 0.0 },
            { 0.0, height + 2*dx_eps, 0.0}
          }
        },-- object
      } -- geometry
    } -- spatial object
  )
end --periodic
--! [Non-periodic]


--! [Use of obstacles]
-- If obstacle is true, we distinguish between 2D (cylinder) and 3D (sphere)
-- case.
if useObstacle == true then
  if case2d == true then
    stlfile  = stlFolder .. 'stl/cylinder.stl'
    stlLabel = 'cylinder'
  else
    if usePeriodic == true then
      stlfile  = stlFolder .. 'stl/cylinder.stl'
      stlLabel = 'cylinder'
    else
      stlfile  = stlFolder .. 'stl/sphere.stl'
      stlLabel = 'sphere'
    end
  end
  -- If we use an obstacle, we can make use of qValues. Here, we also
  -- distinguish between 2D (cylinder) and 3D (sphere). If we don't want to use
  -- them.
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
end -- useObstacle
--! [Use of obstacles]
------------------- End of seeder configurations -------------------------------
