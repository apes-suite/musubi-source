-- Use this file as template. Do not modify this file for running some testcases

--! [common header]

level =  8                -- general level of mesh refinement
refinementLevel  = 0      -- special levels for local refinement
refinementLevel2 = 0
label   = 'channel'       
idLabel = 'channel'
simName = 'channel'
outputname = 'channel_cylinder'
outputpreview = true
folder = 'mesh/' 

--! [common header]

--! [geometry definition]

qValues = true            -- activates qvalues
usePeriodic = true        -- activates periodic boundaries
periodicX = false         -- activates periodic boundaries for x-axis
useObstacle = false       -- activates obstacle in the channel
case2d = false            -- creates a 2D channel
walls = true              -- walls are used
testIntersection = false  -- ?
verbose = true            -- ?
fixHeight = true          -- the channel's height is fixed and is used to compute dx

--! [geometry definition]

--! [geometry relations]

height = 1.               -- will be adapted to be exactly resolved by the discretization (physical)
ratio  = 8.               -- will be adapted to be exactly resolved by the discretization 
length = ratio*height     -- absolute length of the channel (physical)
radius = 0.1*height       -- radius of the channel if a tube is used (physical)

--! [geometry relations]


-- here should be the logic how to fix the height of the channel.
-- probably define a number of elements in the height?
-- such as 2^(level-4) 
-- -> nH(l=7) = 8
-- -> nH(l=8) = 16 ...


--! [dx computed by height]

if fixHeight then

  nElemsHeight = 2^(level-4)                              -- number of elements in the height (lattice)   
  dx = height/nElemsHeight                                
  nElemsLength = nElemsHeight*ratio + 2                   -- need inlet and outlet element
  level = math.ceil( math.log(nElemsLength)/math.log(2))  -- compute the required level
  lengthSeeder = 2^level*dx                               -- length of a bounding cube element


--! [dx computed by height]


-- computation of the discrete sizes.
-- fixing the length and adapting the height to such a ratio
-- that the half-height is always integer multiples of dx


--! [dx computed by length]

else
  lengthSeeder = length/(1.-2./2.^level)                 -- length of a bounding cube element 
  dx  = lengthSeeder/(2^level)                           
  nElemsHeight = 2.*math.floor(length/(dx*2.*ratio))     -- number of elements in the height (lattice)
  height = nElemsHeight*dx                               -- absolute height of the channel (physical)
end 

--! [dx computed by length]


--! [reference values]

minlevel=level
maxLevel = level+math.max(refinementLevel, refinementLevel2) 
nElemsMax = 2^maxLevel
dxMin  = lengthSeeder/(2^maxLevel)
dxDash = 0.5*dxMin
lPhys = length  -- m
dxPhys = dx
heightLB = nElemsHeight
heightPhys = heightLB*dx

ebug = {debugMode = true, debugMesh = './prototree/'}
size_x = 0.1*length
start_x = -0.5*size_x
if useObstacle then
  size_x = 0.68*length
  start_x = -0.45*length
  size_y = 2.*height/5.-5.*dxDash
  if testIntersection then
    start_x = -0.40*length
  end
else
  size_y = height/5-5.*dxDash
end
start_y = -size_y/2
size_z = size_y
start_z = start_y

start2_x = -0.32*length
size2_x  = 0.37*size_x
size2_y = size_y*0.48
start2_y = -size2_y/2
size2_z = size2_y
start2_z = start2_y

--! [reference values]


-- boundingbox: two entries: origin and length in this
-- order, if no keys are used


--! [bounding cube]

bounding_cube = {origin = {-lengthSeeder*0.5, -lengthSeeder*.5, -lengthSeeder*0.5},
               length = lengthSeeder}

--! [bounding cube]


--! [spatial objects]

spatial_object = {
  { attribute = { kind = 'seed', },
    geometry = { kind = 'canoND',
                 object = { origin = { 0.0, 0.0, dx*0.5 },
               }
    } -- geometry
  } -- seed
}
--! [spatial objects]


--! [refinement box]

table.insert( spatial_object,
{
    attribute = {
      kind = 'refinement',
      level = refinementLevel+level,
      label='box1'
    },
    geometry = {
      kind = 'canoND',
      object = {
        origin = {start_x, start_y, start_z },
        vec = {{size_x, 0., 0.},
              {0.,size_y, 0.},
              {0., 0., size_z}}
      }
    }
  })

table.insert( spatial_object,
  {
    attribute = {
      kind = 'refinement',
      level = refinementLevel2+level,
      label='box2'
    },
    geometry = {
      kind = 'canoND',
      object = {
        origin = {start2_x, start2_y, start2_z },
        vec = {{size2_x, 0., 0.},
              {0.,size2_y, 0.},
              {0., 0., size2_z}}
      }
    }
  })

--! [refinement box]

--! [walls]

if walls == true then
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

end

--! [walls]


--! [use periodicX boundaries]

if periodicX == true then
  table.insert(spatial_object, { 
    attribute = { 
      kind = 'periodic', 
      label = 'periodicX'
    },
    geometry = {
      kind = 'periodic',
      object = {
        plane1 = {
          origin = { dx+dxDash, -length/2, -length*0.5},
          vec = { 
            { 0.0, 0.0, length },
            { 0.0, length, 0.0 },
            }
        }, -- plane 1
        plane2 = {
          origin = { -dxDash,  -length/2, -length*0.5},
          vec = { 
            { 0.0, length, 0.0 },
            { 0.0, 0.0, length },
            }
        }, -- plane 2
      } -- object
    } -- geometry

  }) 

--! [use periodicX boundaries]
else

--! [inlet and outlet boundary conditions]
  
table.insert(spatial_object,  { 
    attribute = {
      kind = 'boundary',
      label='east' -- outlet
    },
    geometry = {
      kind = 'canoND',
      object = {
        origin = {length*0.5+dxDash, -length*0.5, -length*0.5},
        vec = {{0.0, length, 0.},
              {0.,0.0, length}}
      }
    }
    }) 

table.insert(spatial_object,  { 
    attribute = {
      kind = 'boundary',
      label='west' -- inlet
    },
    geometry = {
      kind = 'canoND',
      object = {
        origin = {-length*0.5-dxDash, -length*0.5, -length*0.5},
        vec = {{0.0, length, 0.},
              {0.,0.0, length}}
      }
    }
  })

end

--! [inlet and outlet boundary conditions]


--! [use periodic boundaries]

if usePeriodic == true then
  table.insert(spatial_object, { 
    attribute = { 
      kind = 'periodic', 
      label = 'periodic'
    },
    geometry = {
      kind = 'periodic',
      object = {
        plane1 = {
          origin = { -length/2, -length/2, dx+dxDash},
          vec = { { length, 0.0, 0.0},
                  { 0.0, length, 0.0},}
        }, -- plane 1
        plane2 = {
          origin = { -length/2,  -length/2, -dxDash},
          vec = { { 0., length, 0.0, 0.0},
                  { length, 0., 0.0},}
        }, -- plane 2
      } -- object
    } -- geometry

  }) 

--! [use periodic boundaries]

--! [non-periodic]

else
  if case2d then
    depth = 4*dx
  else
    depth=height
  end
  table.insert(spatial_object, { 
    attribute = { 
      kind = 'boundary', 
      label = 'top'
    },
    geometry = {
      kind = 'canoND',
      object = {
          origin = { -length/2, -length/2, depth/2+dxDash},
          vec = { { length, 0.0, 0.0},
                  { 0.0, length, 0.0},} 
        } -- object
      } -- geometry 
    } -- spatial object    
  )

  table.insert(spatial_object, {
    attribute = {
      kind = 'boundary',
      label = 'bottom'
    },
    geometry = {
      kind = 'canoND',
      object = {
          origin = { -length/2, -length/2, -depth*0.5-dxDash},
          vec = { { length, 0.0, 0.0},
                  { 0.0, length, 0.0},}
        },-- object
      } -- geometry
    } -- spatial object
  ) 
end

--! [non-periodic]


--! [use of obstacles]

if useObstacle ==true then

  if case2d == false then
    stlfile='stl/cylinder.stl'
    stlLabel='cylinder'
  else
    stlfile='stl/sphere.stl'
    stlLabel='sphere'
  end 

  if qValues == true then
    if stlLabel == 'cylinder' then
      table.insert(spatial_object,  { 
        attribute = { 
          kind = 'boundary', 
          level = maxLevel,
          label = stlLabel, 
          calc_dist = qValues,
        },
        geometry = {
          kind = 'stl', -- was: sphere
          object = { filename = stlfile
  --          { origin = {-length*0.3,-0.01*height,0.},
  --            radius = radius }
          }
        },
       transformation = {
          deformation =  2.0,
          translation =  {-2., 0., 0. }
        }
        }) 
    elseif stlLabel == 'sphere' then
      table.insert(spatial_object,  { 
        attribute = { 
          kind = 'boundary', 
          level = maxLevel,
          label = stlLabel, 
          calc_dist = qValues,
        },
        geometry = {
          kind = 'stl', -- was: sphere
          object = { filename = stlfile
  --          { origin = {-length*0.3,-0.01*height,0.},
  --            radius = radius }
          }
        },
       transformation = {
          deformation =  2.0,
          translation =  {-2., 0., 0. }
        }
        }) 
      -- elseif stlLabel == 'sphere' then
    end -- stlLabel
  else
--! [use of obstacles]
--! [cylinderObj]
    table.insert(spatial_object,  { 
      attribute = { 
        kind = 'boundary', 
        level = maxLevel,
        label = stlLabel, 
        calc_dist = false,
      },
      geometry = {
        kind = stlLabel, -- was: sphere
        object = { filename = stlfile,
         -- { 
            --origin = {-length*0.40,-0.01*height,-lengthSeeder*0.5},
            --vec = {0.0,.0,lengthSeeder}, -- length and axis of the cylinder
            --radius = radius }
        }
      }
      }) 
--! [cylinderObj]
  end -- qValue == true
end -- useObstacle
