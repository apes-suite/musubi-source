--! [interpolation settings]
interpolation_method = 'linear'
--! [interpolation settings]

-- common header --
--! [minlevel]
-- A minimum level, by which all parts in the computational domain should at
-- least be resolved with. Default is 0.
minlevel = 7 
--! [minlevel]
--
--! [refinement level]
refinementLevel = 1  
--! [refinement level]
label = 'channel'
idLabel = 'channel'
u_in = 0.02
--! [reference length]
length = 2.05
height =  3./16.*length
--! [reference length]
radius = 0.05
usePeriodic = false
useObstacle = true
useCylinder = true
cylOrientation = 'Y'
testIntersection = false

refOmega = 1.9
refLevel = 7

--
m = 2^(refLevel - minlevel)
omega = 1./(1/m*(1/refOmega -0.5)+0.5)
maxlevel = minlevel+refinementLevel

dx     = length/2^minlevel
dxDash = 0.5*dx
viscLB = 1/3*(1/omega-1/2)
--! [physical values]
dxPhys = dx
lPhys = length  -- m
csPhys = 300    -- m/s
csLB = 1/math.sqrt(3)
--acoustic scaling with speeds of sound as reference
-- csPhys = csLB * dxPhys / dtPhys
dtPhys = csLB/csPhys*dxPhys
rho0 = 1  -- kg/m^3
viscPhys = viscLB*dxPhys^2/dtPhys
--! [physical values]

--! [iterations]
nElemsMax = 2^minlevel
nIters = nElemsMax*50
--! [iterations]

-- Use this file as template. Do not modify this file for running some testcases
outputname = 'channel_2D'
folder = 'mesh/' 

seederLength = length+dx

--! [patch size]
size_x = 0.63*length
start_x = -0.42*length
size_y = height/2-2.*dxDash
start_y = -size_y/2
size_z  = 3*size_y
start_z = 3*start_y
--! [patch size]

-- boundingbox: two entries: origin and length in this
-- order, if no keys are used
--! [bounding box]
bounding_cube = {origin = {-seederLength*0.5, -seederLength*.5, -seederLength*0.5},
               length = seederLength}
--! [bounding box]

spatial_object = {
--! [refinement box]
  {
    attribute = {
      kind = 'refinement',
      level = maxlevel,
      label='box'
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
  },
--! [refinement box]
--! [seed definition]
  { attribute = { kind = 'seed', },
    geometry = { kind = 'canoND',
                 object = { origin = { 0.0, 0.0, dx*0.5 },
               }
    } -- geometry
  }, -- seed

}
--! [seed definition]

--! [geometry walls]
    table.insert(spatial_object, {
         attribute = {
           kind = 'boundary',
      	   label='east'
         },
         geometry = {
           kind = 'canoND',
           object = {
             origin = {length*0.5,-length*0.5,-length*0.5},
             vec = {{0.0,length,0.0},
                    {0.0,0.0,length}},
           }
         }
       })
    table.insert(spatial_object, {
         attribute = {
           kind = 'boundary',
      	   label='west'
         },
         geometry = {
           kind = 'canoND',
           object = {
             origin = {-length*0.5,-length*0.5,-length*0.5},
             vec = {{0.0,length,0.0},
                    {0.0,0.0,length}},
           }
         }
       })
    table.insert(spatial_object, {
         attribute = {
           kind = 'boundary',
      	   label='north'
         },
         geometry = {
           kind = 'canoND',
           object = {
          	 origin = {-length*0.5,height*0.5+dxDash,-length*0.5},
             vec = {{length,0.0,0.0},
                    {0.0,0.0,length}},
           }
         }
       })
    table.insert(spatial_object, {
         attribute = {
           kind = 'boundary',
      	   label='south'
         },
         geometry = {
           kind = 'canoND',
           object = {
          	 origin = {-length*0.5,-height*0.5-dxDash,-length*0.5},
             vec = {{length,0.0,0.0},
                    {0.0,0.0,length}},
           }
         }
       })
--! [geometry walls]
--! [geometry cylinder]
    table.insert(spatial_object, {
         attribute = {
           kind = 'boundary',
      	   label='obs'
         },
         geometry = {
           kind = 'stl',
           object = { filename = 'cylinder.stl'}
          }
        })
--! [geometry cylinder]

--! [front back planes]
    table.insert(spatial_object, {
         attribute = {
           kind = 'boundary',
      	   label='front'
         },
         geometry = {
           kind = 'canoND',
           object = {
      	     origin = {-length*0.5,-length*0.5,height*0.5+dxDash},
             vec = {{length,0.0,0.0},
                    {0.0,length,0.0}},
           }
         }
       })
    table.insert(spatial_object, {
         attribute = {
           kind = 'boundary',
      	   label='back'
         },
         geometry = {
           kind = 'canoND',
           object = {
      	     origin = {-length*0.5,-length*0.5,-height*0.5-dxDash},
             vec = {{length,0.0,0.0},
                    {0.0,length,0.0}},
           }
         }
       })
--! [front back planes]
