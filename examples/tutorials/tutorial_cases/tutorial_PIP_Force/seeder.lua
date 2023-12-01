----------------------- PLEASE READ THIS ---------------------------!!!
-- This input file is set up to run for regression check
-- Please make sure you DO NOT MODIFY AND PUSH it to the repository
--------------------------------------------------------------------!!!

-- set true for recheck. We reduce the number of elements for recheck.
recheck = true
if recheck==true then
  stlFolder = '$!stl_path!$'
else
  stlFolder = './'
end
stlFile = stlFolder..'pipe.stl'

--! [Global variables]
-- Diameter of the pipe [m]
diameter =  0.41
-- Radius of the pipe
radius = diameter/2.0
-- diameter to length ratio
d_l_ratio = 0.20
-- length of the pipe [m]
length = d_l_ratio * diameter
-- Number of elements in diameter
if recheck == true then
  nDiameter = 32
else
  nDiameter = 64
end
-- Number of elements in length
nLength = d_l_ratio * nDiameter
-- Element size
dx = diameter/nDiameter
-- Number of elements in length + 1 element for either side of boundary
nLength_bnd = math.max(nLength + 2, nDiameter + 2)
-- Level required to reach defined dx
level = math.ceil(math.log(nLength_bnd)/math.log(2))
-- Length of the bounding cube
length_bnd = (2^level)*dx
-- Radius of the bounding cube
-- Smallest possible element size
dx_eps = length_bnd/2^20

-- Geometry information of pipe.stl
pipe_center = { length/2, radius, radius }
pipe_center = { 0, 0, 0 }
inletX  = pipe_center[1] - length/2 
outletX = pipe_center[1] + length/2
minY    = pipe_center[2] - radius 
maxY    = pipe_center[2] + radius
minZ    = minY
maxZ    = maxY

--! [Global variables]

---------------------------Seeder configurations -------------------------------
-- Directory to write mesh
folder = 'mesh/'

-- file to write measurement time of mesh generation each step
timing_file = 'sdr_timing.res'

-- How detailed should the output be?
-- The higher the level, the more output you'll get from seeder.
logging = { level = 3 }

-- A minimum level, by which all parts in the computational domain should at
-- least be resolved with. Default is 0.
minlevel = level

-- boundingbox: two entries: origin and length in this
-- order, if no keys are used
bounding_cube = {
  origin = { pipe_center[1] - length_bnd/2, 
             pipe_center[2] - length_bnd/2, 
             pipe_center[3] - length_bnd/2  },
  length = length_bnd
}

-- *********************** Table of spatial objects ************************* --
-- Each spatial object is defined by an attribute and some geometric entity
-- attached to this attribute. Attributes might be defined multiple times.
-- Attributes are described by a kind (boundary, periodic, seed or refinement),
-- a level and maybe further kind specific values, like a label for the
-- boundary. Geometric objects might by right now:
-- - canoND (point, line, plane or box)
-- - stl
spatial_object = {
--! [Seed point]
  {
    attribute = {
      kind = 'seed',
    },
    geometry = {
      kind = 'canoND',
      object = { origin = pipe_center }
    }
  },
--! [Seed point]

--! [Pipe STL]
  {
    attribute = {
      kind = 'boundary',
      label = 'pipe',
      calc_dist = true, -- calculate qValues to increase boundary approximation
      flood_diagonal = true,
    },
    geometry = {
      kind = 'stl',
      object = {
        filename = stlFile,
        format = 'binary'
      }
    },
    transformation = {
      translation = {
        pipe_center[1]-1.4,  -- 1,4 is the offstet to move the stl x0 to stl center
        pipe_center[2],
        pipe_center[3]
                    }, -- translate the geometry
     },
  },
--! [Pipe STL]
--! [West-East boundary as periodic]
  {
    attribute = {
      kind = 'periodic',
    },
    geometry = {
      kind = 'periodic', object = {
        plane1 = {
          origin = { pipe_center[1] - length/2 - dx_eps, 
                     pipe_center[2] - radius - dx,
                     pipe_center[3] - radius - dx       },
          vec = {{ 0.0, 0.0, diameter + 2.*dx },
                 { 0.0, diameter + 2.*dx, 0.0 }}
                 },
        plane2 = {
          origin = { pipe_center[1] + length/2 + dx_eps,
                     pipe_center[2] - radius - dx,
                     pipe_center[3] - radius - dx       },
          vec = {{ 0.0, diameter + 2.*dx, 0.0 },
                 { 0.0, 0.0, diameter + 2.*dx }}
                 }
      }
    }
  }
}
  --! [West-East boundary as periodic]
------------------- End of seeder configurations -------------------------------
