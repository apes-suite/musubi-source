--! [Global variables]
-- Path to stl file
shepherd = true
if shepherd then
  stlFolder = '$!stl_path!$'
else
  stlFolder = './'
end
-- Diameter of the cylinder
dia_cyl = 0.00324 --m
-- Radius of the cylinder
rad_cyl = dia_cyl/2.0

-- Height of the channel
height =  50*dia_cyl --m
-- Number of elements in height
nHeight = 50
-- Length to height ratio
l_h = 1.0 --Ascpect ratio of domain
-- Number of elements in length
nLength = nHeight*l_h
-- Element size
dx = height/nHeight
-- Half of element size
dx_half = dx*0.5
-- Length of the channel
length = l_h*height

-- Position of the cylinder
-- pos_cyl = {0.0, 0.2, 0.0}

-- Number of elements in bounding cube
-- = number of elements in channel + inlet + outlet
nLength_bnd = nLength + 2
-- Level required to reach computed dx
level = math.ceil(math.log(nLength_bnd)/math.log(2))
-- Length of the bounding cube
length_bnd = (2^level)*dx
-- Smallest possible element size
dx_eps = length_bnd/2^20
-- refinement level towards the cylinder 
cylLevel = level + 4
-- element size near cylinder
dx_cyl = length_bnd/2^cylLevel

-- Inlet to cylinder distance
-- Origin of the inlet plane
inlet_origin = {-length/2.0-dx_eps, -height/2.0-dx_eps, 0.0025 - dx_half - dx_eps}
-- Outlet to cylinder distance
outlet_origin = {length/2.0+dx_eps, -height/2.0-dx_eps, 0.0025- dx_half -dx_eps} 
-- Bottom origin
bot_origin = inlet_origin
--  Top origin
top_origin = {-length/2.0-dx_eps, height/2.0+dx_eps, 0.0025- dx_half -dx_eps}

--! [Global variables]


---------------------------Seeder configurations -------------------------------
-- Directory to write mesh
folder = 'mesh/'

-- file to write measurement time of mesh generation each step
timing_file = 'sdr_timing.res'

-- How detailed should the output be?
-- The higher the level, the more output you'll get from seeder.
logging = { level = 5 }

-- Debug outputs
NOdebug = {
  debugMode = true,
  debugFiles = true
}

-- Bounding_cube: two entries: origin and length in this
-- order, if no keys are used
bounding_cube = {
  origin = {-length/2.0 - dx, -height/2.0 - dx, 0.0025 -dx -dx_half },
  length = length_bnd
}

-- A minimum level, by which all parts in the computational domain should at
-- least be resolved with. Default is 0.
minlevel = level

-- *********************** Table of spatial objects ************************* --
-- Each spatial object is defined by an attribute and some geometric entity
-- attached to this attribute. Attributes might be defined multiple times.
-- Attributes are described by a kind (boundary, periodic, seed or refinement),
-- a level and maybe further kind specific values, like a label for the
-- boundary. Geometric objects might by right now:
-- - canoND (point, line, plane or box)
-- - stl
-- - sphere
-- - cylinder
-- - Spacer
-- - triangle
-- - ellipsoid
-- - periodic
--
-- Periodic boundaries are special spatial objects.
-- Attribute with kind = 'periodic' can only make use of geometric objects of
-- kind = 'periodic'.
spatial_object = {
--! [Seed point]
  {
    attribute = {
      kind = 'seed',
    },
    geometry = {
      kind = 'canoND',
      object = {
        origin = {
          -length/2.0 + 2*dx ,
          -height/2.0 + height/4.0,
          0.0025
        }, -- origin
        vec = { length - 4*dx, 0.0, 0.0 }
      } --- object
    }
  },
--! [Seed point]

--! [Refinement around cylinder]
--  {
--    attribute = {
--      kind = 'refinement',
--     level = minlevel+3
--    },
--    geometry = {
--      kind = 'canoND',
--     object = {
--        origin = { 0, -1.3*dia_cyl, 0.0025 -dx_half}, 
--        vec = {{ 15*dia_cyl, 0, 0 }, 
--		{ 0, 2.6*dia_cyl, 0 },
--		{ 0, 0, dx }
--		}
--      }
--    }
--  },

--  {
--    attribute = {
--      kind = 'refinement',
--     level = minlevel+2
--    },
--    geometry = {
--      kind = 'canoND',
--     object = {
--        origin = { 0, -4.5*dia_cyl, 0.0025 -dx_half}, 
--        vec = {{ 30*dia_cyl, 0, 0 }, 
--		           { 0, 9*dia_cyl, 0 },
--		           { 0, 0, dx }
--		    }
--      }
--    }
--  },
--! [Refinement around cylinder]

--! [West boundary]
  {
    attribute = {
      kind = 'boundary',
      calc_dist = false,
      label = 'inlet'     -- label for the boundary at x = -x_offset

    },
    geometry = {
      kind = 'canoND',
      object = {
        origin = inlet_origin,
        vec = {
          { 0.0, height + 2*dx_eps, 0.0 },
          { 0.0, 0.0, dx + 2*dx_eps } 
        }
      }
    }
  },
--! [West boundary]

--! [East boundary]
  {
    attribute = {
      kind = 'boundary',
      calc_dist = false,
      label = 'outlet'     -- label for the boundary at x = length + x_offset
    },
    geometry = {
      kind = 'canoND',
      object = {
        origin = outlet_origin,
        vec = {
          { 0.0, height + 2*dx_eps, 0.0 },
          { 0.0, 0.0, dx + 2*dx_eps }
        }
      }
    }
  },
--! [East boundary]


--! [North boundary]
  {
    attribute = {
      kind = 'boundary',
      calc_dist = false, -- to calculate q-values
      label = 'outlet',   -- label for the boundary at y = height + y_offset
    },
    geometry = {
      kind = 'canoND',
      object = {
        origin = top_origin,
        vec = {
          { length + 2*dx_eps, 0.0, 0.0 },
          { 0.0, 0.0, dx + 2*dx_eps }
        }
      }
    }
  },
--! [North boundary]

--! [South boundary]
  {
    attribute = {
      kind = 'boundary',
      calc_dist = false,
      label = 'outlet',   -- label for the boundary at y = -y_offset
    },
    geometry = {
      kind = 'canoND',
      object = {
        origin = bot_origin,
        vec = {
          { length + 2*dx_eps, 0.0, 0.0 },
          { 0.0, 0.0, dx + 2*dx_eps }
        }
      }
    }
  },
--! [South boundary]


--! [Cylinder boundary]
  {
    attribute = {
      kind = 'boundary',
      label = 'cylinder',
      calc_dist = true,
      flood_diagonal = false,
      level = cylLevel,
      distance_refine = {
        {
          radius = dia_cyl/2,--10 * dx_cyl,
          level_offset = 0
        },
        {
          radius = 2*dia_cyl,
          level_offset = -1
        },
	      {
          radius = 4*dia_cyl,
          level_offset = -2
        },
      	{
          radius = 20*dia_cyl,
          level_offset = -3
        }
      }
    },
    geometry = {
      kind = 'stl',
      object = {
        filename = stlFolder..'Acoustic_cylinder_musubi.stl', 
	      format = 'binary'
      }
    },
  },
--! [Top and bottom boundary (periodic)]
  -- No label required for periodic boundary.
  -- Normal direction of each plane should point outward from fluid domain.
  {
    attribute = {
      kind = 'periodic',
    },
    geometry = {
      kind = 'periodic',
      object = {
        -- Top boundary
        plane1 = {
          origin = {inlet_origin[1], inlet_origin[2], 0.0025 + dx_half + dx_eps},
          vec = {
            { length + 2*dx_eps, 0.0, 0.0 },
            { 0.0, height + 2*dx_eps, 0.0}
          }
        },
        -- Bottom boundary
        plane2 = {
          origin = { inlet_origin[1], inlet_origin[2], 0.0025 -dx_half -dx_eps },
          vec = {
            { 0.0, height + 2*dx_eps, 0.0 },
            { length + 2*dx_eps, 0.0, 0.0}
          }
        }
      }
    }
  }
--! [Top and bottom boundary (periodic)]
} -- spatial object
------------------- End of seeder configurations -------------------------------

