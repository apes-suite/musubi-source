require 'params'

folder    = 'mesh/'
comment   = simulation_name

minlevel  = level

bounding_cube = { origin = bnd_origin,
                  length = bnd_length  }

spatial_object = {
    {
    attribute   = {
      kind      = 'boundary',
      label     = 'vessel',
      level     = minlevel,
    },
    geometry  = {
      kind    = 'stl',
      object  = {filename = 'barton_bfs_vessel.stl'}
    }
  },
    {
    attribute   = {
      kind      = 'boundary',
      label     = 'inlet',
      level     = minlevel,
    },
    geometry  = {
      kind    = 'stl',
      object  = {filename = 'barton_bfs_inlet.stl'}
    }
  },
  {
    attribute   = {
      kind      = 'boundary',
      label     = 'outlet',
      level     = minlevel,
    },
    geometry  = {
      kind    = 'stl',
      object  = {filename = 'barton_bfs_outlet.stl'}
    }
  },

--  {
--    attribute   = {
--      kind      = 'boundary',
--      label     = 'inlet', 
--      level     = minlevel,
--    },
--    geometry  = {
--      kind    = 'canoND',
--      object  = {
--        origin = { dx_eps, -dx_eps, -1-dx_eps},
--        vec = { { 0.0, 0.0, 2+2*dx_eps},
--                { 0.0, 2*h+2*dx_eps, 0.0 }
--        }
--      }
--    }
--  },
--  {
--    attribute   = {
--      kind      = 'boundary',
--      label     = 'outlet', 
--      level     = minlevel,
--    },
--    geometry  = {
--      kind    = 'canoND',
--      object  = {
--        origin = { 38*h-dx_eps, -dx_eps, -1-dx_eps },
--        vec = { { 0.0, 0.0, 2+2*dx_eps},
--                { 0.0, 2*h+2*dx_eps, 0.0 }
--        }
--      }
--    }
--  },
  {
    attribute   = {
      kind      = 'periodic',
      label     = 'periodicZ', -- at x = 0
      level     = minlevel,
    },
    geometry  = { 
      kind    = 'periodic',
      object  = {
        plane1 = {
          origin = { -dx_eps, -dx_eps, dx + dx_eps  },
          vec = { { 0.0, 2*h+2*dx_eps, 0.0},
                  { 38*h+2*dx_eps, 0.0, 0.0 }
          },
        },
        plane2 = {
          origin = { -dx_eps, -dx_eps, -dx_eps  },
          vec = { { 0.0, 2*h+2*dx_eps, 0.0},
                  { 38*h + 2*dx_eps, 0.0, 0.0 }
        },
      }
    }
  },
},
{
    attribute = { 
      kind    = 'seed',
      label   = 'seed',
    },
    geometry  = {
      kind    = 'canoND', 
      object  = { origin = seed_orig }
    }                
  }
}
