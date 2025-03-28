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
       label     = 'east', -- at x = 0
       level     = minlevel,
     },
     geometry  = { 
       kind    = 'canoND',
       object  = {
         origin = { -dx_eps, -dx_eps, -dx_eps},
         vec = { { 0.0, 0.0, H+2*dx_eps},
                 { 0.0, W+2*dx_eps, 0.0 } 
         }       
       }
     }
   },
   {
     attribute   = {
       kind      = 'boundary',
       label     = 'west', -- at x = L
       level     = minlevel,
     },
     geometry  = { 
       kind    = 'canoND',
       object  = {
         origin = { L+dx_eps, -dx_eps, -dx_eps },
         vec = { { 0.0, 0.0, H+2*dx_eps},
                 { 0.0, W+2*dx_eps, 0.0 } 
         }       
       }
     }
   },
   {
     attribute   = {
       kind      = 'boundary',
       label     = 'north', -- at y = W
       level     = minlevel,
     },
     geometry  = { 
       kind    = 'canoND',
       object  = {
         origin = { -dx_eps, W+dx_eps, -dx_eps  },
         vec = { { 0.0, 0.0, H+2*dx_eps},
                 { L+2*dx_eps, 0.0, 0.0 } 
         }       
       }
     }
   },
   {
     attribute   = {
       kind      = 'boundary',
       label     = 'south', -- at y = 0
       level     = minlevel,
     },
     geometry  = { 
       kind    = 'canoND',
       object  = {
         origin = { -dx_eps, -dx_eps, -dx_eps  },
         vec = { { 0.0, 0.0, H+2*dx_eps},
                 { L+2*dx_eps, 0.0, 0.0 } 
         }       
       }
     }
   },
   {
     attribute   = {
       kind      = 'boundary',
       label     = 'top',  -- at z = H
       level     = minlevel,
     },
     geometry  = { 
       kind    = 'canoND',
       object  = {
         origin = { -dx_eps, -dx_eps, H+dx_eps  },
         vec = { { 0.0, W+2*dx_eps, 0.0},
                 { L+2*dx_eps, 0.0, 0.0 } 
         }       
       }
     }
   },
   {
     attribute   = {
       kind      = 'boundary',
       label     = 'bottom',    -- at z = 0
       level     = minlevel,
     },
     geometry  = { 
       kind    = 'canoND',
       object  = {
         origin = { -dx_eps, -dx_eps, -dx_eps  },
         vec = { { 0.0, W+2*dx_eps, 0.0},
                 { L+2*dx_eps, 0.0, 0.0 } 
         }       
       }
     }
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
