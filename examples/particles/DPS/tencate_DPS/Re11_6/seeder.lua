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
      label     = 'west', -- at x = 0
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
      label     = 'north',
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
      label     = 'south',
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
      label     = 'top',
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
      label     = 'bottom',
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
