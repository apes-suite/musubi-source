-- Geometry information like length, width, height, dx are loaded from seeder 
-- configuration file because they are required by Seeder to generate the mesh
-- and they are also required for Musubi configuration.
require 'seeder'

--! [Local variables]
-- Flow parameters
Re = 1e7
-- Mach number
Ma = 0.05
-- Physical speed of sound [m/s]
csPhys = 343.2
-- Density of the fluid [kg/m^3]
rho0_phy = 1.0
-- Lattice speed of sound
csLB = 1./math.sqrt(3.)
-- Background mean flow velocity
vel_mean = Ma * csPhys
-- Kinematic viscosity of the fluid [m^2/s]
nuPhys = vel_mean * length / Re
-- Damping factor for absorbing layer
damp_factor = 1.5 -- Maximum limit: 4.0/omega-0.001

------------ Compute physical time step from speed of sound ---- ---------------
dt = csLB/csPhys*dx
--------------------------------------------------------------------------------

----------------------------- Time settings ------------------------------------
-- Time to reach outlet
T0 = length/csPhys
-- Time of wave propagation of one wavelength
tmax = T0
------------------------- End of time settings ---------------------------------
--! [Local variables]

----------------------- Parameters for acoustic source ..........................
amplitude = 1e-3
background = rho0_phy
frequency = 2.0
Tp = 1.0/(frequency*dt)

function bc_acousticLineSrc(x, y, z, t)
  return (background+amplitude*math.sin(2*math.pi*frequency*t*8/T0)) * csPhys^2
end
-------------------------------------------------------------------------------

------------------ Absorbing layer plane as a lua function ---------------------
-- It is not used. Just provided as an exampe
westStart = -length/2.0+abs_thickness+dx/2.0
eastStart = length-abs_thickness-dx/2.0
southStart = -length/2.0+abs_thickness+dx/2.0
northStart = length/2.0-abs_thickness-dx/2.0
function absorbLayer_fun(x,y,z,t)
  fac = 0.0
  if x > eastStart then
    fac = 3125*(abs_thickness+eastStart-x)*(x-eastStart)^4/(256*(abs_thickness)^5)
  end

  sigma_s = damp_factor*fac
  return sigma_s
end
-------------------------------------------------------------------------------

--------------------------- Musubi configuration -------------------------------
-- Simulation name used in tracking and restart outputs
simulation_name = 'acousticLineSrc'
-- Print runtime information like memory usage at the end of the simulation
printRuntimeInfo = false
-- file to write measurement time of simulation for each solver step
timing_file = 'mus_timing.res'
-- Location of mesh files
mesh = './mesh/'
-- Logging output from simulation
logging = {
  level = 3,
  --filename = 'log' -- filename to write logging output
}

-- Debug outputs to write additional information
NOdebug = {
  logging = {
    level = 1,
    filename = 'dbg',
    root_only = false -- all involved MPI processes writes output
  }
}

--! [Simulation control]
sim_control = {
  time_control = {
    max = tmax,
      interval = tmax/50
  }
}
--! [Simulation control]

--! [Physics parameters]
-- Required to convert physical unit to lattice unit
physics = {
  dt = dt,
  rho0 = rho0_phy
}
--! [Physics parameters]

--! [Scheme identifier]
identify = {
  kind = 'fluid',     -- Physics
  relaxation = 'bgk', -- Collision 
  layout = 'd2q9'     -- Stencil
}
--! [Scheme identifier]

--! [Fluid]
fluid = {
  kinematic_viscosity = { 
    predefined = 'combined',
    temporal = 1.0,
    spatial = {
      predefined = 'viscous_spongelayer_plane',
      origin = {length-abs_thickness-dx/2.0, 0.0, dx/2.0},
      normal = {1.0,0.0,0.0},
      thickness = abs_thickness,
      damp_factor = 5.0,
      target_state = {
        viscosity = nuPhys
      }
    }
  }
}
--! [Fluid]


--! [Initial condition]
initial_condition = {
  pressure = background*csPhys^2,
  velocityX = vel_mean,
  velocityY = 0.0,
  velocityZ = 0.0
}
--! [Initial condition]


--! [Boundary condition]
boundary_condition = {
  {
    label = 'west',
    kind = 'pressure_noneq_expol',
    pressure = bc_acousticLineSrc,
  },
  {
    label = 'east',
    kind = 'pressure_noneq_expol',
    pressure = background*csPhys^2,
  },
}
--! [Boundary condition]

--! [Absorb layer as source term]
source = {
  --absorb_layer = absorbLayer_fun
  absorb_layer = 'absorblayer',
  absorb_layer_target = {
    pressure = background*csPhys^2,
    velocity = {vel_mean, 0.,0.},
    --velocity = 'dynamic',
    --nrecord = math.ceil(abs_thickness/dx)
  }
}
--! [Absorb layer as source term]

variable = {
  {
    name = 'absorblayer',
    ncomponents = 1,
    vartype = 'st_fun',
    st_fun = {
      --fun = absorbLayer_fun,
      predefined = 'combined',
      temporal = 1.0,
      spatial = {
        predefined = 'spongelayer_plane',
        origin = {length-abs_thickness-dx/2.0, 0.0, dx/2.0},
        normal = {1.0,0.0,0.0},
        damp_profile = 'polynomial_n6',
        thickness = abs_thickness,
        damp_factor = damp_factor
      },
      shape = {
        inverted = false,
        kind = 'canoND',
        object = {
          origin = {length-abs_thickness-dx/2.0, 0.0, dx/2.0},
          vec = {
            {abs_thickness+dx, 0.0, 0.0},
            {0.0, height, 0.0}
          }
        }
      },
    }
  }
}

tracking = {
--  {
--    label = 'vtk',
--    variable = {'pressure_phy', 'absorblayer', 'density_phy', 'velocity_phy', 
--                'kine_viscosity_phy'
--    },
--    --shape= { kind = 'all'},
--    shape = {
--      kind = 'canoND',
--      object = {
--        origin = {0.0, 0.0, dx/2.0},
--        vec = {
--          {length, 0.0, 0.0},
--          {0.0, height, 0.0}
--        }
--      },
--    },
--    output = {format='vtk'},
--    folder='tracking/',
--    time_control = {interval=tmax/10, min={iter= 0}, max =tmax}
--  },
  {
    label = 'probe',
    variable = {'pressure_phy', 'density_phy'},
    --shape= { kind = 'all'},
    shape = {
      kind = 'canoND',
      object = {
        origin = {dx/2.0, height/2.0, dx/2.0},
      },
    },
    output = {format='ascii'},
    folder='tracking/',
    time_control = {interval=10*dt, min={iter= 0}, max =tmax}
  },
  {
    label = 'line',
    variable = {'pressure_phy', 'density_phy'},
    --shape= { kind = 'all'},
    shape = {
      kind = 'canoND',
      object = {
        origin = {dx/2.0, height/2.0, dx/2.0},
        vec = {length, 0.0, 0.0}
      },
    },
    output = {format='asciiSpatial'},
    folder='tracking/',
    time_control = {interval=tmax, min={iter= 0}, max =tmax}
  }

}
