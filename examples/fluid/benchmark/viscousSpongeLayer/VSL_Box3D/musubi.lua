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

------------ Compute physical time step from speed of sound ---- ---------------
dt = csLB/csPhys*dx
--------------------------------------------------------------------------------
nuLB = nuPhys * dt / dx^2
omega = 1.0 / ( nuLB/csLB^2.0 + 0.5 )
-- Thickness of the absorbing layer
abs_thickness = 0.1 --m
-- Damping factor for absorbing layer
damp_factor = 1.0--4.0/omega-0.001 -- Maximum limit: 4.0/omega-0.001

----------------------------- Time settings ------------------------------------
tmax = 250*dt
------------------------- End of time settings ---------------------------------
--! [Local variables]

----------------------- Parameters for acoustic Pulse ..........................
centerX = 0.0
centerY = 0.
centerZ = 0.0
halfwidth = 1.0/20.
amplitude = 1e-3
background = rho0_phy

-- Function for 1D acoustic pulse
function ic_1Dgauss_pulse(x, y, z, t)
  return background+amplitude*math.exp(-0.5/(halfwidth^2)*( x - centerX )^2)
end
-- Function for 2D acoustic pulse
function ic_2Dgauss_pulse(x, y, z, t)
  r = ( x - centerX )^2+( y - centerY )^2
  return (background + amplitude*math.exp(-0.5/(halfwidth^2)*r)) * csPhys * csPhys
end
-- Function for 3D acoustic pulse
function ic_3Dgauss_pulse(x, y, z, t)
  r = ( x - centerX )^2+( y - centerY )^2 + ( z - centerZ )^2
  return (background + amplitude*math.exp(-0.5/(halfwidth^2)*r)) * csPhys * csPhys
end
-------------------------------------------------------------------------------

--------------------------- Musubi configuration -------------------------------
-- Simulation name used in tracking and restart outputs
simulation_name = 'pulse3D'
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
    interval = tmax/10
  },
  abort_criteria = {
    velocity_lat_max = 0.2 -- Maximum lattice velocity permitted
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
  relaxation = {
    name = 'bgk', -- Collision
    variant = 'improved', -- a variant of collision
  },
  layout = 'd3q27'     -- Stencil
}
--! [Scheme identifier]

--! [Fluid]
fluid = {
  kinematic_viscosity = {
    predefined = 'combined',
    temporal = 1.0,
    spatial = {
      predefined = 'viscous_spongelayer_box',
      origin = {-length/2.0+abs_thickness+dx/2.0,
                -length/2.0+abs_thickness+dx/2.0,
                -length/2.0+abs_thickness+dx/2.0
      },
      extent = { length-2*(abs_thickness+dx/2.0),
                 length-2*(abs_thickness+dx/2.0),
                 length-2*(abs_thickness+dx/2.0)
      },
      thickness = abs_thickness,
      damp_factor = 10.0,
      target_state = {
        viscosity = nuPhys
      }
    }
  },
  bulk_viscosity = 2*nuPhys/3.0
}
--! [Fluid]


--! [Initial condition]
initial_condition = {
  pressure = ic_3Dgauss_pulse,
  velocityX = vel_mean,
  velocityY = 0.0,
  velocityZ = 0.0
}
--! [Initial condition]


tracking = {
  {
    label = 'probe',
    variable = {'pressure_phy', 'velocity_phy'},
    shape = {
      kind = 'canoND',
      object = {
        origin = {-length/4.0, 0.0, dx/2.0},
      },
    },
    output = {format='ascii'},
    folder='tracking/',
    time_control = {interval=10*dt, min={iter= 0}, max = tmax}
  },
--  {
--    label = 'vtk',
--    variable = {'pressure_phy', 'kine_viscosity_phy',
--                'density_phy', 'velocity_phy'},
--    shape= { kind = 'all'},
--    output = {format='vtk'},
--    folder='tracking/',
--    time_control = {interval=tmax/10, min={iter= 0}, max = tmax}
--  }
}
