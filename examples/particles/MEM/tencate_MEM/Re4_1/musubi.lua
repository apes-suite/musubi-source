-- Created: 2021-11-19 9:49AM
-- Test case for a single particle moving around a cubic domain with an initial velocity
-- Size of the domain is L as defined in params.lua
-- Meant to be a simple test case for particles moving across processes and wall interactions

require 'params'
timing_file = 'mus_timing.res'
-------------------------------------------------------------------------------
mesh               = './mesh/'  
-------------------------------------------------------------------------------

tmax               = math.ceil(3.0/dt)
iter_control       = 100

iter_track_main    = 10000
iter_track_slice   = 5000

sim_control        = { 
  time_control     = { 
    min      = { iter = 0        },
    max      = { iter = tmax     },
    interval = { iter = iter_control } 
  }
}
-------------------------------------------------------------------------------
physics   = { dt    = dt,    rho0 = rho_phy }
fluid     = { omega = omega, rho0 = rho_phy, kinematic_viscosity = nu_phy }
identify  = { 
  label      = 'fluid',
  kind       = 'fluid_incompressible',
  relaxation = 'bgk',
  layout     = 'd3q19'
}

particles_start = 1
particles_stop = 1
particles_iter = 1

particles = {
  nParticles = 1,

  kind = 'MEM',

  boundaries = {
    domain_bnd = { 0.0, L, 0.0, W, 0.0, H },
    bnd_kind = {'wall', 'wall', 'wall', 'wall', 'wall', 'wall'}
  },

  nDEMsubcycles = 25,

  particleLogInterval = 1,

  particleBufferSize = 100,

  particle_collision_time = 5*dt,
  particle_collision_tol = 0.0,
  rho0_lat = rho_p0_lat,
  halo_distance = Dia_p_phy,


  position = { { x_p_phy, y_p_phy, z_p_phy, 0.0, 0.0, 0.0} },
 
  velocity = { {0.0, 0.0, 0.0, 0.0, 0.0, 0.0} },

  force    = { {0.0, 0.0, -F_gravity + F_buoyancy, 0.0, 0.0, 0.0} },

  radius = { 0.5*Dia_p_phy },
 
  mass = { m_p_phy }
}

initial_condition = { 
  pressure  = p0, 
  velocityX = 0.0,
  velocityY = 0.0,
  velocityZ = 0.0 
} 

boundary_condition = {  
  { label = 'north', kind = 'wall'},
  { label = 'south', kind = 'wall'},
  { label = 'east', kind = 'wall'},
  { label = 'west', kind = 'wall'},
  { label = 'top', kind = 'wall'},
  { label = 'bottom', kind = 'wall'},
}

-- tracking = {
-- --  {
-- --    label     = 'main',
-- --    folder    = 'tracking/', 
-- --    output    = { format = 'vtk' },
-- --    variable  = {'velocity_phy','density_phy', 'process'},
-- --    shape     = {kind='all'}, 
-- --    time_control = { min = {iter = 5450}, max = {iter = 5550}, interval = {iter = 1} },
-- --  },
--   {
--     label = 'slice',
--     folder = 'tracking/',
--     output = { format = 'vtk' },
--     variable = { 'velocity_phy', 'density_phy', 'process' },
--     shape = { kind = 'canoND',
--               object = {
--                 origin = {0.0, 0.5*W, 0.0 },
--                 vec = {
--                     { 0.0, 0.0, H },
--                     { L, 0.0, 0.0 }
--                 }
--               }
--             },
--    time_control = { min = {iter = 0}, max = {iter = tmax}, interval = {iter = iter_track_slice} },
--   },
-- }

NOrestart = {
  NOread = '',
  NOwrite = 'restart/',
  time_control = { 
    min      = { iter = tmax   }, 
    max      = { iter = tmax    }, 
    interval = { iter = 1    } 
  },
}

