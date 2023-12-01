----------------------- PLEASE READ THIS ---------------------------!!!

-- This input file is set up to run for regression check
-- Please make sure you DO NOT MODIFY AND PUSH it to the repository

--------------------------------------------------------------------!!!
-- Geometry information like length, width, height, dx are loaded from seeder 
-- configuration file because they are required by Seeder to generate the mesh
-- and they are also required for Musubi configuration.
require 'seeder'

--! [Local variables]
-- Flow parameters
-- Reynolds number
Re = 100
-- Mach number
Ma = 0.04
-- speed of sound in air [m/s]
cs_phy = 343
-- Density of the fluid [kg/m^3]
rho0_phy = 1.0
-- Inflow velocity at center of pipe computed from Mach number [m/s]
vel_phy = Ma * cs_phy
-- Kinematic viscosity of the fluid computed from Re number [m^2/s]
nu_phy = vel_phy * diameter / Re

-- Ambient pressure
press_ambient = rho0_phy*cs_phy^2
-- Pressure drop across the channel length for Pipe flow m/s kg/m3 1/s m  = 
press_drop = 4.0 * vel_phy * rho0_phy * nu_phy * length / radius^2.0 
-- Pressure gradient
press_grad = press_drop / length

------------ Compute physical timestep from lattice Mach number ---------------
-- Lattice speed of sound
cs_lat = math.sqrt(1.0/3.0)
-- Lattice maximum velocity
vel_lat = Ma * cs_lat
-- Physical timestep computed from physical and lattice velocity
dt = dx * vel_lat / vel_phy
--------------------------------------------------------------------------------

----------------------------- Time settings ------------------------------------
-- Physical simulation end time [s]
tmax_phy = 1
-- Number of iterations required to reach physical simulation end time.
-- tmax_iter is also number of lattice iterations
tmax_iter =  math.ceil(tmax_phy/dt)
-- Interval to check status of the simulation [s]
interval_phy = 1000*dt
-- Starting time for tracking output [s]
trac_start = 0
-- Starting time for restart output [s]
rest_start = tmax_phy
-- Termination wall clock time [s]
wall_clock = 4*60*60-10*60
------------------------- End of time settings ---------------------------------

---------------------------Output directories-----------------------------------
-- tracking folder
tracking_fol = './tracking/'
-- restart folder
restart_fol = 'restart/'
---------------------End of output directories----------------------------------
--! [Local variables]

---------------------------- Lua functions -------------------------------------
--! [Analytical functions]
-- Analytical solutions are taken from 
-- L.D. LANDAU, E.M. LIFSHITZ, Fluid Mechanics (Second Edition),
-- CHAPTER II - VISCOUS FLUIDS,
-- Editor(s): L.D. LANDAU, E.M. LIFSHITZ,  1987
-- https://doi.org/10.1016/B978-0-08-033933-7.50010-6.
-- Analytical solution for velocity-x profile for Poiseuille flow
function vel_analy(x,y,z)
  rSqr = y^2 + z^2
  return ( 1.0 / (4.0*rho0_phy*nu_phy) ) * press_grad * ( radius^2 - rSqr )
end

-- Analytical solution for wall shear stress across the height
function wss_analy(x, y, z)
  r = math.sqrt(y^2 + z^2)
  return 0.25 * press_grad * 2.0 * r 
end

-- Analytical solution for pressure profile across the length
function press_analy(x,y,z,t)
  return press_ambient + press_drop * ( length/2.0-x ) / length
end
--! [Analytical functions]

--! [Inflow velocity profile]
-- Velocity profile for inflow boundary condition
function vel_inflow(x,y,z,t)
  return {vel_analy(x,y,z), 0.0, 0.0} -- parabolic profile for inflow
  --return {vel_phy, 0.0,0.0}    -- block profile for inflow
end
--! [Inflow velocity profile]
-----------------------End of Lua functions ------------------------------------

--------------------------- Musubi configuration -------------------------------
-- Simulation name used in tracking and restart outputs
simulation_name = 'pipe'
-- Print runtime information like memory usage at the end of the simulation
printRuntimeInfo = false
-- file to write measurement time of simulation for each solver step
timing_file = 'mus_timing.res'
-- Location of mesh files
mesh = './mesh/'
-- Logging output from simulation
logging = {
  level = 5,
  --filename = 'log' -- filename to write logging output
}
-- Scaling for multilevel simulation
scaling = 'acoustic'
-- Interpolation method for multilevel simulation
interpolation_method = 'quadratic'

NOdebug = {logging ={level=1, filename='dbg', root_only=false}}

--! [Simulation control]
sim_control = {
  time_control = { 
    max = { sim=tmax_phy, clock=wall_clock }, 
    interval = interval_phy 
  },
  abort_criteria = {
    stop_file = 'stop',
    steady_state = true,
    convergence = {
      variable = {'pressure_phy','vel_mag_phy'}, 
      shape = { 
        kind = 'canoND',
        object = {
          origin = {pipe_center[1], pipe_center[2], pipe_center[3]}
        }
      }, 
      time_control = { min = 0, max = tmax_phy, interval = 10*dt },
      reduction = {'average','average'},
      norm='average', 
      nvals = 100, 
      absolute = true,
      condition = {
        { threshold = 1.e-6, operator = '<=' },
        { threshold = 1.e-6, operator = '<=' }
      }  
    }
  }  
}

--! [Physics parameters]
-- Required to convert physical unit to lattice unit
physics = {
  dt = dt,
  rho0 = rho0_phy
}
--! [Physics parameters]

--! [Scheme identifier]
identify = {
  label = '3D',
  kind = 'fluid_incompressible', -- Physics
  layout = 'd3q19',              -- Stencil
  relaxation = 'bgk'             -- Collision
}
--! [Scheme identifier]

--! [Fluid]
-- Fluid properties
fluid = {
  kinematic_viscosity = nu_phy
}
--! [Fluid]

--! [Initial condition]
initial_condition = {
  pressure = press_ambient,
  velocityX = 0.0,
  velocityY = 0.0,
  velocityZ = 0.0,
}
--! [Initial condition]

--! [Boundary conditions]
-- Label is a boundary identifier and it should be same as in seeder
-- configuration
boundary_condition = {
  {
    label = 'pipe',
    kind = 'wall_libb', -- wall with q-values
  }
}
--! [Boundary conditions]

--! [Source]
glob_source = {
  force = {press_grad, 0.0, 0.0},
  force_order = 2
}
--! [Source]

--! [User defined variables]
-- Mainly used for tracking.
-- This variable can be refered to as variable in boundary condition and source
variable = {
  { 
    name = 'press_grad', 
    ncomponents = 1, 
    vartype = 'st_fun',
    st_fun = press_grad
  },
  { 
    name = 'vel_analy', 
    ncomponents = 1, 
    vartype = 'st_fun',
    st_fun = vel_analy
  },
  { 
    name = 'wss_analy', 
    ncomponents = 1, 
    vartype = 'st_fun',
    st_fun = wss_analy
  },
  { 
    name = 'press_analy', 
    ncomponents = 1, 
    vartype = 'st_fun',
    st_fun = press_analy
  },
  { name = 'press_diff', 
    ncomponents = 1, 
    vartype = 'operation',
    operation = {
      kind = 'difference',
      input_varname = {'pressure_phy','press_analy'}
    }
  },
  { name = 'wss_diff',
    ncomponents = 1, 
    vartype = 'operation',
    operation = {
      kind = 'difference',
      input_varname = {'wss_phy','wss_analy'},
    }
  },
  { name = 'vel_diff', 
    ncomponents = 1,
    vartype = 'operation',
    operation = {
      kind = 'difference',
      input_varname = {'vel_mag_phy','vel_analy'}
    }
  },
}
--! [User defined variables]

--! [Tracking]
tracking = {
--  {
--    label = 'vtk', 
--    folder = tracking_fol,
--    variable = { 'pressure_phy', 'velocity_phy', 'vel_analy', 'press_analy'}, 
--    shape = {kind = 'all'},
--    time_control = {
--      min = tmax_phy, 
--      max = tmax_phy, 
--      interval = interval_phy
--    },
--    output = {format = 'vtk'}      
--  },
  {
    label = 'probeAtInlet', 
    folder = tracking_fol,
    variable = {'pressure_phy', 'velocity_phy'}, 
    shape = {
      kind = 'canoND', 
      object = {
        origin = { inletX+dx, pipe_center[2], pipe_center[3] }
      }
    },
    time_control = {
      min= trac_start, 
      max = tmax_phy, 
      interval = 10*dt
    },
    output = {format = 'ascii'}
  },
  {
    label = 'probeAtOutlet', 
    folder = tracking_fol,
    variable = {'pressure_phy', 'velocity_phy'}, 
    shape = {
      kind = 'canoND', 
      object = {
        origin = { outletX-dx, pipe_center[2], pipe_center[3] }
      }
    },
    time_control = {
      min= trac_start, 
      max = tmax_phy, 
      interval = 10*dt
    },
    output = {format = 'ascii'}
  },
  {
    label = 'probeAtCenter', 
    folder = tracking_fol,
    variable = {'pressure_phy','velocity_phy'}, 
    shape = {
      kind = 'canoND', 
      object = {
        origin = { pipe_center[1], pipe_center[2], pipe_center[3] }
      }
    },
    time_control = {
      min= trac_start, 
      max = tmax_phy, 
      interval = 10*dt
    },
    output = {format = 'ascii'}
  },
  {
    label = 'pressAlongLength', 
    folder = tracking_fol,
    variable = { 'pressure_phy', 'press_analy', 'press_diff' }, 
    shape = {
      kind = 'canoND', 
      object = {
        origin = { inletX, pipe_center[2], pipe_center[3] },
        vec = { length, 0.0, 0.0 },
      }
    },
    time_control = {
      min = tmax_phy, 
      max = tmax_phy, 
      interval = tmax_phy
    },
    output = { format = 'asciiSpatial', use_get_point = false }
  },
  { 
    label = 'velAlongHeight',
    variable = { 'velocity_phy', 'vel_mag_phy', 'vel_analy' },
    folder = tracking_fol,
    shape = {
      kind = 'canoND', 
      object = {
        origin = { pipe_center[1], minY+dx_eps, pipe_center[3] },
        vec = { 0.0, diameter, 0.0 },
      }
    },
    time_control = {
      min = tmax_phy, 
      max = tmax_phy, 
      interval = tmax_phy
    },
    output = {format = 'asciiSpatial'}
  },
  { 
    label = 'planeCenter',
    variable = { 'velocity_phy', 'vel_mag_phy', 'vel_analy' },
    folder = tracking_fol,
    shape = {
      kind = 'canoND', 
      object = {
        origin = { pipe_center[1], pipe_center[2]-radius, pipe_center[3]-radius },
        vec = {
          { 0.0, diameter, 0.0 },
          { 0.0, 0.0, diameter }
        }  
      }
    },
    time_control = {
      min = tmax_phy, 
      max = tmax_phy, 
      interval = tmax_phy
    },
    output = {format = 'asciiSpatial'}
  },
}
--! [Tracking]

--! [Restart]
restart = {
  NOread = restart_fol..'pipe_lastHeader.lua',
  write = restart_fol,
  -- without timeControl restart will be dumped by default at end
  -- of simulation when write restart is set
  time_control = {
    min = rest_start,
    max = tmax_phy,
    interval = tmax_phy
  }
}
--! [Restart]
--------------------------- Musubi configuration -------------------------------
