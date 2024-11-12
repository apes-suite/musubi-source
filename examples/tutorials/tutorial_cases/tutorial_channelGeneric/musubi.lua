-- Geometry information like length, width, height, dx are loaded from seeder
-- configuration file because they are required by Seeder to generate the mesh
-- and they are also required for Musubi configuration.
require "seeder"

--! [Set options]
-- Initialization of channel, for true --> values>0, false --> initalized with 0
initChannel = false
--! [Set options]


--! [Variables]
-- Collision Model (bgk, trt, mrt)
relaxationModel = 'bgk'
-- Stencil (d1q3, d2q9, d3q6, d3q7, d3q19, d3q27)
if case2d then
  stencil = 'd2q9' -- use 2D stencil
else
  stencil = 'd3q19' --use 3D stencil
end
-- Physics (fluid, fluid_incompressible)
physicsModel = 'fluid_incompressible'
--! [Variables]


--! [Local variables]
-- Flow parameters
-- Density of the fluid [kg/m^3]
rho0_phy = 1.0
-- Reynolds number of the flow
Re = 60
-- speed of sound in air [m/s]
cs_phy = 343

------------ Compute physical timestep from lattice Mach number ---------------
-- Lattice speed of sound
cs_lat = math.sqrt(1.0/3.0)
-- In acoustic scaling, Ma is fixed by user
-- Mach number
Ma = 0.05
-- Inflow velocity computed from Ma number [m/s]
vel_phy = Ma * cs_phy
-- Kinematic viscosity of the fluid calculated from Re [m^2/s]
nu_phy = vel_phy * height / Re
-- Lattice velocity
vel_lat = Ma * cs_lat
-- Physical timestep computed from physical and lattice speed of sound
dt  = cs_lat / cs_phy * dx
-- Lattice viscosity
nu_lat  = nu_phy*dt /dx^2
-- Relaxation parameter
omega   = 1.0/(nu_lat/cs_lat^2 + 0.5)

-- Bulk viscosity
if physicsModel == 'fluid' then
  bulk_visc = 2.0/3.0 * nu_phy
end
--------------------------------------------------------------------------------

--! [Reference LB values]
-- Square of lattice speed of sound
cs2_lat  = 1.0/3.0
-- Lattice density
rho0_lat = 1.0
--! [Reference LB values]


--! [Reference pressure in physical unit]
press_ambient = rho0_phy*cs_phy^2
--! [Reference pressure in physical unit]

----------------------------- Time settings ------------------------------------
-- Physical simulation end time [s]
tmax_phy = 10.0
-- Number of iterations required to reach physical simulation end time.
-- tmax_iter is also number of lattice iterations
tmax_iter =  math.ceil(tmax_phy/dt)
-- Interval to check status of the simulation [s]
interval_phy = tmax_phy/10.0
-- Starting time for tracking output [s]
trac_start = 0.0
-- Starting time for restart output [s]
rest_start = tmax_phy/4.0
------------------------- End of time settings ---------------------------------
--! [Local variables]


---------------------------- Lua functions -------------------------------------
--! [Pressure function]
-- Pressure drop across the channel length for Poiseuille flow
press_drop = 8 * vel_phy * rho0_phy * nu_phy * length / height^2
function pressureRef(x,y,z,t)
  return press_ambient + press_drop * (0.5 - x/length)
end
--! [Pressure function]


--! [Boundary pressure]
-- Pressure at inlet
pressureIn=pressureRef(-length/2,0.0,0.0,0.0)
-- Pressure at outlet
pressureOut=pressureRef(length/2,0.0,0.0,0.0)
--! [Boundary pressure]


-- Reference values for the flow state in the 2d stationary channel at laminar
-- flow state
--! [Velocity function]
function velX(x,y,z,t)
  velX_phy = vel_phy * ( 1.0 - ( 2.0*y/height )^2 )
  return velX_phy
end
--! [Velocity function]


function Sxx(x,y,z)
  return 0.0
end

function Syy(x,y,z)
  return 0.0
end

--! [Shear Stress Function]
function Sxy(x,y,z)
  tauxy= -nu_phy*rho0_phy*8./height^2*vel_phy*y
  S_xy = tauxy/nu_phy/rho0_phy
  return S_xy
end
--! [Shear Stress Function]

-- Reference shear stress
function stressRef(x,y,z)
  return Sxy(x,y,z)*rho0_phy*nu_phy
end

-- Consistent initial conditions for the channel
if initChannel then
  function ic_velX(x,y,z)
    return velX(x,y,z,0.0)
  end
  function ic_pressure(x,y,z)
    return pressureRef(x,y,z)
    end
  function ic_Sxx(x,y,z)
    return Sxx(x,y,z)
    end
  function ic_Syy(x,y,z)
    return Syy(x,y,z)
    end
  function ic_Sxy(x,y,z)
    return Sxy(x,y,z)
    end
else
  -- Initialize with 0
  function ic_velX(x,y,z)
    return 0.
  end
  function ic_pressure(x,y,z)
    return press_ambient
  end
  function ic_Sxx(x,y,z)
    return 0.
  end
  function ic_Syy(x,y,z)
    return 0.
  end
  function ic_Sxy(x,y,z)
    return 0.
  end
end
---------------------------- Lua functions -------------------------------------


--------------------------- Musubi configuration -------------------------------
-- Simulation name used in tracking and restart outputs
simulation_name = 'channel'
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
-- Interpolation method for multilevel simulation (linear, quadratic)
interpolation_method = 'quadratic'
-- Debug outputs to write additional information
NOdebug = {
  logging = {
    level = 1,
    filename = 'dbg',
    root_only = false -- all involved MPI processes writes output
  }
}


--! [Simulation Control]
sim_control = {
  time_control = {
    max = { sim = tmax_phy },
    interval = { sim = interval_phy },
    clock = 3600 --s
   },
  -- Abort conditions for simulation
  abort_criteria = {
    -- Abort if file with the name 'stop' exist
    stop_file = 'stop',
    -- Abort if steady state is reached with condition defined in convergence
    -- table
    steady_state = true,
    -- Convergence condition for steady state check
    convergence = {
      -- Variables to check
      variable = {'pressure_phy'},
      -- Check only point in middle of domain
      shape = {
        kind = 'canoND',
        object = {{origin = {0.,0.,0.} }}
      },
      -- Reduce variables values in space to single average value
      reduction = {'average'},
      -- How often to do convergence?
      time_control = {
        min = tmax_phy/2.0,  -- Start convergence check after half sim. time
        max = tmax_phy,      -- DO convergence until end of simulation
        interval = 10*dt     -- Do convernce check every 10*dt [s]
      },
      norm='average',
      nvals = 50,
      absolute = true,
      -- Condition to statisfy to every variable
      condition = {
        { threshold = 1.e-5, operator = '<=' },
      }
    },
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
  layout = stencil,
  relaxation = relaxationModel,
  kind = physicsModel
}
--! [Scheme identifier]


--! [Fluid]
-- For both, incompressible and compressible kinematic viscosity has to be
-- defined. While for the first one, default values are stored for bulk
-- viscosity,, the user has to explicitly give them for compr. fluid.
if (physicsModel == 'fluid') then
  fluid = {
    kinematic_viscosity = nu_phy,
    bulk_viscosity = bulk_visc
  }
else
  fluid = {
    kinematic_viscosity = nu_phy,
  }
end
--! [Fluid]


--! [Initial conditions]
initial_condition = { pressure  = ic_pressure,
                      velocityX = ic_velX,
                      velocityY = 0.0,
                      velocityZ = 0.0,
                      Sxx = ic_Sxx,
                      Syy = ic_Syy,
                      Sxy = ic_Sxy
                      }
--! [Initial conditions]


--! [Boundary conditions]
-- Label is a boundary identifier and it should be same as in seeder
-- configuration.
boundary_condition = {
  {  label = 'north',
     kind  = 'wall'
  },
  {  label = 'south',
     kind  = 'wall'
  },
  {  label = 'west',
     kind  = 'pressure_expol',
     pressure = 'pressureIn'
  },
  {  label = 'east',
     kind  = 'pressure_expol',
     pressure = 'pressureOut'
  }
}

-- If we deactivated usePeriodic in seeder.lua, we have to add boundaries for
-- top and bottom.
if usePeriodic ==false then
  table.insert( boundary_condition,
  {  label = 'top',
     kind = 'wall'
  }
  )
  table.insert( boundary_condition,
  {
    label = 'bottom',
    kind = 'wall'
  }
  )
end

-- If we activated useObstacle in seeder.lua, we have to add a boundary for it,
-- too. Label depends on 2D (cylinder) or 3D (sphere).
if useObstacle ==true then
  table.insert( boundary_condition,
  {
    label = stlLabel,
    kind = 'wall_libb'
  }
  )
 end
--! [Boundary conditions]


--! [User defined variables]
-- Mainly used for tracking.
-- This variable can be refered to as variable in boundary condition and source
variable = {
  -- Reference pressure dependent on physicsModel
  {
    name='pressureRef',
    ncomponents=1,
    vartype = 'st_fun',
    st_fun = pressureRef
  },
  -- Background pressure
  {
    name='pressure0',
    ncomponents=1,
    vartype = 'st_fun',
    st_fun = press_ambient
  },
  -- Pressure at inlet
  {
    name='pressureIn',
    ncomponents=1,
    vartype = 'st_fun',
    st_fun = pressureIn
  },
  -- Pressure at outlet
  {
    name='pressureOut',
    ncomponents=1,
    vartype = 'st_fun',
    st_fun = pressureOut
  },
  -- Reference velocity
  {
    name='velRef',
    ncomponents=1,
    vartype = 'st_fun',
    st_fun = velX
  },
  -- Reference shear stress
  {
    name='stressRef',
    ncomponents=1,
    vartype = 'st_fun',
    st_fun = stressRef },
  -- Difference between numerical pressure and reference pressure
  {
    name='press_error',
    ncomponents=1,
    vartype = 'operation',
    operation = {
      kind = 'difference',
      input_varname = { 'pressure_phy', 'pressureRef',},
    },
  },
  -- Difference between numerical pressure and background pressure
  {
    name='press_fluc',
    ncomponents=1,
    vartype = 'operation',
    operation = {
      kind = 'difference',
      input_varname = { 'pressure_phy', 'pressure0',},
    },
  },
}
--! [User defined variables]


-- Tracking
--! [Tracking example]
tracking = {
--! [Tracking example]
  {
    label = 'vtk', 
    folder = 'tracking/',
    variable = { 'pressure_phy', 'velocity_phy', 'shear_stress_phy'}, 
    shape = {kind = 'all'},
    time_control = {
      min = 0, 
      max = tmax_phy, 
      interval = interval_phy
    },
    output = {format = 'vtk'}      
  },
  -- Track pressure at the center of the channel over time.
  {
    label = 'probeAtCenter',
    folder = 'tracking/',
    variable = {'pressure_phy', 'press_fluc', 'velocity_phy'},
    shape = {
      kind = 'canoND',
      object = {
        origin = { 0.0, 0.0, 0.0 }
      }
    },
    time_control = { min = {iter=1}, interval = {iter=1} },
    output = { format = 'ascii' },
  },
 {
  -- Track pressure profile along the center axis of the channel.
  -- Write this output only at the end of the simulation.
  --label = 'dpdx',
  label = 'pressAlongLength',
  folder = 'tracking/',
  variable = { 'pressure_phy','pressureRef' },
  shape = {
    kind = 'canoND',
    object = {
      origin = { -length*0.5, -height*0.0, 0.25*dx },
      vec = { length, 0.0, 0.0 },
      },
    },
  time_control = { min = tmax_phy, max = tmax_phy, interval = interval_phy },
  output = { format = 'asciiSpatial' },
 },
 {
  -- Track shear stress along the height of the channel at the middle of the
  -- channel length. Write this output only at the end of the simulation.
  label = 'shearAlongHeight',
  folder = 'tracking/',
  variable = { 'shear_stress_phy', 'stressRef', 'shear_mag_phy' },
  shape = {
    kind = 'canoND',
    object = {
      origin = { -length*0.0, -height*0.5, 0.25*dx },
      vec = { length*0.0,  height, 0.0 },
    },
  },
  time_control = { min = tmax_phy, max = tmax_phy, interval = interval_phy },
  output = { format = 'asciiSpatial' },
 },
 {
  -- Track shear stress along the height of the channel at the middle of the
  -- channel length. Write this output only at the end of the simulation.
  label = 'wssAlongHeight',
  folder = 'tracking/',
  variable = { 'wss_phy'},
  shape = {
    kind = 'canoND',
    object = {
      origin = { -length*0.0, -height*0.5, 0.25*dx },
      vec = { length*0.0,  height, 0.0 },
    },
  },
  time_control = { min = tmax_phy, max = tmax_phy, interval = interval_phy },
  output = { format = 'asciiSpatial' },
 },

 {
  -- Track velocity along the height of the channel at the middle of the channel
  -- length. Write this output only at the end of the simulation.
  label = 'velAlongHeight',
  folder = 'tracking/',
  variable = { 'velocity_phy', 'velRef' },
  shape = {
    kind = 'canoND',
    object = {
      origin = { -length*0.0, -height*0.5, 0.25*dx },
      vec = { length*0.0, height, 0.0 },
    },
  },
  time_control = { min = tmax_phy, max = tmax_phy, interval = interval_phy },
  output = { format = 'asciiSpatial' },
 },
 {
  -- Calculate L2-norm of press_error variable at the end of the simulation.
 -- label = '_Errdpdx',
  label = 'press_l2norm',
  folder = 'tracking/',
  variable = { 'press_error' },
  shape = {
    kind = 'canoND',
    object = {
      origin = { -length*0.5 ,-height*0.0 ,0.25*dx },
      vec = {
        { length, 0.0, 0.0},
        {0.0 ,0.0, 0.0}
      },
    }
  },
  time_control = { min = tmax_phy, max = tmax_phy, interval = interval_phy },
  reduction = {'l2norm'},
  output = { format = 'asciiSpatial' },
 },
}


--! [Restart]
-- Without timeControl restart will be dumped by default at end
-- of simulation when write restart is set.
restart = {
  NOread = 'restart/channel2D_lastHeader.lua',
  write = 'restart/',
  time_control = {
    min = rest_start,
    max = tmax_phy,
    interval = interval_phy
  }
}
--! [Restart]
