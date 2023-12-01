-- Geometry information like length, width, height, dx are loaded from seeder 
-- configuration file because they are required by Seeder to generate the mesh
-- and they are also required for Musubi configuration.
require 'seeder'

------------ Compute physical time step from lattice Mach number ---------------
divider = 1 -- 1 mean M=0.2, 5 M=0.04
-- Lattice Mach number
Ma = 0.2/divider
-- Lattice speed of sound
cs_lat = math.sqrt(1.0/3.0)
-- Lattice maximum velocity
vel_lat = Ma * cs_lat

--! [Local variables]
-- Flow parameters
-- Reynolds number of the flow
Re = 68620 -- not used
-- Density of the fluid [kg/m^3]
rho0_phy = 1.0
-- Speed of sound [m/s]
cs_phy = 343.0
-- Mean inflow velocity computed from Reynolds number [m/s]
vel_phy = Ma * cs_phy
-- Kinematic viscosity of the fluid [m^2/s]
nu_phy = 2e-4 --vel_phy * nLength / length / Re

-- Ambient pressure
press_ambient = rho0_phy * cs_phy^2

-- Physical time step computed from physical and lattice velocity
dt = dx * cs_lat / cs_phy
--------------------------------------------------------------------------------

----------------------------- Time settings ------------------------------------
-- Physical simulation end time [s]
tmax_phy = .04*divider
-- Number of iterations required to reach physical simulation end time.
-- tmax_iter is also number of lattice iterations
tmax_iter =  math.ceil(tmax_phy/dt)
-- Interval to check status of the simulation [s]
interval_phy = tmax_phy/5.
-- Starting time for tracking output [s]
trac_start = 0
-- Starting time for restart output [s]
rest_start = 0
------------------------- End of time settings ---------------------------------
--! [Local variables]

--------------- Vortex parameters ----------
Rc = 0.1 -- [m]
Eps = 0.14 -- nondimensional vortex strength
EpsCs = Eps * cs_phy
U0 = .0
V0 = Ma * cs_phy
-- center of the vortex
P = { 2.5, 1.375 } -- [m]


---------------------------- Lua functions -------------------------------------
--! [Analytical functions] in SI units
function r2(x,y)
  return ( x - P[1] )^2 + ( y - P[2] )^2
end

-- pseudo-isoentropic
--function u0_analy(x,y,z)
--  Ueddy = -EpsCs / Rc * ( y - P[2] ) * math.exp( 0.5 - ( r2(x,y) / ( 2.0 * Rc * Rc ) ) )
--  return U0 + Ueddy 
--end
--
--function v0_analy(x,y,z)
--  Veddy = EpsCs / Rc * ( x - P[1] ) * math.exp( 0.5 - ( r2(x,y) / ( 2.0 * Rc * Rc ) ) )
--  return V0 + Veddy
--end
--
--function rho0_analy(x,y,z)
--  rho1 = - ( ( rho0_phy * Eps * Eps ) / 2.0) * math.exp ( 1.0 - ( r2(x, y) / ( Rc * Rc ) ) )
--  return rho0_phy + rho1 + rho1^2 / 2.8 / rho0_phy
--end

function u0_analy(x,y,z)
  Ueddy = -EpsCs / Rc * ( y - P[2] ) * math.exp( - r2(x,y) / ( 2.0 * Rc * Rc ) )
  return U0 + Ueddy
end

function v0_analy(x,y,z)
  Veddy =  EpsCs / Rc * ( x - P[1] ) * math.exp( - r2(x,y) / ( 2.0 * Rc * Rc ) )
  return V0 + Veddy
end

function rho0_analy(x,y,z)
  rho1 = - Eps * Eps / 2.0 * math.exp ( - r2(x, y) / ( Rc * Rc ) )
  return rho0_phy * math.exp(rho1)
end

function p0_analy(x,y,z)
  rho0_a = rho0_analy(x,y,z)
  return rho0_a * cs_phy^2
end

--------------------------- Musubi configuration -------------------------------
-- Simulation name used in tracking and restart outputs
simulation_name = 'vortex'
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
-- Scaling for multilevel simulation
scaling = 'acoustic'
-- Interpolation method for multilevel simulation
interpolation_method = {
   method = 'quadratic',
   type_of_interpolation = 'pdf', --'moment'
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
    max = { sim = tmax_phy },
    interval = { sim = interval_phy }
   },
  -- Abort conditions for simulation
  abort_criteria = {
    -- Abort if file with the name 'stop' exist
    stop_file = 'stop',
    velocity_lat_max = 0.5
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
  layout = 'd3q27',    -- Stencil
  relaxation = 'hrr_bgk', -- Collision 
  kind = 'fluid'      -- Physics
}
--! [Scheme identifier]

--! [Fluid]
fluid = {
  kinematic_viscosity = nu_phy,
  bulk_viscosity = 2.0*nu_phy/3.0,
  hrr_sigma = 0.98,
  drt_taun = 0.7
}
--! [Fluid]

--! [Initial condition]
initial_condition = {
  pressure = p0_analy,
  velocityX = u0_analy,
  velocityY = v0_analy,
  velocityZ = 0.0
}
--! [Initial condition]


--! [Tracking]
tracking = {
  -- Output file to visualize in Paraview.
  {
    label = 'vtk',
    folder = 'tracking/',
    variable = { 'velocity_phy', 'density_phy', 'grad_velocity' },
    shape = { kind = 'all' },
    time_control = { min= 0, max = tmax_phy, interval = tmax_phy/5. },
    output = { format = 'vtk' }
  }
}
--! [Tracking]

--! [Restart]
-- Without timeControl restart will be dumped by default at end
-- of simulation when write restart is set
restart = {
  NOread = 'restart/'..simulation_name..'_lastHeader.lua',
  write = 'restart/',
}
--! [Restart]

--------------------------- Musubi configuration -------------------------------
