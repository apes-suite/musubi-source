----------------------- PLEASE READ THIS ---------------------------!!!

-- This input file is set up to run for regression check
-- Please make sure you DO NOT MODIFY AND PUSH it to the repository

--------------------------------------------------------------------!!!


-- Musubi configuration file.

require "common"
inlet = 'ubb'
outlet = 'expol'--os.getenv("outlet")
if Re == 20 then
  tracking_fol = 'tracking/Re20/'
  restart_fol = 'restart/Re20/'
else
  tracking_fol = 'tracking/Re100/'
  restart_fol = 'restart/Re100/'
end

-- This is a LUA script.
--dx = getdxFromLevel( {len_bnd=length_bnd, level=level})
--dt = getdtFromVel( {dx = dx, u_p = u_in_phy, u_l = u_in_L } )
--omega = getOmegaFromdt( {dx=dx, dt=dt, nu_p = nu_phy } )

-- Simulation name
simulation_name = 'FlowAroundCyl'
timing_file = 'mus_timing.res'
mesh = 'mesh/' -- Mesh information
printRuntimeInfo = false
scaling = 'acoustic'
control_routine = 'fast'
io_buffer_size = 10 -- default is 80 MB

-- Time step settigs
tmax = 8-- real time in seconds
interval = tmax/10
sim_control = {
  time_control = {
    max = tmax,
    interval = interval
  } -- time control
 ,abort_criteria = {
   steady_state = true
 }
} -- simulation control

-- restart
estart = {
      ead = restart_fol..simulation_name..'_lastHeader.lua',
      write = restart_fol,
      time_control = { min = 0, max = tmax, interval = interval}
 }

logging = {level=10}
--debug = { logging = {level=1, filename='debug'},debugMode = true, debugFiles = true,
--          debugMesh = './debug/mesh_', debugStates = {
--  write = {
--    folder    = './debug/',    -- the folder the restart files are written to
--    interval  = 1,           -- dump restart file interval
--    tmin      = 1,           -- first timestep to output
--    tmax      = tmax+1       -- last timestep to output
--    }
-- }
-- }


-- needed to dump variable in physical unit
physics = { dt = dt, rho0 = rho0_p }

fluid = { omega = omega,
          rho0 = rho0_p }

interpolation_method = 'linear'

-- Initial condition
initial_condition = { pressure = p0_p,
                      velocityX = 0.0,
                      velocityY = 0.0,
                      velocityZ = 0.0 }

identify = {label='2D',layout='d3q19',kind='lbm_incomp', relaxation='bgk'}
-- Boundary conditions
boundary_condition = {
{ label = 'west',
  kind = 'inlet_'..inlet,
  velocityX = u_inflow,
  velocityY = 0.0, velocityZ = 0.0
},

{ label = 'east',
--   kind = 'outlet_zero_prsgrd',
--   kind = 'outlet_eq',
--   kind = 'outlet_pab',
--   kind = 'outlet_dnt',
   kind = 'outlet_'..outlet,
   pressure = p0_p },
{ label = 'north',
   kind = 'wall' },
{ label = 'south',
   kind = 'wall' },
{ label = 'sphere',
   kind = 'wall_linearInterpolation', }
 }

--drag and lift coeff factors
cD = 2 / (rho0_p * u_mean_phy * u_mean_phy * Dia * dx)
cL = cD

-- lua variables
variable = {
  { name = 'coeff_fac',
    ncomponents = 3,
    vartype = 'st_fun',
    st_fun = {cD, cL, 0}
  },
  { name  = 'coeff',
    ncomponents = 3,
    vartype = 'operation',
    operation = {
      kind = 'multiplication',
      input_varname = { 'bnd_force_phy', 'coeff_fac'}
    }
  }
}

-- Tracking
tracking = {
-- global flow
{
  label = 'global',
  folder = tracking_fol,
  variable = {
    'pressure_phy',
    'velocity_phy',
  },
  shape = { kind = 'global' },
  time_control = { min = 6, max = tmax, interval = 1/30 },
  format = 'harvester'
},
-- pressure drop at two points over time: P(0.15,0.2,t) and P(0.25,0.2,t)
{
  label = 'probe_pressure',
  folder = tracking_fol,
  variable = {'pressure_phy'},
  shape = {
          {kind = 'canoND', object = {origin ={0.15-dx,0.2,zpos} }},
          {kind = 'canoND', object = {origin ={0.25+dx,0.2,zpos} }}
          },
  time_control = {min = tmax, max = tmax, interval = tmax},
  format = 'ascii'      -- over time
},
{
  label = 'probe_pressure',
  folder = tracking_fol,
  variable = {'pressure_phy'},
  shape = {
          {kind = 'canoND', object = {origin ={0.15-dx,0.2,zpos} }},
          {kind = 'canoND', object = {origin ={0.25+dx,0.2,zpos} }}
          },
  time_control = {min = 0, max = tmax, interval = 10*dt},
  format = 'ascii'      -- over time
},
-- convergence
-- {
--   label = 'convergence',
--   folder = tracking_fol,
--   variable = {'pressure_phy'},
--   shape = {
--           {kind = 'canoND', object = {origin ={0.15-dx,0.2,zpos} }},
--           {kind = 'canoND', object = {origin ={0.25+dx,0.2,zpos} }}
--           },
--   time_control = {min = 0, max = tmax, interval = 10*dt},
--   format='convergence',
--   reduction = 'average',
--   convergence = {norm='average', nvals = 50, absolute = true,
--   condition = { threshold = 1.e-10, operator = '<=' }}
-- },
-- forces at the obstacle
-- {
-- label = 'force',
-- folder = tracking_fol,
-- variable = {'bnd_force_phy','bnd_force'},
-- shape = { kind = 'boundary', boundary = {'sphere'}},
-- time_control = {min = tmax, max = tmax, interval = tmax},
-- reduction = {'sum','sum'},
-- format = 'ascii'      --},
-- {
--  label = 'hvs_force',
--  folder = tracking_fol,
--  variable = {'bnd_force_phy'},
--  shape = { kind = 'boundary', boundary = {'sphere'}},
--  time_control = {min = tmax, max = tmax, interval = tmax},
--  format = 'harvester'
-- },
-- velocity over time
{
 label = 'probe_velocity',
 folder = tracking_fol,
 variable = {'velocity_phy'},
 shape = {
           {kind = 'canoND', object = {origin = {length*0.5,0.2,zpos}}}},
  time_control = {min = 0, max = tmax, interval = 10*dt},
  format = 'ascii'
},
-- velocity at inlet
{
  label = 'probe_velocity_in',
  folder = tracking_fol,
  variable = {'velocity_phy'},
  shape = {kind = 'canoND', object = {origin = {dx_half,0.0,zpos},
                                     vec = {0.0,height,0.0},
                                     segments = nHeight+2}},
  time_control = {min = 0, max = tmax, interval = {iter=100}},
  reduction = {'average'},
  format = 'ascii'
},
-- velocity over height
-- {
--   label = 'line_velocity_in',
--   folder = tracking_fol,
--   variable = {'velocity_phy'},
--   shape = {kind = 'canoND', object = {origin = {dx_half,0.0,zpos},
--                                      vec = {0.0,height,0.0},
--                                      segments = nHeight+2}},
--   time_control = {min = tmax, max = tmax, interval = tmax},
--   format = 'asciiSpatial'
-- },
-- {
--   label = 'line_velocity_inmean',
--   folder = tracking_fol,
--   variable = {'velocity_phy'},
--   shape = {kind = 'canoND', object = {origin = {0.1/2.0,0.0,zpos},
--                                      vec = {0.0,height,0.0},
--                                      segments = nHeight+2}},
--   reduction = 'average',
--   time_control = {min = 0, max = tmax, interval = tmax},
--   format = 'asciiSpatial'
-- },
--
-- {
--   label = 'line_velocity_middle',
--   folder = tracking_fol,
--   variable = {'velocity_phy'},
--   shape = {kind = 'canoND', object = {origin = {length*0.5,0.0,zpos},
--                                      vec = {0.0,height,0.0},
--                                      segments = nHeight+2 }},
--   time_control = {min = tmax, max = tmax, interval = tmax},
--   format = 'asciiSpatial'
-- },
-- velocity over length
{
  label = 'line',
  folder = tracking_fol,
  variable = {'pressure_phy','velocity_phy'},
  shape = {kind = 'canoND', object = {origin = {0.0,height/2.,zpos},
                                     vec = { length, 0.0,0.0},
                                     segments = nLength+2}
                                         },
  time_control = {min = tmax, max = tmax, interval = tmax},
  format = 'asciiSpatial'
}

}
-- drag and lift coefficient at last time step
if Re == 20 then
  table.insert( tracking, {
    label = 'cyl_coeff',
    folder = tracking_fol,
    variable = {'coeff'},
    shape = { kind = 'boundary', boundary = {'sphere'}},
    time_control = {min = tmax, max = tmax, interval = tmax},
    format = 'ascii',
    reduction = {'sum'}})
-- drag and lift coefficient over time
else
  table.insert( tracking,
  {
    label = 'cyl_coeff',
    folder = tracking_fol,
    variable = {'coeff','bnd_force_phy'},
    shape = { kind = 'boundary', boundary = {'sphere'}},
    time_control = { min = 0.0, max = tmax, interval = {iter=100} },
    format = 'ascii',
    reduction = {'sum','sum'}})
end
