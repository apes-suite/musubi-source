----------------------- PLEASE READ THIS ---------------------------!!!

-- This input file is set up to run for regression check
-- Please make sure you DO NOT MODIFY AND PUSH it to the repository

--------------------------------------------------------------------!!!

-- --------------------------
-- Musubi configuration file.
-- --------------------------
--! [basic info]
require "common"

inlet   = 'ubb'
outlet  = 'expol' --os.getenv("outlet")

if Re == 20 then
  tracking_fol = 'tracking/Re20/'
  restart_fol = 'restart/Re20/'
elseif Re == 100 then
  tracking_fol = 'tracking/Re100/'
  restart_fol = 'restart/Re100/'
end

-- ----------------------
-- Simulation information
-- ----------------------

simulation_name   = 'FlowAroundCyl'
timing_file       = 'mus_timing.res'
mesh              = 'mesh/' -- Mesh information
printRuntimeInfo  = false
control_routine   = 'fast'
io_buffer_size    = 10 -- default is 80 MB
hvs_output        = false

--! [basic info]

-- This is a LUA script.
--dx = getdxFromLevel( {len_bnd=length_bnd, level=level})
--dt = getdtFromVel( {dx = dx, u_p = u_in_phy, u_l = u_in_L } )
--omega = getOmegaFromdt( {dx=dx, dt=dt, nu_p = nu_phy } )

--! [time]
-- ------------------
-- Simulation control
-- ------------------

tmax        = 10      -- real time in seconds
interval    = tmax/10
sim_control = {
  time_control = {
    max = tmax,
    interval = interval
  }, -- time control
  abort_criteria = {
    stop_file     = 'stop',
    steady_state  = true,
    convergence   = {
      variable = {'pressure_phy'},
      shape = {
        {kind = 'canoND', object = {origin ={0.15-dx,0.2,zpos} }},
        {kind = 'canoND', object = {origin ={0.25+dx,0.2,zpos} }}
      },
      time_control = {min = 0, max = tmax, interval = 10*dt},
      format    ='convergence',
      reduction = {'average'},
      norm      = 'average',
      nvals     = 50,
      absolute  = true,
      condition = { threshold = 1.e-10, operator = '<=' }
    }
  },
} -- sim_control
--! [time]

--! [restart]
-- -------
-- restart
-- -------

restart = {
      ead = restart_fol..simulation_name..'_lastHeader.lua',
      write = restart_fol,
      time_control = { min = 0, max = tmax, interval = interval}
 }

--! [restart]


logging = {level=9}--10}
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

--! [physics]
-- needed to dump variable in physical unit
physics = { dt = dt, rho0 = rho0_p }

fluid = { omega = omega }
--! [physics]
--! [interpolation]
interpolation_method = 'linear'
--! [interpolation]

--! [initial]
-- -----------------
-- Initial condition
-- -----------------

initial_condition = { pressure = p0_p,
                      velocityX = 0.0,
                      velocityY = 0.0,
                      velocityZ = 0.0 }
--! [initial]
--! [identify]
identify = {
  label='2D',
  layout='d3q19',
  kind='lbm_incomp'
}
--! [identify]
--! [inlet]
-- -------------------
-- Boundary conditions
-- -------------------

boundary_condition = {
{ label = 'west',
  kind = 'inlet_'..inlet,
  velocity = 'inlet_vel'
},
--! [inlet]
--! [outlet]
{ label = 'east',
   kind = 'outlet_'..outlet,
   pressure = 'p0' },
--! [outlet]
--! [wall]
{ label = 'north',
   kind = 'wall' },
{ label = 'south',
   kind = 'wall' },
{ label = 'sphere',
--   kind = 'wall_linearInterpolation', }
  kind = 'wall'},
}
--! [wall]
--! [coeff]
-- ---------------------------
-- drag and lift coeff factors
-- ---------------------------

cD = 2 / (rho0_p * u_mean_phy * u_mean_phy * Dia * dx)
cL = cD
--! [coeff]
--! [variables]
-- -------------
-- lua variables
-- -------------

variable = {
  { name = 'coeff_fac',
    ncomponents = 3,
    vartype = 'st_fun',
    st_fun = {cD, cL, 0}
  },-- coeff_fac
  { name  = 'coeff',
    ncomponents = 3,
    vartype = 'operation',
    operation = {
      kind = 'multiplication',
      input_varname = { 'bnd_force_phy', 'coeff_fac'}
    }
  },-- coeff
  { name = 'vel_x',
    ncomponents = 1,
    vartype = 'st_fun',
    st_fun = {
      predefined = 'combined',
      temporal = {
        predefined='smooth', min_factor = 0.25,
        max_factor=1.0, from_time=0, to_time=1000*dt
      },
      patial = {
        predefined='parabol',
        shape = {
          kind = 'canoND',
          object = {
            origin  = {0.0,0.0,zpos},
            vec = {0.0,height,0.0}
          },
        }, -- shape
        amplitude = 1.0
      },
      spatial = u_inflow
    },-- st_fun
  },-- vel_x
  { name = 'vel_y',
    ncomponents = 1,
    vartype = 'st_fun',
    st_fun = 0.0
  },-- vel_y
  { name = 'vel_z',
    ncomponents = 1,
    vartype = 'st_fun',
    st_fun = 0.0
  },--vel_z
  { name = 'inlet_vel',
    ncomponents = 3,
    vartype = 'operation',
    operation = {
      kind = 'combine',
      input_varname = {'vel_x','vel_y','vel_z'}
    }
  },--inlet_vel
  { name = 'p0',
    ncomponents = 1,
    vartype = 'st_fun',
    st_fun = p0_p,
  },-- p0
}
--! [variables]
--! [pressure]
-- --------
-- Tracking
-- --------

tracking = {
-- pressure drop at two points over time: P(0.15,0.2,t) and P(0.25,0.2,t)
  {
    label = 'probe_pressure',
    folder = tracking_fol,
    variable = {'pressure_phy'},
    shape = {
            {kind = 'canoND', object = {origin ={0.15-dx,0.2,zpos} }},
            {kind = 'canoND', object = {origin ={0.25+dx,0.2,zpos} }}
            },
    time_control = {min = 0, max = tmax, interval = 10*dt},
    output = {format = 'ascii'}      -- over time
  },
}
--! [pressure]
--! [velocity over time]
table.insert( tracking, {
    label = 'probe_velocity',
    folder = tracking_fol,
    variable = {'velocity_phy'},
    shape = {
              {kind = 'canoND', object = {origin = {length*0.5,height/2.,zpos}}}
            },
    time_control = {min = 0, max = tmax, interval = 10*dt},
    output = {format = 'ascii'}
})
--! [velocity over time]
--! [global flow]
if hvs_output then
table.insert( tracking, {
    label = 'global_shape',
    folder = tracking_fol,
    variable = {
      'pressure_phy',
      'velocity_phy',
    },
    shape = { kind = 'global' },
    time_control = { min = tmax - 3, max = tmax, interval = 1/10 },
    output = {format = 'vtk'}
})
end
--! [global flow]
--! [forces at the obstacle]
table.insert( tracking, {
   label = 'force',
   folder = tracking_fol,
   variable = {'bnd_force_phy','bnd_force'},
   shape = { kind = 'boundary', boundary = {'sphere'}},
   time_control = {min = tmax, max = tmax, interval = tmax},
   reduction = {'sum','sum'},
   output = {format = 'ascii'}
})
if hvs_output then
table.insert( tracking, {
    label = 'hvs_force',
    folder = tracking_fol,
    variable = {'bnd_force_phy'},
    shape = { kind = 'boundary', boundary = {'sphere'}},
    time_control = {min = tmax, max = tmax, interval = tmax},
    output = {format = 'vtk'}
})
end
--! [forces at the obstacle]
--! [velocity over height]
table.insert( tracking, {
     label = 'line_velocity_in',
     folder = tracking_fol,
     variable = {'velocity_phy'},
     shape = {kind = 'canoND', object = {origin = {dx_half,0.0,zpos},
                                        vec = {0.0,height,0.0},
                                        segments = nHeight+2}},
     time_control = {min = tmax, max = tmax, interval = tmax},
     output = {format = 'asciiSpatial'}
})
--! [velocity over height]
--! [mean velocity]
table.insert( tracking, {
     label = 'line_velocity_inmean',
     folder = tracking_fol,
     variable = {'velocity_phy'},
     shape = {kind = 'canoND', object = {origin = {0.1/2.0,0.0,zpos},
                                        vec = {0.0,height,0.0},
                                        segments = nHeight+2}},
     reduction = 'average',
     time_control = {min = 0, max = tmax, interval = tmax},
     output = {format = 'asciiSpatial'}
})
--! [mean velocity]
table.insert( tracking, {
     label = 'line_velocity_middle',
     folder = tracking_fol,
     variable = {'velocity_phy'},
     shape = {kind = 'canoND', object = {origin = {length*0.5,0.0,zpos},
                                        vec = {0.0,height,0.0},
                                        segments = nHeight+2 }},
     time_control = {min = tmax, max = tmax, interval = tmax},
     output = {format = 'asciiSpatial'}
})
--! [velocity over length]
table.insert( tracking, {
    label = 'line',
    folder = tracking_fol,
    variable = {'pressure_phy','velocity_phy'},
    shape = {kind = 'canoND', object = {origin = {0.0,height/2.,zpos},
                                       vec = { length, 0.0,0.0},
                                       segments = nLength+2}
                                           },
    time_control = {min = tmax, max = tmax, interval = tmax},
    output = {format = 'asciiSpatial'}
})
--! [velocity over length]


-- drag and lift coefficient at last time step
if Re == 20 then
--! [Re20]
  table.insert( tracking, {
    label = 'cyl_coeff',
    folder = tracking_fol,
    variable = {'coeff'},
    shape = { kind = 'boundary', boundary = {'sphere'}},
    time_control = {min = tmax, max = tmax, interval = tmax},
    output = {format = 'ascii'},
    reduction = {'sum'}})
--! [Re20]
-- drag and lift coefficient over time
elseif Re == 100 then
--! [Re100]
  table.insert( tracking,
  {
    label = 'cyl_coeff',
    folder = tracking_fol,
    variable = {'coeff'},
    shape = { kind = 'boundary', boundary = {'sphere'}},
    time_control = { min = 0.0, max = tmax, interval = {iter=100} },
    output = {format = 'ascii'},
    reduction = {'sum'},
  --  transient_reduction = {'max'}
  })
--! [Re100]
end
