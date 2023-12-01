require "seeder"

--require "PATH"
tracking_folder = './tracking/'--"TRACKING_FOLDER"
output_folder = './output/'--"OUTPUT_FOLDER"

--! [general settings]
mesh = './mesh/' -- Mesh information
simulation_name = 'channelRefine'
-- file to write measurement time of simulation for each solver step
timing_file = 'mus_timing.res'
--! [general settings]
--! [simulation settings]
identify = {label = idLabel, layout='d3q19'}
fluid = { omega = omega,
          rho0 = 1.0 }
--! [simulation settings]

--! [time settings]
tmax =   nIters    -- total iteration number
interval = tmax/30
tRamping = tmax/20
sim_control = {
  time_control = {
    min = {iter=0},
    max = {iter=tmax},
    interval = {iter=interval}
  }
}
--! [time settings]


--! [initial conditions]
initial_condition = { pressure = 1.0,
                      velocityX = 0.0,
                      velocityY = 0.0,
                      velocityZ = 0.0 }
--! [initial conditions]

--! [wall boundary conditions]
boundary_condition = {
{ label = 'obs',
   kind = 'wall' },
{ label = 'north',
   kind = 'wall' },
{ label = 'south',
   kind = 'wall' },
 { label = 'front',
   kind = 'wall' },
{ label = 'back',
   kind = 'wall' } }
--! [wall boundary conditions]

--! [inlet condition]
table.insert( boundary_condition,
{ label = 'west',
  kind = 'inlet_ubb',
  velocity = {
    predefined = 'combined',
    temporal = {
      predefined='smooth',
      min_factor = 0.0,
      max_factor=1.0,
      from_time=0,
      to_time=tRamping},
    spatial = {
      predefined ='parabol',
      shape = {
        kind = 'canoND',
        object = {
          center = {-0.5*length,0.,0.},
          vec = { {0.0, height, 0.}, {0., 0., height}}
        }
      },
      amplitude = {u_in,0.0,0.0}
    }
  },
})
--! [inlet condition]

--! [outlet condition]
table.insert( boundary_condition,
{ label = 'east',
   kind = 'outlet_pab',
   pressure =  1.0}
 )
--! [outlet condition]


--! [tracking basics]
tracking = {
  {
    label = label..'_probePressure_l'..level,
    variable = {'density'},
    shape = {
      kind = 'canoND',
      object = {
        origin = {0.0,0.0,0.0}
      }
    },
    time_control = {
      min = {iter = 0},
      max = {iter = tmax},
      interval = {iter = 5}
    },
    output = {
      format = 'ascii'
    },
    folder = tracking_folder
  },
  {
  -- tracking object for getting the error in terms of the l2norm
    label = label..'_hvsXY_l'..level..interpolation_method,
    variable = {'density', 'velocity', 'wss'},
    shape = {
      kind = 'canoND',
      object = {
        origin ={ -length*0.5,-height*0.5,0.},
        vec={{length, 0., 0.}, {0.,height,0.}},
        segments = {2*nElemsMax, nElemsMax/2}
      }
    },
    time_control = {
      min = {iter = tmax},
      max = {iter = tmax},
      interval = {iter = interval},
    },
    output = {
      format = 'vtk'
    },
    folder = tracking_folder
  }
}
--! [tracking basics]
