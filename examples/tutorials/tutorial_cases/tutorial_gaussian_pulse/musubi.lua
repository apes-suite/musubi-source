--! [general settings]
simulation_name = 'Gausspulse'
timing_file = 'mus_timing.res'

rho0 = 1.0
cs2=1/3.0
--! [ref pressure]
p0 = rho0*cs2
--! [ref pressure]

nu_phy  = 0.03 -- m^2 / s

fluid = {
  kinematic_viscosity = nu_phy,
  bulk_viscosity = nu_phy
}
--! [general settings]

--! [time_control settings]
sim_control = {
  time_control = {
    max = {iter=50},
    interval = {iter=5}
  }
}
--! [time_control settings]

--! [identify]
identify = {
  label = '',
  kind  = 'fluid',
  layout = 'd3q19',
  relaxation = 'bgk',
}
--! [identify]


--! [geometry]
mesh = {
  predefined = 'cube',
  origin = { 0.0, 0.0, 0.0 },
  length = 10.0,
  refinementLevel = 6
}
--! [geometry]



--! [tracking basics]
tracking = {
  label = 'track_pressure',
--! [tracking basics]
--! [tracking vars]
  variable = { 'density', 'pressure', 'velocity' },
  folder = './tracking/',
--! [tracking vars]
--! [tracking shape]
  shape = {
    kind = 'canoND',
    object = {origin = {1.0, 1.0, 1.0} }
  },
--! [tracking shape]
--! [tracking format]
  output = {format = 'ascii'},
--! [tracking format]
--! [tracking time_control]
  time_control = { min = {iter=1}, max = {iter=50}, interval = {iter=1} },
}
--! [tracking time_control]


--! [ic function]
function gausspulse(x, y, z)
  originX =  5.0
  halfwidth = 1.0
  amplitude = 0.01
  return p0+amplitude*math.exp(-0.5/(halfwidth^2)*( x - originX )^2)
end
--! [ic function]

--! [initial conditions]
initial_condition = {
  velocityX = 0.0,
  velocityY = 0.0,
  velocityZ = 0.0,
  pressure  = gausspulse
}
--! [initial conditions]

--! [restart settings]
restart = {
  write = 'restart/',   -- prefix to write the files to
  time_control = { min = {iter=50}, max = {iter=50}, interval = {iter=50}}  -- timing definitions (either iterations or simulation time)
}
--! [restart settings]
