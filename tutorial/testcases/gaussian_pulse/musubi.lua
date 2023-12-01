--! [general settings]
simulation_name = 'Gausspulse'
timing_file = 'mus_timing.res'

fluid = { omega = 1.7, rho0 = 1.0 }
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
  kind  = 'lbm',
  layout = 'd3q19',
  relaxation = 'bgk',
}
--! [identify]


--! [geometry]
mesh = {
  predefined = 'cube',
  origin = { 0.0, 0.0, 0.0 },
  length = 10.0,
  refinementLevel = 4
}
--! [geometry]



--! [tracking basics]
tracking = {
  label = 'track_pressure',
--! [tracking basics]
--! [tracking vars]
  variable = { 'pressure', 'velocity' },
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


--! [ref pressure]
rho0 = 1.
cs2 = 1./3.
p0 = rho0*cs2
--! [ref pressure]

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
  time_control = { min = 0, max = 10, interval = 10}  -- timing definitions (either iterations or simulation time)
}
--! [restart settings]
