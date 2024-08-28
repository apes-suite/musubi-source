require "seeder"

intp = 'linear'
intp = 'quadratic'
wall_kind = 'wall'
wall_kind = 'wall_libb'
inlet_kind  = 'pressure_expol'
outlet_kind = 'pressure_expol'

scaling='diffusive'
scaling='acoustic'
relaxation = 'bgk'
stencil='d3q19'
model='fluid_incompressible'
tEnd = 10.0
initChannel = false
track = true

interpolation_method = {
  method = intp,
}

--require "PATH"
tracking_folder = './tracking/'--"TRACKING_FOLDER"
--! [Simulation name]
mesh = './mesh/' -- Mesh information
simulation_name = simName
timing_file = 'mus_timing.res'

identify = {
  label = model..relaxation..scaling,
  layout=stencil,
  kind = model,
  relaxation = relaxation
}
--! [Simulation name]
--! [reference physical values]
Ma        = 0.1 -- Mach number
Re        = 60   -- Reynolds number
csPhys    = 300  -- Speed of sound [m/s]
u0Phys    = Ma * csPhys -- Velocity [m/s]
viscPhys  = u0Phys*heightPhys/Re -- Kinematic viscosity [m^2/s]
rho0Phys  = 1. -- Density [kg/m^3]
--! [reference physical values]
--! [reference pressure]
press0Phys = rho0Phys * csPhys^2 -- Ambient pressure [kg/(m s^2)]
--! [reference pressure]
--! [reference LB values]
csLB     = 1/math.sqrt(3) -- Lattice speed of sound
rho0LB   = 1.
press0LB = rho0LB * csLB^2 -- Lattice pressure
omega0   = 1.83 -- relaxation parameter
visc0LB  = 1./6.*(2./omega0 - 1.) -- Lattice viscosity
--! [reference LB values]
-- determine the relationships
--! [time step determination]
if scaling == 'acoustic' then
  uLB     = Ma * csLB
  dt      = csLB / csPhys * dx
  viscLB  = viscPhys * dt / dx / dx
  omega   = 1./(3.*viscLB + 0.5)
else
  omega   = omega0
  viscLB  = 1./6.*(2./omega - 1.)
  dt      = viscLB/viscPhys*dx*dx
  uLB     = u0Phys*dt/dx
end
--! [time step determination]
if verbose then
  print('    Mach number, Ma = '..Ma)
  print('Reynolds number, Re = '..Re)
  print(' ##### Physical values #####')
  print('            viscPhy = '..viscPhys)
  print('          rho0Physy = '..rho0Phys)
  print('        press0Physy = '..press0Phys)
  print('          heightPhy = '..heightPhys)
  print('          lengthPhy = '..length)
  print('               uPhy = '..u0Phys)
  print(' ##### Lattice values  #####')
  print('              omega = '..omega)
  print('             viscLB = '..viscLB)
  print('             rho0LB = '..rho0LB)
  print('           press0LB = '..press0LB)
  print('           heightLB = '..heightLB)
  print('                uLB = '..uLB)
end


amplitude = u0Phys
originX = -length*0.25
originY = height*0.
originZ = 00.0
halfwidth = length*0.015
amplitude = 0.00
background = press0Phys
function ic_2Dgauss_pulse(x, y, z)
  return background+amplitude*math.exp(-0.5/(halfwidth^2)*(( x - originX )^2+( y - originY )^2))
end
-- Reference values for the flow state in the 2d stationary channel at laminar
-- flow state
--! [velocity function]
function uX(x,y,z)
    uPhys = u0Phys*(1-(2.*y/heightPhys)^2)
    return uPhys
end
--! [velocity function]
--! [pressure function]
function pressureRef(x,y,z)
  dp = u0Phys*8.*viscPhys*rho0Phys/heightPhys^2*length
  return press0Phys + dp*0.5 - dp/length*x
end
--! [pressure function]
--! [boundary pressure]
pressureIn=pressureRef(-length/2,0.,0.)
pressureOut=pressureRef(length/2,0.,0.)
--! [boundary pressure]
function Sxx(x,y,z)
  return 0.
end
function Syy(x,y,z)
  return 0.
end
--! [shear stress function]
function Sxy(x,y,z)
  tauxy= -viscPhys*rho0Phys*8./heightPhys^2*u0Phys*y
  S_xy = tauxy/viscPhys/rho0Phys
  return S_xy
end
--! [shear stress function]
function stressRef(x,y,z)
    return Sxy(x,y,z)*rho0Phys*viscPhys
end

-- Consistent initial conditions for the channel
if initChannel then
  function ic_uX(x,y,z)
    return uX(x,y,z)
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
  function ic_uX(x,y,z)
    return 0.
  end
  function ic_pressure(x,y,z)
    return press0Phys
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

init_allElems = false
-- Time step settigs
interval = tEnd/10.
tRamping = tEnd/10.
--! [time settings]
-- simulation will end when simulation wall clock time reaches 1 hr
-- or simulation time reaches tEnd
sim_control = {
  time_control = {
    max = tEnd,
    interval = interval,
    clock = 3600 --s
  },
--! [time settings]
--! [convergence]
  abort_criteria = {
    stop_file    = 'stop',
    steady_state = true,
    convergence = {
      variable = {'pressure_phy'},
      shape = {
        kind = 'canoND',
        object = {{origin = {0.,0.,0.} }}
      },
      time_control = {
        interval = interval/20.,
        min = tEnd/2,
        max = tEnd,
      },
      format = 'convergence',
      reduction = {'average'},
      norm = 'average',
      nvals = 100,
      absolute = false,
      condition = {
        threshold = 1.e-5,
        operator = '<='
      }
    }
  }
}
--! [convergence]



--! [physics table]
physics = { 
  dt = dt,
  rho0 = rho0Phys
}
--! [physics table]

--! [fluid table]
fluid = {
  kinematic_viscosity = viscPhys
}
--! [fluid table]

--! [initial conditions]
initial_condition = {
  pressure  = ic_pressure,
  velocityX = ic_uX,
  velocityY = 0.0,
  velocityZ = 0.0,
  Sxx = 0.,
  Syy = 0.,
  Sxy = ic_Sxy
}
--! [initial conditions]

--! [boundary conditions]
boundary_condition = {
  { 
    label = 'north',
    kind = 'wall'
  },
  {
    label = 'south',
    kind = 'wall'
  },
  {
    label = 'west',
    kind = inlet_kind,
    pressure= 'pressureIn'
  },
  {
    label = 'east',
    kind = outlet_kind,
    pressure= 'pressureOut'
  }
}
--! [boundary conditions]

if usePeriodic ==false then
  table.insert( boundary_condition,
    {
      label = 'top',
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
if useObstacle ==true then
  table.insert(boundary_condition,
    {
      label = stlLabel,
      kind = wall_kind
    }
  )
end


variable = {
  { name='pressureRef', ncomponents=1, vartype = 'st_fun', st_fun = pressureRef },
  { name='pressureIn', ncomponents=1, vartype = 'st_fun', st_fun = pressureIn },
  { name='pressureOut', ncomponents=1, vartype = 'st_fun', st_fun = pressureOut },
  { name='stressRef',   ncomponents=1, vartype = 'st_fun', st_fun = stressRef },
  {
    name='error_p',
    ncomponents=1,
    vartype = 'operation',
    operation = {
      kind = 'difference',
      input_varname = { 'pressure_phy', 'pressureRef',},
    },
  },
}

if track then
-- Tracking
--! [Tracking example]
tracking = {
  {-- tracking object to write pressure at channel center over time
    label = 'probePressure',
    variable = {'pressure_phy'},
    shape = {
      kind = 'canoND',
      object = {
        origin ={0.0, 0., 0.}
      }
    },
    time_control = {min = {iter=1}, interval = {iter=1}},
    output = {format = 'ascii'}, 
    folder = './tracking/'
  }
--! [Tracking example]
 ,
  {
  -- tracking object for Paraview output of macroscopic flow variables on the
  -- channel slice.
    label = 'hvsXY',
    variable = {
      'pressure_phy',
      'velocity_phy',
      'shear_stress_phy',
    },
    shape = {
      kind = 'canoND',
      object = {
        origin ={-length*0.5, -height*0.5, 0.5*dx},
        vec = {
          {length, 0., 0.},
          {0., height, 0.}
        },
        segments = {2*nElemsMax, nElemsMax/2}
      }
    },
    time_control = {min = tEnd, max = tEnd, interval = interval},
    output = {format = 'vtk'},
    folder = tracking_folder
  },
  {
  -- tracking object for Paraview output of pdf variables on the channel slice
    label = 'hvsPdf',
    variable = {'pdf'},
    shape = {
      kind = 'canoND',
      object = {
        origin ={-length*0.5, -height*0.5, 0.5*dx},
        vec = {
          {length, 0., 0.},
          {0., height, 0.}
        },
        segments = {2*nElemsMax, nElemsMax/2}
      }
    },
    time_control = {min = 0., max = tEnd, interval = interval},
    output = {format = 'vtk'},
    folder = tracking_folder
  },
  {
  -- tracking object for getting pressure from simulation and analytic along the
  -- channel center axis at last time step.
    label = 'dpdx',
    variable = {
      'pressure_phy', 'pressureRef',
    },
    shape = {
      kind = 'canoND',
      object = {
        origin ={ -length*0.25, -height*0., 0.5*dx},
        vec = {
          {length*0.6, 0., 0.},
          {0., 0., 0.}
        },
        segments = {2*nElemsMax, nElemsMax/2}
      }
    },
    time_control = {min = tEnd, max = tEnd, interval = interval},
    output = {format = 'asciiSpatial'},
    folder = tracking_folder
  },
  {
  -- tracking object for getting shear stress from simulation and analytic along
  -- channel bottom wall at last time step.
    label = 'shear_stress',
    variable = {
      'shear_stress_phy', 'stressRef',
    },
    shape = {
      kind = 'canoND',
      object = {
        origin ={-length*0., -height*0.5, 0.5*dx},
        vec = {length*0., height, 0.},
        segments = {2*nElemsMax, 1}
      },
    },
    time_control = {min = tEnd, max = tEnd, interval = interval},
    output = {format = 'asciiSpatial'},
    folder = tracking_folder
  },
  {
  -- tracking object for getting velocity from simulation along height
  -- on the channel center at last time step.
    label = 'velocity',
    variable = {'velocity_phy'},
    shape = {
      kind = 'canoND',
      object = {
        origin ={-length*0., -height*0.5, 0.5*dx},
        vec = {length*0., height, 0.},
        segments = {2*nElemsMax}
      },
    },
    time_control = {min = tEnd, max = tEnd, interval = interval},
    output = {format = 'asciiSpatial'},
    folder = tracking_folder
  },
  {
  -- tracking object for getting the error in terms of the l2norm
    label = '_Errdpdx',
    variable = {'error_p'},
    reduction = {'l2norm'},
    shape = {
      kind = 'canoND',
      object = {
      origin = {-length*0.25, -height*0., 0.5*dx},
      vec = {
        {length*0.6, 0., 0.},
        {0., 0., 0.}
      },
      segments = {2*nElemsMax, nElemsMax/2}
      }
    },
    time_control = {min = tEnd, max = tEnd, interval = interval},
    output = {format = 'ascii'},
    folder = tracking_folder
  },
}
else
  tracking = { }
end


--! [restart]
restart = {
  write = 'restart/',
  time_control = {
    min = tEnd/4, -- start to write restart after 1/4 of simulation time
    max = tEnd,
    interval = interval
  }
}
--! [restart]
