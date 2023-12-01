require "seeder"

intp = 'linear'
intp = 'quadratic'
-- inlet_kind = 'wall'
-- outlet_kind = 'wall'
wall_kind = 'wall'
wall_kind = 'wall_linearInterpolation'
inlet_kind  = 'outlet_expol'
outlet_kind = 'outlet_expol'

scaling='acoustic'
relaxation = 'bgk'
stencil='d3q19'
model='lbm_incomp'
tEnd = 10.0
Re = 60.
initChannel = false
testIntp = true
control_routine = true
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
  relaxation = relaxation}
--! [Simulation name]
--! [reference physical values]
u0Phys    = 1.
viscPhys  = u0Phys*heightPhys/Re
csPhys    = 300    -- m/s
rho0Phys  = 1.
--! [reference physical values]
--! [reference LB values]
csLB    = 1/math.sqrt(3)
cs2     = 1./3.
u0LB    = 0.05
rho0LB  = 1.
rho0LB0 = 0.
omega0  = 1.7
--! [reference LB values]
-- determine the relationships
--! [time step determination]
if scaling == 'acoustic' then
  uLB     = u0LB
  dt      = uLB/u0Phys*dx
  viscLB  = viscPhys*dt/dx/dx
  omega   = 1./(3.*viscLB + 0.5)
else
  omega   = omega0
  viscLB  = 1./6.*(2./omega - 1.)
  dt      = viscLB/viscPhys*dx*dx
  uLB     = u0Phys*dt/dx
end
--! [time step determination]
-- Reynolds number
Re = heightPhys*u0Phys/viscPhys
ReLB = heightLB*uLB/viscLB
if verbose then
  print('           omega   = '..omega)
  print('          viscLB   = '..viscLB)
  print('        heightLB   = '..heightLB)
  print('             uLB   = '..uLB)
  print('Reynolds number Re = '..ReLB)
  print('          viscPh   = '..viscPhys)
  print('        heightPh   = '..heightPhys)
  print('        lengthPh   = '..length)
  print('             uPh   = '..u0Phys)
  print('Reynolds number Re = '..Re)
end


amplitude = u_in
--! [reference pressure]
if model == 'lbm_incomp' then
  p0 = 0.
else
  p0 = rho0LB*cs2*dx*dx/dt/dt
end
--! [reference pressure]

originX = -length*0.25
originY = height*0.
originZ = 00.0
halfwidth = length*0.015
amplitude = 0.00
background = p0
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
  return p0 + dp*0.5 - dp/length*x
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
    if( model=='lbm') then
      return pressureRef(x,y,z) + rho0Phys*cs2*dx/dt
    else
      return pressureRef(x,y,z)
    end
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
    return p0
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
physics = { dt = dt, rho0 = rho0Phys }
--! [physics table]

--! [fluid table]
fluid = {
       omega = omega,
       rho0 = rho0Phys,
  }
--! [fluid table]

--! [initial conditions]
initial_condition = { pressure  = ic_pressure,
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
{  label = 'north',
   kind = 'wall' },
{  label = 'south',
   kind = 'wall' },
{  label = 'west',
   kind = inlet_kind,
   pressure= 'pressureIn'},
{  label = 'east',
   kind = outlet_kind,
   pressure= 'pressureOut'}}
--! [boundary conditions]

if usePeriodic ==false then
  table.insert( boundary_condition,
{ label = 'top',
   kind = 'wall' } )
  table.insert( boundary_condition,
{ label = 'bottom',
   kind = 'wall' } )
end
if useObstacle ==true then
  table.insert( boundary_condition,
{ label = 'sphere',
   kind = wall_kind } )
 end


variable = {
  { name='pressureRef', ncomponents=1, vartype = 'st_fun', st_fun = pressureRef },
  { name='pressureIn', ncomponents=1, vartype = 'st_fun', st_fun = pressureIn },
  { name='pressureOut', ncomponents=1, vartype = 'st_fun', st_fun = pressureOut },
  { name='stressRef',   ncomponents=1, vartype = 'st_fun', st_fun = stressRef },
  { name='error_p',   ncomponents=1, vartype = 'operation',
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
  {
  label = 'probePressure',
 variable = {'pressure_phy'},
 shape = {kind = 'canoND', object = {origin ={0.0,0.,0.} } },
 time_control = {min = {iter=1}, interval = {iter=1}},
  output = {format = 'ascii'}, folder = './tracking/'
 }
--! [Tracking example]
 ,
  {
  -- tracking object for getting the error in terms of the l2norm
  label = 'hvsXY',
  variable = {
    'pressure_phy',
    'velocity_phy',
    'shear_stress_phy',
  },
  shape = {kind = 'canoND', object = {origin ={ -length*0.5,-height*0.5,0.5*dx},
            vec={{length, 0., 0.}, {0.,height,0.}},
            segments = {2*nElemsMax, nElemsMax/2} } },
  time_control = {min = 0., max = tEnd, interval = interval},
  output = {format = 'vtk'},
  folder = tracking_folder
 },
  {
  -- tracking object for getting the error in terms of the l2norm
  label = 'hvsPdf',
  variable = {
    'pdf',
  },
  shape = {kind = 'canoND', object = {origin ={ -length*0.5,-height*0.5,0.5*dx},
            vec={{length, 0., 0.}, {0.,height,0.}},
            segments = {2*nElemsMax, nElemsMax/2} } },
  time_control = {min = 0., max = tEnd, interval = interval},
  output = {format = 'vtk'},
  folder = tracking_folder
 },
 {
  -- tracking object for getting the error in terms of the l2norm
  label = 'dpdx',
  variable = {
    'pressure_phy',
    'pressureRef',
  },
  shape = {kind = 'canoND', object = {origin ={ -length*0.25,-height*0.,0.5*dx},
            vec={{length*0.6, 0., 0.}, {0.,0.,0.}},
            segments = {2*nElemsMax, nElemsMax/2} } },
  time_control = {min = tEnd, max = tEnd, interval = interval},
  output = {format = 'asciiSpatial'},
  folder = tracking_folder
 },
 {
  -- tracking object for getting the error in terms of the l2norm
  label = 'shear_stress',
  variable = {
     'shear_stress_phy', 'stressRef',
  },
  shape = {
    kind = 'canoND',
    object = {
      origin ={ -length*0.,-height*0.5,0.5*dx},
      vec={length*0., height, 0.},
      segments = {2*nElemsMax, 1}
    },
  },
  time_control = {min = tEnd, max = tEnd, interval = interval},
  output = {format = 'asciiSpatial'},
  folder = tracking_folder
 },
 {
  -- tracking object for getting the error in terms of the l2norm
  label = 'velocity',
  variable = {
     'velocity_phy'
  },
  shape = {
    kind = 'canoND',
    object = {
      origin ={ -length*0.,-height*0.5,0.5*dx},
      vec={length*0., height, 0.},
      segments = {2*nElemsMax, 1}
    },
  },
  time_control = {min = tEnd, max = tEnd, interval = interval},
  output = {format = 'asciiSpatial'},
  folder = tracking_folder
 },
 {
  -- tracking object for getting the error in terms of the l2norm
  label = '_Errdpdx',
  variable = {
    'error_p',
  },
  reduction = {'l2norm'},
  shape = {kind = 'canoND', object = {origin ={ -length*0.25,-height*0.,0.5*dx},
            vec={{length*0.6, 0., 0.}, {0.,0.,0.}},
            segments = {2*nElemsMax, nElemsMax/2} } },
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
    min = {iter = tEnd/4},
    max = {iter = tEnd},
    interval = {iter = interval}
  }
}
--! [restart]
