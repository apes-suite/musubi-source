require 'params'
-------------------------------------------------------------------------------
mesh               = './mesh/'  
-------------------------------------------------------------------------------

tmax               = 25*t_ref_lat
iter_control       = 1000

iter_track_main    = 10000
iter_track_slice   = 10000

sim_control        = { 
  time_control     = { 
    min      = { iter = 0        },
    max      = { iter = tmax     },
    interval = { iter = iter_control } 
  },
  abort_criteria = {
    stop_file = 'stop'
  }
}
-------------------------------------------------------------------------------
  physics   = { dt    = dt,    rho0 = rho_phy }
fluid     = { omega = omega, rho0 = rho_phy, kinematic_viscosity = nu_phy }
identify  = { 
  label      = 'fluid',
  kind       = 'fluid_GNS',
  relaxation = 'bgk',
  layout     = 'd2q9'
}



initial_condition = { 
  pressure  = 0.0, 
  velocityX = 0.0,
  velocityY = 0.0,
  velocityZ = 0.0 
} 

glob_source = {
  force = { 0.0, -0.0*rho_phy*g, 0.0 },
  force_order = 2
}

function inlet_vel(x,y,z,t)
  ys = y - 1.5*h  -- shifted y-coordinate = 0 at center of parabola
  ux = 1.5*um*(1-ys^2)

  if(ux < 0.0) then
    ux = 0.0
  end

  -- If time is less than ramp time, scale velocity linearly
  if(t < t_ramp_phy) then
    fac = 1.0 - (t_ramp_phy-t)/t_ramp_phy
    ux = ux*fac
  end

  -- y and z components of velocity are 0
  uy = 0.0
  uz = 0.0
  return { ux, uy, uz }
end

boundary_condition = {  
  { label = 'vessel', kind = 'wall'},
  { label     = 'inlet',
    kind      = 'velocity_bounceback',
    velocity  = inlet_vel,
  },
  { label     = 'outlet',
    kind      = 'outlet_dnt',
    pressure  = 0.0,
  },
}

variable = {
} -- variable

tracking = {
  {
    label     = 'main',
    folder    = 'tracking/', 
    output    = { format = 'vtk' },
    variable  = {'velocity_phy','density_phy', 'vol_frac'},
    shape     = {kind='all'}, 
    time_control = { min = {iter = 0}, max = {iter = tmax}, interval = {iter = iter_track_main} },
  },

}
restart = {
  NOread = 'restart/barton_bfs_h64_g_header_37.924E+06.lua',
  write = 'restart/',
  time_control = { 
    min      = { iter = 5*t_ref_lat   }, 
    max      = { iter = tmax          }, 
    interval = { iter = 5*t_ref_lat   } 
  },
}

