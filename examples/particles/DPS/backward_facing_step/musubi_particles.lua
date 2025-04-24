require 'params'
-------------------------------------------------------------------------------
mesh               = './mesh/'  
-------------------------------------------------------------------------------
timing_file        = 'mus_timing.res'

tmax               = 25*t_ref_lat + 2*t_ref_lat
iter_control       = math.floor(t_ref_lat/20)

iter_track_main    = math.floor(t_ref_lat/10)
iter_track_slice   = math.floor(t_ref_lat/10)

iter_track_zoom_start = math.floor(25.4*t_ref_lat)
iter_track_zoom_end   = iter_track_zoom_start + 100


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

nP = 10         -- Number of particles
particles = {

  kind = 'DPS', 

  nParticles = nP,

  -- No boundaries table present means particles are removed from the global 
  -- domain once they hit a domain boundary
  
  interpolation_kind = 'delta',
  particle_collision_time = 0.1*dt,
  particle_collision_tol = 0.0,

  particleLogInterval = 1,

  particleBufferSize = 100,

  position = { 
  },
 
  velocity = {
             },

  force = { 
          },

  radius = { 
           },
 
  mass = { 
         },
}

xp = 3*dx      -- particle starting x-coordinate, same for all particles
yp = 1.05*h    -- starting y-coordinate of first particle 
zp = 0.5*dx    -- starting z-coordinate of all particles

for iP = 1, nP do
  -- Calculate velocity of particle
  up = inlet_vel(xp, yp, zp, 2*t_ramp_lat*dt)

  table.insert(particles['position'], {xp, yp, zp, 0.0, 0.0, 0.0})
  table.insert(particles['velocity'], { up[1], up[2], up[3], 0.0, 0.0, 0.0})
  table.insert(particles['force'], {0.0, Fbuoyancy, 0.0, 0.0, 0.0, 0.0})
  table.insert(particles['radius'], 0.5*dp)
  table.insert(particles['mass'], mp)

  -- increment y-position for next particle
  yp = yp + 0.1*h
end



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
  {
    label     = 'zoom',
    folder    = 'tracking/', 
    output    = { format = 'vtk' },
    variable  = {'velocity_phy','density_phy', 'vol_frac'},
    shape     = {kind='all'}, 
    time_control = { min = {iter = iter_track_zoom_start}, 
                     max = {iter = iter_track_zoom_end}, 
                     interval = {iter = 1} },
  },

}
restart = {
  read = 'restart/barton_bfs_h64_nog_header_18.962E+06.lua',
  NOwrite = 'restart/',
  time_control = { 
    min      = { iter = 5*t_ref_lat   }, 
    max      = { iter = tmax          }, 
    interval = { iter = 5*t_ref_lat   } 
  },
}

