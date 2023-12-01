-- This is a sample Musubi configuration file in lua format
-- 
-- WARNING!!!
-- LEAVE THIS FILE UNCHANGED IN THE REPOSITORY
-- 
-- Thank you!
--
--
io_buffer_size = 1

originX =  -1.3
originY =  0.8
originZ =  0.1
halfwidth = 0.50
amplitude = 0.01
background = 1.0

function ic_1Dgauss_pulse(x, y, z)
  return background+amplitude*math.exp(-0.5/(halfwidth^2)*( x - originX )^2)
end
function ic_2Dgauss_pulse(x, y, z)
  return background +amplitude*math.exp(-0.5/(halfwidth^2)*(( x - originX )^2+( y - originY )^2))
end
-- Simulation name
simulation_name = 'Gausspulse'
--mesh = 'testsuite/physics/lidcavity/mesh/'
mesh = { predefined='cube', 
         origin = {0.,0.,0.}, 
         length = 10., 
         refinementLevel = 4 }

-- Interpolation method
-- average, copyfirst, eqneq, debug
interpolation_method = 'average'                   

-- Time step settings
tmax           =  50     -- total iteration number
time = {useIterations = true,
        min = 1, max = tmax, interval = tmax }

-- scheme model for single fluid simulation
scheme = {
  identify = {
      label = 'default',
      kind = 'lbm', 
      relaxation = 'bgk', 
-- scheme layout
      layout = 'd3q19'},
  -- field which defines fluid or specie
  -- Single fluid simulation
  field = {
--     label = 'fluid',
-- Fluid properties
     fluid = { omega = 1.7, rho0 = 1.0 },
-- Initial condition for each field
     initial_condition = { 
-- initial with lua function
              density = ic_1Dgauss_pulse,

-- possible velocity definitions
              -- const value for each
              velocityX = 0.0,
              velocityY = 0.0,
              velocityZ = 0.0}
              --,
           -- vectorial lua function
              --velocity = lua_func_vec
           -- lua function for each
              --velocity = {lua_fun1, lua_fun2, lua_fun3}
           -- lua function for one and constant for others
              --velocity = {lua_fun, const, const}
           -- Single Transient and spatial. vector val is given by array of ref_value
              --velocity = {ref_value = {1.0,2.0,3.0}, 
              --            transient ,
              --            spatial }
           -- Individual Transient and spatial 
              --velocity = { {transient, spatial, ref_value}, const, lua_fun}
              }
-- boundary condition for each field              
--    boundary_condition = {{}}
, tracking = {
   label = 'point', 
   variable = {'density'},
   folder = 'tracking/',
   shape = {kind = 'canoND', object = {origin = {0.0,0.0,0.0}} },
   format = 'ascii',
   time = {min = 0, max = 10, interval = 1}, 
  }
    }
--}

-- scheme model for multispecies simulation
--scheme = {
--  identifier = {
--      label = 'default',
--      kind = 'LBM_multispecies', 
--      relaxation = 'BGK', 
---- scheme layout
--      layout = 'd3q19'},
--      depend = { vay_sys = {variable={'velocity'}}, usescheme='fluid'},
--  field = {
--   {
--   label = 'spc1'
--   species = { 
--        molecular_weight = 1.0,
--        diff_coeff = 1.0
--       }
--    -- Initial condition 
--     initial_condition = { 
---- initial with lua function
--              density = ic_2Dgauss_pulse
--              velocity = {0.0,0.0,0.0}}
---- boundary_condition = {{}}
--   
--   },
--   {
--   label = 'spc2'
--   species = { 
--        molecular_weight = 2.0,
--        diff_coeff = 1.0
--       }
---- Initial condition 
--     initial_condition = { 
---- initial with lua function
--              density = ic_2Dgauss_pulse
--              velocityX = 0.0,
--              velocityY = 0.0,
--              velocityZ = 0.0 }
---- boundary_condition = {{}}
--
--   }
--  }
--}
--
---- Output settings
--output = { active = true, -- VTK output activated?
--  folder = 'output/',     -- Output location
--  vtk = true,             -- VTK output activated?
--  time = {min = 0, max = -1, interval = 1}
--  } 
--
---- Tracking              
--tracking = {
--  label = 'probe_press', 
--   variable = {'density'},   -- options: density, velocity
--   shape={kind = 'canoND', object = {origin ={3.0,0.0,0.0}} },
--   time = {min = 0, max = -1, interval = 1}, 
--   format = 'ascii',
--   folder = './', 
-- }
--

