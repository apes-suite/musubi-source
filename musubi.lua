-- This is a sample Musubi configuration file in lua format
-- Its main purpose is to show the options you might want to set for your
-- simulation. It thus serves as a little documentation for the configuration.
--
-- WARNING!!!
-- LEAVE THIS FILE UNCHANGED IN THE REPOSITORY
-- unless you change configurable options, in this case you should update this
-- script accordingly.
--
-- This script should always serve as an example which runs "out of the box".
--
-- Thank you!
--------------------------------------------------------------------------------

-- With the simulation name you can put a label on your simulation run, to have
-- some concise reference. The name is arbitrary and will be provided in the
-- timing.res output.
simulation_name = 'Gausspulse'

-- Global parameters
nu_phy = 3e-2
rho0   = 1.
cs2    = 1./3.
p0     = 1

-- Include general treelm settings
require('treelm/config')

-- Some parameters for the functions given below
originX =  -1.3
originY =  0.8
originZ =  0.1
halfwidth = 0.50
amplitude = 0.01

-- Some spatial lua functions, they take three coordinates and return one value.
function ic_1Dgauss_pulse(x, y, z)
  return p0+amplitude*math.exp(-0.5/(halfwidth^2)*( x - originX )^2)
end
function ic_2Dgauss_pulse(x, y, z)
  return p0 +amplitude*math.exp(-0.5/(halfwidth^2)*(( x - originX )^2+( y - originY )^2))
end
function ic_1Dgauss_pulse2(x, y, z,t)
  return p0+amplitude*math.exp(-0.5/(halfwidth^2)*( x - originX )^2)
end

-- Define the default interpolation method to use between different refinement
-- levels, if there are insufficient elements for this kind it will fall back to
-- a lower order interpolation methode:
-- quadratic > linear > weighted_average
interpolation_method = 'quadratic' -- default: quadratic

-- Time step settings
tmax           =  20     -- some helping variable to use below

-- simulation control settings
sim_control = {
-- Actual time definitions
-- times can be given in terms of iterations or simulation time.
  time_control = {
    max = { iter = tmax }, -- maximum number of iterations to reach in the finest level
    interval = { iter = tmax } -- interval to perform checks total density
  },
  abort_criteria = {
    stop_file = '', -- provide some filename, if exist simultion will stop
    steady_state = false -- Set true to stop simulation if simulation reach steady state
                         -- use tracking format='convergence' to check for steady state
                         -- based on single variable on domain with specified
                         -- convergence conditionsis done only when convergence
  }
}

-- file to write measurement time of simulation for each solver step
timing_file = 'mus_timing.res'

io_buffer_size = 8    -- in terms of MB, default value is 8

printRuntimeInfo = true   -- print runtime info after simulation

control_routine = 'multischeme'  -- another option is 'benchmark'

scaling = 'acoustic' -- or diffusive

comm_reduced = true   -- communicate reduced set of links?

init_allElems = false -- besides fluids, initialize ghost and halo elements?
                      -- only use this feature for debugging!

init_allTimeLayers = true -- initialize all time layers

-- mpi communication type, other options are:
-- isend_irecv_overlap
-- typed_isend_irecv
-- gathered_type
commpattern = 'isend_irecv'


-- Debug options
debug = {
  debugMode = false,  -- activate global debug modus
  debugFiles = false, -- open the debugFiles for each process dbgOut0000**.out

  dumpTreeIDs    = false,
  dumpPropBits   = false,
  dumpAuxLists   = false, -- write connectivity
  dumpBoundaries = false, -- write boundary elements information to dbgUnit

   dumpDependencies = false, -- write Ghost and Source elements to dbgUnit
  debugDependencies = false, -- write interpolation matrix
  checkDependencies = false, -- write neighbor informaiton

  dumpState     = false,
  dumpHaloState = false,

  checkNans = false,
  checkSteps = false,

  debugMesh = false,
  debugSource = false,
  debugRestart = false,

  traceMemory = false,

  logging = {level=1, filename = 'dbg', root_only = false },
}
-- Physical reference values, used for LB to physical unit convertion
-- if this table is present, all input paramters to musubi like
-- fluid/species property, initial_condition, boundary condition (amplitude
-- and transient times )must defined in physical units.
-- physics = {
--             dt = dt_phy, -- dt of coarsest simulation time
--             -- reference kg can be defined by mass0/molWeight0/rho0
--             rho0 = rho0_phy, -- reference density
--             -- reference mole can be defined by mole0/moleDens0
--             -- if neither defined default is set to inverse of Avogadro's constant
--             moleDens0 = moleDens0_phy, -- reference mole density
--             temp0 = t_phy, -- reference temperature
-- }


-- scheme model for single fluid simulation
identify = {
    kind       = 'fluid',     -- simulation type of this scheme
                            -- ( fluid, fluid_incomp, passive_scalar, ...)
    relaxation = 'bgk',     -- relaxation type (bgk, mrt, ...)
    -- Scheme layout
    -- This describes the stencil to use in the simulation.
   layout = 'd3q19'
}

  -- field which defines fluid or species
  -- Single fluid simulation
  --  field = {
  --     label = 'fluid',

-- Define the Lattice Boltzmann fluid properties
fluid = { kinematic_viscosity = nu_phy,
          bulk_viscosity      = 2/3*nu_phy }

-- For single field, source can be added via both key word:
-- "glob_source" or "source.
-- For multi field, field "source" must be defined inside field table
-- and "glob_source" must be defined in scheme table.
-- Source variable is defined as space-time function
-- if shape is not provided then source is applied on all fluid elements
-- For detail example on how to define spacetime function refer to
-- treelm/utest/tem_spacetime_fun_test.f90
-- Possible source variables are:
-- For fluid and fluid_incomp: force (unit: Newton)
-- For passive_scalar: injection, equal_injection
-- For multi-species_liquid: electric_field (unit: volt/meter),
--                           gravity_field (unit: m/s2)
source = {
  force = 'zeroforce'
}

-- Define nonNewtonian model
-- add nonNewtonian table to fluid table when nonNewtonian feature is needed
--  fluid = {
--    kinematic_viscosity = nu_phy,
--    nonNewtonian = {
--      model = 'power_law', or 'Casson', 'Carreau-Yasuda'
--      n = n, -- parameter for power law
--      k = viscPhy, -- parameter for power law
--    },
--  }
-- Also, identify table needs to set
--  identify = {
--      kind = 'fluid', or 'fluid_incompressible'
--      relaxation = 'bgk_pl', or 'bgk_cy', 'bgk_cs', or 'bgk_pl_explicit'
--      layout = 'd3q19'},

-- Initial condition for each field
-- Define initial conditions per variable.
-- They might be constants, predefined functions (described with a table), or
-- lua functions.
initial_condition = {
--     density =  {
--       predefined='gausspulse',
--       center={5.0,5.0, 5.0},
--       halfwidth=1.,
--       amplitude=1.20,
--       background=1.000},
-- initial with lua function
     pressure  = ic_1Dgauss_pulse, -- see above for its definition / computation.
     velocityX = 0.0,  -- constant IC
     velocityY = 0.0,  -- constant IC
     velocityZ = 0.0,  -- constant IC
     Sxx = 0.0, Syy = 0.0, Szz = 0.0, Sxy = 0.0, Syz = 0.0, Sxz = 0.0,
} -- constant IC

-- boundary condition for each field
-- For this example (Gaussian pulse), there is no boundary condition,
-- thus no need to define this table.
--    boundary_condition = {}
-- For other case, boundary conditions can be defined as:
-- boundary_condition = {
--  { label = 'inlet',  kind = 'inlet_ubb',  velocity = 'inlet_vel', }
--  { label = 'outlet', kind = 'outlet_pab', pressure = 'outlet_p',  },
--  { label = 'top',    kind = 'wall', },
--  { label = 'bottom', kind = 'wall', },
-- }
-- In musubi, there are following BCs provided:
-- wall, wall_libb, slip_wall
-- inlet_ubb,
-- outlet_expol, outlet_pab, outlet_eq, outlet_zero_prsgrd

variable = {
  {
    name = 'inlet_vel',
    ncomponents = 3,
    vartype = 'st_fun',
    st_fun = {
      predefined = 'combined',
      temporal  = {
        predefined = 'smooth',
        min_factor = 0, max_factor = u_in_L,
        from_time = 0, to_time = tmax/4,
      },
      spatial = { const = {1.0,0.0,0.0} },
    },
  },
  {
    name = 'outlet_p',
    ncomponents = 1,
    vartype = 'st_fun',
    st_fun = 0,
  },
  { name = 'spacetime',
    ncomponents = 1,
    vartype = 'st_fun',
    st_fun = ic_1Dgauss_pulse2
  },
  {
    name = 'difference',
    ncomponents =1,
    vartype = 'operation',
    operation = {kind='difference',
                input_varname = {'density','spacetime'}}
  },
  { name = 'zeroforce',
    ncomponents = 3,
    vartype = 'st_fun',
    st_fun = {0.0,0.0,0.0}
  },
}

-- See treelm/config.lua for explanations on tracking settings.
tracking = {
   label = 'point',
   variable = { 'density', 'velocity'
    ,'spacetime'
--    ,'difference',
-- possible names are:
--    pdf, density, pressure, velocity, velMag,
--    temperature, molfrac, kinetic_energy
--    equilFromState, equilFromMacro, nonEquilFromState,
--    shearstress, wss, shearMag
--    strainRate, shearRate,
--    strainRatenNwtn, when nonNewtonian is on, use this instead of strainRate
--    spacetime, difference
--    {'velocity',sub_index={1}}
              }, -- variable table
   folder = 'tracking_',
   shape = {kind = 'canoND', object = {origin = {3.0,3.1,3.0} } },
--   shape = { kind = 'canoND', object = {origin = {0.0, 0.0, 0.0},
--                                        vec = {0.0, 10.0, 0.0},
--                                        segments = 100,
--                                        distribution = 'equal'}},
--    track boundary elements
--    shape = {kind = 'property', property = {'boundary'}},
   output = { format = 'ascii' }, -- 'asciiSpatial', 'harvester', 'convergence'
   time_control = { min = {iter = 10},
                    max = {iter=tmax},
                    interval = {iter = 1}
   }
-- convergence = {
--   norm = 'average', 'simple'
--   nvals = 100,
--   condition = { threshold = 2.0e-10, operator = '<='}
--   }
}   -- end of tracking table
--}


-- scheme model for multispecies simulation
--  identifier = {
--      label = 'default',
--      kind = 'multispecies_liquid',
--      relaxation = 'BGK',
---- scheme layout
--      layout = 'd3q19'},
--      depend = { vay_sys = {variable={'velocity'}}, usescheme='fluid'}
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
--
