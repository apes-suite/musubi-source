title: Convected vortex in a domain with a coarser grid at the center

# Vortex convection # {#eg_VC}

In this example, we will simulate a Vortex convected in a uniform flow.

The objectives of this example is to introduce how to:
* How to use the new interpolation based on the transfer of both
  equilibrium and non equilibrium pdfs at multi level interfaces.
  
The choice is available in the musubi.lua file.
interpolation_method = {
  method = 'quadratic', -- or any other: linear, weighted_average
  type_of_interpolation = 'moment'  -- 'pdf' 
}

type_of_interpolation = pdf => feq_neq (new version)
type_of_interpolation = moment => moment (old version)

* how to set the new collision models: 
    Recursive Regularized bgk (rr_bgk), 
    Projected RR_bgk (prr_bgk), 
    Hybrid RR_bgk (hrr_bgk), 
    Dual Relaxation Time rr_bgk (drt_bgk). 
  For HRR_bgk we need to give in input the value of sigma in the fluid 
  table, see below. The default value is 0.98. For DRT_bgk we need to 
  give in input the tauN value. The dafault value is 0.70.

--! [Scheme identifier]
identify = {
  layout = 'd2q9',    -- Stencil
  relaxation = 'prr_bgk', -- Collision: rr_bgk, prr_bgk, hrr_bgk
  kind = 'fluid'      -- Physics
}
--! [Scheme identifier]

--! [Fluid]
fluid = {
  kinematic_viscosity = nu_phy,
  hrr_sigma = 0.00, -- default: 0.98
  drt_taun = 0.7, -- default: 0.70
}
--! [Fluid]


## ToDo ##
This file has to be filled with content like:
* Problem description
* Formulas
* Results
* Comparison for different resoultions
* How to run the simulation
