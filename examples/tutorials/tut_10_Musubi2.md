title: Prerequisites
@warning WORK IN PROGRESS @endwarning

Navigate: [Overview](index.html)

# Changes in Musubi 2.0

Musubi 2.0 comes with several changes. Old config files will not work any more 
since some layout and names have been changes.

Here is a list of important changes:
* Added lattice Mach number and lattice velocity check
* Changed layout names
* Renamed boundary routine names
* fluid table: we don't specify omega anymore, now we need kinematic_viscosity 
and bulk_viscosity

# Bc Names

Boundary routine names changed:

> **old name**                --> **new name**
> inlet_bfl                   --> velocity_bfl
> inlet_bfl_incomp            --> velocity_bfl_incomp
> inlet_eq                    --> velocity_eq
> inlet_mfr                   --> mfr_bounceback
> inlet_mfr_eq                --> mfr_eq
> inlet_ubb                   --> velocity_bounceback
> inlet_ubb_incomp            --> velocity_bounceback_incomp
> 
> moleDens_dirichlet_curved   --> moleDens_nonEqExpol_curved
> moleDens_dirichlet_straight --> moleDens_nonEqExpol
> moleDens_neumann_straight   --> moleDens_neumann
> 
> moments_press               --> pressure_momentsbased
> moments_press_incomp        --> pressure_momentsbased_incomp
> moments_vel                 --> velocity_momentsbased
> moments_vel_incomp          --> velocity_momentsbased_incomp
> 
> outlet_eq                   --> pressure_eq
> outlet_expol                --> pressure_expol
> outlet_expol_slow           --> pressure_expol_slow
> outlet_pab                  --> pressure_antiBounceBack
> 
> potential_dirichlet         --> potential_nonEqExpol
> potential_dirichlet_curved  --> potential_nonEqExpol_curved

# Layout Names

> lbm  --> fluid
> lbm_incomp --> fluid_incompressible

# Lattice velocity check

When running musubi we check if Lattice velocity exceed the default stability
threshold of 0.15. If this happens we abort the simulation and end the code in
a propper way. This theshold can be modified inside the config file via:

```lua
latv_max=0.2
```

