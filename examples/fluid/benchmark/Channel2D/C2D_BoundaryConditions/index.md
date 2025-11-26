title: Channel 2D Boundary Conditions
@warning WORK IN PROGRESS @endwarning

## Channel 2D Boundary Conditions ## {#ex_channel2D_BC}

The testcases of channel 2D boundary present some
combination of boundary conditions using the Poisueille testcase:

* Boundary conditions to define pressure, velocity or mass flow rate at the 
boundaries, using the analytical solution as input.
* Some of the cases use symmetry conditions.
* Q-Values at boundaries to improve the accuracy of boundary conditions.

Here is the list of test cases:

* [inlet: mass flow rate bounce back, outlet: pressure extrapolation](C2D_BC_MfrBB_PressExpol/index.md)
* [inlet: mass flow rate equilibrium, outlet: pressure equilibrium](C2D_BC_MfrEq_PressEq/index.md)
* [inlet: pressure extrapolation, outlet: pressure extrapolation](C2D_BC_PressExpol_PressExpol/index.md)
* [inlet: velocity bounce back, outlet: pressure extrapolation](C2D_BC_VelBB_PressExpol/index.md)
* [inlet: velocity linear extrapolation , outlet: pressure extrapolation](C2D_BC_VelBFL_PressExpol/index.md)
* [inlet: velocity equilibrium, outlet: pressure equilibrium](C2D_BC_VelEq_PressEq/index.md)
* [inlet: velocity non-equilibrium etrapolation, outlet: pressure non-equilibrium etrapolation](C2D_BC_VelNonEqExpol_PressNonEqExpol/index.md)
