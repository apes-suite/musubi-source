require 'seeder'
require 'musubi'

print('-----------------------------------------------------------------------')
print('Simulation name: ', simulation_name)
print('-----------------------------------------------------------------------')
if case2d then
  print('You are using a 2D case,')
else
  print('You are using a 3D case,')
end
if usePeriodic then
  print('with periodicity in z-direction')
end
if useObstacle then
  print('and the following obstacle:', stlLabel)
--  print(stlLabel)
end
if qValues then
  print('with qValues.')
end
if useRefine then
  print('You are also using refinement (multi-level mesh).')
else
  print('You are not using refinement (single-level mesh).')
end
print('-----------------------------------------------------------------------')
print('---Mesh parameters---')
print('length            =', length)
print('height            =', height)
print('width             =', depth)
print('nLength           =', nLength)
print('nHeight           =', nHeight)
print('length_bnd        =', length_bnd)
print('level             =', level)
print('level of refine 1 =', level+refinementLevel)
print('level of refine 2 =', level+refinementLevel2)
print('-----------------------------------------------------------------------')
print('---Flow parameters---')
print('Re                =', Re)
print('Mach              =', Ma)
print('---In physical units---')
print('Velocity          =', vel_phy, '[m/s]')
print('Kinematic visc.   =', nu_phy, '[m^s/2]')
print('Density           =', rho0_phy, '[kg/m^3]')
print('Element size (dx) =', dx, '[m]')
print('Timestep (dt)     =', dt, '[s]')
print('Press. drop       =', press_drop, '[N/m^2]')
print('Press. west       =', pressureIn, '[N/m^2]')
print('Press. east       =', pressureOut, '[N/m^2]')
print('--- In lattice units---')
print('You are using', scaling, 'scaling')
print('Lattice vel.         =', vel_lat)
print('Latt. viscosity      =', nu_lat)
print('Latt. speed of sound =', cs_lat)
print('Relaxation param.    =', omega)
print('-----------------------------------------------------------------------')
