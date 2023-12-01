! See copyright notice in the COPYRIGHT file.
! This utest program tests the the non-equilibirium calcuation method by
! diffusive scaling (Asymptotic analysis)
program mus_fNeq_diffusive_test
  use env_module,         only: rk, eps
  use tem_general_module, only: tem_start, tem_finalize, tem_general_type

  use mus_scheme_layout_module, only: mus_scheme_layout_type, &
    &                                 mus_define_d3q19, &
    &                                 mus_define_d2q9
  use mus_derivedQuantities_module2, only: getNEq_diffusive, &
    &                                      secondMom_3D, convPrePost
  use mus_physics_module, only: mus_physics_type, &
    &                           set_values_by_levels, mus_set_convFac

  implicit none

  logical :: error
  real(kind=rk) :: tolerance, max_error

  type( tem_general_type ) :: general
  type( mus_scheme_layout_type ) :: layout
  real(kind=rk) :: stress_phy(3,3)
  real(kind=rk) :: strain_phy(3,3)
  real(kind=rk) :: strain_minLevel_lb(3,3)
  real(kind=rk) :: strain_maxLevel_lb(3,3)
  real(kind=rk), allocatable :: fNeq_minLevel(:)
  real(kind=rk), allocatable :: fNeq_maxLevel(:)
  real(kind=rk) :: rho0  = 1.0_rk
  real(kind=rk) :: omega = 1.3333_rk
  real(kind=rk) :: nu, nu_phy
  real(kind=rk) :: stress_ini(6), diffs(6)
  real(kind=rk) :: stress_res_minLevel_phy(6)
  real(kind=rk) :: stress_res_maxLevel_phy(6)
  integer :: nDims = 3
  type( mus_physics_type ) :: physics

  integer :: iDir, iVal, minLevel, maxLevel

  call tem_start( 'fNeq diffusive utest', 'utest', general )
  error = .false.
  tolerance = eps * 2500._rk
  max_error = 0.0_rk
  write(*,"(A, F14.9)") 'tolerance = ', tolerance

  write(*,"(A)")    'Flow parameters:'
  write(*,"(A, F10.6)") ' rho0  = ', rho0
  write(*,"(A, F10.6)") ' omega = ', omega
  nu = (1.0_rk / omega - 0.5_rk ) / 3.0_rk
  write(*,"(A, F10.6)") ' nu    = ', nu

  minLevel = 2
  maxLevel = 3
  ! Set up physisc convertion factor
  physics%dx   = 1.0_rk
  physics%dt   = 1.0_rk
  physics%rho0 = rho0
  nu_phy = 1.0_rk
  ! Assign dx and dt for each level
  allocate(physics%dxLvl( minLevel:maxLevel ))
  allocate(physics%dtLvl( minLevel:maxLevel ))
  allocate(physics%fac(   minLevel:maxLevel ))
  ! dx and dt is set according scaling type
  physics%dxLvl( minLevel:maxLevel ) =  &
    & set_values_by_levels( physics%dx, minLevel, maxLevel, 2 )
  physics%dtLvl( minLevel:maxLevel ) =  &
    & set_values_by_levels( physics%dt, minLevel, maxLevel, 4 )
  ! compute and store conversion factors in converstionFac type
  call mus_set_convFac( me = physics, minLevel = minLevel, maxLevel = maxLevel )

  ! Physical stress
  ! initialize stress(3,3), it has to be symmetric
  stress_phy = transpose(reshape([  0.01_rk,   0.04_rk,  0.06_rk, &
    &                           0.04_rk,   0.02_rk,  0.05_rk, &
    &                           0.06_rk,   0.05_rk,  0.03_rk ], shape(stress_phy)))
  strain_phy = stress_phy / 2.0_rk / nu_phy / rho0
  stress_ini = [  stress_phy(1,1), stress_phy(2,2), stress_phy(3,3), &
    &             stress_phy(1,2), stress_phy(2,3), stress_phy(1,3)]
  ! Physical to LB unit
  strain_minLevel_lb = strain_phy / physics%fac(minLevel)%strainRate
  strain_maxLevel_lb = strain_phy / physics%fac(maxLevel)%strainRate


  if ( nDims == 2 ) then
    call mus_define_d2q9( layout = layout, nElems = 1 )
  else
    call mus_define_d3q19( layout = layout, nElems = 1 )
  end if

  ! calculate fNeq at minLevel and maxLevel
  allocate( fNeq_minLevel(layout%fStencil%QQ) )
  allocate( fNeq_maxLevel(layout%fStencil%QQ) )
  fNeq_minLevel(:) = getnEq_diffusive( layout = layout, &
    &                         omega  = omega,  &
    &                         Sxx    = strain_minLevel_lb  )

  fNeq_maxLevel(:) = getnEq_diffusive( layout = layout, &
    &                         omega  = omega,  &
    &                         Sxx    = strain_maxLevel_lb  )

  ! calculate shear stress
  stress_res_minLevel_phy =   &
    &    secondMom_3D( layout%fStencil%cxcx, fNeq_minLevel, layout%fStencil%QQ ) &
    &  * ( omega * .5_rk - 1._rk ) * convPrePost( omega ) &
    &  * physics%fac(minLevel)%press

  stress_res_maxLevel_phy =   &
    &    secondMom_3D( layout%fStencil%cxcx, fNeq_maxLevel, layout%fStencil%QQ ) &
    &  * ( omega * .5_rk - 1._rk ) * convPrePost( omega ) &
    &  * physics%fac(maxLevel)%press


  ! print out fNeq values
    write(*, "(A4, 2A14)") "iDir", "fNeq_minL", "fNeq_maxL"
  do iDir = 1, layout%fStencil%QQ
    write(*, "(I4, 2F14.8)") iDir, fNeq_minLevel(iDir), fNeq_maxLevel(iDir)
  end do
  write(*,*) ''

  ! check stress_ini and stress_res

  diffs = abs(stress_res_minLevel_phy - stress_res_maxLevel_phy)
  if ( maxval( diffs ) > tolerance ) then
    error = .true.
    write(*, "(A)")    "Max Error exceeds tolerance!!!"
    write(*, "(3A14)") "stress_minL", "stress_maxL", "difference"
    do iVal = 1, 6
      write(*, "(3F14.8)") stress_res_minLevel_phy(iVal), stress_res_maxLevel_phy(iVal), diffs(iVal)
    end do
  else
    write(*, "(A)") "Max Error within tolerance!"
    error = .false.
  end if

  call tem_finalize(general)
  if (.not. error) then
    write(*,'(A)') 'PASSED'
  else
    write(*,'(A)') 'FAILED'
  end if

end program mus_fNeq_diffusive_test
!******************************************************************************
