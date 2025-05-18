! See copyright notice in the COPYRIGHT file.
! This utest program tests the the non-equilibirium calcuation method by
! acoustic scaling (CE expansion)
program mus_fNeq_acoustic_test
  use env_module,         only: rk, eps
  use tem_general_module, only: tem_start, tem_finalize, tem_general_type

  use mus_scheme_layout_module, only: mus_scheme_layout_type, &
    &                                 mus_define_d3q19, &
    &                                 mus_define_d2q9
  use mus_derivedQuantities_module2, only: getNEq_acoustic, &
    &                                      secondMom_3D, convPrePost

  implicit none

  logical :: error
  real(kind=rk) :: tolerance, max_error

  type( tem_general_type ) :: general
  type( mus_scheme_layout_type ) :: layout
  real(kind=rk) :: stress(3,3)
  real(kind=rk) :: strain(3,3)
  real(kind=rk), allocatable :: fNeq(:)
  real(kind=rk) :: rho0  = 1.0_rk
  real(kind=rk) :: omega = 1.8_rk
  real(kind=rk) :: nu
  real(kind=rk) :: stress_ini(6), stress_res(6), diffs(6)
  integer :: nDims = 3

  integer :: iDir, iVal

  call tem_start( 'fNeq acoustic utest', general )
  error = .false.
  tolerance = eps * 2500._rk
  max_error = 0.0_rk
  write(*,"(A, F14.9)") 'tolerance = ', tolerance

  write(*,"(A)")    'Flow parameters:'
  write(*,"(A, F10.6)") ' rho0  = ', rho0
  write(*,"(A, F10.6)") ' omega = ', omega
  nu = (1.0_rk / omega - 0.5_rk ) / 3.0_rk
  write(*,"(A, F10.6)") ' nu    = ', nu

  ! initialize stress(3,3), it has to be symmetric
  stress = transpose(reshape([  0.01_rk,   0.04_rk,  0.06_rk, &
    &                           0.04_rk,   0.02_rk,  0.05_rk, &
    &                           0.06_rk,   0.05_rk,  0.03_rk ], shape(stress)))
  strain = stress / 2.0_rk / nu / rho0
  stress_ini = [  stress(1,1), stress(2,2), stress(3,3), &
    &             stress(1,2), stress(2,3), stress(1,3)]

  if ( nDims == 2 ) then
    call mus_define_d2q9( layout = layout, nElems = 1 )
  else
    call mus_define_d3q19( layout = layout, nElems = 1 )
  end if

  ! calculate fNeq
  allocate( fNeq(layout%fStencil%QQ) )
  fNeq(:) = getnEq_acoustic( layout = layout, &
    &                        omega  = omega,  &
    &                        Sxx    = strain  )

  ! calculate shear stress
  stress_res =   secondMom_3D( layout%fStencil%cxcx, fNeq, layout%fStencil%QQ ) &
    &          * ( omega * .5_rk - 1._rk ) * convPrePost( omega )

  ! print out fNeq values
    write(*, "(A4, A14)") "iDir", "fNeq"
  do iDir = 1, layout%fStencil%QQ
    write(*, "(I4, F14.8)") iDir, fNeq(iDir)
  end do
  write(*,*) ''

  ! print out stress and stress_res

  diffs = abs(stress_ini - stress_res)
  if ( maxval( diffs ) > tolerance ) then
    error = .true.
    write(*, "(A)")    "Max Error exceeds tolerance!!!"
    write(*, "(3A14)") "stress_ini", "stress_res", "difference"
    do iVal = 1, 6
      write(*, "(3F14.8)") stress_ini(iVal), stress_res(iVal), diffs(iVal)
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

end program mus_fNeq_acoustic_test
!******************************************************************************
