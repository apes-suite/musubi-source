! See copyright notice in the COPYRIGHT file.
! This utest program checks if the bgk d3Q19 optimized(Op) and explicit compute
! kernel output the same results. Input pdf is initilized by equilibrium with
! randomly generated density and velocity values.
program mus_fNeq_fromPDF_test
  use env_module, only: rk, eps
  use tem_param_module,      only: div1_4, div1_8, div1_36,&
    &                              cs2inv, t2cs2inv, t2cs4inv
  ! use tem_logging_module, only: logUnit
  use tem_general_module, only: tem_general_type, tem_start, tem_finalize

  use mus_scheme_layout_module, only: mus_scheme_layout_type, &
    &                                 mus_define_d3q19

  implicit none

  integer, parameter :: QQ       = 19   !> number of pdf directions

  integer, parameter :: q__W     = 1     !< west             x-
  integer, parameter :: q__S     = 2     !< south            y-
  integer, parameter :: q__B     = 3     !< bottom           z-
  integer, parameter :: q__E     = 4     !< east             x+
  integer, parameter :: q__N     = 5     !< north            y+
  integer, parameter :: q__T     = 6     !< top              z+
  integer, parameter :: q_BS     = 7     !<                  z-,y-
  integer, parameter :: q_TS     = 8     !<                  z+,y-
  integer, parameter :: q_BN     = 9     !<                  z-,y+
  integer, parameter :: q_TN     = 10    !<                  z+,y+
  integer, parameter :: q_BW     = 11    !<                  x-,z-
  integer, parameter :: q_BE     = 12    !<                  x+,z-
  integer, parameter :: q_TW     = 13    !<                  x-,z+
  integer, parameter :: q_TE     = 14    !<                  x+,z+
  integer, parameter :: q_SW     = 15    !<                  y-,x-
  integer, parameter :: q_NW     = 16    !<                  y+,x-
  integer, parameter :: q_SE     = 17    !<                  y-,x+
  integer, parameter :: q_NE     = 18    !<                  y+,x+
  ! integer, parameter :: q__0     = 19    !< rest density is last

  logical :: error
  real(kind=rk) :: tolerance
  integer :: iDir

  real(kind=rk) :: rho, u_x, u_y, u_z, usq, ucx !, omega
  real(kind=rk) :: f(QQ), fEq(QQ), fNeq(QQ)

  real(kind=rk) :: usqn, usqn_o2
  real(kind=rk) :: coeff_1, coeff_2
  real(kind=rk) :: ui1, fac_1, sum1_2
  real(kind=rk) :: fac_2, sum2_2, fac_4, sum4_2, fac_9, sum9_2
  real(kind=rk) :: ui3, fac_3, sum3_2
  real(kind=rk) :: ui10, fac_10, sum10_2
  real(kind=rk) :: ui11, fac_11, sum11_2
  real(kind=rk) :: ui12, fac_12, sum12_2
  real(kind=rk) :: ui13, fac_13, sum13_2
  real(kind=rk) :: fNeq_NE_SW
  real(kind=rk) :: fNeq_NW_SE
  real(kind=rk) :: fNeq_BW_TE
  real(kind=rk) :: fNeq_BE_TW
  real(kind=rk) :: fNeq_BS_TN
  real(kind=rk) :: fNeq_BN_TS
  real(kind=rk) :: fNeq_N_S
  real(kind=rk) :: fNeq_W_E
  real(kind=rk) :: fNeq_T_B

  type( mus_scheme_layout_type ) :: layout
  type( tem_general_type ) :: general

  call tem_start('fNeq from PDF calculation utest', 'utest', general)
  error = .false.
  tolerance = eps * 25._rk
  write(*,*) 'tolerance = ', tolerance

  call mus_define_d3q19( layout = layout, nElems = 1 )

  ! initialize PDF
  do iDir = 1, QQ
    f(iDir) = 0.01_rk * iDir
  end do

  ! calculate rho and u_x, u_y, u_z, usq
  rho = sum( f(:) )
  u_x = sum( f(:) * layout%fStencil%cxDir(1,:) )
  u_y = sum( f(:) * layout%fStencil%cxDir(2,:) )
  u_z = sum( f(:) * layout%fStencil%cxDir(3,:) )
  usq = u_x * u_x + u_y * u_y + u_z * u_z

  ! calculate fEq
  do iDir = 1, QQ

    ! velocity times lattice unit velocity
    ucx =   dble( layout%fStencil%cxDir( 1, iDir ))*u_x  &
      &   + dble( layout%fStencil%cxDir( 2, iDir ))*u_y  &
      &   + dble( layout%fStencil%cxDir( 3, iDir ))*u_z

    ! calculate equilibrium density
    fEq( iDir ) = layout%weight( iDir ) * rho * ( 1.d0 + ucx*cs2inv &
      &         + ucx*ucx*t2cs4inv - usq*t2cs2inv )

  enddo

  ! Calculate the non-equilibrium part
  fNeq(:) = f(:) - fEq(:)

  ! Calculate according to BGK optimized kernel
  ! omega = 1.0_rk
  usqn = div1_36 * (1._rk - 1.5_rk * usq) * rho

  coeff_1 = div1_8 * rho

  ui1     = u_x + u_y
  fac_1   = coeff_1 * ui1
  sum1_2  = fac_1 * ui1 + usqn
  fNeq_NE_SW = f(Q_NE) + f(Q_SW) - sum1_2*2.0_rk
  call check_error( fNeq_NE_SW, fNeq(Q_NE) + fNeq(Q_SW), tolerance, error )

  ui3     = -u_x + u_y
  fac_3   = coeff_1 * ui3
  sum3_2  = fac_3 * ui3 + usqn
  fNeq_NW_SE = f(Q_NW) + f(Q_SE) - sum3_2*2.0_rk
  call check_error( fNeq_NW_SE, fNeq(Q_NW) + fNeq(Q_SE), tolerance, error )

  ui10    =  u_x + u_z
  fac_10  = coeff_1 * ui10
  sum10_2 = fac_10 * ui10 + usqn
  fNeq_BW_TE = f(Q_BW) + f(Q_TE) - sum10_2*2.0_rk
  call check_error( fNeq_BW_TE, fNeq(Q_BW) + fNeq(Q_TE), tolerance, error )

  ui12    = -u_x + u_z
  fac_12  = coeff_1 * ui12
  sum12_2 = fac_12 * ui12 + usqn
  fNeq_BE_TW = f(Q_BE) + f(Q_TW) - sum12_2*2.0_rk
  call check_error( fNeq_BE_TW, fNeq(Q_BE) + fNeq(Q_TW), tolerance, error )

  ui11    =  u_y + u_z
  fac_11  = coeff_1 * ui11
  sum11_2 = fac_11 * ui11 + usqn
  fNeq_BS_TN = f(Q_BS) + f(Q_TN) - sum11_2*2.0_rk
  call check_error( fNeq_BS_TN, fNeq(Q_BS) + fNeq(Q_TN), tolerance, error )

  ui13    = -u_y + u_z
  fac_13  = coeff_1 * ui13
  sum13_2 = fac_13 * ui13 + usqn
  fNeq_BN_TS = f(Q_BN) + f(Q_TS) - sum13_2*2.0_rk
  call check_error( fNeq_BN_TS, fNeq(Q_BN) + fNeq(Q_TS), tolerance, error )

  coeff_2 = div1_4 * rho
  usqn_o2 = 2.0_rk * usqn

  fac_2   = coeff_2 * u_y
  sum2_2  = fac_2 * u_y + usqn_o2
  fNeq_N_S = f(Q__N) + f(Q__S) - sum2_2*2.0_rk
  call check_error( fNeq_N_S, fNeq(Q__N) + fNeq(Q__S), tolerance, error )

  fac_4   = coeff_2 * u_x
  sum4_2  = fac_4 * u_x + usqn_o2
  fNeq_W_E = f(Q__W) + f(Q__E) - sum4_2*2.0_rk
  call check_error( fNeq_W_E, fNeq(Q__W) + fNeq(Q__E), tolerance, error )

  fac_9   = coeff_2 * u_z
  sum9_2  = fac_9 * u_z + usqn_o2
  fNeq_T_B = f(Q__T) + f(Q__B) - sum9_2*2.0_rk
  call check_error( fNeq_T_B, fNeq(Q__T) + fNeq(Q__B), tolerance, error )

  if ( error ) then
    write(*,*) 'Error exceeds tolerance!'
  end if
  write(*,*) 'optimized,           normal'
  write(*,*) fNeq_NE_SW, fNeq(Q_NE) + fNeq(Q_SW)
  write(*,*) fNeq_NW_SE, fNeq(Q_NW) + fNeq(Q_SE)
  write(*,*) fNeq_BW_TE, fNeq(Q_BW) + fNeq(Q_TE)
  write(*,*) fNeq_BE_TW, fNeq(Q_BE) + fNeq(Q_TW)
  write(*,*) fNeq_BS_TN, fNeq(Q_BS) + fNeq(Q_TN)
  write(*,*) fNeq_BN_TS, fNeq(Q_BN) + fNeq(Q_TS)
  write(*,*) fNeq_N_S,   fNeq(Q__N) + fNeq(Q__S)
  write(*,*) fNeq_W_E,   fNeq(Q__W) + fNeq(Q__E)
  write(*,*) fNeq_T_B,   fNeq(Q__T) + fNeq(Q__B)

  call tem_finalize(general)

  if (.not. error) then
    write(*,*) 'PASSED'
  else
    write(*,*) 'FAILED'
  end if

  !*****************************************************************************

contains

  subroutine check_error( a, b, tolerance, error )
    real(kind=rk), intent(in) :: a, b, tolerance
    logical :: error
    if ( tolerance < abs( a-b ) ) then
      error = .true.
    end if
  end subroutine check_error

end program mus_fNeq_fromPDF_test
!******************************************************************************
