! Copyright (c) 2023 Gregorio Gerardo Spinelli <gregoriogerardo.spinelli@dlr.de>
! Copyright (c) 2023 Harald Klimach <harald.klimach@dlr.de>
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are met:
!
! 1. Redistributions of source code must retain the above copyright notice,
! this list of conditions and the following disclaimer.
!
! 2. Redistributions in binary form must reproduce the above copyright notice,
! this list of conditions and the following disclaimer in the documentation
! and/or other materials provided with the distribution.
!
! THIS SOFTWARE IS PROVIDED BY THE UNIVERSITY OF SIEGEN “AS IS” AND ANY EXPRESS
! OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
! OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
! IN NO EVENT SHALL UNIVERSITY OF SIEGEN OR CONTRIBUTORS BE LIABLE FOR ANY
! DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
! (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
! LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
! ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
! (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
! SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
! **************************************************************************** !
program mus_wall_function_schmitt_test
  use env_module,         only: rk
  use tem_general_module, only: tem_start, tem_finalize, tem_general_type

  use mus_wall_function_abstract_module, only: mus_wall_function_type
  use mus_wall_function_schmitt_module,  only: mus_wall_function_schmitt_type, &
    &                                          get_uTau_logLayer,              &
    &                                          get_uTau_subVisousLayer

  implicit none

  !> tolerance value for the test
  real(kind=rk), parameter :: eps = 1e-10_rk

  integer, parameter :: uPlus_sub_test = 1
  integer, parameter :: df_du_sub_test = 2
  integer, parameter :: uPlus_log_test = 3
  integer, parameter :: df_du_log_test = 4
  integer, parameter :: uPlus_buf_test = 5
  integer, parameter :: df_du_buf_test = 6

  logical :: error(6)

  type( tem_general_type ) :: general

  class(mus_wall_function_type), allocatable :: wall_function
  !> yPlus value for the test
  real(kind=rk) :: yPlus
  !> uPlus value for the test
  real(kind=rk) :: uPlus
  !> friction velocity, dynamic viscosity and wall distance
  real(kind=rk) :: uTau, nu, y, visc_div_dist, result, velSW

  call tem_start( 'Schmitt wall profile test', general )

  ! Assume all tests to fail
  error = .true.

  allocate(mus_wall_function_schmitt_type :: wall_function)

  ! sub-layer
  yPlus = 2._rk
  ! Reference obtained with python sympy
  uPlus = 2._rk

  ! test uPlus
  error(uPlus_sub_test) = (abs(wall_function%get_uPlus(yPlus) - uPlus) > eps)

  ! buffer-layer
  yPlus = 25._rk
  ! Reference obtained with python sympy
  uPlus = 8.49668862456533_rk

  ! test uPlus
  error(uPlus_buf_test) = (abs(wall_function%get_uPlus(yPlus) - uPlus) > eps)

  ! log-layer
  yPlus = 100._rk
  ! Reference obtained with python sympy
  uPlus = 16.0247911497310_rk

  ! test uPlus
  error(uPlus_log_test) = (abs(wall_function%get_uPlus(yPlus) - uPlus) > eps)

  ! test derivative uPlus w.r.t. uTau
  ! sub-layer
  velSW = 25._rk
  nu = 1e-5_rk
  y = 1._rk
  visc_div_dist = nu / y
  uTau = 0.0158113883008419_rk

  ! copied from line 369 of mus_turb_wallFunc_module.f90
  result = get_uTau_subVisousLayer( visc_div_dist = visc_div_dist, &
    &                               velSW = velSW                  )

  error(df_du_sub_test) = (abs(result - uTau) > eps)

  ! buffer-layer
  uPlus = 0.0869045919087855_rk
  uTau = 25._rk

  error(df_du_buf_test) = ( abs(wall_function%get_d_uPlus_d_uTau(y, uTau, nu) - uPlus) > eps )

  ! log-layer
  uTau = 0.622306270073191_rk
  velSW = 25._rk

  ! copied from line 366-367 of mus_turb_wallFunc_module.f90
  result = get_uTau_logLayer( visc_div_dist = visc_div_dist, &
    &                         velSW = velSW                  )

  error(df_du_log_test) = (abs(result - uTau) > eps)

  call tem_finalize(general)

  if (.not. any(error)) then
    write(*,*) 'PASSED'
  else
    write(*,*) 'FAILED'
  end if

end program mus_wall_function_schmitt_test
! **************************************************************************** !
