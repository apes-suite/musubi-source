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
program mus_wall_function_musker_test
  use env_module,         only: rk
  use tem_general_module, only: tem_start, tem_finalize, tem_general_type

  use mus_wall_function_abstract_module, only: mus_wall_function_type
  use mus_wall_function_musker_module,   only: mus_wall_function_musker_type

  implicit none

  !> tolerance value for the test
  real(kind=rk), parameter :: eps = 1e-10_rk

  integer, parameter :: uPlus_test = 1
  integer, parameter :: df_du_test = 2

  logical :: error(2)

  type( tem_general_type ) :: general

  class(mus_wall_function_type), allocatable :: wall_function
  !> yPlus value for the test
  real(kind=rk) :: yPlus
  !> uPlus value for the test
  real(kind=rk) :: uPlus
  !> friction velocity, dynamic viscosity and wall distance
  real(kind=rk) :: uTau, nu, y

  call tem_start( 'Musker wall profile test', general )

  ! Assume all tests to fail
  error = .true.

  allocate(mus_wall_function_musker_type :: wall_function)
  yPlus = 25._rk
  ! Reference obtained with python sympy
  uPlus = 12.4368001725141_rk

  ! test uPlus
  error(uPlus_test) = (abs(wall_function%get_uPlus(yPlus) - uPlus) > eps)

  ! test derivative uPlus w.r.t. uTau
  uTau = 25._rk
  nu = 1.e-5
  y = 1._rk
  uPlus = 0.0972818682337712_rk

  error(df_du_test) = ( abs(wall_function%get_d_uPlus_d_uTau(y, uTau, nu) - uPlus) > eps )

  call tem_finalize(general)

  if (.not. any(error)) then
    write(*,*) 'PASSED'
  else
    write(*,*) 'FAILED'
  end if

end program mus_wall_function_musker_test
! **************************************************************************** !
