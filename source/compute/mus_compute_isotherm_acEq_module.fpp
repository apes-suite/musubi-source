! Copyright (c) 2016 Philipp Otte <otte@mathcces.rwth-aachen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2016, 2019-2020 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2016-2017 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2017 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2018 Raphael Haupt <raphael.haupt@uni-siegen.de>
! Copyright (c) 2020 Peter Vitt <peter.vitt2@uni-siegen.de>
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
! ****************************************************************************** !
!> author: Philipp Otte
!! Routines and parameter definitions for the isothermal acoustic Eq D3Q19 model
?? include 'header/lbm_macros.inc'
module mus_isotherm_acEq_module
  use iso_c_binding, only: c_f_pointer

  ! include treelm modules
  use env_module,            only: rk
  use tem_param_module,      only: rho0, rho0Inv
  use tem_varSys_module,     only: tem_varSys_type, tem_varSys_op_type
  use tem_dyn_array_module,  only: PositionOfVal
  use tem_aux_module,        only: tem_abort

  ! include musubi modules
  use mus_field_prop_module,    only: mus_field_prop_type
  use mus_scheme_layout_module, only: mus_scheme_layout_type
  use mus_param_module,         only: mus_param_type
  use mus_varSys_module,        only: mus_varSys_data_type
  use mus_derVarPos_module,     only: mus_derVarPos_type

  ! Just for debugging
  ! use tem_logging_module, only: logUnit

  implicit none

  private

  public :: bgk_advRel_isotherm_acEq_d3q19

  ! =============================================================================
  ! D3Q19 flow model
  ! =============================================================================
  !> Definition of the discrete velocity set

  ! integer,parameter :: block = 32
  integer,parameter :: QQ   = 19  !< number of pdf directions

  integer,parameter :: qN00 = 1   !< west             x-
  integer,parameter :: q0N0 = 2   !< south            y-
  integer,parameter :: q00N = 3   !< bottom           z-
  integer,parameter :: q100 = 4   !< east             x+
  integer,parameter :: q010 = 5   !< north            y+
  integer,parameter :: q001 = 6   !< top              z+
  integer,parameter :: q0NN = 7   !<                  z-,y-
  integer,parameter :: q0N1 = 8   !<                  z+,y-
  integer,parameter :: q01N = 9   !<                  z-,y+
  integer,parameter :: q011 = 10  !<                  z+,y+
  integer,parameter :: qN0N = 11  !<                  x-,z-
  integer,parameter :: q10N = 12  !<                  x+,z-
  integer,parameter :: qN01 = 13  !<                  x-,z+
  integer,parameter :: q101 = 14  !<                  x+,z+
  integer,parameter :: qNN0 = 15  !<                  y-,x-
  integer,parameter :: qN10 = 16  !<                  y+,x-
  integer,parameter :: q1N0 = 17  !<                  y-,x+
  integer,parameter :: q110 = 18  !<                  y+,x+
  integer,parameter :: q000 = 19  !< rest density is last

  real(kind=rk), parameter :: f1 = 2.0_rk / 5.0_rk
  real(kind=rk), parameter :: f2 = 1.0_rk / 30.0_rk
  real(kind=rk), parameter :: f8 = 1.0_rk / 30.0_rk


contains

! ****************************************************************************** !
  !> Advection relaxation routine for the D3Q19 model with BGK for the
  !> isothermal acoustic equation.
  !!
  !! This subroutine interface must match the abstract interface definition
  !! [[kernel]] in scheme/[[mus_scheme_type_module]].f90 in order to be callable
  !! via [[mus_scheme_type:compute]] function pointer.
  subroutine bgk_advRel_isotherm_acEq_d3q19( fieldProp, inState, outState,    &
    &                                        auxField, neigh, nElems, nSolve, &
    &                                        level, layout, params, varSys,   &
    &                                        derVarPos                        )
    ! -------------------------------------------------------------------- !
    !> Array of field properties (fluid or species)
    type(mus_field_prop_type), intent(in) :: fieldProp(:)
    !> variable system definition
    type(tem_varSys_type), intent(in) :: varSys
    !> current layout
    type(mus_scheme_layout_type), intent(in) :: layout
    !> number of elements in state Array
    integer, intent(in) :: nElems
    !> input  pdf vector
    real(kind=rk), intent(in)  ::  inState(nElems * varSys%nScalars)
    !> output pdf vector
    real(kind=rk), intent(out) :: outState(nElems * varSys%nScalars)
    !> Auxiliary field computed from pre-collision state
    !! Is updated with correct velocity field for multicomponent models
    real(kind=rk), intent(inout) :: auxField(nElems * varSys%nAuxScalars)
    !> connectivity vector
    integer, intent(in) :: neigh(nElems * layout%fStencil%QQ)
    !> number of elements solved in kernel
    integer, intent(in) :: nSolve
    !> current level
    integer,intent(in) :: level
    !> global parameters
    type(mus_param_type),intent(in) :: params
    !> position of derived quantities in varsys for all fields
    type( mus_derVarPos_type ), intent(in) :: derVarPos(:)
    ! -------------------------------------------------------------------- !
?? IF (SOA) THEN
?? ELSE
    ! ---------------------------------------------------------------------------
    integer :: iElem, nScalars
    real(kind=rk) :: fN00, f0N0, f00N, f100, f010, f001, f0NN, f0N1, f01N, &
      &              f011, fN0N, f10N, fN01, f101, fNN0, fN10, f1N0, f110, &
      &              f000
    real(kind=rk) :: rho     ! local density
    real(kind=rk) :: u_x     ! local x-velocity
    real(kind=rk) :: u_y     ! local y-velocity
    real(kind=rk) :: u_z     ! local z-velocity
    real(kind=rk) :: omega, cmpl_o
    integer :: dens_pos, vel_pos(3), elemOff
    ! ---------------------------------------------------------------------------
    dens_pos = varSys%method%val(derVarPos(1)%density)%auxField_varPos(1)
    vel_pos = varSys%method%val(derVarPos(1)%velocity)%auxField_varPos(1:3)

    nScalars = varSys%nScalars

!$omp do schedule(static)

?? copy :: dir_novec
    nodeloop: do iElem = 1, nSolve
      ! First load all local values into temp array
      fN00 = inState(?FETCH?( qN00, 1, iElem, QQ, nScalars, nElems,neigh ))
      f0N0 = inState(?FETCH?( q0N0, 1, iElem, QQ, nScalars, nElems,neigh ))
      f00N = inState(?FETCH?( q00N, 1, iElem, QQ, nScalars, nElems,neigh ))
      f100 = inState(?FETCH?( q100, 1, iElem, QQ, nScalars, nElems,neigh ))
      f010 = inState(?FETCH?( q010, 1, iElem, QQ, nScalars, nElems,neigh ))
      f001 = inState(?FETCH?( q001, 1, iElem, QQ, nScalars, nElems,neigh ))
      f0NN = inState(?FETCH?( q0NN, 1, iElem, QQ, nScalars, nElems,neigh ))
      f0N1 = inState(?FETCH?( q0N1, 1, iElem, QQ, nScalars, nElems,neigh ))
      f01N = inState(?FETCH?( q01N, 1, iElem, QQ, nScalars, nElems,neigh ))
      f011 = inState(?FETCH?( q011, 1, iElem, QQ, nScalars, nElems,neigh ))
      fN0N = inState(?FETCH?( qN0N, 1, iElem, QQ, nScalars, nElems,neigh ))
      f10N = inState(?FETCH?( q10N, 1, iElem, QQ, nScalars, nElems,neigh ))
      fN01 = inState(?FETCH?( qN01, 1, iElem, QQ, nScalars, nElems,neigh ))
      f101 = inState(?FETCH?( q101, 1, iElem, QQ, nScalars, nElems,neigh ))
      fNN0 = inState(?FETCH?( qNN0, 1, iElem, QQ, nScalars, nElems,neigh ))
      fN10 = inState(?FETCH?( qN10, 1, iElem, QQ, nScalars, nElems,neigh ))
      f1N0 = inState(?FETCH?( q1N0, 1, iElem, QQ, nScalars, nElems,neigh ))
      f110 = inState(?FETCH?( q110, 1, iElem, QQ, nScalars, nElems,neigh ))
      f000 = inState(?FETCH?( q000, 1, iElem, QQ, nScalars, nElems,neigh ))

      ! compute density and velocity
      ! element offset for auxField array
      elemOff = (iElem-1)*varSys%nAuxScalars
      ! local density
      rho = auxField(elemOff + dens_pos)
      ! local x-, y- and z-velocity
      u_x = auxField(elemOff + vel_pos(1))
      u_y = auxField(elemOff + vel_pos(2))
      u_z = auxField(elemOff + vel_pos(3))

      ! read the relaxation parameter omega for the current level
      omega = fieldProp(1)%fluid%viscKine%omLvl(level)%val(iElem)
      cmpl_o  = 1._rk - omega

      ! set output
      outState( ?SAVE?( qN00, 1, iElem, QQ, nScalars, nElems, neigh )) = &
        & cmpl_o * fN00 + omega * ( rho - 3 * rho0 * u_x ) * f2
      outState( ?SAVE?( q0N0, 1, iElem, QQ, nScalars, nElems, neigh )) = &
        & cmpl_o * f0N0 + omega * ( rho - 3 * rho0 * u_y ) * f2
      outState( ?SAVE?( q00N, 1, iElem, QQ, nScalars, nElems, neigh )) = &
        & cmpl_o * f00N + omega * ( rho - 3 * rho0 * u_z ) * f2
      outState( ?SAVE?( q100, 1, iElem, QQ, nScalars, nElems, neigh )) = &
        & cmpl_o * f100 + omega * ( rho + 3 * rho0 * u_x ) * f2
      outState( ?SAVE?( q010, 1, iElem, QQ, nScalars, nElems, neigh )) = &
        & cmpl_o * f010 + omega * ( rho + 3 * rho0 * u_y ) * f2
      outState( ?SAVE?( q001, 1, iElem, QQ, nScalars, nElems, neigh )) = &
        & cmpl_o * f001 + omega * ( rho + 3 * rho0 * u_z ) * f2
      outState( ?SAVE?( q0NN, 1, iElem, QQ, nScalars, nElems, neigh )) = &
        & cmpl_o * f0NN + omega * ( rho - 3 * rho0 * ( u_y + u_z ) ) * f8
      outState( ?SAVE?( q0N1, 1, iElem, QQ, nScalars, nElems, neigh )) = &
        & cmpl_o * f0N1 + omega * ( rho - 3 * rho0 * ( u_y - u_z ) ) * f8
      outState( ?SAVE?( q01N, 1, iElem, QQ, nScalars, nElems, neigh )) = &
        & cmpl_o * f01N + omega * ( rho + 3 * rho0 * ( u_y - u_z ) ) * f8
      outState( ?SAVE?( q011, 1, iElem, QQ, nScalars, nElems, neigh )) = &
        & cmpl_o * f011 + omega * ( rho + 3 * rho0 * ( u_y + u_z ) ) * f8
      outState( ?SAVE?( qN0N, 1, iElem, QQ, nScalars, nElems, neigh )) = &
        & cmpl_o * fN0N + omega * ( rho - 3 * rho0 * ( u_x + u_z ) ) * f8
      outState( ?SAVE?( qN01, 1, iElem, QQ, nScalars, nElems, neigh )) = &
        & cmpl_o * fN01 + omega * ( rho - 3 * rho0 * ( u_x - u_z ) ) * f8
      outState( ?SAVE?( q10N, 1, iElem, QQ, nScalars, nElems, neigh )) = &
        & cmpl_o * f10N + omega * ( rho + 3 * rho0 * ( u_x - u_z ) ) * f8
      outState( ?SAVE?( q101, 1, iElem, QQ, nScalars, nElems, neigh )) = &
        & cmpl_o * f101 + omega * ( rho + 3 * rho0 * ( u_x + u_z ) ) * f8
      outState( ?SAVE?( qNN0, 1, iElem, QQ, nScalars, nElems, neigh )) = &
        & cmpl_o * fNN0 + omega * ( rho - 3 * rho0 * ( u_x + u_y ) ) * f8
      outState( ?SAVE?( qN10, 1, iElem, QQ, nScalars, nElems, neigh )) = &
        & cmpl_o * fN10 + omega * ( rho - 3 * rho0 * ( u_x - u_y ) ) * f8
      outState( ?SAVE?( q1N0, 1, iElem, QQ, nScalars, nElems, neigh )) = &
        & cmpl_o * f1N0 + omega * ( rho + 3 * rho0 * ( u_x - u_y ) ) * f8
      outState( ?SAVE?( q110, 1, iElem, QQ, nScalars, nElems, neigh )) = &
        & cmpl_o * f110 + omega * ( rho + 3 * rho0 * ( u_x + u_y ) ) * f8
      outState( ?SAVE?( q000, 1, iElem, QQ, nScalars, nElems, neigh )) = &
        & cmpl_o * f000 + omega * ( rho ) * f1
    enddo nodeloop
!$omp end do nowait
?? END IF

  end subroutine bgk_advRel_isotherm_acEq_d3q19
! ****************************************************************************** !

end module mus_isotherm_acEq_module
! ****************************************************************************** !
