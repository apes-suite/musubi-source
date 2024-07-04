! Copyright (c) 2011-2013 Manuel Hasert <m.hasert@grs-sim.de>
! Copyright (c) 2011-2014 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2011, 2013 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2011-2016, 2019-2020 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2011 Jan Hueckelheim <j.hueckelheim@grs-sim.de>
! Copyright (c) 2012-2016 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2012, 2014 Kartik Jain <kartik.jain@uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
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
?? include 'header/lbm_macros.inc'
?? include 'tem/source/logMacros.inc'
! ****************************************************************************** !
!> This module provides the definition and methods for
!! BGK advection relaxation scheme.
module mus_bgk_module
  use iso_c_binding, only: c_f_pointer

  ! include treelm modules
  use env_module,               only: rk
  use tem_varSys_module,        only: tem_varSys_type !, tem_varSys_op_type
  use tem_param_module,         only: div1_3, div1_36, div1_8, div3_4h, div1_4,&
    &                                 div3_8, div9_16, div3_16, cs2inv, cs4inv,&
    &                                 rho0
  ! use tem_spacetime_fun_module, only: tem_st_fun_listElem_type,                &
  !   &                                 tem_spacetime_for
  ! use tem_logging_module,       only: logUnit
  ! use tem_dyn_array_module,     only: PositionOfVal
  ! use tem_aux_module,           only: tem_abort

  ! include musubi modules
  use mus_field_prop_module,    only: mus_field_prop_type
  use mus_scheme_layout_module, only: mus_scheme_layout_type
  use mus_scheme_type_module,   only: mus_scheme_type
  use mus_derVarPos_module,     only: mus_derVarPos_type
  use mus_param_module,         only: mus_param_type
  use mus_varSys_module,        only: mus_varSys_data_type

  implicit none

  private

  public :: mus_advRel_kCFD_rBGK_vStdNoOpt_l
  public :: mus_advRel_kPS_rBGK_v1st
  public :: bgk_advRel_flekkoy_noFluid
  public :: mus_advRel_kPS_rBGK_v2nd
  public :: mus_advRel_kPS_rTRT_vStdNoOpt

contains

! **************************************************************************** !
  !> Advection relaxation routine for the flekkoy diffusion model.
  !!
  !! This subroutine interface must match the abstract interface definition
  !! [[kernel]] in scheme/[[mus_scheme_type_module]].f90 in order to be callable
  !! via [[mus_scheme_type:compute]] function pointer.
  subroutine mus_advRel_kPS_rBGK_v1st( fieldProp, inState, outState, auxField, &
    &                            neigh, nElems, nSolve, level, layout,   &
    &                            params, varSys, derVarPos               )
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
    integer :: iElem, iDir
    type(mus_varSys_data_type), pointer :: fPtr
    type(mus_scheme_type), pointer :: scheme
    real(kind=rk) :: pdfTmp( layout%fStencil%QQ ) ! temporary local pdf values
    real(kind=rk) :: rho, feq
    real(kind=rk) :: d_omega
    real(kind=rk) :: transVel( nSolve*3 ) ! velocity from the transport field
    real(kind=rk) :: uc ! u_i,fluid * c_i
    integer :: vel_varPos ! position of transport velocity variable in varSys
    real(kind=rk) :: inv_vel, u_fluid(3)
    ! --------------------------------------------------------------------------
    ! access scheme via 1st variable method data which is a state variable
    call C_F_POINTER( varSys%method%val(derVarPos(1)%pdf)%method_Data, fPtr )
    scheme => fPtr%solverData%scheme

    ! passive scalar has only one transport Variable
    vel_varPos = scheme%transVar%method(1)%data_varPos
    ! Get velocity field
    call varSys%method%val(vel_varPos)%get_valOfIndex( &
      & varSys  = varSys,                              &
      & time    = params%general%simControl%now,       &
      & iLevel  = level,                               &
      & idx     = scheme%transVar%method(1)            &
      &           %pntIndex%indexLvl(level)            &
      &           %val(1:nSolve),                      &
      & nVals   = nSolve,                              &
      & res     = transVel                             )

    ! convert physical velocity into LB velocity
    inv_vel = 1.0_rk / params%physics%fac( level )%vel
    transVel = transVel * inv_vel

    ! initialize and define some field wise parameters
    d_omega = 2._rk / ( 1._rk + 6._rk                               &
      &                       * fieldProp(1)%species%diff_coeff( 1 ))

    nodeloop: do iElem = 1, nSolve
      ! x-, y- and z-velocity from transport field
      u_fluid = transVel( (iElem-1)*3+1 : iElem*3 )

      do iDir = 1, layout%fStencil%QQ
        pdfTmp( iDir ) = &
& instate( ?FETCH?( iDir, 1, iElem, layout%fStencil%QQ, varSys%nScalars, nElems,neigh ) )
      end do
      rho = sum( pdfTmp )

      do iDir = 1, layout%fStencil%QQ
        ! compute c_i * u
        uc = dble( layout%fStencil%cxDir(1, iDir)) * u_fluid(1) + &
          &  dble( layout%fStencil%cxDir(2, iDir)) * u_fluid(2) + &
          &  dble( layout%fStencil%cxDir(3, iDir)) * u_fluid(3)

        ! compute the equilibrium (fi_eq = weight_i * rho * ( 1+c_i*u / cs^2))
        feq = rho * layout%weight( iDir ) * (1._rk+3._rk*uc)

        outstate(                                                            &
& ?SAVE?( iDir, 1, iElem, layout%fStencil%QQ, varSys%nScalars, nElems,neigh ) ) =     &
          &                pdfTmp( iDir ) + d_omega * ( feq - pdfTmp( idir ))
      end do

    end do nodeloop

  end subroutine mus_advRel_kPS_rBGK_v1st
! ****************************************************************************** !


! ****************************************************************************** !
  !> Advection relaxation routine for the flekkoy diffusion model without prior
  !! definition of a dependent fluid scheme. Velocity is read from the depend
  !! table.
  !!
  !! This subroutine interface must match the abstract interface definition
  !! [[kernel]] in scheme/[[mus_scheme_type_module]].f90 in order to be callable
  !! via [[mus_scheme_type:compute]] function pointer.
  subroutine bgk_advRel_flekkoy_noFluid( fieldProp, inState, outState,    &
    &                                    auxField, neigh, nElems, nSolve, &
    &                                    level, layout, params, varSys,   &
    &                                    derVarPos                        )
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
    integer :: iElem
    integer :: iDir
    integer :: QQ, nScalars
    real(kind=rk) :: pdfTmp( layout%fStencil%QQ ) ! temporary local pdf values
    real(kind=rk) :: rho
    real(kind=rk) :: feq
    real(kind=rk) :: d_omega
    ! ---------------------------------------------------------------------------
    QQ       = layout%fStencil%QQ
    nScalars = varSys%nScalars

    ! initialize and define some field wise parameters
    d_omega = 2._rk / ( 1._rk + 6._rk * fieldProp(1)%species%diff_coeff(1) )

    nodeloop: do iElem = 1, nSolve

      rho = 0._rk
      do iDir = 1, QQ
        ! store the pdfs for this element and compute the density locally
        pdfTmp( iDir ) = &
          & instate( ?FETCH?( iDir, 1, iElem, QQ, nScalars, nElems,neigh ))
        rho = rho + pdfTmp( iDir )
      end do

      directionloop: do iDir = 1, QQ
        ! compute the equilibrium
        ! fi_eq = weight_i * rho * ( 1+c_i*u / cs^2))
        !       = weight_i * rho
        feq = rho * layout%weight( iDir )

        outstate( ?SAVE?( iDir, 1, iElem, QQ, nScalars, nElems,neigh ) ) =  &
          &                pdfTmp( iDir ) + d_omega * ( feq - pdfTmp( idir ))
      end do directionloop
    end do nodeloop

  end subroutine bgk_advRel_flekkoy_noFluid
! ****************************************************************************** !


! ****************************************************************************** !
  !> Advection relaxation routine for the
  !! BGK model with an explicit calculation of all equilibrium
  !! quantities. Slow and simple. This routine should only be
  !! used for testing purposes
  !!
  !! This subroutine interface must match the abstract interface definition
  !! [[kernel]] in scheme/[[mus_scheme_type_module]].f90 in order to be callable
  !! via [[mus_scheme_type:compute]] function pointer.
  subroutine mus_advRel_kCFD_rBGK_vStdNoOpt_l( fieldProp, inState, outState,   &
    &                                          auxField, neigh, nElems, nSolve,&
    &                                          level, layout, params, varSys,  &
    &                                          derVarPos )
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
    integer :: iElem, iDir                       ! voxel element counter
    integer :: QQ, nScalars
    ! temporary distribution variables
    real(kind=rk) pdfTmp(layout%fStencil%QQ)
    real(kind=rk) ux(3)   ! local velocity
    real(kind=rk) rho     ! local density
    ! derived constants
    ! equilibrium calculation variables
    real(kind=rk) fEq(layout%fStencil%QQ)
    real(kind=rk) omega
    integer :: dens_pos, vel_pos(3), elemOff
    ! ---------------------------------------------------------------------------
    dens_pos = varSys%method%val(derVarPos(1)%density)%auxField_varPos(1)
    vel_pos = varSys%method%val(derVarPos(1)%velocity)%auxField_varPos(1:3)

    QQ = layout%fStencil%QQ
    ! nElems = size(neigh)/QQ
    nScalars = varSys%nScalars

    nodeloop: do iElem=1,nSolve
      !> Generic fetching step:
      !! Streaming for pull
      !! Local copy for push
      ux = 0._rk
      rho = 0._rk
      do iDir = 1, QQ
        pdfTmp( iDir ) = inState( ?FETCH?( iDir, 1, iElem, QQ, nScalars, nElems,neigh))
      end do

      ! element offset for auxField array
      elemOff = (iElem-1)*varSys%nAuxScalars
      ! local density
      rho = auxField(elemOff + dens_pos)
      ! local x-, y- and z-velocity
      ux(1) = auxField(elemOff + vel_pos(1))
      ux(2) = auxField(elemOff + vel_pos(2))
      ux(3) = auxField(elemOff + vel_pos(3))

      !> relaxation parameter
      omega = fieldProp(1)%fluid%viscKine%omLvl(level)%val(iElem)

      fEq = layout%quantities%pdfEq_ptr(rho = rho, vel = ux, QQ = QQ)

      do iDir = 1, QQ
        !> Relaxation
        outState( ?SAVE?( iDir, 1, iElem, QQ, nScalars, nElems,neigh )) =     &
          &            pdfTmp( iDir ) - omega*( pdfTmp(iDir ) - fEq( iDir ))
      end do  !< iDir

    end do nodeloop

  end subroutine mus_advRel_kCFD_rBGK_vStdNoOpt_l
! ****************************************************************************** !


! **************************************************************************** !
!> Advection relaxation routine for the 2nd order diffusion model.

!! A comparison to the previous flekkoy model can be found in 
!! Chopard, B., Falcone, J. & Latt, J. "The lattice Boltzmann advection
!! -diffusion model revisited." Eur. Phys. J. Spec. Top. 171, 245–249 (2009). 
!! https://doi.org/10.1140/epjst/e2009-01035-5
!!
!! This subroutine interface must match the abstract interface definition
!! [[kernel]] in scheme/[[mus_scheme_type_module]].f90 in order to be callable
!! via [[mus_scheme_type:compute]] function pointer.
  subroutine mus_advRel_kPS_rBGK_v2nd( fieldProp, inState, outState, auxField, &
    &                            neigh, nElems, nSolve, level, layout,   &
    &                            params, varSys, derVarPos               )
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
    integer :: iElem, iDir
    type(mus_varSys_data_type), pointer :: fPtr
    type(mus_scheme_type), pointer :: scheme
    real(kind=rk) :: pdfTmp( layout%fStencil%QQ ) ! temporary local pdf values
    real(kind=rk) :: rho, feq
    real(kind=rk) :: d_omega
    real(kind=rk) :: transVel( nSolve*3 ) ! velocity from the transport field
    real(kind=rk) :: uc, usq ! u_i,fluid * c_i
    integer :: vel_varPos ! position of transport velocity variable in varSys
    real(kind=rk) :: inv_vel, u_fluid(3)
    ! --------------------------------------------------------------------------
    ! access scheme via 1st variable method data which is a state variable
    call C_F_POINTER( varSys%method%val(derVarPos(1)%pdf)%method_Data, fPtr )
    scheme => fPtr%solverData%scheme
  
    ! passive scalar has only one transport Variable
    vel_varPos = scheme%transVar%method(1)%data_varPos
    ! Get velocity field
    call varSys%method%val(vel_varPos)%get_valOfIndex( &
      & varSys  = varSys,                              &
      & time    = params%general%simControl%now,       &
      & iLevel  = level,                               &
      & idx     = scheme%transVar%method(1)            &
      &           %pntIndex%indexLvl(level)            &
      &           %val(1:nSolve),                      &
      & nVals   = nSolve,                              &
      & res     = transVel                             )
  
    ! convert physical velocity into LB velocity
    inv_vel = 1.0_rk / params%physics%fac( level )%vel
    transVel = transVel * inv_vel
  
    ! initialize and define some field wise parameters
    d_omega = 2._rk / ( 1._rk + 6._rk                               &
      &                       * fieldProp(1)%species%diff_coeff( 1 ))
  
    nodeloop: do iElem = 1, nSolve
      ! x-, y- and z-velocity from transport field
      u_fluid = transVel( (iElem-1)*3+1 : iElem*3 )
  
      do iDir = 1, layout%fStencil%QQ
        pdfTmp( iDir ) = &
  & instate( ?FETCH?( iDir, 1, iElem, layout%fStencil%QQ, varSys%nScalars, nElems,neigh ) )
      end do
      rho = sum( pdfTmp )
  
      do iDir = 1, layout%fStencil%QQ
        ! compute c_i * u
        uc = dble( layout%fStencil%cxDir(1, iDir)) * u_fluid(1) + &
          &  dble( layout%fStencil%cxDir(2, iDir)) * u_fluid(2) + &
          &  dble( layout%fStencil%cxDir(3, iDir)) * u_fluid(3)
  
        usq = u_fluid(1)*u_fluid(1) + u_fluid(2)*u_fluid(2) + u_fluid(3)*u_fluid(3)
  
        ! compute the equilibrium (fi_eq = weight_i * rho * ( 1+c_i*u / cs^2))
        feq = rho * layout%weight( iDir ) * (1._rk + 3._rk*uc    &
        &           +  9._rk*uc*uc*0.5_rk     &
        &           -  usq*0.5_rk*3._rk )
  
        outstate(                                                            &
  & ?SAVE?( iDir, 1, iElem, layout%fStencil%QQ, varSys%nScalars, nElems,neigh ) ) =     &
          &                pdfTmp( iDir ) + d_omega * ( feq - pdfTmp( idir ))
      end do
  
    end do nodeloop
  
    end subroutine mus_advRel_kPS_rBGK_v2nd
  ! ****************************************************************************** !
  

! **************************************************************************** !
!> Advection relaxation routine for the TRT diffusion model.
!!
!!   Irina Ginzburg (2005), "Equilibrium-type and link-type lattice Boltzmann 
!!   models for generic advection and anisotropic-dispersion equation",
!!   Advances in Water Resources, Volume 28, Issue 11
!!
!! This subroutine interface must match the abstract interface definition
!! [[kernel]] in scheme/[[mus_scheme_type_module]].f90 in order to be callable
!! via [[mus_scheme_type:compute]] function pointer.
    subroutine mus_advRel_kPS_rTRT_vStdNoOpt( fieldProp, inState, outState, auxField, &
      &                            neigh, nElems, nSolve, level, layout,   &
      &                            params, varSys, derVarPos               )
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
      integer :: iElem, iDir, invDir
      type(mus_varSys_data_type), pointer :: fPtr
      type(mus_scheme_type), pointer :: scheme
      real(kind=rk) :: pdfTmp( layout%fStencil%QQ ) ! temporary local pdf values
      real(kind=rk) :: rho, feqPlus, feqMinus, fPlus, fMinus, magicParam
      real(kind=rk) :: d_omega, aux_omega
      real(kind=rk) :: transVel( nSolve*3 ) ! velocity from the transport field
      real(kind=rk) :: uc, usq ! u_i,fluid * c_i
      integer :: vel_varPos ! position of transport velocity variable in varSys
      real(kind=rk) :: inv_vel, u_fluid(3)
      ! --------------------------------------------------------------------------
      ! access scheme via 1st variable method data which is a state variable
      call C_F_POINTER( varSys%method%val(derVarPos(1)%pdf)%method_Data, fPtr )
      scheme => fPtr%solverData%scheme
    
      ! passive scalar has only one transport Variable
      vel_varPos = scheme%transVar%method(1)%data_varPos
      ! Get velocity field
      call varSys%method%val(vel_varPos)%get_valOfIndex( &
        & varSys  = varSys,                              &
        & time    = params%general%simControl%now,       &
        & iLevel  = level,                               &
        & idx     = scheme%transVar%method(1)            &
        &           %pntIndex%indexLvl(level)            &
        &           %val(1:nSolve),                      &
        & nVals   = nSolve,                              &
        & res     = transVel                             )
    
      ! convert physical velocity into LB velocity
      inv_vel = 1.0_rk / params%physics%fac( level )%vel
      transVel = transVel * inv_vel
    
      ! initialize and define some field wise parameters
      d_omega = 2._rk / ( 1._rk + 6._rk                               &
        &                       * fieldProp(1)%species%diff_coeff( 1 ))
      magicParam = fieldProp(1)%species%lambda
      aux_omega = 1.0_rk / (magicParam / (1.0_rk / d_omega - 0.5_rk) + 0.5_rk)
    
      nodeloop: do iElem = 1, nSolve
        ! x-, y- and z-velocity from transport field
        u_fluid = transVel( (iElem-1)*3+1 : iElem*3 )
  
        do iDir = 1, layout%fStencil%QQ
          pdfTmp( iDir ) = &
    & instate( ?FETCH?( iDir, 1, iElem, layout%fStencil%QQ, varSys%nScalars, nElems,neigh ) )
        end do
        rho = sum( pdfTmp )
    
        do iDir = 1, layout%fStencil%QQ
          ! compute c_i * u
          uc = dble( layout%fStencil%cxDir(1, iDir)) * u_fluid(1) + &
            &  dble( layout%fStencil%cxDir(2, iDir)) * u_fluid(2) + &
            &  dble( layout%fStencil%cxDir(3, iDir)) * u_fluid(3)
    
          usq = u_fluid(1)*u_fluid(1) + u_fluid(2)*u_fluid(2) + u_fluid(3)*u_fluid(3)
  
          ! compute the equilibrium (fi_eq = weight_i * rho * ( 1+c_i*u / cs^2))
          feqPlus = rho * layout%weight( iDir ) * (1._rk    &
          &           +  9._rk*uc*uc*0.5_rk     &
          &           -  usq*0.5_rk*3._rk )
  
          feqMinus = rho * layout%weight( iDir ) * 3._rk * uc
  
          invDir = layout%fStencil%cxDirInv(iDir)
          fPlus = 0.5_rk * (pdfTmp(iDir) + pdfTmp(invDir))
          fMinus = 0.5_rk * (pdfTmp(iDir) - pdfTmp(invDir))
          
          outstate(                                                            &
    & ?SAVE?( iDir, 1, iElem, layout%fStencil%QQ, varSys%nScalars, nElems,neigh ) ) =     &
            &                pdfTmp( iDir ) + d_omega * ( feqMinus - fMinus)              &
            &                + aux_omega * (feqPlus - fPlus)
        end do
    
      end do nodeloop
    
      end subroutine mus_advRel_kPS_rtrt_vStdNoOpt
    ! ****************************************************************************** !
    

end module mus_bgk_module
! ****************************************************************************** !
