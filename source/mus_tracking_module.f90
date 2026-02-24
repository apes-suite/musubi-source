! Copyright (c) 2015-2016 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2015-2016 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2016 Raphael Haupt <raphael.haupt@uni-siegen.de>
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
!> author: Kannan Masilamani
!! This module provides the MUSUBI subroutines needed for the tracking
!! functionality.
!!
module mus_tracking_module

  ! include treelm modules
  use tem_logging_module,    only: logUnit
  use tem_tracking_module,   only: tem_init_tracker, tem_tracker,              &
    &                              tem_init_tracker_subTree
  use tem_debug_module,      only: dbgUnit

  ! include musubi modules
  use mus_scheme_type_module, only: mus_scheme_type
  use mus_param_module,       only: mus_param_type
  use mus_geom_module,        only: mus_geom_type
  use mus_tools_module,       only: mus_writeSolverSpecInfo
  use mus_time_module,        only: mus_timeControl_homogenize

  ! include libharvester modules
  use hvs_output_module, only: hvs_Internal

  implicit none
  private

  public :: mus_init_tracker


contains


! **************************************************************************** !
  !> This routine initialize tracking subTree to remove empty tracking objects.
  !! On active tracking objects: Homogenize time control, write solver speific
  !! info for harvester output format and initialize output using
  !! tem_init_tracker
  subroutine mus_init_tracker( scheme, geometry, params )
    ! --------------------------------------------------------------------------
    !> scheme type
    type( mus_scheme_type ), intent(inout) :: scheme
    !> Treelmesh data
    type( mus_geom_type ), intent(in) :: geometry
    !> Global parameters
    type( mus_param_type ), intent(in) :: params
    ! --------------------------------------------------------------------------
    integer :: iTrack, iConfig
    ! --------------------------------------------------------------------------

    select case (trim(scheme%header%kind))
    case ('fluid', 'fluid_incompressible')
      scheme%header%needAuxHaloComm = scheme%field(1)%fieldProp%fluid%turbulence%active
    end select

    write(dbgUnit(1),*) 'Enter mus_init_tracker'
    write(dbgUnit(1),*) 'Tracking control active is: ', &
      &                 scheme%track%control%active
    flush(dbgUnit(1))
    if ( .not. scheme%track%control%active ) return

    ! tracking objects
    write(logUnit(1),*) 'Initializing tracker...'

    write(dbgUnit(1),*) 'init tracker subTree'
    flush(dbgUnit(1))
    ! Initialize tracker subTree to remove empty tracking objects
    call tem_init_tracker_subTree(                           &
      &                    me      = scheme%track,           &
      &                    tree    = geometry%tree,          &
      &                    bc_prop = geometry%boundary,      &
      &                    stencil = scheme%layout%fStencil, &
      &                    solver  = params%general%solver )

    write(dbgUnit(1),*) 'homogenize track time'
    flush(dbgUnit(1))
    do iTrack = 1, scheme%track%control%nActive
      iConfig = scheme%track%instance(iTrack)%pntConfig
      call mus_timeControl_homogenize(                                &
        &        me     = scheme%track%config( iConfig )%timeControl, &
        &        dt     = params%physics                              &
        &                 %dtLvl( geometry%tree%global%maxLevel ),    &
        &        reqInt = params%reqInterval                          )

    end do !iTrack

    write(dbgUnit(1),*) 'to call tem init tracker'
    flush(dbgUnit(1))
    ! Creating tracking varMap with variable position in varSys and init
    ! output
    call tem_init_tracker( me       = scheme%track,           &
      &                    tree     = geometry%tree,          &
      &                    globProc = params%general%proc,    &
      &                    solver   = params%general%solver,  &
      &                    varSys   = scheme%varSys           )

    scheme%header%needAuxHaloComm = scheme%header%needAuxHaloComm      &
      &                     .or. mus_tracking_needs_velocity_gradient( &
      &                           scheme = scheme                      )

    if (scheme%header%needAuxHaloComm) then
      write(logUnit(1),*) 'Auxiliary field halo communication activated.'
    else
      write(logUnit(1),*) 'Auxiliary field halo communication deactivated.'
    end if

  end subroutine mus_init_tracker
! **************************************************************************** !


! **************************************************************************** !
  !> Check whether any active tracking object requests velocity-gradient based
  !! variables that require auxField halo values.
  !!
  !! TODO rather than setting this flag according to the tracking variables,
  !!      change the computation in the tracking variables based on this
  !!      flag and compute the values on the fly if they are not available
  !!      in the aux field.
  logical function mus_tracking_needs_velocity_gradient(scheme)
    ! --------------------------------------------------------------------------
    type(mus_scheme_type), intent(in) :: scheme
    ! --------------------------------------------------------------------------
    integer :: iTrack, iConfig, iVar
    ! --------------------------------------------------------------------------

    mus_tracking_needs_velocity_gradient = .false.

    if (.not. scheme%track%control%active) return

    select case (trim(scheme%header%kind))
    case ('fluid', 'fluid_incompressible')
      continue
    case default
      return
    end select

    do iTrack = 1, scheme%track%control%nActive
      iConfig = scheme%track%instance(iTrack)%pntConfig
      do iVar = 1, size(scheme%track%config(iConfig)%varName)
        select case(trim(scheme%track%config(iConfig)%varName(iVar)))
        case ('grad_velocity', 'vorticity', 'q_criterion',             &
              & 'grad_velocity_phy', 'vorticity_phy', 'q_criterion_phy')
          mus_tracking_needs_velocity_gradient = .true.
          return
        end select
      end do
    end do
  end function mus_tracking_needs_velocity_gradient
! **************************************************************************** !

end module mus_tracking_module
! **************************************************************************** !
