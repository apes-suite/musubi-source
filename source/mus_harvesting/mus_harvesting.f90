! Copyright (c) 2015-2019 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2016 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2016 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2016-2017 Raphael Haupt <raphael.haupt@uni-siegen.de>
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
!> Musubi Harvesting Tool
!! Visualization of restart file or tracking harvester formar, generated by
!! musubi
!! (c) 2015 University of Siegen
!!
!! For a documentation, run ./waf gendoxy and find the documentation at
!! ./Documentation/html/index.html
program mus_harvesting
  use iso_c_binding,            only: c_loc

  ! treelm modules
  use mpi
  use env_module,               only: pathLen
  use tem_general_module,       only: tem_start, tem_finalize
  use tem_timeControl_module,   only: tem_timeControl_start_at_sim
  use tem_tracking_module,      only: tem_tracking_finalize, &
    &                                 tem_tracking_print_last_VTK_files
  use tem_aux_module,           only: tem_abort
  use tem_logging_module,       only: logUnit

  ! musubi modules
  use mus_scheme_type_module,        only: mus_scheme_type
  use mus_scheme_module,             only: mus_init_scheme
  use mus_param_module,              only: mus_param_type
  use mus_timer_module,              only: mus_init_mainTimer, &
    &                                      mus_init_levelTimer
  use mus_geom_module,               only: mus_geom_type
  use mus_varSys_module,             only: mus_varSys_solverData_type
  use mus_restart_module,            only: mus_readRestart
  use mus_flow_module,               only: fillHelperElementsCoarseToFine, &
    &                                      fillHelperElementsFineToCoarse, &
    &                                      mus_initAuxField
  use mus_construction_module,       only: mus_construct

  ! libharvesting
  use hvs_output_module,        only: hvs_output_file_type, hvs_output_open,   &
    &                                 hvs_output_close, hvs_output_write,      &
    &                                 hvs_output_init, hvs_output_finalize
  use hvs_aux_module,           only: hvs_banner
  use mus_bndForce_module,      only: mus_calcBndForce

  ! mus_harvesting
  use mus_hvs_config_module,       only: mus_hvs_config_type,                  &
    &                                    mus_hvs_config_load
  use mus_hvs_aux_module,          only: mus_hvs_init_aux
  use mus_hvs_construction_module, only: mus_hvs_construct

  ! aotus
  use aotus_module,             only: flu_State

  implicit none

  ! -----------------------------------------------------------------------------
  type(mus_scheme_type), target :: scheme
  type(mus_geom_type), target :: geometry
  type(mus_param_type), target :: params
  type(mus_varSys_solverData_type), target :: solverData
  type(mus_hvs_config_type)  :: config
  type(hvs_output_file_type) :: out_file
  integer :: ierr
  integer :: iTrack, iConfig
  integer :: nVars
  integer :: minLevel, maxLevel
  ! basename for output if tracking table is not defined
  ! trim(config%prefix)//trim(general%solver%simName)//trim(varSys%systemName)
  character(len=pathLen) :: basename
  ! -----------------------------------------------------------------------------

  ! Initialize environment.
  call tem_start( codeName = 'mus_harvesting', &
    &             general  = params%general    )

  if (params%general%proc%rank == 0) then
    call hvs_banner( solveHead      = params%general%solver, &
      &              supportSolName = 'Musubi'               )
  end if

  call mus_init_mainTimer()

  ! Load mus_harvesting config file
  ! From this config file, try to load restart file from restart table key "read"
  ! If no restart read defined then abort.
  !
  ! If restart read is defined then load basic scheme info
  ! like identify table, field, physics table from restart header file
  !
  ! Load tracking and variable table from config file
  call mus_hvs_config_load( me         = config,     &
    &                       scheme     = scheme,     &
    &                       solverData = solverData, &
    &                       geometry   = geometry,   &
    &                       params     = params      )

  minLevel = geometry%tree%global%minLevel
  maxLevel = geometry%tree%global%maxLevel

  call mus_init_levelTimer( minLevel, maxLevel )

  ! initialize scheme
  ! append variables loaded from restart header
  ! file as state variables
  call mus_init_scheme( me         = scheme,        &
    &                   tree       = geometry%tree, &
    &                   solverData = solverData     )

  !! If restart read variable is not pdf then set interpolation
  !! method to none to deactivate interpolation for
  !! derived variable
  if ( .not. scheme%readVarIsPdf ) then
    scheme%intp%config%method = 'none'
    if ( minLevel /= maxLevel ) then
      write(logUnit(1), *) 'WARNING: Multi-level interpolation is not supported'
      write(logUnit(1), *) 'for derived variables read from restart file'
    end if
  end if

  ! construct levelDescriptor, connectivity array
  if ( scheme%readVarIsPdf ) then
    call mus_construct( scheme    = scheme,   &
      &                 geometry  = geometry, &
      &                 params    = params    )
  else
    ! Initialize only fluid element list. Ghost and halo elements are
    ! not required if variable read from restart is not PDF.
    call mus_hvs_construct( scheme    = scheme,   &
      &                      geometry = geometry, &
      &                      params   = params    )
  end if


  ! scheme%state%val array is allocated with single buffer so set nNow and nNext to 1
  scheme%pdf(minLevel:maxLevel)%nNow  = 1
  scheme%pdf(minLevel:maxLevel)%nNext = 1

  ! Init auxiliary features such as interpolation, boundaries, restart
  ! and the tracker
  call mus_hvs_init_aux( scheme    = scheme,   &
    &                    geometry  = geometry, &
    &                    params    = params    )

  ! Read restart file
  call mus_readRestart( levelPointer = geometry%levelPointer,  &
    &                   restart      = params%general%restart, &
    &                   scheme       = scheme,                 &
    &                   tree         = geometry%tree           )

  if ( scheme%readVarIsPdf ) then

    ! init auxiliary field variable from state for fluid elements in state and
    ! interpolate for ghost elements in do_intpArbiVal routine.
    call mus_initAuxField( scheme, params%general, minLevel, maxLevel )

    ! Fill all elements (ghost, halo) with valid values from fluid elements
    call fillHelperElementsFineToCoarse( scheme     = scheme,         &
      &                                  general    = params%general, &
      &                                  physics    = params%physics, &
      &                                  iLevel     = minLevel,       &
      &                                  maxLevel   = maxLevel        )

    call fillHelperElementsCoarseToFine( scheme     = scheme,         &
      &                                  general    = params%general, &
      &                                  physics    = params%physics, &
      &                                  iLevel     = minLevel,       &
      &                                  minLevel   = minLevel,       &
      &                                  maxLevel   = maxLevel        )

    ! Force on boundary elements are computed from post-collision so calculate
    ! here after apply source and before tracking output
    call mus_calcBndForce( bndForce   = geometry%bndForce,          &
      &                    bndMoment  = geometry%bndMoment,         &
      &                    posInBndID = geometry%posInBndID,        &
      &                    nBCs       = geometry%boundary%nBCtypes, &
      &                    field      = scheme%field,               &
      &                    globBC     = scheme%globBC,              &
      &                    minLevel   = minLevel,                   &
      &                    maxLevel   = maxLevel,                   &
      &                    state      = scheme%state,               &
      &                    pdf        = scheme%pdf,                 &
      &                    levelDesc  = scheme%levelDesc,           &
      &                    layout     = scheme%layout,              &
      &                    varSys     = scheme%varSys,              &
      &                    physics    = params%physics              )
  end if


  call mpi_barrier( MPI_COMM_WORLD, ierr )

  ! For active tracking output is initialized in init_tracker in ini_aux
  if ( scheme%track%control%active ) then

   do iTrack = 1, scheme%track%control%nActive
     iConfig = scheme%track%instance(iTrack)%pntConfig
     if( scheme%track%instance(iTrack)%subTree%useGlobalMesh ) then
       ! Open the output files, this also generates the vertices for the mesh,
       ! and writes the mesh data to disk. Also writes header file depends
       ! on output vis_kind
       call hvs_output_open(                                             &
         &         out_file = scheme%track%instance(iTrack)%output_file, &
         &         mesh     = geometry%tree,                             &
         &         varSys   = scheme%varSys,                             &
         &         time     = params%general%simControl%now              )

       ! Evaluate and write results to disk
       call hvs_output_write(                                            &
         &         out_file = scheme%track%instance(iTrack)%output_file, &
         &         varSys   = scheme%varSys,                             &
         &         mesh     = geometry%tree                              )

       ! Close opened files
       call hvs_output_close(                                            &
         &         out_file = scheme%track%instance(iTrack)%output_file, &
         &         varSys   = scheme%varSys,                             &
         &         mesh     = geometry%tree                              )
     else
       ! Open the output files, this also generates the vertices for the mesh,
       ! and writes the mesh data to disk. Also writes header file depends
       ! on output vis_kind
       call hvs_output_open(                                             &
         &         out_file = scheme%track%instance(iTrack)%output_file, &
         &         mesh     = geometry%tree,                             &
         &         varSys   = scheme%varSys,                             &
         &         subTree  = scheme%track%instance(iTrack)%subTree,     &
         &         time     = params%general%simControl%now              )

       ! Evaluate and write results to disk
       call hvs_output_write(                                            &
         &         out_file = scheme%track%instance(iTrack)%output_file, &
         &         varSys   = scheme%varSys,                             &
         &         mesh     = geometry%tree,                             &
         &         subTree  = scheme%track%instance(iTrack)%subTree      )

       ! Close opened files
       call hvs_output_close(                                            &
         &         out_file = scheme%track%instance(iTrack)%output_file, &
         &         varSys   = scheme%varSys,                             &
         &         mesh     = geometry%tree,                             &
         &         subTree  = scheme%track%instance(iTrack)%subTree      )
     end if !Global mesh
    end do !iTrack
    ! Finialize output
    call tem_tracking_finalize( scheme%track )
  else
    ! number of variables loaded from restart file
    nVars = scheme%stateVarMap%varName%nVals

    ! filename for output
    basename = trim(config%prefix) // trim(params%general%solver%simName)

    ! De-activate transient and spatial reduction
    !out_file%isTransientReduce = .false.
    out_file%ascii%isReduce = .false.

    ! No tracking defined, just dump the restart file
    ! Open the output files, this also generates the vertices for the mesh,
    ! and writes the mesh data to disk.
    call hvs_output_init(out_file    = out_file,                 &
      &                  out_config  = config%output,            &
      &                  tree        = geometry%tree,            &
      &                  varSys      = scheme%varSys,            &
      &                  varPos      = scheme%stateVarMap%varPos &
      &                                %val(:nVars),             &
      &                  basename    = trim(basename),           &
      &                  globProc    = params%general%proc,      &
      &                  solver      = params%general%solver     )

    ! Open output file handle
    call hvs_output_open( out_file   = out_file,                     &
      &                   mesh       = geometry%tree,                &
      &                   varsys     = scheme%varsys,                &
      &                   time       = params%general%simControl%now )

    ! Fill output files with data.
    call hvs_output_write( out_file = out_file,      &
      &                    varsys   = scheme%varsys, &
      &                    mesh     = geometry%tree  )

    ! Close output files again.
    call hvs_output_close( out_file = out_file,      &
      &                    varsys   = scheme%varsys, &
      &                    mesh     = geometry%tree  )

    ! Finialize output
    call hvs_output_finalize( out_file = out_file )
  end if

  ! Finalize environment.
  call tem_finalize(params%general)

end program mus_harvesting
