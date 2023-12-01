! See copyright notice in the COPYRIGHT file.
!> This program test state variable access routines
?? include 'header/lbm_macros.inc'
program mus_varSys_accessStateVar_test
  use, intrinsic :: iso_c_binding,  only: C_NEW_LINE, c_loc

  use env_module,                   only: rk, eps, globalMaxLevels

  use treelmesh_module,             only: treelmesh_type
  use tem_geometry_module,          only: tem_CoordOfReal, &
    &                                     tem_PosofId, tem_BaryOfId, &
    &                                     tem_ElemSize
  use tem_topology_module,          only: tem_IdOfCoord, &
    &                                     tem_levelOf
  use tem_tools_module,             only: tem_horizontalSpacer
  use tem_logging_module,            only: logUnit
  use tem_varSys_module,             only: tem_varSys_init
  use tem_general_module,      only: tem_finalize
  use tem_construction_module, only: tem_init_elemLevels, tem_find_allElements,&
    &                                tem_build_verticalDependencies,           &
    &                                tem_build_horizontalDependencies,         &
    &                                tem_cleanup_arrays

  use mus_connectivity_module,  only: mus_construct_connectivity
  use mus_scheme_layout_module, only: mus_define_d3q19
  use mus_param_module,         only: mus_param_type
  use mus_geom_module,          only: mus_geom_type
  use mus_scheme_type_module,   only: mus_scheme_type
  use mus_variable_module,      only: mus_append_stateVar
  use mus_varSys_module,        only: mus_varSys_data_type,       &
    &                                 mus_varSys_solverData_type, &
    &                                 mus_init_varSys_solverData

  use mus_utestEnv_module, only: load_env

  implicit none

  !****************************************************************************!
  ! PAREMETER
  ! write output to screen
  ! number of elements to track
  integer, parameter :: nElems_track = 5
  ! position in global treeID list to track
  integer, dimension(nElems_track), parameter :: elemPos = (/ 1, 3, 5, 7, 8 /)
  !*****************************************************************************
  type(mus_geom_type), target :: geometry
  type(mus_scheme_type), target :: scheme
  type(mus_param_type), target :: params

  type(mus_varSys_solverData_type), target :: solverData

  integer :: iElem, iComp, iLevel
  integer :: level, nElems
  logical :: res_correct(3)
  ! Number of state var components
  integer :: QQ, nScalars, varPos
  real(kind=rk), allocatable :: ref_state(:)
  real(kind=rk), allocatable :: res(:)
  real(kind=rk), allocatable :: diff(:)
  real(kind=rk) :: point(nElems_track,3)
  integer :: idx(nElems_track)
  integer :: statePos, nSize
  !*****************************************************************************
  write(*,*) 'UTEST: varSys access state variable'

  ! Load mesh with 8 elements in refinement level 1
  ! and init treElm environment
  call load_env( tree     = geometry%tree,     &
    &            boundary = geometry%boundary, &
    &            general  = params%general     ) 

  ! minlevel must be one
  level    = geometry%tree%global%minLevel
  if (level/=1) then
    write(logUnit(1),*) 'Error: Cannot continue for level /=1'
    write(logUnit(1),*) 'FAILED'
    stop
  end if
  nElems = geometry%tree%nElems

  write(logUnit(1),*)'Initializing musubi scheme'
  scheme%nFields = 1
  allocate(scheme%field(1))
  allocate(scheme%field(1)%bc(0))
  scheme%field(1)%label = 'test_'

  scheme%header%kind       = 'lbm'
  scheme%header%layout     = 'd3q19'
  scheme%header%relaxation = 'bgk'

  call mus_define_d3q19( layout = scheme%layout, nElems = nElems )
  QQ = scheme%layout%fStencil%QQ

  ! Set solverData pointer for variable method data
  call mus_init_varSys_solverData( me        = solverData,     &
    &                              scheme    = scheme,         &
    &                              physics   = params%physics, &
    &                              geometry  = geometry        )

  write(logUnit(1),*)'Initializing musubi variables'
  call tem_varSys_init( me         = scheme%varSys,            &
    &                   systemName = trim(scheme%header%kind), &
    &                   length     = 8                         )

  ! append state variable depends on scheme kind
  call mus_append_stateVar( varSys       = scheme%varSys,          &
    &                       stateVarMap  = scheme%stateVarMap,     &
    &                       solverData   = solverData,             &
    &                       schemeHeader = scheme%header,          &
    &                       stencil      = scheme%layout%fstencil, &
    &                       nFields      = scheme%nFields,         &
    &                       fldLabel     = scheme%field(:)%label   )
  nScalars = scheme%varSys%nScalars  
  statePos = scheme%stateVarMap%varPos%val(1)

  allocate( scheme%pdf( globalMaxLevels ))
  allocate( scheme%state( globalMaxLevels ))
  allocate( scheme%layout%stencil(1) )
  scheme%layout%stencil(1) = scheme%layout%fStencil
  ! initialize level descriptor
  call tem_init_elemLevels(                             &
    &                me        = scheme%levelDesc,      &
    &                boundary  = geometry%boundary,     &
    &                tree      = geometry%tree,         &
    &                stencils  = scheme%layout%stencil  )

  ! Create level Descriptor
  ! the neigh array is created using the LD and communication buffers 
  ! are filled up. 
  call tem_find_allElements( tree            = geometry%tree,              &
    &                        levelDesc       = scheme%levelDesc,           &
    &                        levelPointer    = geometry%levelPointer,      &
    &                        computeStencil  = scheme%layout%stencil,      &
    &                        commPattern     = params%general%commPattern, &
    &                        cleanup         = .true.,                     &
    &                        reqNesting      = params%nNesting,            &
    &                        proc            = params%general%proc         )

  call tem_build_horizontalDependencies(            &
    &       iStencil       = 1,                     &
    &       levelDesc      = scheme%levelDesc,      &
    &       tree           = geometry%tree,         &
    &       computeStencil = scheme%layout%fStencil )
  call tem_cleanup_arrays( levelDesc = scheme%levelDesc )

  write(logUnit(1),*) 'Initializing state'
  scheme%pdf(level)%nElems_local = geometry%tree%nElems
  scheme%pdf(level)%nSize = geometry%tree%nElems
  ! allocate state array
  allocate( scheme%state( level )%val(         &
    & scheme%pdf( level )%nSize * nScalars, 1 ))
  ! allocate the connectivity array
  allocate( scheme%pdf( level )%neigh( &
    & scheme%pdf( level )%nSize * QQ ) )
  scheme%pdf(level)%nNow = 1
  scheme%pdf(level)%nNext = 1

  do iLevel = geometry%tree%global%minLevel, geometry%tree%global%maxLevel
    ! construct connectivity vector pdf( iLevel )%neigh including
    ! bounce back rules
    call mus_construct_connectivity(                   &
      & neigh       = scheme%pdf(iLevel)%neigh,        &
      & nSize       = scheme%pdf(iLevel)%nSize,        &
      & nElems      = scheme%pdf(iLevel)%nElems_local, &
      & levelDesc   = scheme%levelDesc(iLevel),        &
      & stencil     = scheme%layout%fStencil,          &
      & varSys      = scheme%varSys,                   &
      & stateVarMap = scheme%stateVarMap               )
  end do   

  ! fill pdf array
  call fill_state( scheme   = scheme,                        &
    &              minLevel = geometry%tree%global%minLevel, &
    &              maxLevel = geometry%tree%global%maxLevel  )

  write(*,*) 'Setup indices '
  do iElem = 1, nElems_track
    point(iElem,:) = tem_BaryOfId( tree   = geometry%tree,            &
      &                            treeID = geometry%tree             &
      &                                     %treeID(elemPos(iElem)) )
  end do

  call scheme%varSys%method%val(statePos)     &
    & %setup_indices( varSys = scheme%varSys, &
    &                 point  = point,         &
    &                 iLevel = level,         &
    &                 tree   = geometry%tree, &
    &                 nPnts  = nElems_track,  &
    &                 idx    = idx            )

  ! store reference result
  allocate(ref_state(nElems_track*QQ))
  nSize = scheme%pdf(level)%nSize
  do iElem = 1, nElems_track
    do iComp = 1, QQ
      varPos = scheme%varSys%method%val(statePos)%state_varPos(iComp)
      ref_state( (iElem-1)*QQ + varPos ) &
        & = scheme%state(level)%val( &
& ?SAVE?(iComp,1,elemPos(iElem),QQ,nScalars,nSize,scheme%pdf(Level)%neigh ), 1)
    end do
  end do

  allocate( res(nElems_track*QQ))
  allocate(diff(nElems_track*QQ))
  write(*,*) 'Checking get_element '
  call scheme%varSys%method%val(statePos)%get_element(                &
    &                                varSys  = scheme%varSys,         &
    &                                elemPos = elemPos,               &
    &                                time    = params%general         &
    &                                                %simControl%now, &
    &                                tree    = geometry%tree,         &
    &                                nElems  = nElems_track,          &
    &                                nDofs   = 1,                     &
    &                                res     = res                    )

  diff = res - ref_state
  if (any(diff > eps)) then
    res_correct(1) = .false.
    write(*,*) 'Failed access state for element'
    do iElem = 1, nElems_track
      do iComp = 1, QQ
        varPos = scheme%varSys%method%val(statePos)%state_varPos(iComp)
        write(logUnit(1), "(2F18.13)") ref_state( (iElem-1)*QQ + varPos ), &
                                     &       res( (iElem-1)*QQ + varPos )
      end do
    end do
  else
    res_correct(1) = .true.
  end if  

  write(*,*) 'Checking get_point '
  call scheme%varSys%method%val(statePos)%get_point(                 &
    &                                varSys = scheme%varSys,         &
    &                                point  = point,                 &
    &                                time   = params%general         &
    &                                               %simControl%now, &
    &                                tree   = geometry%tree,         &
    &                                nPnts  = nElems_track,          &
    &                                res    = res                    )

  diff = res - ref_state
  if (any(diff > eps)) then
    res_correct(2) = .false.
    write(*,*) 'Failed access state for point'
  else
    res_correct(2) = .true.
  end if  

  write(*,*) 'Checking get_valOfIndex '
  call scheme%varSys%method%val(statePos)%get_valOfIndex(            &
    &                                varSys = scheme%varSys,         &
    &                                time   = params%general         &
    &                                               %simControl%now, &
    &                                iLevel = level,                 &
    &                                idx    = idx,                   &
    &                                nVals  = nElems_track,          &
    &                                res    = res                    )

  diff = res - ref_state
  if (any(diff > eps)) then
    res_correct(3) = .false.
    write(*,*) 'Failed access state for getValOfIndex'
  else
    res_correct(3) = .true.
  end if  

  call tem_finalize(params%general)
  if (all(res_correct)) then
    write(*,*) 'PASSED'
  else
    write(*,*) 'FAILED'
  end if

contains

  ! Fill state array
  subroutine fill_state( scheme, minLevel, maxLevel )
    !---------------------------------------------------------------------------
    !> Current scheme
    type( mus_scheme_type ) :: scheme
    !> Minlevel and maxLevel
    integer, intent(in) :: minLevel
    integer, intent(in) :: maxLevel
    !---------------------------------------------------------------------------
    integer :: iLevel, iElem, iState, iDir, varPos, nScalars
    integer :: statePos, QQ, nSize
    character(len=120) :: buffer
    !---------------------------------------------------------------------------
    nScalars = scheme%varSys%nScalars
    QQ = scheme%layout%fStencil%QQ

    write(logUnit(1),*) ' STATE:'
    levelLoop: do iLevel = minLevel, maxLevel
      nSize = scheme%pdf(iLevel)%nSize
      write(logUnit(1),*) 'iLevel: ', iLevel
      write(logUnit(1),*) ' nSize: ', nSize
      elemLoop: do iElem = 1, nSize
        write(logUnit(1),*) 'iElem', iElem
        stateLoop: do iState = 1, scheme%varSys%nStateVars
          write(logUnit(1),*) 'iState', iState, &
            &                 trim(scheme%varSys%varname%val(iState))
          ! Loop over every link
          buffer = ''
          dirLoop: do iDir  = 1, scheme%varSys%method%val(iState)%nComponents
            varPos = scheme%varSys%method%val(                     &
              &            scheme%stateVarMap%varPos%val(iState) ) &
              &                  %state_varPos(iDir)
            !statePos = (iElem-1)*nScalars + varPos
            statePos = &
& ?SAVE?(iDir, 1, iElem, QQ, nScalars, nSize, scheme%pdf(iLevel)%neigh )
            scheme%state(iLevel)%val( statePos, 1) = statePos
            write(buffer, '(a,i4)') trim(buffer), statePos
          end do dirLoop
          write(logUnit(1),*) trim(buffer)
        end do stateLoop
      end do elemLoop
    end do levelLoop

  end subroutine fill_state


end program mus_varSys_accessStateVar_test
