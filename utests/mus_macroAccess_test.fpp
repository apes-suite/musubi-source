! See copyright notice in the COPYRIGHT file.
!> This UTEST is to test lbm macros FETCH, SAVE, IDX and NGPOS for
!! PULL and Array of structures
! main program
?? include 'header/lbm_macros.inc'
program mus_macroAccess_test
  use iso_c_binding, only: c_loc

  use env_module,          only: labelLen, globalMaxLevels
  use tem_general_module,  only: tem_start, tem_finalize
  use tem_tools_module,    only: tem_horizontalSpacer
  use tem_logging_module,  only: logUnit
  use treelmesh_module,    only: treelmesh_type
  use tem_varMap_module,   only: tem_varMap_type
  use tem_grow_array_module,    only: init, append, truncate  
  use tem_varSys_module,   only: tem_varSys_type, tem_varSys_init, &
    &                            tem_varSys_append_stateVar,       &
    &                            tem_varSys_proc_point,            &
    &                            tem_varSys_proc_element,          &
    &                            tem_varSys_proc_setParams,        &
    &                            tem_varSys_proc_getParams,        &
    &                            tem_varSys_proc_setupIndices,     &
    &                            tem_varSys_proc_getValOfIndex             
  use tem_construction_module, only: tem_init_elemLevels, tem_find_allElements,&
    &                                tem_build_verticalDependencies,           &
    &                                tem_build_horizontalDependencies,         &
    &                                tem_cleanup_arrays

  use mus_scheme_type_module,   only: mus_scheme_type
  use mus_scheme_layout_module, only: mus_define_d3q19
  use mus_param_module,         only: mus_param_type
  use mus_varSys_module,        only: mus_varSys_solverData_type,  &
    &                                 mus_deriveVar_ForPoint
  use mus_stateVar_module,      only: mus_access_state_ForElement
  use mus_geom_module,          only: mus_geom_type
  use mus_connectivity_module,  only: mus_construct_connectivity
  use mus_utestEnv_module,      only: load_env

  use mus_macroAccess_module

  implicit none

  type( mus_scheme_type ), target   :: scheme
  type( mus_param_type ), target    :: params
  type( mus_geom_type ), target     :: geometry
  type( mus_varSys_solverData_type ), target   :: solverData

  character(len=labelLen) :: sysName
  integer :: level, minLevel, maxLevel, iState, iElem, iDir, iLevel
  ! integer :: neighPos, elemPos, pos
  integer :: nScalars
  character(len=120) :: buffer
  integer, allocatable :: fetchRef(:,:,:), fetchIndex(:,:,:)
  logical :: success

  write(*,*) 'UTEST: MACRO Access for state and neigh array'
  !write(*,*) 'IDX ', IDX(1, 2, 9)
  ! load utest mesh
  call load_env( tree     = geometry%tree,     &
    &            boundary = geometry%boundary, &
    &            general  = params%general     ) 

  level    = geometry%tree%global%minLevel
  minlevel = geometry%tree%global%minLevel
  maxLevel = geometry%tree%global%maxLevel 

  ! if level is not 1 and minlevel is not equal to maxLevel stop proceeding
  ! This program works only for 8 elements i.e octree with single level
  if (level /= 1 .or. minLevel /= maxLevel) STOP

  sysName = 'var_system'
  solverData%scheme    => scheme
  solverData%geometry => geometry

  ! set nFields
  scheme%nFields = 2

  call mus_define_d3q19( layout = scheme%layout, nElems = geometry%tree%nElems )

  ! initialize variable system
  call init_varSys( scheme%varSys, scheme%stateVarMap, sysName, solverData )

  allocate( scheme%pdf( globalMaxLevels ))
  allocate( scheme%state( globalMaxLevels ))
  allocate( scheme%layout%stencil(1) )
  scheme%layout%stencil(1) = scheme%layout%fStencil
  ! initialize level descriptor
  call tem_init_elemLevels(                             &
    &                me        = scheme%levelDesc, &
    &                boundary  = geometry%boundary,     &
    &                tree      = geometry%tree,         &
    &                stencils  = scheme%layout%stencil  )

   ! Create level Descriptor
   ! the neigh array is created using the LD and communication buffers are filled up.
   call tem_find_allElements( tree            = geometry%tree,              &
     &                        levelDesc       = scheme%levelDesc,      &
     &                        levelPointer    = geometry%levelPointer,     &
     &                        computeStencil  = scheme%layout%stencil,      &
     &                        commPattern     = params%general%commPattern, &
     &                        cleanup         = .true.,                     &
     &                        reqNesting      = params%nNesting,            &
     &                        proc            = params%general%proc         )

  call tem_build_horizontalDependencies(            &
    &       iStencil       = 1,                     &
    &       levelDesc      = scheme%levelDesc, &
    &       tree           = geometry%tree,         &
    &       computeStencil = scheme%layout%fStencil )
  call tem_cleanup_arrays( levelDesc = scheme%levelDesc )

  scheme%pdf(level)%nElems_local = geometry%tree%nElems
  scheme%pdf(level)%nSize = geometry%tree%nElems
  ! allocate state array
  allocate( scheme%state( level )%val(                       &
    & scheme%pdf( level )%nSize * scheme%varSys%nScalars, 1 ))
  ! allocate the connectivity array
  allocate( scheme%pdf( level )%neigh(                       &
    & scheme%pdf( level )%nSize * scheme%layout%fStencil%QQ ))

  do iLevel = minLevel, maxLevel
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
  call fill_state( scheme   = scheme,   &
    &              minLevel = minLevel, &
    &              maxLevel = maxLevel  )

  nScalars = scheme%varSys%nScalars
!  ! dump neigh array
!  write(logUnit(1),*) 
!  write(logUnit(1),*) 'NEIGH:'
!  do iElem = 1, scheme%pdf( level )%nElems_local
!    write(logUnit(1),*) 'iElem', iElem
!    ! Loop over every link
!    buffer = ''
!    do iDir  = 1, scheme%layout%fStencil%QQ
!      neighPos = NGPOS(iDir, iElem, scheme%layout%fStencil%QQ)
!      write(buffer, '(a,i4)') trim(buffer), scheme%pdf(level)%neigh( neighPos ) 
!    end do 
!    write(logUnit(1),*) trim(buffer)
!  end do
!
!  ! SAVE INDEX
!  write(logUnit(1),*) 
!  write(logUnit(1),*) 'PDF SAVE:'
!  do iElem = 1, scheme%pdf( level )%nElems_local
!    write(logUnit(1),*) 'iElem', iElem
!    do iState = 1, scheme%varSys%nStateVars
!      if ( scheme%layout%fStencil%QQ ==                 &
!        &  scheme%varSys%method%val(iState)%nComponents ) then
!        write(logUnit(1),*) 'iState', iState, trim(scheme%varSys%varname%val(iState))
!        ! Loop over every link
!        buffer = ''
!        !do iDir  = 1, scheme%varSys%method%val(iState)%nComponents
!        do iDir  = 1, scheme%layout%fStencil%QQ 
!          write(buffer, '(a,i4)') trim(buffer), &
! & SAVE(iDir, iState, iElem, scheme%layout%fStencil%QQ, nScalars, scheme%pdf(level)%neigh)
!        end do 
!        write(logUnit(1),*) trim(buffer)
!      end if  
!    end do 
!  end do
!
  ! FETCH INDEX
  write(logUnit(1),*) 
  write(logUnit(1),*) 'PDF FETCH Index:'
  do iElem = 1, scheme%pdf( level )%nElems_local
    write(logUnit(1),*) 'iElem', iElem
    do iState = 1, scheme%varSys%nStateVars
      if ( scheme%layout%fStencil%QQ ==                 &
        &  scheme%varSys%method%val(iState)%nComponents ) then
        write(logUnit(1),*) 'iState', iState, trim(scheme%varSys%varname%val(iState))
        ! Loop over every link
        buffer = ''
        !do iDir  = 1, scheme%varSys%method%val(iState)%nComponents
        do iDir  = 1, scheme%layout%fStencil%QQ 
?? if ( PUSH ) then
          write(buffer, '(a,i4)') trim(buffer), &
& SAVE(iDir, iState, iElem, scheme%layout%fStencil%QQ, nScalars, scheme%pdf(level)%nSize, scheme%pdf(level)%neigh)
?? else
          write(buffer, '(a,i4)') trim(buffer), &
& FETCH(iDir, iState, iElem, scheme%layout%fStencil%QQ, nScalars, scheme%pdf(level)%nSize, scheme%pdf(level)%neigh)
?? end if
        end do 
        write(logUnit(1),*) trim(buffer)
      end if  
    end do 
  end do
!
!  write(logUnit(1),*) 
!  write(logUnit(1),*) 'NEIGH ELEM:'
!  do iElem = 1, scheme%pdf( level )%nElems_local
!    write(logUnit(1),*) 'iElem', iElem
!    do iState = 1, scheme%varSys%nStateVars
!      if ( scheme%layout%fStencil%QQ ==                 &
!        &  scheme%varSys%method%val(iState)%nComponents ) then
!        write(logUnit(1),*) 'iState', iState, trim(scheme%varSys%varname%val(iState))
!        ! Loop over every link
!        buffer = ''
!        !do iDir  = 1, scheme%varSys%method%val(iState)%nComponents
!        do iDir  = 1, scheme%layout%fStencil%QQ 
!          write(buffer, '(a,i4)') trim(buffer), &
! & ELEMIDX(FETCH(iDir, iState, iElem, scheme%layout%fStencil%QQ, scheme%pdf(level)%neigh), &
! & scheme%varSys%nScalars)
!        end do 
!        write(logUnit(1),*) trim(buffer)
!      end if  
!    end do 
!  end do

  allocate(fetchRef(geometry%tree%nElems, scheme%nFields, scheme%layout%fStencil%QQ)) 
  allocate(fetchIndex(geometry%tree%nElems, scheme%nFields, scheme%layout%fStencil%QQ)) 

  ! neighbor reference index for cube with 8 elements for stencil d3q19 
  ! with 2 pdf fields and one omega variable in state array
  ! field 1
  fetchRef(1, 1, :) = (/  40,  80, 159,  43,  83, 162, 241, 242, 243, 244, 206, 207, 208, 209, 132, 133, 134, 135,  19 /)
  fetchRef(2, 1, :) = (/   1, 119, 198,   4, 122, 201, 280, 281, 282, 283, 167, 168, 169, 170,  93,  94,  95,  96,  58 /)
  fetchRef(3, 1, :) = (/ 118,   2, 237, 121,   5, 240, 163, 164, 165, 166, 284, 285, 286, 287,  54,  55,  56,  57,  97 /)
  fetchRef(4, 1, :) = (/  79,  41, 276,  82,  44, 279, 202, 203, 204, 205, 245, 246, 247, 248,  15,  16,  17,  18, 136 /)
  fetchRef(5, 1, :) = (/ 196, 236,   3, 199, 239,   6,  85,  86,  87,  88,  50,  51,  52,  53, 288, 289, 290, 291, 175 /)
  fetchRef(6, 1, :) = (/ 157, 275,  42, 160, 278,  45, 124, 125, 126, 127,  11,  12,  13,  14, 249, 250, 251, 252, 214 /)
  fetchRef(7, 1, :) = (/ 274, 158,  81, 277, 161,  84,   7,   8,   9,  10, 128, 129, 130, 131, 210, 211, 212, 213, 253 /)
  fetchRef(8, 1, :) = (/ 235, 197, 120, 238, 200, 123,  46,  47,  48,  49,  89,  90,  91,  92, 171, 172, 173, 174, 292 /)

  ! field 2
  fetchRef(1, 2, :) = (/  59,  99, 178,  62, 102, 181, 260, 261, 262, 263, 225, 226, 227, 228, 151, 152, 153, 154,  38 /)
  fetchRef(2, 2, :) = (/  20, 138, 217,  23, 141, 220, 299, 300, 301, 302, 186, 187, 188, 189, 112, 113, 114, 115,  77 /)
  fetchRef(3, 2, :) = (/ 137,  21, 256, 140,  24, 259, 182, 183, 184, 185, 303, 304, 305, 306,  73,  74,  75,  76, 116 /)
  fetchRef(4, 2, :) = (/  98,  60, 295, 101,  63, 298, 221, 222, 223, 224, 264, 265, 266, 267,  34,  35,  36,  37, 155 /)
  fetchRef(5, 2, :) = (/ 215, 255,  22, 218, 258,  25, 104, 105, 106, 107,  69,  70,  71,  72, 307, 308, 309, 310, 194 /)
  fetchRef(6, 2, :) = (/ 176, 294,  61, 179, 297,  64, 143, 144, 145, 146,  30,  31,  32,  33, 268, 269, 270, 271, 233 /)
  fetchRef(7, 2, :) = (/ 293, 177, 100, 296, 180, 103,  26,  27,  28,  29, 147, 148, 149, 150, 229, 230, 231, 232, 272 /)
  fetchRef(8, 2, :) = (/ 254, 216, 139, 257, 219, 142,  65,  66,  67,  68, 108, 109, 110, 111, 190, 191, 192, 193, 311 /)

  success = .true.
  ! check only neighbor for pdfs since omega does not require neighbor
  do iElem = 1, scheme%pdf( level )%nElems_local
    do iState = 1, scheme%varSys%nStateVars
      if ( scheme%layout%fStencil%QQ ==                 &
        &  scheme%varSys%method%val(iState)%nComponents ) then
        do iDir  = 1, scheme%layout%fStencil%QQ 
?? if (PUSH) then        
          fetchIndex(iElem, iState, iDir) = &
 & SAVE(iDir, iState, iElem, scheme%layout%fStencil%QQ, nScalars, scheme%pdf(level)%nElems_local, scheme%pdf(level)%neigh)
?? else
          fetchIndex(iElem, iState, iDir) = &
 & FETCH(iDir, iState, iElem, scheme%layout%fStencil%QQ, nScalars, scheme%pdf(level)%nElems_local, scheme%pdf(level)%neigh)
?? end if
        end do
        if (any(fetchIndex(iElem, iState,:) /= fetchRef(iElem,iState,:))) success = .false.
      end if  
    end do
  end do

  call tem_finalize(params%general)
  if (success) then
    write(*,*) 'PASSED'
  else  
    write(logUnit(1),*) 
    write(logUnit(1),*) 'Reference FETCH Index:'
    do iElem = 1, scheme%pdf( level )%nElems_local
      write(logUnit(1),*) 'iElem', iElem
      do iState = 1, scheme%varSys%nStateVars
        if ( scheme%layout%fStencil%QQ ==                 &
          &  scheme%varSys%method%val(iState)%nComponents ) then
          write(logUnit(1),*) 'iState', iState, trim(scheme%varSys%varname%val(iState))
          ! Loop over every link
          buffer = ''
          !do iDir  = 1, scheme%varSys%method%val(iState)%nComponents
          do iDir  = 1, scheme%layout%fStencil%QQ
            write(buffer, '(a,i4)') trim(buffer), fetchRef(iElem, iState, iDir)
          end do
          write(logUnit(1),*) trim(buffer)
        end if
      end do
    end do

    write(*,*) 'FAILED'
  end if

contains

  ! append pdf to variable system
  subroutine init_varSys( varSys, stateVarMap, sysName, solverData )
    type( tem_varSys_type ) :: varSys
    type(tem_varMap_type), intent(out) :: stateVarMap
    character(len=labelLen) :: sysName
    type(mus_varSys_solverData_type), target :: solverData

    integer :: pdfPos, omegaPos
    integer :: iField
    logical :: wasAdded
    character(len=2) :: fldID
    character(len=labelLen) :: buffer
    procedure(tem_varSys_proc_point), pointer :: get_point => NULL()
    procedure(tem_varSys_proc_element), pointer :: get_element => NULL()
    procedure(tem_varSys_proc_setParams), pointer :: set_params => null()
    procedure(tem_varSys_proc_getParams), pointer :: get_params => null()
    procedure(tem_varSys_proc_setupIndices), pointer :: &
      &                                      setup_indices => null()
    procedure(tem_varSys_proc_getValOfIndex), pointer :: &
      &                                       get_valOfIndex => null()

    call tem_varSys_init( varSys, sysName )
    ! initialize state varMap
    call init(stateVarMap%varName)
    call init(stateVarMap%varPos)
    get_element => mus_access_state_ForElement
    get_point => mus_deriveVar_ForPoint

    do iField = 1, solverData%scheme%nFields
      write(fldID,'(i2)') iField
      write(buffer,'(a)') 'pdf_'//trim(fldID)
      call tem_varSys_append_stateVar( me             = varSys,                &
        &                              varName        = buffer,                &
        &                              nComponents    =                        &
        &                                solverData%scheme%layout%fStencil%QQ, &
        &                              method_data    = c_loc(solverData),     &
        &                              get_point      = get_point,             &
        &                              get_element    = get_element,           &
        &                              set_params     = set_params,            &
        &                              get_params     = get_params,            &
        &                              setup_indices  = setup_indices,         &
        &                              get_valOfIndex = get_valOfIndex,        &
        &                              pos            = pdfPos,                &
        &                              wasAdded       = wasAdded               )
      call append( me = stateVarMap%varPos, val = pdfPos )
      call append( me = stateVarMap%varname, val = buffer )
    end do    

    write(buffer,'(a)') 'omega'
    call tem_varSys_append_stateVar( me             = varSys,            &
      &                              varName        = buffer,            &
      &                              nComponents    = 1,                 &
      &                              method_data    = c_loc(solverData), &
      &                              get_point      = get_point,         &
      &                              get_element    = get_element,       &
      &                              set_params     = set_params,        &
      &                              get_params     = get_params,        &
      &                              setup_indices  = setup_indices,     &
      &                              get_valOfIndex = get_valOfIndex,    &
      &                              pos            = omegaPos,          &
      &                              wasAdded       = wasAdded           )

    call append( me = stateVarMap%varPos, val = omegaPos )
    call append( me = stateVarMap%varname, val = buffer ) 

    stateVarMap%nScalars = varSys%nScalars
    call truncate(me = stateVarMap%varPos)
    call truncate(me = stateVarMap%varname)

  end subroutine init_varsys
 

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
    integer :: statePos
    character(len=120) :: buffer
    !---------------------------------------------------------------------------
    nScalars = scheme%varSys%nScalars

    write(logUnit(1),*) 'STATE:'
    levelLoop: do iLevel = minLevel, maxLevel
      write(logUnit(1),*) 'iLevel', iLevel
      elemLoop: do iElem = 1, scheme%pdf( iLevel )%nElems_local
        write(logUnit(1),*) 'iElem', iElem
        stateLoop: do iState = 1, scheme%varSys%nStateVars
          write(logUnit(1),*) 'iState', iState, trim(scheme%varSys%varname%val(iState))
          ! Loop over every link
          buffer = ''
          dirLoop: do iDir  = 1, scheme%varSys%method%val(iState)%nComponents
            varPos = scheme%varSys%method%val(                     &
              &            scheme%stateVarMap%varPos%val(iState) ) &
              &                  %state_varPos(iDir)
            statePos = IDX(varPos, iElem, nScalars, scheme%pdf( iLevel )%nElems_local)
            scheme%state(iLevel)%val( statePos, 1) = statePos
            write(buffer, '(a,i4)') trim(buffer), statePos
          end do dirLoop
          write(logUnit(1),*) trim(buffer)
        end do stateLoop
      end do elemLoop
    end do levelLoop

  end subroutine fill_state

end program mus_macroAccess_test
