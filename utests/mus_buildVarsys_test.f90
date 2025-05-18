     ! See copyright notice in the COPYRIGHT file.
     !> This program build musubi variable system mainly for cray compiler bug
program mus_buildVarsys_test
  use tem_general_module, only: tem_general_type, tem_start, tem_finalize
  use tem_logging_module, only: logUnit
  use tem_stringKeyValuePair_module, only: init
  use tem_varSys_module,  only: tem_varSys_init

  use mus_scheme_type_module, only: mus_scheme_type
  use mus_scheme_layout_module, only: mus_define_d3q19
  use mus_variable_module,    only: mus_build_varSys, mus_append_stateVar
  use mus_varSys_module, only: mus_varSys_solverData_type

  implicit none

  type( tem_general_type ) :: general
  type( mus_scheme_type ),target :: scheme
  type( mus_varSys_solverData_type ), target :: solverData

  call tem_start('Build varSys unit test', general)

  solverData%scheme => scheme
  scheme%nFields = 1
  allocate(scheme%field(1))
  allocate(scheme%field(1)%bc(0))
  scheme%field(1)%label = 'test_'
  call init(scheme%field(1)%source%varDict) 
  call init(scheme%globSrc%varDict) 

  scheme%header%kind       = 'fluid'
  scheme%header%layout     = 'd3q19'
  scheme%header%relaxation = 'bgk'

  call mus_define_d3q19( layout = scheme%layout, nElems = 1 )

  allocate(scheme%luaVar(0))

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

  ! initialize all musubi derived variables
  scheme%field(1)%fieldProp%fluid%turbulence%active = .false.
  call mus_build_varSys( varSys       = scheme%varSys,          &
    &                    solverData   = solverData,             &
    &                    schemeHeader = scheme%header,          &
    &                    stencil      = scheme%layout%fStencil, &
    &                    nFields      = scheme%nFields,         &
    &                    derVarPos    = scheme%derVarPos,       &
    &                    luaVar       = scheme%luaVar,          &
    &                    field        = scheme%field(:),        &
    &                    globSrc      = scheme%globSrc,         &
    &                    poss_srcVar  = scheme%poss_srcVar,     &
    &                    st_funList   = scheme%st_funList       )

  call tem_finalize(general)
  write(logUnit(1),*)'PASSED'

end program mus_buildVarsys_test
