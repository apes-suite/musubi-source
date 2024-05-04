! See copyright notice in the COPYRIGHT file.
! This utest program check the output of bgk D3Q19 optimized compute kernel with
! input value of stencil weights. Output should be the same as input.
program mus_nNwtn_CY_weights_test
  use iso_c_binding, only: c_loc

  use env_module,         only: rk, eps, labelLen, init_random_seed
  use tem_tools_module,   only: tem_horizontalSpacer
  use treelmesh_module,   only: treelmesh_type
  use tem_bc_prop_module, only: tem_bc_prop_type
  use tem_logging_module, only: logUnit
  use tem_general_module, only: tem_start, tem_finalize
  use tem_varSys_module,  only: tem_varSys_type, tem_varSys_op_type, &
    &                           tem_varSys_init,                     &
    &                           tem_varSys_append_stateVar,          &
    &                           tem_varSys_proc_point,               &
    &                           tem_varSys_proc_element,             &
    &                           tem_varSys_proc_setParams,           &
    &                           tem_varSys_proc_getParams,           &
    &                           tem_varSys_proc_setupIndices,        &
    &                           tem_varSys_proc_getValOfIndex

  use mus_varSys_module,        only: mus_varSys_solverData_type
  use mus_variable_module,      only: mus_append_stateVar
  use mus_param_module,         only: mus_param_type
  use mus_scheme_type_module,   only: mus_scheme_type
  use mus_scheme_layout_module, only: mus_scheme_layout_type,    &
    &                                 mus_define_d3q19,          &
    &                                 mus_set_weights_d3q19
  use mus_field_prop_module,    only: mus_field_prop_type
  use mus_d3q19_module,         only: mus_advRel_kFluid_rBGK_vStd_lD3Q19
  use mus_physics_module,       only: mus_physics_type,          &
    &                                 mus_set_convFac,           &
    &                                 set_values_by_levels,      &
    &                                 mus_set_scaleFac,          &
    &                                 mus_physics_dump2outUnit
  use mus_nonNewtonian_module,  only: carreauYasuda, calcVisc_CY
  use mus_scheme_derived_quantities_module, only: mus_assign_derived_functions_ptr

  use mus_utestEnv_module,      only: init_fluid, load_env, init_varSys

  implicit none

  character(len=labelLen), parameter :: scaling = 'diffusive'
  integer, parameter :: scaleFactor = 4

  logical :: error = .true.
  real(kind=rk) :: tolerance, max_error
  integer :: iDir

  integer :: QQ
  real(kind=rk) :: omega, auxField(4)
  real(kind=rk), allocatable :: instate(:), outState(:), diff(:)
  integer,       allocatable :: neigh(:)
  integer :: level  = 1
  character(len=labelLen) :: sysName
  integer :: clock

  ! dummy variables
  type( mus_scheme_type ), target :: scheme
  type( mus_param_type ), target     :: param
  type( mus_varSys_solverData_type ), target    :: solverData
  type( treelmesh_type )   :: tree
  type( tem_bc_prop_type ) :: boundary

  call load_env( tree, boundary, param%general )

  tolerance = eps * 2500._rk
  write(logUnit(1), "('tolerance = ', ES11.4 )" ) tolerance

  ! define scheme header
  scheme%header%kind = 'fluid'
  scheme%header%relaxation = 'bgk'
  scheme%header%layout = 'd3q19'
  scheme%header%relaxHeader%variant = 'standard'

  ! generate random omega
  CALL SYSTEM_CLOCK( COUNT = clock )
  call init_random_seed( clock )
  call random_number( omega )
  omega = omega * 2.0_rk ! omega range [0,2)
  write( logUnit(1), "('omega = ', F5.3)") omega

  ! fill param
  param%scaling     = scaling
  param%scaleFactor = scaleFactor

  scheme%nFields = 1
  allocate( scheme%field(1) )

  ! fill fieldProp
  write( logUnit(1), *) 'fill fieldProp'
  call init_fluid( scheme%field(1)%fieldProp%fluid, omega, 1, scheme%header )

  ! fill nonNewtonian model parameters
  scheme%field(1)%label = ''
  scheme%field(1)%fieldProp%fluid%nNwtn%active = .true.
  scheme%field(1)%fieldProp%fluid%nNwtn%label = 'carreau_yasuda'
  scheme%field(1)%fieldProp%fluid%nNwtn%model  = carreauYasuda

  ! init layout
  write( logUnit(1), *) 'fill layout'
  call mus_define_d3q19( layout = scheme%layout, nElems = 1 )
  QQ = scheme%layout%fStencil%QQ

  ! set up dx, dt, physics table
  write( logUnit(1), *) 'set physics table'
  call set_physics( param%physics, omega, level, level )

  ! init var system
  write( logUnit(1), *) 'fill variable system'
  sysName = 'var_system'
  solverData%scheme    => scheme
  solverData%physics   => param%physics

  ! Initialize quantities pointers
  solverData%scheme%layout%quantities = mus_assign_derived_functions_ptr( &
    & label_stencil = scheme%header%layout,  &
    & label_fluid = scheme%header%kind       )

  call init_varSys( scheme%varSys, sysName, QQ, solverData, scheme%header, &
    &               scheme%derVarPos )

  ! test pdf with standard weights values
  call tem_horizontalSpacer(fUnit = logUnit(1))
  allocate( inState( scheme%varSys%nScalars ) )
  allocate(outState( scheme%varSys%nScalars ) )
  allocate(   neigh( scheme%varSys%nScalars ) )
  allocate(    diff( scheme%varSys%nScalars ) )
  write( logUnit(1), *) 'Input:  inState with standard weights'
  write( logUnit(1), *) 'Expect: outState be the same as inState'
  call mus_set_weights_d3q19( inState(1:QQ) )
  Outstate = -1.0_rk

  ! fill neigh array
  write( logUnit(1), *) 'fill neigh array'
  neigh = [ (iDir, iDir = 1, scheme%varSys%nScalars) ]
  allocate(scheme%pdf(level))
  scheme%pdf(level)%nElems_fluid = 1

  ! initialize auxField
  ! density
  auxField(1) = sum(inState(1:QQ))
  ! velocity
  auxField(2) = sum(inState(1:QQ) * scheme%layout%fStencil%cxDirRK(1,:))
  auxField(3) = sum(inState(1:QQ) * scheme%layout%fStencil%cxDirRK(2,:))
  auxField(4) = sum(inState(1:QQ) * scheme%layout%fStencil%cxDirRK(3,:))

  ! Calculate non-Newtonian viscosity
  call calcVisc_CY(                                                            &
    & nNwtn    = scheme%field(1)%fieldProp%fluid%nNwtn,                        &
    & viscKine = scheme%field(1)%fieldProp%fluid%viscKine%dataOnLvl(1)%val(:), &
    & omega    = scheme%field(1)%fieldProp%fluid%viscKine%omLvl(1)%val(:),     &
    & state    = instate,                                                      &
    & neigh    = neigh,                                                        &
    & densPos  = scheme%varSys%method%val(scheme%derVarPos(1)%density)         &
    &                  %auxField_varPos(1),                                    &
    & velPos   = scheme%varSys%method%val(scheme%derVarPos(1)%velocity)        &
    &                  %auxField_varPos(:),                                    &
    & auxField = auxField,                                                     &
    & nSize    = 1,                                                            &
    & nSolve   = 1,                                                            &
    & nScalars = scheme%varSys%nScalars,                                       &
    & nAuxScalars = scheme%varSys%nAuxScalars,                                 &
    & layout   = scheme%layout,                                                &
    & convFac  = param%physics%fac(1)                                          )

  ! call compute kernel
  write( logUnit(1), *) 'Calling compute kernel routine.'
  call mus_advRel_kFluid_rBGK_vStd_lD3Q19(      &
    &    fieldProp = scheme%field(:)%fieldProp, &
    &    inState   = inState,                   &
    &    outState  = outState,                  &
    &    auxField  = auxField,                  &
    &    neigh     = neigh,                     &
    &    nElems    = 1,                         &
    &    nSolve    = 1,                         &
    &    level     = level,                     &
    &    layout    = scheme%layout,             &
    &    params    = param,                     &
    &    derVarPos = scheme%derVarPos,          &
    &    varSys    = scheme%varSys              )

  write( logUnit(1), *) 'Calculating errors.'
  diff = outState - inState
  max_error = maxval( diff(:) )
  if ( max_error > tolerance ) then
    write(logUnit(1), *) "Max error exceeds tolerance!"
    write(logUnit(1), *) "iDir   inState  outState       diff"
    do iDir = 1, QQ
      if ( diff(idir) > tolerance ) &
        & write(logUnit(1),'( I4, 2F10.3, ES11.4)') iDir, inState(iDir), outState(iDir), diff(iDir)
    end do
  else
    error = .false.
    write(logUnit(1), *) "Max error within tolerance!"
  end if

  ! deallocate( fieldProp )
  call tem_finalize(param%general)

  if (.not. error) then
    write(*,*) 'PASSED'
  else
    write(*,*) 'FAILED'
  end if

  !*****************************************************************************

contains

  subroutine set_physics( physics, omega, minLevel, maxLevel )
    type( mus_physics_type )  :: physics
    real(kind=rk), intent(in) :: omega
    integer, intent(in) :: minLevel
    integer, intent(in) :: maxLevel

    real(kind=rk) :: dx, nuLB

    dx = 1.0_rk
    nuLB = ( 1.0 / omega - 0.5 ) / 3.0

    physics%dx   = 1.0_rk
    physics%dt   = nuLB * dx * dx / ( 0.0035_rk / 1050._rk )
    physics%rho0 = 1050._rk

    ! Assign dx and dt for each level
    allocate(physics%dxLvl( minLevel:maxLevel ))
    allocate(physics%dtLvl( minLevel:maxLevel ))
    allocate(physics%fac(   minLevel:maxLevel ))

    physics%dxLvl( minLevel:maxLevel ) = &
      &       set_values_by_levels( physics%dx, minLevel, maxLevel, 2 )
    physics%dtLvl( minLevel:maxLevel ) = &
      &       set_values_by_levels( physics%dt, minLevel, maxLevel, scaleFactor)

    call mus_set_convFac( physics, minLevel, maxLevel )
    call mus_set_scaleFac( physics, minLevel, maxLevel )
    call mus_physics_dump2outUnit( physics, logUnit(1), minLevel, maxLevel )

  end subroutine set_physics

end program mus_nNwtn_CY_weights_test
!******************************************************************************
