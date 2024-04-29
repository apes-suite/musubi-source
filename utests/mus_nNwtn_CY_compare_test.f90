! See copyright notice in the COPYRIGHT file.
! This utest program compare the outputs of nonNewtonian CY optimized kernel and
! explicit kernel.
! input pdf is initialized by equilibrium function, where density and velocity
! is randomly generated.
program mus_nNwtn_CY_compare_test
  use iso_c_binding, only: c_loc

  use env_module,         only: rk, eps, labelLen, init_random_seed
  use tem_tools_module,   only: tem_horizontalSpacer
  use tem_logging_module, only: logUnit
  use treelmesh_module,   only: treelmesh_type
  use tem_bc_prop_module, only: tem_bc_prop_type
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
  use mus_scheme_layout_module, only: mus_define_d3q19
  use mus_field_prop_module,    only: mus_field_prop_type
  use mus_d3q19_module,         only: bgk_advRel_d3q19
  use mus_bgk_module,           only: bgk_advRel_generic
  use mus_physics_module,       only: mus_physics_type,         &
    &                                 mus_set_convFac,          &
    &                                 mus_set_scaleFac,         &
    &                                 set_values_by_levels,     &
    &                                 mus_physics_dump2outUnit
  use mus_nonNewtonian_module,  only: carreauYasuda, calcVisc_CY
  use mus_derivedQuantities_module2, only: getEquilibrium
  use mus_scheme_derived_quantities_module, only: mus_assign_derived_functions_ptr

  use mus_utestEnv_module, only: init_fluid, load_env, init_varSys

  implicit none

  integer, parameter :: QQ = 19   !> number of pdf directions
  integer, parameter :: nScalars = QQ
  real(kind=rk), parameter :: densPhy = 1050._rk
  real(kind=rk), parameter :: nuPhy   = 0.0035_rk / densPhy
  character(len=labelLen), parameter :: scaling = 'diffusive'
  integer, parameter :: scaleFactor = 2

  logical :: error = .true.
  real(kind=rk) :: tolerance, max_error
  integer :: iDir

  real(kind=rk) :: omega, auxField(4)
  real(kind=rk) :: instate(nScalars), diff(nScalars)
  real(kind=rk) :: outOP(nScalars), outEx(nScalars)
  real(kind=rk) :: density, velocity(3)
  integer :: neigh(nScalars)
  integer :: level  = 1
  character(len=labelLen) :: sysName
  integer :: clock ! used for generate random number

  ! dummy variables
  type( mus_scheme_type ), target :: scheme
  type( mus_param_type ), target     :: param
  type( mus_varSys_solverData_type ), target    :: solverData
  type( treelmesh_type )   :: tree
  type( tem_bc_prop_type ) :: boundary

  ! start program, set up error tolerance
  call load_env( tree, boundary, param%general )

  tolerance = eps * 10000._rk
  write(logUnit(1), "('tolerance = ', ES11.4 )" ) tolerance

  ! define scheme header
  scheme%header%kind = 'fluid'
  scheme%header%relaxation = 'bgk'
  scheme%header%layout = 'd3q19'

  ! generate random omega
  CALL SYSTEM_CLOCK( COUNT = clock )
  call init_random_seed( clock )
  call random_number( omega )
  omega = omega * 2.0_rk ! omega range [0,2)
  write( logUnit(1), "('omega = ', F12.7)") omega

  ! generate density within range [0.9, 1.1]
  call random_number( density )
  density = density * 0.2 + 0.9_rk
  write( logUnit(1), "('density = ', F5.3)" ) density

  ! generate velocity within range [-0.1, 0.1]
  call random_number( velocity(1) )
  call random_number( velocity(2) )
  call random_number( velocity(3) )
  velocity = velocity * 0.2_rk - 0.1_rk
  write( logUnit(1), "('velocity = ', 3F7.3) " ) velocity

  ! fill param
  param%scaling = scaling
  param%scaleFactor = scaleFactor

  scheme%nFields = 1
  allocate( scheme%field(1) )

  ! fill fieldProp
  write( logUnit(1), '(A)') 'fill fieldProp'
  call init_fluid( scheme%field(1)%fieldProp%fluid, omega, 1, scheme%header )
  ! fill nonNewtonian model parameters
  scheme%field(1)%label = ''
  scheme%field(1)%fieldProp%fluid%nNwtn%active = .true.
  scheme%field(1)%fieldProp%fluid%nNwtn%label = 'carreau_yasuda'
  scheme%field(1)%fieldProp%fluid%nNwtn%model  = carreauYasuda

  ! fill neigh array
  write( logUnit(1), '(A)') 'fill neigh array'
  neigh = [ (iDir, iDir = 1, nScalars) ]
  allocate(scheme%pdf(level))
  scheme%pdf(level)%nElems_fluid = 1

  ! init layout
  write( logUnit(1), '(A)') 'fill layout'
  call mus_define_d3q19( layout = scheme%layout, nElems = 1 )

  ! set up dx, dt, physics table, write to solSpec string
  write( logUnit(1), '(A)') 'set physics table'
  call set_physics( param%physics, omega, level, level )

  write( logUnit(1), '(A)') 'fill variable system'
  sysName = 'var_system'
  solverData%scheme    => scheme
  solverData%physics   => param%physics

  ! Initialize quantities pointers
  solverData%scheme%layout%quantities = mus_assign_derived_functions_ptr(       &
    &                                     label_stencil = scheme%header%layout, &
    &                                     scheme_kind = scheme%header%kind      )

  call init_varSys( scheme%varSys, sysName, QQ, solverData, scheme%header, &
    &               scheme%derVarPos )

  call tem_horizontalSpacer(fUnit = logUnit(1))
  write( logUnit(1), '(A)') 'Input:  inState initialized by equilibrium'
  write( logUnit(1), '(A)') 'Expect: outState from two kernels should be the same.'
  ! initialize pdf by calculating equilibirum
  inState(1:QQ) = getEquilibrium( density, velocity, scheme%layout )

  ! initialize auxField
  ! density
  auxField(1) = density 
  ! velocity
  auxField(2) = velocity(1) 
  auxField(3) = velocity(2) 
  auxField(4) = velocity(3) 

  ! Calculate non-Newtonian viscosity
  call calcVisc_CY(                                                            &
    & nNwtn    = scheme%field(1)%fieldProp%fluid%nNwtn,                        &
    & viscKine = scheme%field(1)%fieldProp%fluid%viscKine%dataOnLvl(1)%val(:), &
    & omega    = scheme%field(1)%fieldProp%fluid%viscKine%omLvl(1)%val(:),     &
    & state    = instate,                                                      &
    & neigh    = neigh,                                                        &
    & auxField = auxField,                                                     &
    & densPos  = scheme%varSys%method%val(scheme%derVarPos(1)%density)         &
    &                  %auxField_varPos(1),                                    &
    & velPos   = scheme%varSys%method%val(scheme%derVarPos(1)%velocity)        &
    &                  %auxField_varPos(:),                                    &
    & nSize    = 1,                                                            &
    & nSolve   = 1,                                                            &
    & nScalars = scheme%varSys%nScalars,                                       &
    & nAuxScalars = scheme%varSys%nAuxScalars,                                 &
    & layout   = scheme%layout,                                                &
    & convFac  = param%physics%fac(1)                                          )

  outOP = -1.0_rk
  ! call Optimized compute kernel
  write( logUnit(1), '(A)') 'Calling compute kernel routine.'
  call bgk_advRel_d3q19( fieldProp = scheme%field(:)%fieldProp, &
    &                    inState   = inState,                   &
    &                    outState  = outOP,                     &
    &                    auxField  = auxField,                  &
    &                    neigh     = neigh,                     &
    &                    nElems    = 1,                         &
    &                    nSolve    = 1,                         &
    &                    level     = level,                     &
    &                    layout    = scheme%layout,             &
    &                    params    = param,                     &
    &                    derVarPos = scheme%derVarPos,          &
    &                    varSys    = scheme%varSys              )

  ! call Explicit compute kernel
  outEx = -1.0_rk
  write( logUnit(1), '(A)') 'Calling compute kernel routine.'
  call bgk_advRel_generic( fieldProp = scheme%field(:)%fieldProp, &
    &                       inState   = inState,                  &
    &                       outState  = outEx,                    &
    &                       auxField  = auxField,                 &
    &                       neigh     = neigh,                    &
    &                       nElems    = 1,                        &
    &                       nSolve    = 1,                        &
    &                       level     = level,                    &
    &                       layout    = scheme%layout,            &
    &                       params    = param,                    &
    &                       derVarPos = scheme%derVarPos,         &
    &                       varSys    = scheme%varSys             )

  write( logUnit(1), '(A)') 'Calculating errors.'
  diff = outOp - outEx
  max_error = maxval( abs(diff(:)) )
  if ( max_error > tolerance ) then
    write(logUnit(1), *) "Max error exceeds tolerance: "
    write(logUnit(1), *) " Max error = ", max_error
    write(logUnit(1), *) " Tolerance = ", tolerance
    write(logUnit(1), '(A4,4A10)') "iDir", "inState", "outOp", "outEx", "diff"
    do iDir = 1, nScalars
      write(logUnit(1),*) iDir, inState(iDir), &
        &                             outOp(iDir), outEx(iDir), diff(iDir)
    end do
  else
    error = .false.
    write(logUnit(1), '(A)') "Max error within tolerance!"
    write(logUnit(1), *) " Max error = ", max_error
    write(logUnit(1), *) " Tolerance = ", tolerance
  end if

  call tem_finalize(param%general)

  if (.not. error) then
    write(*,'(A)') 'PASSED'
  else
    write(*,'(A)') 'FAILED'
  end if

  !*****************************************************************************

contains

  subroutine set_physics( physics, omega, minLevel, maxLevel )
    type( mus_physics_type ) :: physics
    real(kind=rk), intent(in) :: omega
    integer, intent(in) :: minLevel
    integer, intent(in) :: maxLevel

    real(kind=rk) :: dx, dt, nuLB

    dx = 1.0_rk
    nuLB = ( 1.0 / omega - 0.5 ) / 3.0
    dt = nuLB * dx * dx / nuPhy

    physics%dt = dt
    physics%dx = dx
    physics%rho0 = densPhy

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

end program mus_nNwtn_CY_compare_test
!******************************************************************************
