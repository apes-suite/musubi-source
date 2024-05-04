! See copyright notice in the COPYRIGHT file.
! This utest program checks if the bgk d3Q19 optimized(Op) and explicit compute
! kernel output the same results. Input pdf is initilized by equilibrium with
! randomly generated density and velocity values.
program mus_bgk_d3q19_incomp_compare_test

  use env_module,         only: rk, eps, labelLen, init_random_seed
  use tem_tools_module,   only: tem_horizontalSpacer
  use tem_logging_module, only: logUnit
  use tem_general_module, only: tem_start, tem_finalize
  use tem_param_module,   only: rho0

  use mus_param_module,              only: mus_param_type
  use mus_scheme_type_module,        only: mus_scheme_type
  use mus_scheme_layout_module,      only: mus_scheme_layout_type, &
    &                                      mus_define_d3q19
  use mus_varSys_module,             only: mus_varSys_solverData_type
  use mus_field_prop_module,         only: mus_field_prop_type
  use mus_d3q19_module,              only: mus_advRel_kFluidIncomp_rBGK_vStd_lD3Q19
  use mus_bgk_module,                only: mus_advRel_kCFD_rBGK_vStdNoOpt_l
  use mus_derivedQuantities_module2, only: getEquilibriumIncomp

  use mus_utestEnv_module,           only: init_fluid, init_varSys
  use mus_scheme_derived_quantities_module, only: mus_assign_derived_functions_ptr

  implicit none

  integer,  parameter :: QQ  = 19   !> number of pdf directions

  logical :: error = .true.
  real(kind=rk) :: tolerance, max_error
  integer :: iDir

  real(kind=rk) :: omega, velocity(3), density, auxField(4)
  real(kind=rk) :: instate(QQ), outEx(QQ), outOp(QQ), diff(QQ)
  integer :: neigh(QQ)
  type( mus_scheme_layout_type ) :: layout
  character(len=labelLen) :: sysName
  integer :: clock

  ! dummy variables
  integer :: level
  type( mus_param_type ), target               :: params
  type( mus_scheme_type ), target :: scheme
  type( mus_varSys_solverData_type ), target              :: solverData

  call tem_start('BGK D3Q19 incompressible kernels comparison utest', &
    &            'utest', params%general)
  tolerance = eps * 2500._rk
  write(*,*) 'tolerance = ', tolerance

  ! define scheme header
  scheme%header%kind = 'fluid_incompressible'
  scheme%header%relaxation = 'bgk'
  scheme%header%layout = 'd3q19'
  scheme%header%relaxHeader%variant = 'standard'

  ! use system clock as seed for random number generator
  call system_clock( count = clock )
  call init_random_seed( clock )
  ! generate omega within range [0,2)
  call random_number( omega )
  write( logUnit(1), *) ' omega = ', omega

  ! generate density within range [0.9, 1.1]
  call random_number( density )
  density = density * 0.2 + 0.9_rk
  write( logUnit(1), *) ' density = ', density

  ! generate velocity within range [-0.1, 0.1]
  call random_number( velocity(1) )
  call random_number( velocity(2) )
  call random_number( velocity(3) )
  velocity = velocity * 0.2_rk - 0.1_rk
  write( logUnit(1), *) ' velocity = ', velocity

  ! fill fieldProp with omega and rho0
  allocate( scheme%field(1) )
  call init_fluid( scheme%field(1)%fieldProp%fluid, omega, 1, scheme%header )

  neigh = [ (iDir, iDir = 1, QQ) ]

  ! init layout
  call mus_define_d3q19( layout = layout, nElems = 1 )

  sysName = 'var_system'
  solverData%scheme    => scheme
  solverData%physics   => params%physics

  ! Initialize quantities pointers
  layout%quantities = mus_assign_derived_functions_ptr( &
    & label_stencil = scheme%header%layout,  &
    & label_fluid = scheme%header%kind       )

  call init_varSys( scheme%varSys, sysName, QQ, solverData, scheme%header, &
    &               scheme%derVarPos )

  level = 1

  call tem_horizontalSpacer(fUnit = logUnit(1))
  write( logUnit(1), *) 'Input:  inState initilized by Eq from random den and vel'
  write( logUnit(1), *) 'Expect: outState from optimized kernel and explicity'
  write( logUnit(1), *) '        kernel be the same as inState'
  ! initialize pdf by calculating equilibirum
  inState =  getEquilibriumIncomp( density, velocity, layout, rho0 )

  ! initialize auxField
  ! density
  auxField(1) = density 
  ! velocity
  auxField(2) = velocity(1)
  auxField(3) = velocity(2)
  auxField(4) = velocity(3)

  outOP = -1.0_rk
  ! call optimized compute kernel
  write( logUnit(1), *) 'Calling Optimized compute kernel routine.'
  call mus_advRel_kFluidIncomp_rBGK_vStd_lD3Q19( &
    &    fieldProp = scheme%field(:)%fieldProp,  &
    &    inState   = inState,                    &
    &    outState  = outOP,                      &
    &    auxField  = auxField,                   &
    &    neigh     = neigh,                      &
    &    nElems    = 1,                          &
    &    nSolve    = 1,                          &
    &    level     = level,                      &
    &    layout    = layout,                     &
    &    params    = params,                     &
    &    derVarPos = scheme%derVarPos,           &
    &    varSys    = scheme%varSys               )

  outEx = -1.0_rk
  ! call explicit compute kernel
  write( logUnit(1), *) 'Calling Explicit compute kernel routine.'
  call mus_advRel_kCFD_rBGK_vStdNoOpt_l(        &
    &    fieldProp = scheme%field(:)%fieldProp, &
    &    inState   = inState,                   &
    &    outState  = outEx,                     &
    &    auxField  = auxField,                  &
    &    neigh     = neigh,                     &
    &    nElems    = 1,                         &
    &    nSolve    = 1,                         &
    &    level     = level,                     &
    &    layout    = layout,                    &
    &    params    = params,                    &
    &    derVarPos = scheme%derVarPos,          &
    &    varSys    = scheme%varSys              )

  write( logUnit(1), *) 'Calculating errors between results from two kernels.'
  diff = abs(outEx - outOp)
  max_error = maxval( diff(:) )
  if ( max_error > tolerance ) then
    write(logUnit(1), *) "Max error exceeds tolerance!"
    write(logUnit(1), "(A2,4A22)") " ", "inState", "outEx", "outOp", "diff"
    do iDir = 1, QQ
      if ( diff(idir) > tolerance ) then
        write(logUnit(1),"(I2,4ES22.14)") iDir, inState(iDir), outEx(iDir), outOp(iDir), &
                            diff(iDir)
      end if
    end do
  else
    error = .false.
    write(logUnit(1), *) "Max error within tolerance!"
    write(logUnit(1), "(A2,4A22)") " ", "inState", "outEx", "outOp", "diff"
    do iDir = 1, QQ
      write(logUnit(1),"(I2,4E22.14)") iDir, inState(iDir), outEx(iDir), outOp(iDir), &
                          diff(iDir)
    end do
  end if

  call tem_finalize(params%general)

  if (.not. error) then
    write(*,*) 'PASSED'
  else
    write(*,*) 'FAILED'
  end if

  !*****************************************************************************

end program mus_bgk_d3q19_incomp_compare_test
!******************************************************************************
