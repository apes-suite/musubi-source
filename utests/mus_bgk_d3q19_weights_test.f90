! See copyright notice in the COPYRIGHT file.
! This utest program check the output of bgk D3Q19 optimized compute kernel with
! input value of stencil weights. Output should be the same as input.
program mus_bgk_d3q19_weights_test
  use iso_c_binding,      only: c_loc
  use env_module,         only: rk, eps, labelLen, init_random_seed
  use tem_tools_module,   only: tem_horizontalSpacer
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

  use mus_param_module,         only: mus_param_type
  use mus_scheme_type_module,   only: mus_scheme_type
  use mus_varSys_module,        only: mus_varSys_solverData_type
  use mus_scheme_layout_module, only: mus_scheme_layout_type,    &
    &                                 mus_define_d3q19,          &
    &                                 mus_set_weights_d3q19
  use mus_field_prop_module,    only: mus_field_prop_type
  use mus_d3q19_module,         only: mus_advRel_kFluid_rBGK_vStd_lD3Q19

  use mus_utestEnv_module,      only: init_fluid, init_varSys

  implicit none

  integer, parameter :: QQ  = 19   !> number of pdf directions

  logical :: error
  real(kind=rk) :: tolerance, max_error
  integer :: iDir

  real(kind=rk) :: omega
  real(kind=rk) :: instate(QQ), outState(QQ), diff(QQ), auxField(4)
  integer :: neigh(QQ)
  integer :: nStart, nEnd
  type( mus_scheme_layout_type ) :: layout
  character(len=labelLen) :: sysName
  integer :: clock

  ! dummy variables
  integer :: level
  type( mus_param_type ), target :: params
  type( mus_scheme_type ), target :: scheme
  type(mus_varSys_solverData_type), target :: solverData

  call tem_start('BGK D3Q19 optimized kernel weights utest', params%general)
  error = .true.
  tolerance = eps * 2500._rk
  write(*,*) 'tolerance', tolerance

  ! define scheme header
  scheme%header%kind = 'fluid'
  scheme%header%relaxation = 'bgk'
  scheme%header%layout = 'd3q19'
  scheme%header%relaxHeader%variant = 'standard'

  CALL SYSTEM_CLOCK( COUNT = clock )
  call init_random_seed( clock )
  call random_number( omega )
  write( logUnit(1), "(A,F8.4)") 'omega = ', omega
  omega = omega * 2.0_rk ! omega range [0,2)
  allocate( scheme%field(1) )
  level = 1
  call init_fluid( scheme%field(1)%fieldProp%fluid, omega, 1, scheme%header )

  neigh = [ (iDir, iDir = 1, QQ) ]
  allocate(scheme%pdf(level))
  scheme%pdf(level)%nElems_fluid = 1

  nStart = 1
  nEnd   = 1

  ! init layout
  call mus_define_d3q19( layout = layout, nElems = 1 )

  sysName = 'var_system'
  solverData%scheme    => scheme
  solverData%physics   => params%physics

  call init_varSys( scheme%varSys, sysName, QQ, solverData, scheme%header, &
    &               scheme%derVarPos )

  ! test pdf with standard weights values
  call tem_horizontalSpacer(fUnit = logUnit(1))
  write( logUnit(1), *) 'Input:  inState with standard weights'
  write( logUnit(1), *) 'Expect: outState be the same as inState'
  call mus_set_weights_d3q19( inState )
  Outstate = -1.0_rk

  ! initialize auxField
  ! density
  auxField(1) = sum(inState)
  ! velocity
  auxField(2) = sum(inState * layout%fStencil%cxDirRK(1,:))
  auxField(3) = sum(inState * layout%fStencil%cxDirRK(2,:))
  auxField(4) = sum(inState * layout%fStencil%cxDirRK(3,:))

  ! call compute kernel
  write( logUnit(1), *) 'Calling compute kernel routine.'
!$omp single
  call mus_advRel_kFluid_rBGK_vStd_lD3Q19(      &
    &    fieldProp = scheme%field(:)%fieldProp, &
    &    inState   = inState,                   &
    &    outState  = outState,                  &
    &    auxField  = auxField,                  &
    &    neigh     = neigh,                     &
    &    nElems    = 1,                         &
    &    nSolve    = 1,                         &
    &    level     = level,                     &
    &    layout    = layout,                    &
    &    params    = params,                    &
    &    derVarPos = scheme%derVarPos,          &
    &    varSys    = scheme%varSys              )
!$omp end single

  write( logUnit(1), *) 'Calculating errors.'
  diff = outState - inState
  max_error = maxval( diff(:) )
  if ( max_error > tolerance ) then
    write(logUnit(1), *) "Max error exceeds tolerance!"
    write(logUnit(1), *) "  inState                 outState               diff"
    do iDir = 1, QQ
      if ( diff(idir) > tolerance ) &
        & write(logUnit(1),*) iDir, inState(iDir), outState(iDir), diff(iDir)
    end do
  else
    error = .false.
    write(logUnit(1), *) "Max error within tolerance!"
  end if

  call tem_finalize(params%general)

  if (.not. error) then
    write(*,*) 'PASSED'
  else
    write(*,*) 'FAILED'
  end if

end program mus_bgk_d3q19_weights_test
!******************************************************************************
