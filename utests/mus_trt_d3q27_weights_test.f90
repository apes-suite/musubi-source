! See copyright notice in the COPYRIGHT file.
! This utest program check the output of bgk D3Q19 optimized compute kernel with
! input value of stencil weights. Output should be the same as input.
program mus_trt_d3q27_weights_test
  use iso_c_binding,      only: c_loc
  use env_module,         only: rk, eps, labelLen, init_random_seed
  use tem_tools_module,   only: tem_horizontalSpacer
  use tem_logging_module, only: logUnit
  use tem_general_module, only: tem_start, tem_finalize

  use mus_param_module,         only: mus_param_type
  use mus_scheme_type_module,   only: mus_scheme_type
  use mus_varSys_module,        only: mus_varSys_solverData_type
  use mus_scheme_layout_module, only: mus_scheme_layout_type,    &
    &                                 mus_define_d3q27,          &
    &                                 mus_set_weights_d3q27
  use mus_field_prop_module,    only: mus_field_prop_type
  use mus_d3q27_module,         only: trt_advRel_d3q27 !, bgk_d3q27

  use mus_utestEnv_module,      only: init_fluid, init_varSys

  implicit none

  integer, parameter :: QQ  = 27

  logical :: error
  real(kind=rk) :: tolerance, max_error
  integer :: iDir

  type( mus_field_prop_type ), allocatable :: fieldProp(:)
  real(kind=rk) :: omega, rho0
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

  call tem_start('TRT D3Q27 optimized kernel weights utest', 'utest', params%general)
  error = .true.
  tolerance = eps * 2500._rk
  write(*,*) 'tolerance', tolerance

  ! define scheme header
  scheme%header%kind = 'fluid'
  scheme%header%relaxation = 'trt'
  scheme%header%layout = 'd3q27'

  CALL SYSTEM_CLOCK( COUNT = clock )
  call init_random_seed( clock )
  omega = 1.9_rk
  write( logUnit(1), *) 'omega = ', omega
  rho0  = 1.0_rk
  allocate( fieldProp(1) )
  level = 1
  call init_fluid( fieldProp(1)%fluid, omega, 1, scheme%header )

  neigh = [ (iDir, iDir = 1, QQ) ]

  nStart = 1
  nEnd   = 1

  ! init layout
  call mus_define_d3q27( layout = layout, nElems = 1 )

  sysName = 'var_system'
  solverData%scheme => scheme
  solverData%physics => params%physics

  call init_varSys( scheme%varSys, sysName, QQ, solverData, scheme%header, &
    &               scheme%derVarPos )

  ! test pdf with standard weights values
  call tem_horizontalSpacer(fUnit = logUnit(1))
  write( logUnit(1), *) 'Input:  inState with standard weights'
  write( logUnit(1), *) 'Expect: outState be the same as inState'
  call mus_set_weights_d3q27( inState )
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
  call trt_advRel_d3q27( fieldProp = fieldProp,        &
    &                    inState   = inState,          &
    &                    outState  = outState,         &
    &                    auxField  = auxField,         &
    &                    neigh     = neigh,            &
    &                    nElems    = 1,                &
    &                    nSolve    = 1,                &
    &                    level     = level,            &
    &                    layout    = layout,           &
    &                    params    = params,           &
    &                    derVarPos = scheme%derVarPos, &
    &                    varSys    = scheme%varSys     )

  write( logUnit(1), *) 'Calculating errors.'
  diff = outState - inState
  max_error = maxval( abs(diff(:)) )
  if ( max_error > tolerance ) then
    write(logUnit(1), *) "Max error exceeds tolerance!"
    write(logUnit(1), "(3A18)") "inState", "outState", "diff"
    do iDir = 1, QQ
      if ( diff(idir) > tolerance ) &
        & write(logUnit(1),"(I2, 3F18.13)") iDir, inState(iDir), outState(iDir), diff(iDir)
    end do
  else
    error = .false.
    do iDir = 1, QQ
      write(logUnit(1),"(I2, 3F18.13)") iDir, inState(iDir), outState(iDir), diff(iDir)
    end do
    write(logUnit(1), *) "Max error within tolerance!"
  end if

  deallocate( fieldProp )
  call tem_finalize(params%general)

  if (.not. error) then
    write(*,*) 'PASSED'
  else
    write(*,*) 'FAILED'
  end if

  !*****************************************************************************

end program mus_trt_d3q27_weights_test
!******************************************************************************
