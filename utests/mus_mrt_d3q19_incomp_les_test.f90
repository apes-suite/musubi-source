! See copyright notice in the COPYRIGHT file.
! This utest program checks if the bgk d3Q19 optimized(Op) and explicit compute
! kernel output the same results. Input pdf is initilized by equilibrium with
! randomly generated density and velocity values.
program mus_mrt_d3q19_incomp_les_test
  use iso_c_binding, only: c_loc

  use env_module,         only: rk, labelLen, init_random_seed
  use tem_param_module,   only: rho0
  use tem_tools_module,   only: tem_horizontalSpacer
  use tem_logging_module, only: logUnit
  use tem_general_module, only: tem_start, tem_finalize

  use mus_param_module,              only: mus_param_type
  use mus_scheme_type_module,        only: mus_scheme_type
  use mus_scheme_layout_module,      only: mus_scheme_layout_type, &
    &                                      mus_define_d3q19
  use mus_varSys_module,             only: mus_varSys_solverData_type
  use mus_field_prop_module,         only: mus_field_prop_type
  use mus_mrt_d3q19_module,          only: mus_advRel_kFluidIncomp_rMRT_vStd_lD3Q19, &
    &                                      mus_advRel_kFluidIncomp_rMRT_vStdNoOpt_lD3Q19
  use mus_derivedQuantities_module2, only: getEquilibriumIncomp
  use mus_relaxationParam_module,    only: mus_calcOmegaFromVisc
  use mus_Smagorinsky_module,        only: &
    & mus_turbVisc_Smagorinsky_fromPreColPDF_incomp_3D

  use mus_utestEnv_module,           only: init_fluid, init_varSys

  implicit none

  integer,  parameter :: QQ  = 19   !> number of pdf directions

  logical :: error
  real(kind=rk) :: tolerance, max_error
  integer :: iDir

  real(kind=rk) :: omega, velocity(3), density, auxField(4)
  real(kind=rk) :: instate(QQ), outEx(QQ), outOp(QQ), diff(QQ)
  integer :: neigh(QQ)
  integer :: nStart, nEnd
  type( mus_scheme_layout_type ) :: layout
  character(len=labelLen) :: sysName
  integer :: clock

  ! dummy variables
  integer :: level
  type( mus_param_type ), target               :: params
  type( mus_scheme_type ), target              :: scheme
  type( mus_varSys_solverData_type ), target   :: solverData
  real(kind=rk) :: tot_visc, turbVisc(1) 

  call tem_start('MRT D3Q19 incomp les kernels comparison utest', params%general)
  error = .true.
  tolerance = 1.0E-11
  write(*,*) 'tolerance = ', tolerance

  ! define scheme header
  scheme%header%kind = 'fluid_incompressible'
  scheme%header%relaxation = 'mrt'
  scheme%header%layout = 'd3q19'
  scheme%header%relaxHeader%variant = 'standard'

  ! generate random omega
  CALL SYSTEM_CLOCK( COUNT = clock )
  call init_random_seed( clock )
  call random_number( omega )
  omega = omega + 1.0_rk ! omega range [1,2)
  write( logUnit(1), "(' omega = ', F5.3)" ) omega

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

  neigh = [ (iDir, iDir = 1, QQ) ]

  nStart = 1
  nEnd   = 1
  level = 1

  ! init layout
  call mus_define_d3q19( layout = layout, nElems = 1 )

  sysName = 'var_system'
  solverData%scheme    => scheme
  solverData%physics   => params%physics

  call init_varSys( scheme%varSys, sysName, QQ, solverData, scheme%header, &
    &               scheme%derVarPos )

  ! fill fieldProp with omega
  write(*,*) 'Fill fluid property: omega, s_mrt ...'

  allocate( scheme%field(1) )
  call init_fluid( scheme%field(1)%fieldProp%fluid, omega, 1, scheme%header )
  ! define turbulence parameters
  scheme%field(1)%fieldProp%fluid%turbulence%config%coeff%C_s = 0.17_rk

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

  ! compute turbulent viscosity
  call mus_turbVisc_Smagorinsky_fromPreColPDF_incomp_3D(                   &
    & turbVisc    = turbVisc,                                              &
    & turbConfig  = scheme%field(1)%fieldProp%fluid                        &
    &                     %turbulence%config,                              &
    & state       = inState,                                               &
    & neigh       = neigh,                                                 &
    & auxField    = auxField,                                              &
    & densPos     = scheme%varSys%method%val(scheme%derVarPos(1)%density)  &
    &                     %auxField_varPos(1),                             &
    & velPos      = scheme%varSys%method%val(scheme%derVarPos(1)%velocity) &
    &                     %auxField_varPos(:),                             &
    & nSize       = 1,                                                     &
    & nSolve      = 1,                                                     &
    & nScalars    = scheme%varSys%nScalars,                                &
    & nAuxScalars = 4,                                                     &
    & layout      = layout,                                                &
    & dxL         = 1.0_rk,                                                &
    & dtL         = 1.0_rk,                                                &
    & viscKine    = scheme%field(1)%fieldProp%fluid%viscKine               &
    &                              %dataOnLvl(1)%val(:)                    )

  ! compute total viscosity and update omega
  tot_visc = scheme%field(1)%fieldProp%fluid%viscKine%dataOnLvl(1)%val(1) &
    &      + turbVisc(1)
  scheme%field(1)%fieldProp%fluid%viscKine%omLvl(1)%val(1) &
    & = mus_calcOmegaFromVisc(tot_visc)
  

  ! call optimized compute kernel
  outOP = -1.0_rk
  write( logUnit(1), *) 'Calling Optimized compute kernel routine.'
  call mus_advRel_kFluidIncomp_rMRT_vStd_lD3Q19( &
    &    fieldProp = scheme%field(:)%fieldProp,  &
    &    inState   = inState,                    &
    &    auxField  = auxField,                   &
    &    outState  = outOP,                      &
    &    neigh     = neigh,                      &
    &    nElems    = 1,                          &
    &    nSolve    = 1,                          &
    &    level     = level,                      &
    &    layout    = layout,                     &
    &    params    = params,                     &
    &    derVarPos = scheme%derVarPos,           &
    &    varSys    = scheme%varSys               )

  ! call explicit compute kernel
  outEx = -1.0_rk
  write( logUnit(1), *) 'Calling Explicit compute kernel routine.'
  call mus_advRel_kFluidIncomp_rMRT_vStdNoOpt_lD3Q19( &
    &    fieldProp = scheme%field(:)%fieldProp,       &
    &    inState   = inState,                         &
    &    outState  = outEx,                           &
    &    auxField  = auxField,                        &
    &    neigh     = neigh,                           &
    &    nElems    = 1,                               &
    &    nSolve    = 1,                               &
    &    level     = level,                           &
    &    layout    = layout,                          &
    &    params    = params,                          &
    &    derVarPos = scheme%derVarPos,                &
    &    varSys    = scheme%varSys                    )

  write( logUnit(1), *) 'Calculating errors between results from two kernels.'
  diff = abs(outEx - outOp)
  max_error = maxval( diff(:) )
  if ( max_error > tolerance ) then
    error = .true.
    write(logUnit(1), *) "Max error exceeds tolerance!"
    write(logUnit(1), *) "  inState        outEx              outOp         diff"
    do iDir = 1, QQ
      if ( diff(idir) > tolerance ) &
        & write(logUnit(1),*) iDir, inState(iDir), outEx(iDir), outOp(iDir), &
        &                     diff(iDir)
    end do
  else
    error = .false.
    write(logUnit(1), *) "Max error within tolerance!"
    write(logUnit(1), *) "  inState        outEx              outOp         diff"
    do iDir = 1, QQ
      write(logUnit(1),*) iDir, inState(iDir), outEx(iDir), outOp(iDir), &
                          diff(iDir)
    end do
  end if

  call tem_finalize(params%general)

  if (.not. error) then
    write(*,*) 'PASSED'
  else
    write(*,*) 'FAILED'
  end if

end program mus_mrt_d3q19_incomp_les_test
!******************************************************************************
