! See copyright notice in the COPYRIGHT file.
program mus_intpAcoustic_test
  use aotus_module, only: flu_State, close_config

  use env_module,               only: rk, pathLen, eps, labelLen, newUnit
  use treelmesh_module,         only: treelmesh_type
  use tem_construction_module,  only: tem_levelDesc_type
  use tem_grow_array_module,    only: grw_intArray_type, init, append
  use tem_dyn_array_module,     only: dyn_intArray_type, init, append
  use tem_tools_module,         only: tem_horizontalSpacer
  use tem_logging_module,       only: logUnit
  use tem_general_module,       only: tem_general_type, tem_start,             &
    &                                 tem_load_general, tem_finalize
  use tem_param_module,         only: childPosition, cs2inv, cs2, rho0
  use tem_var_system_module,    only: tem_var_sys_init
  use tem_variable_module,      only: tem_variable_type, tem_variable_out,     &
    &                                 tem_variable_from, tem_variable_dump,    &
    &                                 assignment(=), dyn_varArray_type, init,  &
    &                                 append, PositionOfVal
  use tem_polynomial_module,    only: polynomial, set_divU, set_coefficients

  use mus_fluid_module,         only: set_omegas_acoustic
  use mus_field_prop_module,    only: mus_field_prop_type
  use mus_param_module,         only: mus_param_type, mus_load_param
  use mus_interpolate_header_module, only: interpolation_type,                 &
                                           determine_timeLayer
  use mus_interpolate_linear_module, only:                            &
          fillMyGhostsFromFiner_linearIncomp,                         &
          fillFinerGhostsFromMe_linearIncomp,                         &
          fillMyGhostsFromFiner_linear,                               &
          fillFinerGhostsFromMe_linear
  use mus_interpolate_quadratic_module, only:                         &
          fillFinerGhostsFromMe_quadratic2dIncomp,                    &
          fillFinerGhostsFromMe_quadratic3dIncomp,                    &
          fillFinerGhostsFromMe_quadratic2d,                          &
          fillFinerGhostsFromMe_quadratic3d,                          &
          mus_init_intpMatrices
  use mus_pdf_module,          only: pdf_data_type, pdf_global_type, &
  use mus_moments_module,      only: set_momentIndices, &
                                     mus_init_moments
  use mus_variable_module,     only: mus_init_varSys, mus_build_derVarList, &
    &                                mus_kill_derVarList, mus_init_variables
  use mus_scheme_module,       only: mus_init_stateVars, &
    &                                mus_kill_stateVars
  use mus_scheme_type_module,  only: mus_scheme_type
  use mus_scheme_layout_module,only: mus_scheme_layout_type,            &
    &                                mus_define_d3q19, mus_define_d2q9, &
    &                                mus_init_stencil, mus_finish_stencil
  use mus_derivedQuantities_module2, only: set_pdfAcoustic, &
    &                                      set_pdfDiffusive, &
    &                                      getShearRateTensor_acoustic, &
    &                                      getShearRateTensor_diffusive, &
    &                                      getEquilibriumIncomp, &
    &                                      getEquilibrium, &
    &                                      convPrePost, secondMom
  use mus_config_module,       only: mus_open_config
  use mus_physics_module,      only: mus_load_physics

  implicit none

  !****************************************************************************
  integer :: returnFromFiner, nUnit
  integer, parameter :: iField = 1
  integer, parameter :: iFromFiner   = 1
  integer, parameter :: iFromCoarser = 2
  logical :: error, incompressible
  character(len=pathLen) :: filename
  character(len=labelLen) :: intpLabel(2), label
  type( flu_State ), allocatable :: conf(:) !< flu state
  real(kind=rk) :: error_limit
  type(mus_param_type) :: params
  integer :: iRun

  intpLabel(iFromFiner) = 'fromFiner'
  intpLabel(iFromCoarser) = 'fromCoarser'
  params%scaling = 'acoustic'

  call tem_start('unit test', params%general)
  error_limit = eps * 250000._rk
  write(*,*) 'error limit', error_limit

  ! create a empty dummy lua file
  filename = 'test.lua'
  nUnit = newUnit()
  open( file = filename, unit = nUnit )
  write(nUnit,*) '--dummy lua file'
  close(nUnit)

  ! open musubi config file and solver specific lua functions as chunk
  call mus_open_config( conf = conf, filename = filename, &
    &                   proc = params%general%proc )

  ! load and initialize solver general information
  call tem_load_general( me = params%general, conf = conf(1))
  call mus_load_param( params = params, conf = conf(1))
  call mus_init_variables()

  error = .false.
  write(*,*) 'Running interpolation test for '//trim(params%scaling)//' scaling...'
  call tem_horizontalSpacer(fUnit = logUnit(1))

  call tem_horizontalSpacer( fUnit = logUnit(1), before = 2)
  write(*,*) 'INTERPOLATION FROM FINER...'
  write(*,*)

  ! Run for two cases: compressible and incompressible
  do iRun = 1, 2
    if ( iRun == 1 ) then
      incompressible = .false.
      label = 'compressible'
    else
      incompressible = .true.
      label = 'incompressible'
    end if

    call check_interpolation( params = params, errCode = returnFromFiner, &
      &                       incompressible = incompressible )

    write(*,*) 'Result for '//trim(label), returnFromFiner
    if( returnFromFiner .ne. 0 ) error = .true.
    call tem_horizontalSpacer(fUnit = logUnit(1))
    if ( error ) exit
  end do

  call tem_horizontalSpacer(fUnit = logUnit(1))
  call close_config( conf(1) )
  call tem_finalize(params%general)
  if (.not. error) then
    write(*,*) 'PASSED'
  else
    write(*,*) 'FAILED'
  end if
  !****************************************************************************

contains

  !****************************************************************************
  ! Check the fromFiner interpolation routines
  ! 2d linear
  ! 3d linear
  subroutine check_interpolation( params, errCode, incompressible )
    !---------------------------------------------------------------------------
    type( mus_param_type ), intent(inout) :: params
    integer, intent(out) :: errCode
    logical, intent(in ) :: incompressible
    !---------------------------------------------------------------------------
    integer :: minLevel, maxLevel
    type( mus_scheme_type ) :: scheme
    type( interpolation_type ) :: intp, refIntp
    type( mus_field_prop_type ), allocatable :: fieldProp(:)
    type( mus_scheme_layout_type ) :: refLayout
    type( pdf_data_type ), allocatable :: pdf(:)
    type( pdf_global_type ) :: glob
    integer :: iVal, nSourceElems, iElem, iPos
    type( grw_intArray_type )   :: ind
    type( treelmesh_type )    :: tree !< treelm mesh
    real(kind=rk), allocatable :: momentCoefficients(:,:)
    real(kind=rk), allocatable :: srcMom(:,:)
    real(kind=rk), allocatable :: srcPdf(:)
    real(kind=rk), allocatable :: tgtMom(:), sourceCoord(:,:)
    real(kind=rk), allocatable :: tgtPdf(:), refMom(:), targetEq(:)
    real(kind=rk) :: omSrc, omTgt, targetCoord(3), omega
    integer :: nDOFs, offset, iChild, iDimn, sourceLevel, targetLevel, iEntry
    real(kind=rk) :: error, Sxx(3,3), Sxx_utest(3,3)
    !> moment index(min and max) for pressure, velocity, strain rate
    integer :: iPress, iVelMin, iVelMax, iSxxMin, iSxxMax
    type( tem_variable_type ) :: tVar(0)
    integer :: order, iIntp, iMin, iMax, iLevel, iDir, QQ
    integer, allocatable :: children(:), varPos(:) !< Which child configurations to check
    integer :: refOrder !< reference order
    logical :: setDivergence0

    real(kind=rk) :: rho, vel(3), strain_res(6)
    real(kind=rk), allocatable :: fEq(:), fNeq(:)
    !---------------------------------------------------------------------------

    if( incompressible ) then
      scheme%header%kind = 'lbm_incomp'
      setDivergence0 = .false.
    else
      scheme%header%kind = 'lbm'
      setDivergence0 = .false.
    end if

    refOrder = 1
    errCode = 0
    minLevel = 3
    maxLevel = 4
    omega = 1.8_rk
    ! Set the time layer to 1, so the interpolation routine always uses
    refIntp%timeLayer = 1

    allocate( glob%levelDesc( 1:maxLevel ))
    allocate( fieldProp(1))

    ! set level, omega, viscosity
    tree%global%minLevel = minLevel
    tree%global%maxLevel = maxLevel
    tree%global%boundingCubeLength = 1._rk
    fieldProp(iField)%fluid%kine_viscosityLB = (1._rk/omega - 0.5_rk) / 3._rk

    ! Initialize the LB time and space steps (for multilevel)
    scheme%nFields = 1
    allocate( scheme%field(scheme%nFields ))
    scheme%field(1)%label = 'lbm'
    call mus_load_physics( me= params%physics, conf = conf(1),            &
      &                    tree = tree, scaleFactor = params%scaleFactor, &
      &                    dtRef = 0.007_rk, dxRef = 0.03_rk )


    ! set physical viscosity and omega on levels
    fieldProp(iField)%fluid%kine_viscosity = &
      &             fieldProp(iField)%fluid%kine_viscosityLB &
      &           * params%physics%fac(minLevel)%visc
    call set_omegas_acoustic( fieldProp(iField)%fluid, params%physics, &
    do iLevel = minlevel, maxlevel
      write(*,"(A,I0,A,F6.3)") 'level: ', iLevel, &
        &        ' omega: ', fieldProp(iField)%fluid%omLvl(iLevel)
    end do

    ! order = 1: Linear
    ! order = 2: Quadratic
    do order = 1, 2
      ! Set the reference polynomial order to be the same as the interpolation
      ! order. This has to yield results exact up to numerical accuracy
      refOrder = order
      if ( order == 1 ) then
        ! Test fromFiner and fromCoarser for linear
        iMin = iFromFiner
        iMax = iFromCoarser
      else ! order == 2
        ! For quadratic, test only the fromCoarser
        ! (fromFiner does not exist)
        iMin = iFromCoarser
        iMax = iFromCoarser
      end if

      do iIntp = iMin, iMax
        write(*,"(A)")    'Interpolation: '//trim(intpLabel( iIntp ))
        write(*,"(A,I0)") 'Order:         ', order

        ! set sourceLevel and targetLevel
        call set_intpLevel( minLevel = minLevel, maxLevel = maxLevel,  &
          &                 sourceLevel = sourceLevel, targetLevel = targetLevel, &
          &                 fromFiner = (iIntp==iFromFiner ))

        ! set omega on sourceLevel and targetLevel
        omSrc = fieldProp( iField )%fluid%omLvl( sourceLevel )
        omTgt = fieldProp( iField )%fluid%omLvl( targetLevel )

        if ( allocated( children )) deallocate(children)
        if ( iIntp == iFromFiner ) then
          allocate(children(1))
          children(1) = 1
        else
          allocate(children(8))
          children(:) = [(iChild, iChild = 1,8)]
        end if

        ! Test 2D and 3D
        do iDimn = 2,3
          do iChild = 1, 2**iDimn
            write(*,"(A, I0)")   'Dimension: ', iDimn
            write(*,"(A, I0)")   'Children:  ', iChild
            write(*,"(A)")       'intp: '//trim(intpLabel( iIntp ))
            write(*,"(A, F6.3)") 'omega on source level: ', omSrc

            call set_momentIndices( iDimn, iPress, iVelMin, iVelMax, iSxxMin, iSxxMax)

            intp = refIntp
            scheme%layout = refLayout
            call mus_init_stencil( stencil = scheme%layout%stencil, nStencils = 1 )
            tree%nElems = 1

            if (iDimn == 3) then ! 3D
              if( order == 1 ) then
                nDOfs = 4  ! 3D linear polynomial with 4 coefficients
              else
                nDofs = 10 ! 3D quadratic polynomial with 10 coefficients
              end if
              call mus_define_d3q19( layout = scheme%layout, nElems=tree%nElems)
            else ! 2D
              if( order == 1) then
                nDOfs = 3 ! 2D linear polynomial with 3 coefficients
              else
                nDofs = 6 ! 2D quadratic polynomial with 6 coefficients
              end if
              call mus_define_d2q9( layout = scheme%layout, nElems=tree%nElems)
            end if ! iDimn
            call mus_finish_stencil( scheme%layout )
            QQ = scheme%layout%stencil(1)%QQ

            if( iIntp == iFromFiner .or. order == 1 ) then
              nSourceElems = 2**iDimn
            else
              nSourceElems = QQ
            end if

            call tem_var_sys_init( scheme%varSys )

            call mus_build_derVarList( scheme%derVarList, scheme%header,  &
              &                        scheme%layout%stencil(1))
            call mus_init_stateVars( scheme, size(tVar) )

            ! initialize dependencies
            call init_deps( fromFiner = (iIntp == iFromFiner ),    &
              &             levelDesc = glob%levelDesc( targetLevel ), glob = glob,&
              &             nSourceElems = nSourceElems, layout = scheme%layout,  &
              &             ind = ind, order = order, sourceCoord = sourceCoord,  &
              &             minLevel = minLevel, maxLevel = maxLevel, &
              &             fieldProp = fieldProp, intp = intp )

            ! Initialize the polynomial coefficients for each moment
            call set_Coefficients( me = momentCoefficients, nDOFs = nDOFs, &
              &                    QQ = QQ, &
              &                    valMin = 0.3_rk, valMax = 1.6_rk )

            allocate( srcMom(QQ,QQ ))
            allocate( srcPdf(QQ ))
            allocate( refMom(QQ ))
            allocate( tgtMom(QQ ))
            allocate( tgtPdf(QQ ))
            allocate( targetEq(QQ ))
            allocate( varPos(QQ ))
            varPos(:) = [( iDir, iDir = 1,QQ )]

            write(*,*) "Initialize PDF array"
            call init_stateArrays( pdf = pdf, &
              &                    fromFiner = ( iIntp == iFromFiner ), &
              &                    levelDesc = glob%levelDesc, &
              &                    layout = scheme%layout,&
              &                    minLevel = minLevel,&
              &                    maxLevel = maxLevel,&
              &                    nSourceElems = nSourceElems )

            ! Set states for the source elements
            write(*,*) "Set srcPDF of source elements"
            do iElem = 1, nSourceElems
              write(*,'(A,I2,A,3F5.1)') 'sourceElem: ', iElem, ' crd: ', sourceCoord(:,iElem)

              ! Set each quantity in terms of the polynomial.
              ! They are PHYISCAL states which have to be converted to PDF
              do iDir = 1, QQ
                ! The order is the same as in Z-curve (sources are the children)
                srcMom( iDir, iElem ) = polynomial(                               &
                  &                               a = momentCoefficients(:,iDir), &
                  &                               coord = sourceCoord(:, iElem ), &
                  &                               order = refOrder, nDims = iDimn )
              end do
              write(*,*) 'srcMom physical rho: ', srcMom(1, iElem )
              write(*,*) 'srcMom physical vel: ', srcMom(iVelMin:iVelMax, iElem )
              write(*,*) 'srcMom physical Sxx: ', srcMom(iSxxMin:iSxxMax, iElem )
              write(*,*) ''

              ! Set the moment coefficients to satisfy div(u) = 0
              if( setDivergence0 ) then
                write(*,*) "Set Divergence"
                call set_divU( uOut = srcMom( iVelMax, iElem ), nDims = iDimn,  &
                  &            coeff = momentcoefficients(:, iVelMin:iVelMax ), &
                  &            nStresses = iSxxMax - iSxxMin + 1,               &
                  &         SxxCoeff = momentcoefficients(:, iSxxMin:iSxxMax ), &
                  &            SxxOut = srcMom( iSxxMin: iSxxMax, iElem ),      &
                  &            nDofs = nDofs,  coord = sourceCoord(:, iElem ), order = reforder )
              end if

              ! Now each quantity is being set in terms of the polynomial.
              ! Convert moment from physical unit to LB unit
              srcMom(              1, iElem) =   srcMom(1,iElem) * cs2inv &
                &                              / params%physics%fac(sourceLevel)%press
              srcMom(iVelMin:iVelMax, iElem) =   srcMom( iVelMin:iVelMax, iElem ) &
                &                              / params%physics%fac(sourceLevel)%vel
              srcMom(iSxxMin:iSxxMax, iElem) =   srcMom( iSxxMin:iSxxMax, iElem ) &
                &                              / params%physics%fac(sourceLevel)%strainRate
              write(*,*) 'srcMom LB rho: ', srcMom(1, iElem )
              write(*,*) 'srcMom LB vel: ', srcMom(iVelMin:iVelMax, iElem )
              write(*,*) 'srcMom LB Sxx: ', srcMom(iSxxMin:iSxxMax, iElem )
              write(*,*) ''

              ! fill Sxx tensor
              ! It has to be symmetric
              Sxx = 0._rk
              if ( iDimn == 3 ) then
                Sxx(1,1) = srcMom( iSxxMin+0, iElem )
                Sxx(2,2) = srcMom( iSxxMin+1, iElem )
                Sxx(3,3) = srcMom( iSxxMin+2, iElem )
                Sxx(1,2) = srcMom( iSxxMin+3, iElem )
                Sxx(2,3) = srcMom( iSxxMin+4, iElem )
                Sxx(1,3) = srcMom( iSxxMin+5, iElem )
                Sxx(2,1) = Sxx(1,2)
                Sxx(3,2) = Sxx(2,3)
                Sxx(3,1) = Sxx(1,3)
              else
                Sxx(1,1) = srcMom( iSxxMin+0, iElem )
                Sxx(2,2) = srcMom( iSxxMin+1, iElem )
                Sxx(1,2) = srcMom( iSxxMin+2, iElem )
                Sxx(2,1) = Sxx(1,2)
              end if

              ! set pdf of source element by acoustic scaling
              write(*,*) 'Convert moment to pdf by set_PDFAcoustic'
              ! srcPdf(:) = set_pdfDiffusive( layout = scheme%layout, &
              !   &                          omega = omSrc, &
              !   &                          rho0 = rho0, &
              !   &                          mom = srcMom(:,iElem) )
              srcPdf(:) = set_pdfAcoustic( layout = scheme%layout, &
                &                          omega = omSrc, &
                &                          rho0 = rho0, &
                &                          mom = srcMom(:,iElem),&
                &                          incompressible = incompressible )
              do iDir = 1, QQ
                write(*,"(A,I2,A,F20.15)") "iDir: ", iDir, "  srcPDF: ", srcPDF(iDir)
              end do
              write(*,*) ''

              offset = QQ * ( iElem-1 )
              pdf(sourceLevel)%state( offset+1:offset+QQ, pdf(sourceLevel)%nNext ) = srcPdf

              ! check strain rate calculation
              write(*,*) 'Calculate Sxx from srcPDF'
              write(*,*) 'Sxx init: ', Sxx(1,1), Sxx(2,2), Sxx(3,3),Sxx(1,2), Sxx(2,3), Sxx(1,3)

              Sxx = 0._rk
              if ( incompressible ) then
                Sxx = getShearRateTensor_acoustic( subset = srcPdf, &
                  &                                layout = scheme%layout,   &
                  &                                omega = omSrc, &
                  &                                rho0 = rho0 )
              else
                Sxx = getShearRateTensor_acoustic( subset = srcPdf, &
                  &                                layout = scheme%layout,   &
                  &                                omega = omSrc )
              end if
              write(*,*) 'in deriv: ', Sxx(1,1), Sxx(2,2), Sxx(3,3),Sxx(1,2), Sxx(2,3), Sxx(1,3)
Sxx_utest = 0._rk
strain_res = 0._rk
allocate( fEq(QQ) )
allocate( fNeq(QQ) )
fEq  = 0._rk
fNeq = 0._rk
rho = sum( srcPDF )
vel = 0._rk
do iDir = 1, QQ
  Vel(:) = Vel(:) + srcPDF(iDir) * real(scheme%layout%stencil(1)%cxDir(:,iDir),rk)
end do
if ( incompressible ) then
  vel = vel / rho0
  fEq(:)  = getEquilibriumIncomp( rho, vel, scheme%layout, rho0 )
  fNeq(:) = ( srcPDF(:) - fEq(:) ) * convPrePost( omSrc )
  strain_res(:) = secondMom( scheme%layout%stencil(1)%cxcx, fNeq, QQ )
  strain_res(:) =  -1.5_rk * omSrc / rho0 * strain_res(:)
else
  vel = vel / rho
  fEq(:)  = getEquilibrium( rho, vel, scheme%layout )
  fNeq(:) = ( srcPDF(:) - fEq(:) ) * convPrePost( omSrc )
  strain_res(:) = secondMom( scheme%layout%stencil(1)%cxcx, fNeq, QQ )
  strain_res(:) =  -1.5_rk * omSrc / rho * strain_res(:)
end if
if ( iDimn == 3 ) then
  Sxx_utest(1,1) = strain_res(1)
  Sxx_utest(2,2) = strain_res(2)
  Sxx_utest(3,3) = strain_res(3)
  Sxx_utest(1,2) = strain_res(4)
  Sxx_utest(2,3) = strain_res(5)
  Sxx_utest(1,3) = strain_res(6)
else
  Sxx_utest(1,1) = strain_res(1)
  Sxx_utest(2,2) = strain_res(2)
  Sxx_utest(1,2) = strain_res(3)
end if
deallocate( fEq )
deallocate( fNeq)
              write(*,*) 'in utest: ', Sxx_utest(1,1), Sxx_utest(2,2), Sxx_utest(3,3),Sxx_utest(1,2), Sxx_utest(2,3), Sxx_utest(1,3)
              write(*,*) ''

              ! Sxx = getShearRateTensor_diffusive( subset = srcPdf, &
              !   &                                layout = scheme%layout,   &
              !   &                                omega = omSrc )
              ! write(*,*) 'Sxx Diffusiv: ', Sxx(1,1), Sxx(2,2), Sxx(3,3),Sxx(1,2), Sxx(2,3), Sxx(1,3)
              ! write(*,*) ''


              !write(*,*) 'Sxx calc: ', Sxx(1,1), Sxx(2,2), Sxx(3,3),Sxx(1,2), Sxx(2,3), Sxx(1,3)
              ! write(*,*) ''
            end do ! iElem

            !---------------------------------------------------------------------------
            ! Do the interpolation
            if( iIntp == iFromFiner ) then
              targetCoord(:) = [ 0.5_rk, 0.5_rk, 0.5_rk  ]
              if( incompressible ) then
                call fillMyGhostsFromFiner_linearIncomp(                       &
                  &                          pdf       = pdf,                  &
                  &                          scheme    = scheme,               &
                  &                          glob      = glob,                 &
                  &                          fieldProp = fieldProp,            &
                  &                          level     = targetLevel,          &
                  &                          layout    = scheme%layout,        &
                  &                          intp      = intp,                 &
                  &                          ind       = ind,                  &
                  &                          varSys    = scheme%varSys,        &
                  &                          physics   = params%physics,       &
                  &                          derVarPos = scheme%derVarPos  )
              else
                call fillMyGhostsFromFiner_linear(                             &
                  &                          pdf       = pdf,                  &
                  &                          scheme    = scheme,               &
                  &                          glob      = glob,                 &
                  &                          fieldProp = fieldProp,            &
                  &                          level     = targetLevel,          &
                  &                          layout    = scheme%layout,        &
                  &                          intp      = intp,                 &
                  &                          ind       = ind,                  &
                  &                          varSys    = scheme%varSys,        &
                  &                          physics   = params%physics,       &
                  &                          derVarPos = scheme%derVarPos  )
              end if
            else
              if( order == 1 ) then
                targetCoord(:) = [ 0.25_rk, 0.25_rk, 0.25_rk  ]
                if( incompressible ) then
                  call fillFinerGhostsFromMe_linearIncomp(                     &
                    &                          pdf       = pdf,                &
                    &                          scheme    = scheme,             &
                    &                          glob      = glob,               &
                    &                          fieldProp = fieldProp,          &
                    &                          level     = sourceLevel,        &
                    &                          layout    = scheme%layout,      &
                    &                          intp      = intp,               &
                    &                          ind       = ind,                &
                    &                          varSys    = scheme%varSys,      &
                    &                          physics   = params%physics,     &
                    &                          derVarPos = scheme%derVarPos  )
                else
                  call fillFinerGhostsFromMe_linear(                           &
                    &                          pdf       = pdf,                &
                    &                          scheme    = scheme,             &
                    &                          glob      = glob,               &
                    &                          fieldProp = fieldProp,          &
                    &                          level     = sourceLevel,        &
                    &                          layout    = scheme%layout,      &
                    &                          intp      = intp,               &
                    &                          ind       = ind,                &
                    &                          varSys    = scheme%varSys,      &
                    &                          physics   = params%physics,     &
                    &                          derVarPos = scheme%derVarPos  )
                end if
              else
                targetCoord(:) = real( childPosition(iChild,:), rk) * 0.25_rk
                write(*,*) 'Testing iChild: ', iChild, 'targetCoord', real( targetCoord )
                glob%levelDesc( targetLevel )%depFromCoarser(1)%coord = targetCoord
                if( iDimn == 2 ) then
                  if( incompressible ) then
                    call fillFinerGhostsFromMe_quadratic2dincomp(              &
                      &                          pdf       = pdf,              &
                      &                          scheme    = scheme,           &
                      &                          glob      = glob,             &
                      &                          fieldProp = fieldProp,        &
                      &                          level     = sourceLevel,      &
                      &                          layout    = scheme%layout,    &
                      &                          intp      = intp,             &
                      &                          ind       = ind,              &
                      &                          varSys    = scheme%varSys,    &
                      &                          physics   = params%physics,   &
                      &                          derVarPos = scheme%derVarPos  )
                  else
                    call fillFinerGhostsFromMe_quadratic2d(                    &
                      &                          pdf       = pdf,              &
                      &                          scheme    = scheme,           &
                      &                          glob      = glob,             &
                      &                          fieldProp = fieldProp,        &
                      &                          level     = sourceLevel,      &
                      &                          layout    = scheme%layout,    &
                      &                          intp      = intp,             &
                      &                          ind       = ind,              &
                      &                          varSys    = scheme%varSys,    &
                      &                          physics   = params%physics,   &
                      &                          derVarPos = scheme%derVarPos  )
                  end if
                else
                  if( incompressible ) then
                    call fillFinerGhostsFromMe_quadratic3dIncomp(              &
                      &                          pdf       = pdf,              &
                      &                          scheme    = scheme,           &
                      &                          glob      = glob,             &
                      &                          fieldProp = fieldProp,        &
                      &                          level     = sourceLevel,      &
                      &                          layout    = scheme%layout,    &
                      &                          intp      = intp,             &
                      &                          ind       = ind,              &
                      &                          varSys    = scheme%varSys,    &
                      &                          physics   = params%physics,   &
                      &                          derVarPos = scheme%derVarPos  )
                  else
                    call fillFinerGhostsFromMe_quadratic3d(                    &
                      &                          pdf       = pdf,              &
                      &                          scheme    = scheme,           &
                      &                          glob      = glob,             &
                      &                          fieldProp = fieldProp,        &
                      &                          level     = sourceLevel,      &
                      &                          layout    = scheme%layout,    &
                      &                          intp      = intp,             &
                      &                          ind       = ind,              &
                      &                          varSys    = scheme%varSys,    &
                      &                          physics   = params%physics,   &
                      &                          derVarPos = scheme%derVarPos  )
                  end if
                end if
              end if
            end if
            !---------------------------------------------------------------------------

            tgtPdf = pdf( targetLevel )%state( 1:QQ, pdf(targetLevel)%nNext )

            ! calculate target moments from target pdf
            write(*,*) "Calculate tgtMom from tgtPDF"
            ! rho
            tgtMom(1) = (sum(tgtPDF)) * cs2
            tgtMom(1) = tgtMom(1) * params%physics%fac( targetLevel )%press

            ! velocity
            do iVal = iVelMin, iVelMax
              iPos = iVal - iVelMin + 1
              tgtMom(iVal) =   sum( tgtPdf(:) &
                &            * real( scheme%layout%stencil(1)%cxDir(iPos,:), rk))
            end do
            if( .not. incompressible ) then
              tgtMom( iVelMin:iVelMax ) = tgtMom(iVelMin:iVelMax) / sum(tgtPDF)
            end if
            ! convert LB to physical
            tgtMom( iVelMin:iVelMax ) =   tgtMom(iVelMin:iVelMax)  &
              &                         * params%physics%fac( targetLevel )%vel

            ! strain rate
            if ( incompressible ) then
              Sxx = getShearRateTensor_acoustic(  subset = tgtPdf, &
                &                                 layout = scheme%layout, &
                &                                 rho0 = rho0,  &
                &                                 omega = omTgt ) &
                &   * params%physics%fac( targetLevel )%strainRate
            else
              Sxx = getShearRateTensor_acoustic( subset = tgtPdf, &
                &                                layout = scheme%layout, &
                &                                omega = omTgt ) &
                &   * params%physics%fac( targetLevel )%strainRate
            end if

            if ( scheme%layout%stencil(1)%nDims == 3 ) then
              tgtMom( iSxxMin+0 ) = Sxx(1,1)
              tgtMom( iSxxMin+1 ) = Sxx(2,2)
              tgtMom( iSxxMin+2 ) = Sxx(3,3)
              tgtMom( iSxxMin+3 ) = Sxx(1,2)
              tgtMom( iSxxMin+4 ) = Sxx(2,3)
              tgtMom( iSxxMax )   = Sxx(1,3)
            else
              tgtMom( iSxxMin+0 ) = Sxx(1,1)
              tgtMom( iSxxMin+1 ) = Sxx(2,2)
              tgtMom( iSxxMax )   = Sxx(1,2)
            end if
            ! Setting higher moments to 0.
            tgtMom( iSxxMax+1: ) = 0._rk
            write(*,*) 'tgtMom LB rho: ', tgtMom(1)
            write(*,*) 'tgtMom LB vel: ', tgtMom(iVelMin:iVelMax)
            write(*,*) 'tgtMom LB Sxx: ', tgtMom(iSxxMin:iSxxMax)
            write(*,*) ''

            ! Evaluate the reference polynomials at the targetCoord to obtain
            ! the reference solution (refMom) for each quantity
            do iElem = 1, QQ
              refMom(iElem) = polynomial( a = momentCoefficients(:,iElem), &
              &                       coord = targetCoord, &
              &                       order = refOrder, &
              &                       nDims = iDimn )
            end do
            ! Set the higher order moments to 0
            refMom(iSxxMax+1:) = 0._rk

            ! calculate error by comparing refMom and tgtMom
            call evaluate_error( error, refMom, tgtMom, scheme%layout%stencil(1)%QQ )

            errCode = check_success( error, trim(intpLabel( iIntp ))//'_'//trim(params%scaling ))
            if( errCode .ne. 0 ) return errCode

            deallocate( srcMom )
            deallocate( srcPdf )
            deallocate( refMom )
            deallocate( tgtMom )
            deallocate( tgtPdf )
            deallocate( targetEq )
            deallocate( varPos )
            call mus_kill_derVarList( scheme%derVarList )
            call mus_kill_stateVars( scheme )
            call tem_horizontalSpacer(fUnit = logUnit(1))
            call kill_deps( levelDesc = glob%levelDesc( targetLevel ), &
              &             fromFiner = (iIntp .eq. iFromFiner))

          end do ! iChild
        end do ! iDim
      end do ! iFromFiner, iFromCoarser
    end do ! order

    return errCode

  end subroutine check_interpolation
  !****************************************************************************


!****************************************************************************
  ! initialize dependencies
  subroutine evaluate_error( error, refVal, val, nVals  )
    !---------------------------------------------------------------------------
    real(kind=rk), intent(out) :: error
    integer, intent(in) :: nVals
    real(kind=rk), intent(in) :: refVal( nVals )
    real(kind=rk), intent(in) :: Val( nVals )
    !---------------------------------------------------------------------------
    integer :: iElem
    real(kind=rk) :: errRel
    !---------------------------------------------------------------------------
    error = 0._rk
    write(*,*)  ' ERROR EVALUATION '
    ! Evaluate the error for each moment
    do iElem = 1, nVals
      ! Relative error if feasible,
      if( abs(refVal(iElem )) .gt. eps ) then
        errRel = abs((refVal(iElem) - val( iElem ))/refVal(iElem))
      else !... if not, absolute error
        errRel = abs(refVal(iElem) - val( iElem ))
      end if
      write(*,'(3(a,e18.8))') '  ref  ', refVal(iElem), ' val ', val( iElem ), &
      &  '  E ', errRel
      error = error + errRel
    end do
  end subroutine evaluate_error
  !****************************************************************************


!****************************************************************************
  ! determine if the error is below a certain threshold
  function check_success( error, text ) result( errCode )
    !---------------------------------------------------------------------------
    real(kind=rk), intent(in) :: error
    integer :: errCode
    character(len=*), intent(in), optional :: text
    !---------------------------------------------------------------------------
    character(len=labelLen) :: buffer
    !---------------------------------------------------------------------------
    if( present( text )) then
      buffer = ' for '//trim(text)
    else
      buffer = ''
    endif
    write(*,*) 'Error '//trim(buffer)//':', error
    if( abs(error) .gt. error_limit ) then
      write(*,*) 'ERROR! ', abs(error), ' > ', error_limit
      errCode = 10
    else
      write(*,*) 'finished with SUCCESS.'
      errCode = 0
    end if
  end function check_success
  !****************************************************************************

  !****************************************************************************
  ! initialize dependencies
  ! set either fromFiner dependencies
  ! or fromCoarser dependencies
  subroutine set_intpLevel( minLevel, maxLevel, sourceLevel, targetLevel, &
    &                       fromFiner )
    !---------------------------------------------------------------------------
    integer, intent(in)  ::  minLevel, maxLevel
    integer, intent(out) ::  sourceLevel, targetLevel
    logical, intent(in)  ::  fromFiner
    !---------------------------------------------------------------------------
    if( fromFiner ) then
      sourceLevel = maxLevel
      targetLevel = minLevel 
    else
      sourceLevel = minLevel
      targetLevel = maxLevel 
    end if
  end subroutine set_intpLevel
  !****************************************************************************

  !****************************************************************************
  ! initialize dependencies
  ! set either fromFiner dependencies 
  ! or fromCoarser dependencies
  subroutine init_stateArrays( levelDesc, fromFiner, layout, pdf, minLevel,  &
    &                          maxLevel, nSourceElems  )
    !---------------------------------------------------------------------------
    integer, intent(in) ::  minLevel, maxLevel, nSourceElems
    logical, intent(in) ::  fromFiner
    type( tem_levelDesc_type )   :: levelDesc(minLevel:)
    type( pdf_data_type ), allocatable, intent(out)  :: pdf(:)
    type( mus_scheme_layout_type ), intent(in) :: layout
    !---------------------------------------------------------------------------
    integer :: iLevel, nElems
    integer :: targetLevel, sourceLevel
    !---------------------------------------------------------------------------

    call set_intpLevel( minLevel = minLevel, maxLevel = maxLevel,  &
      & sourceLevel = sourceLevel, targetLevel = targetLevel, fromFiner = fromFiner )

    allocate( pdf( 1:maxLevel ))
    pdf(:)%nNext = 1
    pdf(:)%nNow  = 1

    ! Assign initial values
    do iLevel = minLevel, maxLevel
      if( iLevel == targetLevel ) then ! assign QQ source elements
        nElems = 1
      else ! .. and 1 target element
        nElems = nSourceElems
      endif
      ! Allocate and reset the state array
      if( allocated( pdf( iLevel )%state )) deallocate( pdf( iLevel )%state )
      allocate( pdf( iLevel )%state( nElems*layout%stencil(1)%QQ, 1 ))
      pdf( iLevel )%state(:,:) = -1._rk
      ! Allocate the neighbor array. Is it needed?
      if( allocated( pdf( iLevel )%neigh )) deallocate( pdf( iLevel )%neigh )
      allocate( pdf( iLevel )%neigh( nElems*layout%stencil(1)%QQ ))
      ! Allocate moment buffer with size ( QQ:nSourceELems )
      if( allocated( levelDesc( iLevel )%intpBufFromFiner )) &
        deallocate( levelDesc( iLevel )%intpBufFromFiner )
      allocate( levelDesc( iLevel )%intpBufFromFiner(   &
        &          layout%stencil(1)%QQ , nSourceElems ))
      ! Offset for all element types ( actually only needed for gfC here)
      levelDesc( iLevel )%offset(:,:) = 0
    enddo

  end subroutine init_stateArrays
  !****************************************************************************

  !****************************************************************************
  ! deallocate dependencies
  subroutine kill_deps( levelDesc, fromFiner )
    !---------------------------------------------------------------------------
    type( tem_levelDesc_type )   :: levelDesc
    logical, intent(in) ::  fromFiner
    !---------------------------------------------------------------------------
    if( fromFiner ) then
      deallocate( levelDesc%depFromFiner )
    else
      deallocate( levelDesc%depFromCoarser )
    end if

  end subroutine kill_deps
  !****************************************************************************

  !****************************************************************************
  ! initialize dependencies
  ! set either fromFiner dependencies
  ! or fromCoarser dependencies
  subroutine init_deps( levelDesc, nSourceElems, fromFiner, ind, layout,       &
    &                   glob, order, sourceCoord, intp, minLevel, maxLevel,  &
    &                   fieldProp )
    !---------------------------------------------------------------------------
    integer, intent(in) ::  nSourceElems, minLevel, maxLevel
    integer, intent(in) ::  order !< interpolation order: 1 linear, 2 quadratic
    logical, intent(in) ::  fromFiner
    type( pdf_global_type ) :: glob
    type( interpolation_type ), intent(out) :: intp
    type( grw_intArray_type ), intent(out)   :: ind
    type( tem_levelDesc_type )   :: levelDesc
    type( mus_scheme_layout_type ), intent(in) :: layout
    real(kind=rk), intent(out), allocatable :: sourceCoord(:,:)
    type( mus_field_prop_type ) :: fieldProp(:)
    !---------------------------------------------------------------------------
    integer :: iElem, dPos
    real(kind=rk) :: targetCoord(3), distance(3), sumWeight
    real(kind=rk), allocatable :: weight(:)
    !---------------------------------------------------------------------------
    intp%nMaxSources = nSourceElems        ! Initialize the variable system
    if( allocated( sourceCoord )) deallocate( sourceCoord )
    allocate( sourceCoord(3, nSourceElems ))

    if( fromFiner ) then
      ! Initialize the source from coarser for the dependencies
      call init( me = levelDesc%sourceFromFiner,           &
        &        unique = .true. ) 
      do iElem = 1, nSourceElems
        call append( me = levelDesc%sourceFromFiner,       &
          &          val = iElem, pos = dPos ) 
      end do
      levelDesc%intpFromFiner_high%nVals = 1
      if( allocated( levelDesc%intpFromFiner_high%val) ) &
         deallocate( levelDesc%intpFromFiner_high%val) 
      allocate( levelDesc%intpFromFiner_high%val( 1))
      levelDesc%intpFromFiner_high%val( 1) = 1
      if( allocated( levelDesc%depFromFiner)) &
         deallocate( levelDesc%depFromFiner)
      allocate( levelDesc%depFromFiner(1))
      allocate(levelDesc%depFromFiner(1)                   &
        &          %stencil(3, nSourceElems))
      do iElem = 1, nSourceElems
        levelDesc%depFromFiner(1)%stencil(:, iElem) = &
        & layout%stencil(1)%cxDir(:, iElem )
      end do

      call init( me = levelDesc%depFromFiner(1)%elem )
      call init( me = levelDesc%depFromFiner(1)%elemBuffer )
      do iElem = 1, nSourceElems
        call append( me = levelDesc%depFromFiner(1)%elem, &
          &          val = iElem )
        call append( me = levelDesc%depFromFiner(1)%elemBuffer, &
          &          val = iElem )
      enddo
      ! Source elements are always sorted according to the morton space-filling
      ! curve
      do iElem = 1, nSourceElems
        sourceCoord(:, iElem ) = 0.25_rk*real( childPosition( iElem, :), kind=rk) + 0.5_rk
      end do
    else
      ! Initialize the source from coarser for the dependencies
      call init( me = levelDesc%sourceFromCoarser,           &
        &        unique = .true. ) 
      do iElem = 1, nSourceElems
        call append( me = levelDesc%sourceFromCoarser,       &
          &          val = iElem, pos = dPos ) 
      end do
      if( allocated( levelDesc%depFromCoarser)) &
         deallocate( levelDesc%depFromCoarser )
      allocate( levelDesc%depFromCoarser(1))
      levelDesc%intpFromCoarser_high%nVals = 1
      if( allocated( levelDesc%intpFromCoarser_high%val) ) &
         deallocate( levelDesc%intpFromCoarser_high%val) 
      allocate( levelDesc%intpFromCoarser_high%val( 1))
      allocate(levelDesc%depFromCoarser(1)                   &
        &          %stencil(3, nSourceElems))
      levelDesc%intpFromCoarser_high%val( 1) = 1
      do iElem = 1, nSourceElems
        levelDesc%depFromCoarser(1)%stencil(:, iElem) = &
        & layout%stencil(1)%cxDir(:, iElem )
      end do

      call init( me = levelDesc%depFromCoarser(1)%elem )
      call init( me = levelDesc%depFromCoarser(1)%elemBuffer )
      do iElem = 1, nSourceElems
        call append( me = levelDesc%depFromCoarser(1)%elem, &
          &          val = iElem )
        call append( me = levelDesc%depFromCoarser(1)%elemBuffer, &
          &          val = iElem )
      enddo 

      select case( order ) 
      case( 1 )
        ! Assume for the linear from coarser the following neighborhood
        ! source elements are aligned according to the morton curve. 
        ! The parent element of the target element has its barycenter at 
        ! (0,0,0)
        targetCoord = [ 0.25_rk, 0.25_rk, 0.25_rk ]
        ! The target element barycenter is for all configurations (0.25, 0.25, 0.25)
        ! For the last child, this means there is no rotation.
        ! For all other child configurations, the parent elements must be 
        ! rotated intrinsically for the interpolation according to the 
        ! tem_rotationMatrix to always get the same problem as for the 
        ! last child configuration
        allocate( levelDesc%depFromCoarser(1)%weight( nSourceElems ))
        allocate( weight( nSourceElems ))
        do iElem = 1, nSourceElems
          sourceCoord(:, iElem ) = real(( childPosition( iElem, :) + 1 )/2, kind=rk)
          distance(:) = abs(sourceCoord(:, iElem ) - targetCoord(:))
  !        if( layout%stencil(1)%nDims .eq. 2 ) distance(3) = 0._rk
          write(*,*) iElem, 'dist', distance
          ! Initialize the weighs and coordinates
          weight( iElem ) =    &
             (1._rk-abs(distance(1)))*(1._rk-abs(distance(2)))*(1._rk-abs(distance(3)))
        end do
        sumWeight = sum( weight( 1:nSourceElems ))
        levelDesc%depFromCoarser(1)%weight( 1:nSourceElems ) = weight( 1:nSourceElems ) / sumWeight
      case( 2)
        call mus_init_moments( me = layout%moment, layout = layout )
        do iElem = 1, nSourceElems
          sourceCoord(:, iElem ) = layout%stencil(1)%cxDirRK( :, iElem )
        end do

        ! Initialize the matrix
        if( allocated(intp%fromFiner)) deallocate(intp%fromFiner )
        if( allocated(intp%fromCoarser)) deallocate(intp%fromCoarser )
        if( allocated(levelDesc%intpBufFromCoarser))                           &
          &    deallocate(levelDesc%intpBufFromCoarser )
        allocate( intp%fromFiner(   minLevel:maxLevel ))
        allocate( intp%fromCoarser( minLevel:maxLevel ))
        allocate( levelDesc%intpBufFromCoarser( layout%stencil(1)%QQ, nSourceElems ))
        call mus_init_intpMatrices( me = intp%fromCoarser( maxLevel ),         &
          &                         intp = intp,                               &
          &                         layout = layout,                           &
          &                         dep  = levelDesc%depFromCoarser,           &
          &                         glob  = glob,                              &
          &                         fieldProp = fieldProp(1),                  &
          &                         compact   = .false.,                       &
          &                         debug = .true.,                            &
          &                         targetLevel = maxLevel,                    &
          &                         sourceLevel = minLevel )

      end select
    end if
    ! Indexing array for pointing to a position where the target ghost is stored
    call init( me = ind )
    call append( me = ind, val = 1 )

  end subroutine init_deps
  !****************************************************************************


end program mus_intpAcoustic_test
!******************************************************************************
