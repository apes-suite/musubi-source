!> This program serves as a unit test for diffusive interpolation.
!! Key steps includes:
!! Set up polynomial coefficients for each moments by random number (set Coefficients)
!! calculate source element moments by polynomial
!! convert physical moments into LB units
!! calculate pdf from macroscopic quantities (rho, vel, strain)
!! computer terget pdf by interpolation routine
!! calculate target moments from target pdf
!! calculate ref moment from polynomial
!! compare the difference between ref moment and target moment
program mus_intpDiffusive_test
  use aotus_module ,            only: flu_State, &
    &                                 close_config

  use env_module,               only: rk, pathLen, eps, labelLen, newUnit
  use treelmesh_module,         only: treelmesh_type
  use tem_logging_module,       only: logUnit
  use tem_construction_module,  only: tem_levelDesc_type
  use tem_grow_array_module,    only: grw_intArray_type, init, append
  use tem_dyn_array_module,     only: dyn_intArray_type, init, append
  use tem_tools_module,         only: tem_horizontalSpacer
  use tem_param_module,         only: childPosition, cs2inv, cs2, rho0
  ! use tem_var_system_module,   only: tem_var_sys_init
  use tem_variable_module,      only: tem_variable_type, &
    &                                 assignment(=), init, &
    &                                 append
  use tem_polynomial_module,    only: polynomial, set_divU, set_coefficients
  use tem_general_module,       only: tem_general_type, tem_load_general, &
    &                              tem_start, tem_finalize

  use mus_field_prop_module,         only: mus_field_prop_type
  use mus_param_module,              only: mus_param_type, mus_load_param
  use mus_interpolate_header_module, only: interpolation_type
  use mus_interpolate_linear_module, only:                            &
    &                                 fillMyGhostsFromFiner_linearDiffusive,  &
    &                                 fillFinerGhostsFromMe_linearDiffusive
  use mus_interpolate_quadratic_module, only: mus_init_intpMatrices,          &
    &                            fillFinerGhostsFromMe_quadratic2dDiffusive,  &
    &                            fillFinerGhostsFromMe_quadratic3dDiffusive
  use mus_pdf_module,            only: pdf_data_type, pdf_global_type
  use mus_moments_module,        only: set_momentIndices, &
    &                                  mus_init_moments
  use mus_variable_module,       only: mus_build_varSys
  ! use mus_scheme_module,        only: mus_init_stateVars, mus_writeSolSpecChar, &
  !   &                                 mus_kill_stateVars
  use mus_scheme_type_module,    only: mus_scheme_type
  use mus_scheme_layout_module,  only: mus_scheme_layout_type, &
    &                                 mus_define_d3q19,       &
    &                                 mus_define_d2q9,        &
    &                                 mus_init_stencil
  use mus_derivedQuantities_module2, only: set_pdfDiffusive, &
    &                                      getShearRateTensor_diffusive
  use mus_config_module,         only: mus_open_config
  use mus_physics_module,        only: mus_load_physics
  use mus_solver_type_module,    only: mus_solver_type

  implicit none

  !****************************************************************************
  integer :: returnFromFiner, nUnit
  integer, parameter :: iField = 1
  integer, parameter :: iFromFiner   = 1
  integer, parameter :: iFromCoarser = 2
  logical :: error
  character(len=pathLen) :: filename
  character(len=labelLen) :: intpLabel(2)
  type( flu_State ), allocatable :: conf(:) !< flu state
  real(kind=rk) :: error_limit
  type(mus_param_type) :: params

  intpLabel(iFromFiner)   = 'fromFiner'
  intpLabel(iFromCoarser) = 'fromCoarser'

  call tem_start('unit test', 'utest', params%general)
  error_limit = eps * 100000._rk
  write(*,*) 'error limit', error_limit

  ! create dummy lua file
  filename = 'test.lua'
  nUnit = newUnit()
  open( file = filename, unit = nUnit )
  write(nUnit,*) '--dummy lua file'
  close(nUnit)

  ! open musubi config file and solver specific lua functions as chunk
  call mus_open_config( conf = conf, filename = filename, &
    &                   proc = params%general%proc )
  ! load and initialize the main_debug and primary logging
  call tem_load_general( me = params%general, conf = conf(1) )
  call mus_load_param( params=params, conf = conf(1))
  call mus_init_variables

  error = .false.
  write(*,*) 'Running interpolation test for DIFFUSIVE scaling...'
  call tem_horizontalSpacer(fUnit = logUnit(1))

  call tem_horizontalSpacer(fUnit = logUnit(1), before = 2)
  write(*,*) 'INTERPOLATION FROM FINER...'
  write(*,*)
  call check_intpDiffusive( params = params, errCode = returnFromFiner )
  write(*,*) 'Result', returnFromFiner
  if( returnFromFiner .ne. 0 ) error = .true.
  call tem_horizontalSpacer(fUnit = logUnit(1))

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
  subroutine check_intpDiffusive( params, errCode )
    !---------------------------------------------------------------------------
    type( mus_param_type ), intent(inout) :: params
    integer, intent(out) :: errCode
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
    real(kind=rk) :: omega, targetCoord(3)
    integer :: nDOFs, offset, iChild, iDimn, sourceLevel, targetLevel
    real(kind=rk) :: error, Sxx(3,3)
    integer :: iPress, iVelMin, iVelMax, iSxxMin, iSxxMax
    type( tem_variable_type ) :: tVar(0)
    integer :: order, iIntp, iMin, iMax
    integer, allocatable :: children(:) !< Which child configurations to check
    integer :: refOrder !< reference order
    integer :: QQ
    type( mus_solver_type ), target :: solver
    !---------------------------------------------------------------------------
    scheme%header%kind = 'lbm_incomp'
    params%scaling = 'diffusive'
    refOrder = 1
    errCode = 0
    minLevel = 3
    maxLevel = 4
    omega = 1.5_rk
    ! Set the time layer to 1, so the interpolation routine always uses
    refIntp%timeLayer = 1

    allocate( glob%levelDesc( 1:maxLevel ))
    allocate( fieldProp(1))
    ! Diffusive scaling: omega constant throughout all levels
    fieldProp(iField)%fluid%omLvl( : ) = omega
    tree%global%minLevel = minLevel
    tree%global%maxLevel = maxLevel
    tree%global%boundingCubeLength = 1._rk
    scheme%nFields = 1
    allocate( scheme%field(scheme%nFields ))
    scheme%field(1)%label = 'lbm'
    call mus_load_physics( me= params%physics, conf = conf(1),            &
      &                    tree = tree, scaleFactor = params%scaleFactor, &
      &                    dtRef = 0.007_rk, dxRef = 0.03_rk )

    tree%nElems = 1
    intp = refIntp

    do order = 1, 2
      ! Set the reference polynomial order to be the same as the interpolation
      ! order. This has to yield results exact up to numerical accuracy
      refOrder = order
      if ( order == 1 ) then ! Test fromFiner and fromCoarser for linear
        iMin = iFromFiner
        iMax = iFromCoarser
      else ! For quadratic, test only the fromCoarser
           ! (fromFiner does not exist)
        iMin = iFromCoarser
        iMax = iFromCoarser
      end if

      do iIntp = iMin, iMax

        write(logUnit(1),*)'Checking '//trim(intpLabel( iIntp ))
        write(logUnit(1),*)'ORDER    ', order

        ! set up sourceLevel and targetLevel
        call set_intpLevel( minLevel = minLevel, &
          &                 maxLevel = maxLevel, &
          &                 sourceLevel = sourceLevel, &
          &                 targetLevel = targetLevel,    &
          &                 fromFiner = (iIntp==iFromFiner ) )

        if( allocated( children )) deallocate(children)
        if( iIntp == iFromFiner   ) then
          allocate(children(1))
          children(1) = 1
        else
          allocate(children(8))
          children(1:8) = [1,2,3,4,5,6,7,8]
        end if

        ! Test 2D and 3D dimension
        do iDimn = 2,3

          scheme%layout = refLayout
          call mus_init_stencil( stencil = scheme%layout%stencil, nStencils = 1 )
          select case( iDimn )
          case( 3 )
            if( order == 1 ) then
              nDOfs = 4  ! 3d linear polynomial with 4 coefficients
            else if ( order == 2 ) then
              nDofs = 10 ! 3d quadratic polynomial with 10 coefficients
            end if
            call mus_define_d3q19( layout = scheme%layout, nElems=tree%nElems)
          case( 2 )
            if ( order == 1 ) then
              nDOfs = 3 ! 2d linear polynomial with 3 coefficients
            else if ( order == 2 ) then
              nDofs = 6 ! 2d quadratic polynomial with 6 coefficients
            end if
            call mus_define_d2q9( layout = scheme%layout, nElems=tree%nElems)
          end select
          QQ = scheme%layout%stencil(1)%QQ
          call mus_init_moments( me             = scheme%layout%moment, &
            &                    layout         = scheme%layout         )

          do iChild = 1, 2**iDimn
            write(*,"(a, I1)") 'Dimension: ', iDimn
            write(*,"(a, I1)") 'iChild:    ', iChild
            write(*,"(a    )") 'interpolation: '//trim(intpLabel( iIntp ))
            write(*,"(a, F5.2)") 'omega: ', omega

            ! set up quantity index according to dimension
            call set_momentIndices( iDimn, iPress, iVelMin, iVelMax, iSxxMin, iSxxMax)

            if( iIntp == iFromFiner .or. order == 1 ) then
              ! if linear or from finer
              nSourceElems = 2**iDimn
            else
              nSourceElems = QQ
            end if

            ! prepare variable system
            ! call tem_var_sys_init( scheme%varSys )

            ! call mus_build_derVarList( scheme%derVarList, scheme%header,  &
            !   &                          scheme%layout%stencil(1))
            ! call mus_init_stateVars( scheme, size(tVar) )
            call mus_build_varSys( varSys = scheme%varSys, &
              &                    solver = solver, &
              &                    schemeHeader = scheme%header, &
              &                    stencil = scheme%layout%stencil(1), &
              &                    nFields = scheme%nFields, &
              &                    derVarPos = scheme%derVarPos, &
              &                    luaVar = scheme%luaVar, &
              &                    stateVarMap = scheme%stateVarMap,&
              &                    st_funList = scheme%st_funList  )

            call init_deps( fromFiner = (iIntp == iFromFiner ),           &
             &      levelDesc = glob%levelDesc( targetLevel ),  glob=glob,&
             &      nSourceElems = nSourceElems, layout = scheme%layout,  &
             &      ind = ind, order = order, sourceCoord = sourceCoord,  &
             &      minLevel = minLevel, maxLevel = maxLevel, &
             &      fieldProp = fieldProp, intp = intp )

            ! Initialize the polynomial coefficients for each moment
            call set_Coefficients( me = momentCoefficients, &
              &                    nDOFs = nDOFs, &
              &                    QQ = QQ, &
              &                    valMin = 0.3_rk, valMax = 1.6_rk )
            allocate( srcMom(QQ,QQ) )
            allocate( srcPdf(QQ) )
            allocate( refMom(QQ) )
            allocate( tgtMom(QQ) )
            allocate( tgtPdf(QQ) )
            allocate( targetEq(QQ) )

            call init_stateArrays( pdf = pdf, &
              &                    fromFiner = ( iIntp == iFromFiner ), &
              &                    levelDesc = glob%levelDesc, &
              &                    layout = scheme%layout, &
              &                    minLevel = minLevel,&
              &                    maxLevel = maxLevel, &
              &                    nSourceElems = nSourceElems )

            ! Set states for the source elements
            do iElem = 1, nSourceElems
              write(*,'(a,i2,a,3f5.1)') 'sourceElem: ', iElem, ' crd ', sourceCoord(:, iElem )
              ! Set each quantity in terms of the polynomial.
              ! These quantities are PHYISCAL states which have to be converted to the
              ! pdfs
              do iVal = 1, QQ
                ! Set the source moments one after another.
                ! The order is the same as in Z-curve (sources are the children)
                srcMom( iVal, iElem ) = polynomial(                             &
                  &                             a = momentCoefficients(:,iVal), &
                  &                             coord = sourceCoord(:, iElem ), &
                  &                             order = refOrder, nDims = iDimn )
              end do

              ! Set the moment coefficients to satisfy div(u) = 0
!write(*,*) 'Sxx: ', srcMom( iSxxMin: iSxxMax, iElem )
!write(*,*) 'set divU'
              call set_divU( uOut = srcMom( iVelMax, iElem ), nDims = iDimn,  &
                &            coeff = momentcoefficients(:, iVelMin:iVelMax ), &
                &            nStresses = iSxxMax - iSxxMin + 1,               &
                &         SxxCoeff = momentcoefficients(:, iSxxMin:iSxxMax ), &
                &            SxxOut = srcMom( iSxxMin: iSxxMax, iElem ),      &
                &            nDofs = nDofs,  coord = sourceCoord(:, iElem ), order = reforder )
!write(*,*) 'Sxx: ', srcMom( iSxxMin: iSxxMax, iElem )

              ! convert physical moments into LB units
              srcMom( 1, iElem ) =   srcMom( 1, iElem )  &
                &                  * cs2inv/params%physics%fac( sourceLevel )%press
              srcMom(iVelMin:iVelMax, iElem ) = srcMom( iVelMin:iVelMax, iElem )     &
                &                           / params%physics%fac( sourceLevel )%vel
              !write(*,*) 'shear rate', srcMom(iSxxMin:iSxxMax, iElem )
              srcMom(iSxxMin:iSxxMax, iElem ) = srcMom( iSxxMin:iSxxMax, iElem )     &
                &                     / params%physics%fac( sourceLevel )%strainRate

              ! calculate pdf from macroscopic quantities
              ! NOTICE: srcMom(iSxxMin:iSxxMax, iElem ) /= Sxx
              srcPdf(:) = set_pdfDiffusive( scheme%layout, omega, rho0, srcMom(:,iElem))

              ! fill pdf state vector
              offset = QQ * ( iElem - 1 )
              pdf( sourceLevel )%state( offset+1:offset+QQ, &
                                        pdf(sourceLevel)%nNext ) = srcPdf

              ! calculate Sxx from source pdf.
              ! They should be the same as srcMom( iSxxMin:iSxxMax, iElem )
              ! Sxx = getShearRateTensor_diffusive( f = srcPdf, &
              !   &                                 layout = scheme%layout,   &
              !   &                                 omega = omega ) &
              !   &     * params%physics%fac( sourceLevel )%strainRate
!write(*,*) 'getShearRateTensor_diffusive'
!write(*,*) 'Sxx: ', Sxx(1,1), Sxx(2,2), Sxx(3,3), Sxx(1,2), Sxx(2,3), Sxx(1,3)
            end do

            !---------------------------------------------------------------------------
            ! Do the interpolation
            if ( iIntp == iFromFiner ) then
              ! fill coarse from finer
              targetCoord(:) = [ 0.5_rk, 0.5_rk, 0.5_rk  ]
              call fillMyGhostsFromFiner_linearDiffusive(          &
                &              pdf       = pdf,                    &
                &              scheme    = scheme,                 &
                &              glob      = glob,                   &
                &              fieldProp = fieldProp,              &
                &              level     = targetLevel,            &
                &              layout    = scheme%layout,          &
                &              intp      = intp,                   &
                &              ind       = ind,                    &
                &              varSys    = scheme%varSys,          &
                &              physics   = params%physics,         &
                &              derVarPos = scheme%derVarPos,       &
                &              time      = 0.0_rk  )
            else
              ! fill finer from coarser
              if( order == 1 ) then
                targetCoord(:) = [ 0.25_rk, 0.25_rk, 0.25_rk  ]
                call fillFinerGhostsFromMe_linearDiffusive(        &
                  &              pdf       = pdf,                  &
                  &              scheme    = scheme,               &
                  &              glob      = glob,                 &
                  &              fieldProp = fieldProp,            &
                  &              level     = sourceLevel,          &
                  &              layout    = scheme%layout,        &
                  &              intp      = intp,                 &
                  &              ind       = ind,                  &
                  &              varSys    = scheme%varSys,        &
                  &              physics   = params%physics,       &
                  &              derVarPos = scheme%derVarPos,     &
                  &              time      = 0.0_rk  )
              else
                targetCoord(:) = real( childPosition(iChild,:), kind=rk)*0.25_rk
                write(*,*) 'Testing child', iChild, 'targetCoord', real( targetCoord )
                glob%levelDesc( targetLevel )%depFromCoarser(1)%coord = targetCoord
                if( iDimn == 2 ) then
                  call fillFinerGhostsFromMe_quadratic2dDiffusive(             &
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
                    &                          derVarPos = scheme%derVarPos,   &
                    &                          time      = 0.0_rk              )
                else
                  call fillFinerGhostsFromMe_quadratic3dDiffusive(             &
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
                    &                          derVarPos = scheme%derVarPos,   &
                    &                          time      = 0.0_rk  )
                end if
              end if
            end if
            !---------------------------------------------------------------------------

            ! get target pdf
            tgtPdf = pdf( targetLevel )%state( 1:QQ, &
              &      pdf(targetLevel )%nNext )

            ! Get the physical quantities from the pdf
            ! Pressure
            tgtMom(1) = (sum(tgtPdf))*cs2*params%physics%fac(targetLevel)%press
            ! Velocity? momentum?
            do iVal = iVelMin, iVelMax
              iPos = iVal - iVelMin + 1
              tgtMom( iVal ) = sum( tgtPdf(:) &
                &          * real( scheme%layout%stencil(1)%cxDir(iPos,:), rk))
            enddo
            tgtMom( iVelMin:iVelMax ) = tgtMom( iVelMin:iVelMax ) &
              &                    * params%physics%fac( targetLevel )%vel

            ! strain rate
            Sxx = getShearRateTensor_diffusive( f = tgtPdf, &
              &                                 layout = scheme%layout,   &
              &                                 omega  = omega ) &
              &       * params%physics%fac( targetLevel )%strainRate

            if( scheme%layout%stencil(1)%nDims == 3 ) then
              tgtMom( iSxxMin )   = Sxx(1,1)
              tgtMom( iSxxMin+1 ) = Sxx(2,2)
              tgtMom( iSxxMin+2 ) = Sxx(3,3)
              tgtMom( iSxxMin+3 ) = Sxx(1,2)
              tgtMom( iSxxMin+4 ) = Sxx(2,3)
              tgtMom( iSxxMax )   = Sxx(1,3)
            else
              tgtMom( iSxxMin )   = Sxx(1,1)
              tgtMom( iSxxMin+1 ) = Sxx(2,2)
              tgtMom( iSxxMax )   = Sxx(1,2)
            end if

            ! Setting higher moments to 0
            tgtMom( iSxxMax+1: ) = 0._rk

            ! Evaluate the reference polynomials at the targetCoord to obtain
            ! the reference solution for each quantity
            do iElem = 1, QQ
              refMom(iElem) = polynomial( a = momentCoefficients(:,iElem), &
                &                         coord = targetCoord, &
                &                         order = refOrder, nDims = iDimn )
            end do
            ! Set the higher order moments to 0
            refMom(iSxxMax+1:) = 0._rk

            ! and now compare refMom <-> tgtMom
            call evaluate_error( error, refMom, tgtMom, QQ )

            errCode = check_success( error, trim(intpLabel( iIntp ))//&
              &                             '_'//trim(params%scaling ))
            if ( errCode /= 0 ) return errCode

            deallocate( srcMom)
            deallocate( srcPdf)
            deallocate( refMom)
            deallocate( tgtMom)
            deallocate( tgtPdf)
            deallocate( targetEq)
            ! call mus_kill_derVarList( scheme%derVarList )
            ! call mus_kill_stateVars( scheme )
            call tem_horizontalSpacer(fUnit = logUnit(1))
            call kill_deps( levelDesc = glob%levelDesc( targetLevel ), &
              &             fromFiner = (iIntp .eq. iFromFiner))

          end do ! iChild
        end do ! iDim
      end do ! iFromFiner, iFromCoarser
    end do ! order

    return errCode

  end subroutine check_intpDiffusive
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
  function check_success( error, text ) result( errCode)
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
    &          glob,    order, sourceCoord, intp, minLevel, maxLevel,  &
    &          fieldProp )
    !---------------------------------------------------------------------------
    integer, intent(in) ::  nSourceElems, minLevel, maxLevel
    integer, intent(in) ::  order !< interpolation order: 1 linear, 2 quadratic
    logical, intent(in) ::  fromFiner
    type( pdf_global_type ) :: glob
    type( interpolation_type ), intent(out) :: intp
    type( grw_intArray_type ),  intent(out) :: ind
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
        do iElem = 1, nSourceElems
          sourceCoord(:, iElem ) = real( layout%stencil(1)%cxDir( :, iElem ), kind=rk)
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


end program mus_intpDiffusive_test
!******************************************************************************
