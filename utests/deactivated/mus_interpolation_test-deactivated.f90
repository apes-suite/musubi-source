program mus_interpolation_test
  use aotus_module, only: flu_State, aot_get_val
  use aot_table_module, only: aot_table_open, aot_table_close, aot_table_length

  use env_module, only: rk
  use treelmesh_module, only: treelmesh_type
  use tem_aux_module, only: tem_abort, tem_start
  use tem_element_module, only: eT_ghostFromCoarserElem,                 &
    &                                eT_ghostFromFinerElem
  use tem_geometry_module, only: tem_baryOfId, tem_intp_bilinear
  use tem_grow_array_module, only: grw_intArray_type, init, append
  use tem_dyn_array_module, only: dyn_intArray_type, init, append
  use tem_tools_module, only: tem_write, tem_horizontalSpacer
  use tem_logging_module, only: logUnit
  use tem_comm_env_module, only: process, tem_comm_env_init_empty
  use tem_param_module, only: childPosition
  use tem_var_system_module, only: tem_var_system_type, tem_var_append, &
  & tem_var_sys_init, tem_var_sys_mark
  use tem_variable_module, only: tem_variable_type, tem_variable_out,  &
    &                            tem_variable_from, tem_variable_dump, &
    &                            assignment(=), dyn_varArray_type, init, &
    &                            append, PositionOfVal
  use tem_polynomial_module, only: polynomial, set_divU, set_coefficients

  use mus_field_prop_module, only: mus_field_prop_type 
  use mus_param_module,  only: mus_param_type
  use mus_interpolate_header_module, only: interpolation_type,                 &
                                           determine_timeLayer
  use mus_interpolate_linear_module, only:                                  &
          fillMyGhostsFromFiner_linear, &
          fillMyGhostsFromFiner_linMoment3d, &
          fillMyGhostsFromFiner_linMoment2d
  use mus_interpolate_quadratic_module, only:                                  &
          fillFinerGhostsFromMe_quadratic2d,                           &
          fillFinerGhostsFromMe_quadratic3d,                           &
          mus_init_intpMatrices
  use mus_moments_module, only: mus_init_moments, get_moment
  use mus_pdf_module,    only: pdf_data_type, pdf_global_type, pdf_globInfo_type
  use mus_scheme_layout_module, only: mus_scheme_layout_type
  use mus_variable_module, only: mus_derVarPos_type
  use mus_scheme_type_module, only: mus_scheme_type
  use mus_scheme_layout_module, only: mus_scheme_layout_type, mus_load_newLayout,&
    &                          mus_weights_out, mus_define_d3q19, mus_define_d3q7, &
    &                          mus_define_d2q9, mus_define_d3q6,         &
    &                          mus_destroy_stencil
  use mus_derivedQuantities_module2, only: getEquilibrium, set_pdfAcoustic, convPrePost

  implicit none
  integer :: returnFromCoarser
  integer :: returnFromFiner
  logical :: error
  type(mus_param_type) :: params

  call tem_start('unit test', 'utest', params%general%solver, params%general%proc)
  
  error = .false.
  write(logUnit(1),*) 'Running interpolation_test...'
  call tem_horizontalSpacer(fUnit=logUnit(1))
  write(logUnit(1),*) 'INTERPOLATION FROM COARSER...'
  call check_fromCoarser3D( errCode = returnFromCoarser )
  write(logUnit(1),*) 'Result', returnFromCoarser
  if( returnFromCoarser .ne. 0 ) error = .true.
  call tem_horizontalSpacer(fUnit=logUnit(1))

  call tem_horizontalSpacer( before = 2)
  !write(logUnit(1),*) 'INTERPOLATION FROM FINER...'
  !write(logUnit(1),*) 
  !call check_fromFiner( errCode = returnFromFiner )
  !write(logUnit(1),*) 'Result', returnFromCoarser
  !if( returnFromFiner .ne. 0 ) error = .true.
  !call tem_horizontalSpacer(fUnit=logUnit(1))

  call tem_horizontalSpacer(fUnit=logUnit(1))
  if (.not. error) then
    write(logUnit(1),*) 'PASSED'
  else
    write(logUnit(1),*) 'FAILED'
  end if

contains

  !****************************************************************************
  subroutine check_fromCoarser3D( errCode )
    !---------------------------------------------------------------------------
    integer, intent(out) :: errCode
    !---------------------------------------------------------------------------
    integer :: minLevel, maxLevel
    type( mus_scheme_type ) :: scheme
    type( interpolation_type ) :: intp, refIntp
    type( mus_field_prop_type ), allocatable :: fieldProp(:)
    type( mus_scheme_layout_type ) :: layout, refLayout
    type( pdf_data_type ), allocatable :: pdf(:)
    type( pdf_global_type ) :: glob
    type( tem_var_system_type ) :: varSys
    integer :: dPos, iVal, nElems, nSourceElems, iElem
    type( grw_intArray_type )   :: ind  
    type( mus_derVarPos_type ), allocatable :: derVarPos(:)
    type( treelmesh_type )    :: tree !< treelm mesh
    real(kind=rk), allocatable :: momentCoefficients(:,:)
    real(kind=rk), allocatable :: srcMom(:,:)
    real(kind=rk), allocatable :: srcPdf(:)
    real(kind=rk), allocatable :: tgtMom(:)
    real(kind=rk), allocatable :: tgtPdf(:), refMom(:), targetEq(:)
    real(kind=rk) :: omega, targetCoord(3), refTgtCoord(3)
    integer :: nDOFs, offset, iChild, iDimn
    real(kind=rk) :: error
    !---------------------------------------------------------------------------

    errCode = 0
    minLevel = 3
    maxLevel = 4
    omega = 1.8_rk
    ! Set the coordinate of the reference targetElement
    refTgtCoord = (/0.25_rk, 0.25_rk, 0.25_rk/)
    ! Set the time layer to 1, so the interpolation routine always uses
    refIntp%timeLayer = 1
   
    allocate( derVarPos(1))
    allocate( pdf( 1:maxLevel ))
    pdf(:)%nNext = 1
    pdf(:)%nNow  = 1
    allocate( glob%levelDesc( 1:maxLevel ))
    allocate( fieldProp(1))
    fieldProp( 1)%fluid%omLvl( minLevel ) = omega
    ! See relationship between omegas in Hasert 2012, ECCOMAS Eq. 25 (p.5)
    fieldProp( 1)%fluid%omLvl( maxLevel ) = 2._rk/(2._rk*(2._rk/omega - 1._rk ) + 1._rk )
!    write(logUnit(1),*) 'warning!! setting om Coarse = om Fine!'
!    fieldProp( 1)%fluid%omLvl( maxLevel ) = omega
    write(logUnit(1),*) 'omega coarse', fieldProp(1)%fluid%omLvl( minLevel ), &
               'om fine ', fieldProp(1)%fluid%omLvl( maxLevel )

    allocate( glob%levelDesc( maxLevel )%depFromCoarser(1))
    glob%levelDesc( maxLevel   )%intpFromCoarser_high%nVals = 1
    allocate( glob%levelDesc( maxLevel  )%intpFromCoarser_high%val( 1))

    do iDimn = 2, 2
      intp = refIntp
      layout = refLayout
      select case( iDimn)
      case( 3 )
        nSourceElems = 19
        nDOfs = 10      ! 3d polynomial with 10 coefficients
        intp%nMaxSources = 19
        tree%nElems = nSourceElems+1
        ! Initialize the stencil
        call mus_define_d3q19( layout = layout, tree = tree ) 
      case(2)
        nSourceElems = 9
        nDOfs = 6 ! 2d polynomial with 6 coefficients
        intp%nMaxSources = 9
        tree%nElems = nSourceElems+1
        ! Initialize the stencil
        call mus_define_d2q9( layout = layout, tree = tree ) 
      end select
      ! Initialize the source from coarser for the dependencies
      call init( me = glob%levelDesc( maxLevel )%sourceFromCoarser, unique = .true. ) 
      do iElem = 1, nSourceElems
        call append( me = glob%levelDesc( maxLevel )%sourceFromCoarser, val = iElem, pos = dPos ) 
      end do
      glob%levelDesc( maxLevel   )%intpFromCoarser_high%val( 1) = 1
      allocate(glob%levelDesc( maxLevel )%depFromCoarser(1)%stencil(3, nSourceElems))
      do iElem = 1, nSourceElems
        glob%levelDesc( maxLevel )%depFromCoarser(1)%stencil(:, iElem) = &
        & layout%fStencil%cxDir(:, iElem )
      end do

      call init( me = glob%levelDesc( maxlevel )%depFromCoarser(1)%elem )
      call init( me = glob%levelDesc( maxlevel )%depFromCoarser(1)%elemBuffer )
      do iElem = 1, nSourceElems
        call append( me = glob%levelDesc( maxlevel )%depFromCoarser(1)%elem, &
          &          val = iElem )
        call append( me = glob%levelDesc( maxlevel )%depFromCoarser(1)%elemBuffer, &
          &          val = iElem )
      enddo
      ! Indexing array for pointing to a position where the target ghost is stored
      call init( me = ind )
      call append( me = ind, val = 1 ) 
      ! initialize the transformation from and into moment space
      call mus_init_moments( me = layout%moment, layout = layout )
      ! Initialize the polynomial coefficients for each moment
      call set_coefficients( me = momentCoefficients, nDOFs = nDOFs, &
      &  QQ = layout%fStencil%QQ, valMin = 0.3_rk, valmax = 1.5_rk )
      allocate( srcMom(layout%fStencil%QQ ,layout%fStencil%QQ ))
      allocate( srcPdf(layout%fStencil%QQ ))
      allocate( refMom(layout%fStencil%QQ ))
      allocate( tgtMom(layout%fStencil%QQ ))
      allocate( tgtPdf(layout%fStencil%QQ ))
      allocate( targetEq(layout%fStencil%QQ ))
      ! Initialize the matrix
      allocate( intp%fromFiner(   minLevel:maxLevel ))
      allocate( intp%fromCoarser( minLevel:maxLevel ))
      call mus_init_intpMatrices( me = intp%fromCoarser( maxLevel ),        &
        &                            intp = intp,                           &
        &                            layout = layout,                       &
        &                 dep  = glob%levelDesc( maxLevel )%depFromCoarser, &
        &                            glob  = glob,                          &
        &                            fieldProp = fieldProp(1),       &
        &                            compact   = .false.,              &
        &                            debug = .true.,&
        &                            targetLevel = maxLevel,                &
        &                            sourceLevel = minLevel )

      ! Assign initial values
      do iVal = minLevel, maxLevel
        if( iVal == minLevel ) then ! assign QQ source elements
          nElems = layout%fStencil%QQ
        else ! .. and 1 target element
          nElems = 1
        endif
        ! Allocate and reset the state array
        if( allocated( pdf( iVal )%state )) deallocate( pdf( iVal )%state )
        allocate( pdf( iVal )%state( nElems*layout%fStencil%QQ, 1 ))
        pdf( iVal )%state(:,:) = -1._rk
        ! Allocate the neighbor array. Is it needed?
        if( allocated( pdf( iVal )%neigh )) deallocate( pdf( iVal )%neigh )
        allocate( pdf( iVal )%neigh( nElems*layout%fStencil%QQ ))
        ! Allocate moment buffer with size ( QQ:nSourceELems )
        if( allocated( glob%levelDesc( iVal )%intpBufFromCoarser )) &
          deallocate( glob%levelDesc( iVal )%intpBufFromCoarser )
        allocate( glob%levelDesc( iVal )%intpBufFromCoarser(   &
          &          layout%fStencil%QQ , nSourceElems ))
        ! Offset for all element types ( actually only needed for gfC here)
        glob%levelDesc( iVal )%offset(:,:) = 0
      enddo

      ! Set states for the source elements
      do iElem = 1, nSourceElems
        ! Set each moment
        do iVal = 1, layout%fStencil%QQ
          ! Set the source moments one after another.
          ! The order is the same as in the stencil ordering
          if( layout%fStencil%nDims .eq. 3 ) then
            srcMom( iVal, iElem ) = polynomial( a = momentCoefficients(:,iVal), &
            & coord = real( layout%fStencil%cxDir(:, iElem ), kind=rk),  &
            & order = 2, nDims =iDimn )
          else
            srcMom( iVal, iElem ) = polynomial( a = momentCoefficients(:,iVal), &
            & coord = real( layout%fStencil%cxDir(:, iElem ), kind=rk), &
            & order = 2, nDims = iDimn )
          end if
        end do
        ! shear stress:
        ! set the shear stress sigma = -(1-omega/2)*Pi1
        ! But the resulting moment Pi = Pi0 + Pi1 (+ Pi2...) 
        !                  PiA, A>1 is not available, so neglect.
        ! Pi = Pi0 + sigma/(-(1-omega/2))
        if( layout%fStencil%QQ .eq. 19 ) then
          srcMom(5:10, iElem) = srcMom(5:10, iElem)/&
            & (-(1._rk - fieldProp( 1)%fluid%omLvl( minLevel)  &
            & *0.5_rk)*convPrePost(fieldProp( 1)%fluid%omLvl( minLevel )))
          ! Add the Pi(0) part to the shear stress
          ! pressure to the diagonal components ...
          srcMom( 5:7, iElem )  = srcMom(5:7, iElem ) + srcMom(1, iElem)/3._rk
          ! ... and momentum flux to all components
          srcMom( 5, iELem ) = srcMom( 5, iElem ) + srcMom(2, iElem)*srcMom(2, iElem)/srcMom(1, iElem)
          srcMom( 6, iELem ) = srcMom( 6, iElem ) + srcMom(3, iElem)*srcMom(3, iElem)/srcMom(1, iElem)
          srcMom( 7, iELem ) = srcMom( 7, iElem ) + srcMom(4, iElem)*srcMom(4, iElem)/srcMom(1, iElem)
          srcMom( 8, iELem ) = srcMom( 8, iElem ) + srcMom(2, iElem)*srcMom(3, iElem)/srcMom(1, iElem)
          srcMom( 9, iELem ) = srcMom( 9, iElem ) + srcMom(3, iElem)*srcMom(4, iElem)/srcMom(1, iElem)
          srcMom(10, iELem ) = srcMom(10, iElem ) + srcMom(2, iElem)*srcMom(4, iElem)/srcMom(1, iElem)
        else if( layout%fStencil%QQ .eq. 9 ) then
          srcMom(4:6, iElem) = srcMom(4:6, iElem)/   &
          &  (-(1._rk - fieldProp( 1)%fluid%omLvl( minLevel )*0.5_rk)  &
          &    *convPrePost(fieldProp( 1)%fluid%omLvl( minLevel )))
          ! Add the Pi(0) part to the shear stress
          ! pressure to the diagonal components ...
          srcMom( 4:5, iElem )  = srcMom(4:5, iElem ) + srcMom(1, iElem)/3._rk
          ! ... and momentum flux to all components
          srcMom( 4, iELem ) = srcMom( 4, iElem ) + srcMom(2, iElem)*srcMom(2, iElem)/srcMom(1, iElem)
          srcMom( 5, iELem ) = srcMom( 5, iElem ) + srcMom(3, iElem)*srcMom(3, iElem)/srcMom(1, iElem)
          srcMom( 6, iELem ) = srcMom( 6, iElem ) + srcMom(2, iElem)*srcMom(3, iElem)/srcMom(1, iElem)
        end if
        write(logUnit(1),'(i4,a,19f8.4)') iElem, 'generated srcMom', srcMom(:, iElem )
        ! Transform to pdf and store to the state vector
        srcPdf( : ) = matmul( layout%moment%toPdf%A, srcMom(:, iElem ))
        ! Set the pdfs by 
!        srcPdf(:) = set_pdfAcoustic( layout, omSrc, rho0, srcMom(:, iElem ), .false. )
        offset = layout%fStencil%QQ*( iElem -1)
        pdf( minLevel )%state( offset+1:offset+layout%fStencil%QQ, &
        pdf(minLevel)%nNext ) = srcPdf
      end do

      !------------append variables to variable system------------------------------
      call tem_var_sys_init( sys = varSys )
      call tem_var_append( sys = varSys, varname = 'pdf',     &
        &                  nComponents=layout%fStencil%QQ, pos = dPos )
      call tem_var_sys_mark( sys = varSys, systemName = 'lbm')
      !-----------------------------------------------------------------------------

      do iChild = 1, 8
        targetCoord(:) = real( childPosition(iChild,:), kind=rk)*refTgtCoord(:)
        write(logUnit(1),*) 'Testing child', iChild, 'targetCoord', real( targetCoord )
        glob%levelDesc( maxLevel )%depFromCoarser(1)%coord = targetCoord
        !-----------------------------------------------------------------------------
        ! Do the interpolation
        if( layout%fStencil%nDims .eq. 3 ) then
          call fillFinerGhostsFromMe_quadratic3d(                         &
            &                          pdf       = pdf,                   &
            &                          scheme    = scheme,                &
            &                          glob      = glob,                  &
            &                          fieldProp = fieldProp,             &
            &                          level     = minLevel,              &
            &                          layout    = layout,                &
            &                          intp      = intp,                  &
            &                          ind       = ind,                   &
            &                          varSys    = varSys,                &
            &                          derVarPos = derVarPos  )  
        else
          call fillFinerGhostsFromMe_quadratic2d(                         &
            &                          pdf       = pdf,                   &
            &                          scheme    = scheme,                &
            &                          glob      = glob,                  &
            &                          fieldProp = fieldProp,             &
            &                          level     = minLevel,              &
            &                          layout    = layout,                &
            &                          intp      = intp,                  &
            &                          ind       = ind,                   &
            &                          varSys    = varSys,                &
            &                          derVarPos = derVarPos  )  
        end if
        !-----------------------------------------------------------------------------

        tgtPdf = pdf( maxLevel )%state( 1:layout%fStencil%QQ, pdf(maxLevel )%nNext )
        ! and now compare the results
        tgtMom( : ) = matmul( layout%moment%toMoments%A, tgtPdf(: ))
        write(logUnit(1),*)  ' Moment results'
        error = 0._rk
        if( layout%fStencil%QQ .eq. 19 ) then
          ! Loop over all the moments
          do iElem = 1, layout%fStencil%QQ
            ! Determine the analytical solution from the given polynomial
            ! Be aware that the higher moments have to be re-scaled!
            refMom(iElem) = polynomial( a = momentCoefficients(:,iElem), &
            & coord = targetCoord, order = 2, nDims = iDimn )
          end do

          ! Extract the shear stress
          tgtMom(5:10)  = get_shearstress3d( tgtMom, fieldProp( 1)%fluid%omLvl( maxLevel ))
          ! And calculate the equilibrium distribution from the macroscopic 
          ! quantities in tMom
          targetEq = getEquilibrium( tgtMom( 1),    & ! density
                          (/ tgtMom( 2)/tgtMom(1), & ! velX
                             tgtMom( 3)/tgtMom(1), & ! velY
                             tgtMom( 4)/tgtMom(1)/), & ! velZ
                             layout )

          ! Set the further higher moments to their equilibrium distribution
          refMom(11 ) = get_moment( targetEq, layout, (/2,1,0/) ) !mxxy
          refMom(12 ) = get_moment( targetEq, layout, (/2,0,1/) ) !mxxz
          refMom(13 ) = get_moment( targetEq, layout, (/1,2,0/) ) !myyx
          refMom(14 ) = get_moment( targetEq, layout, (/0,2,1/) ) !myyz
          refMom(15 ) = get_moment( targetEq, layout, (/1,0,2/) ) !mzzx
          refMom(16 ) = get_moment( targetEq, layout, (/0,1,2/) ) !mzzy
          refMom(17 ) = get_moment( targetEq, layout, (/2,2,0/) ) !mxxyy
          refMom(18 ) = get_moment( targetEq, layout, (/0,2,2/) ) !myyzz
          refMom(19 ) = get_moment( targetEq, layout, (/2,0,2/) ) !mzzxx      
        elseif( layout%fStencil%QQ .eq. 9 ) then
          ! Loop over all the moments
          do iElem = 1, layout%fStencil%QQ
            ! Determine the analytical solution from the given polynomial
            ! Be aware that the higher moments have to be re-scaled!
            refMom(iElem) = polynomial( a = momentCoefficients(:,iElem), &
            & coord = targetCoord, order = 2, nDims = iDimn)
          end do

          ! Extract the shear stress
          tgtMom(4:6)  = get_shearstress2d( tgtMom, fieldProp( 1)%fluid%omLvl( maxLevel ))
          ! And calculate the equilibrium distribution from the macroscopic 
          ! quantities in tMom
          targetEq = getEquilibrium( tgtMom( 1),    & ! density
                          (/ tgtMom( 2)/tgtMom(1), & ! velX
                             tgtMom( 3)/tgtMom(1), & ! velY
                             0._rk/), & ! velZ
                             layout )

          ! Set the further higher moments to their equilibrium distribution
          refMom(7) = get_moment( targetEq, layout, (/2,1,0/) ) !mxxy
          refMom(8) = get_moment( targetEq, layout, (/1,2,0/) ) !mxxz
          refMom(9) = get_moment( targetEq, layout, (/2,2,0/) ) !myyx
        else
          write(logUnit(1),*) 'Interpolation not completely implemted yet. Abort.'
          stop
        end if

        ! Evaluate the error for each moment
        do iElem = 1,layout%fStencil%QQ
          write(logUnit(1),*) 'ref moment', refMom(iElem), 'intp mom', tgtMom( iElem ), &
          &  'E', abs((refMom(iElem) - tgtMom( iElem ))/refMom(iElem))
          error = error + abs(refMom(iELem) - tgtMom( iElem ))
          write(logUnit(1),*) 'error', error
        end do
        if( error .gt. 1.E-13 ) then
          write(logUnit(1),*) 'ERROR!'
          write(logUnit(1),*) 'Error from the fillFinerGhostsFromMe routine:', error
          errCode = max( errCode, 10)
        else
          write(logUnit(1),*) 'finished with SUCCESS.'
          errCode = max( errCode, 0)
        end if
        call tem_horizontalSpacer(fUnit=logUnit(1))
      end do
      deallocate( srcMom)
      deallocate( srcPdf)
      deallocate( refMom)
      deallocate( tgtMom)
      deallocate( tgtPdf)
      deallocate( targetEq)
      ! Initialize the matrix
      deallocate( intp%fromFiner)
      deallocate( intp%fromCoarser)
      deallocate(glob%levelDesc( maxLevel )%depFromCoarser(1)%stencil )
    end do ! iDim

    return errCode

  end subroutine check_fromCoarser3D
  !****************************************************************************
  
  !****************************************************************************
  ! Check the fromFiner interpolation routines 
  ! 2d linear
  ! 3d linear
  subroutine check_fromFiner( errCode )
    !---------------------------------------------------------------------------
    integer, intent(out) :: errCode
    !---------------------------------------------------------------------------
    integer :: minLevel, maxLevel
    type( mus_scheme_type ) :: scheme
    type( interpolation_type ) :: intp, refIntp
    type( mus_field_prop_type ), allocatable :: fieldProp(:)
    type( mus_scheme_layout_type ) :: layout, refLayout
    type( pdf_data_type ), allocatable :: pdf(:)
    type( pdf_global_type ) :: glob
    type( tem_var_system_type ) :: varSys
    integer :: dPos, iVal, nElems, nSourceElems, iElem
    type( grw_intArray_type )   :: ind  
    type( mus_derVarPos_type ), allocatable :: derVarPos(:)
    type( treelmesh_type )    :: tree !< treelm mesh
    real(kind=rk), allocatable :: momentCoefficients(:,:)
    real(kind=rk), allocatable :: srcMom(:,:)
    real(kind=rk), allocatable :: srcPdf(:)
    real(kind=rk), allocatable :: tgtMom(:)
    real(kind=rk), allocatable :: tgtPdf(:), refMom(:), targetEq(:)
    real(kind=rk) :: omega, targetCoord(3), refTgtCoord(3)
    integer :: nDOFs, offset, iChild, iDimn, sourceLevel, targetLevel
    real(kind=rk) :: error, sourceCoord(3)
    !---------------------------------------------------------------------------

    errCode = 0
    minLevel = 3
    maxLevel = 4
    targetLevel = minLevel
    sourceLevel = maxLevel
    omega = 1.8_rk
    ! Set the coordinate of the reference targetElement
    refTgtCoord = (/0.25_rk, 0.25_rk, 0.25_rk/)
    ! Set the time layer to 1, so the interpolation routine always uses
    refIntp%timeLayer = 1
   
    allocate( derVarPos(1))
    allocate( pdf( 1:maxLevel ))
    pdf(:)%nNext = 1
    pdf(:)%nNow  = 1
    allocate( glob%levelDesc( 1:maxLevel ))
    allocate( fieldProp(1))
    fieldProp( 1)%fluid%omLvl( targetLevel ) = omega
    ! See relationship between omegas in Hasert 2012, ECCOMAS Eq. 25 (p.5)
    fieldProp( 1)%fluid%omLvl( sourceLevel ) = 2._rk/(2._rk*(2._rk/omega - 1._rk ) + 1._rk )
    write(logUnit(1),*) 'omega coarse', fieldProp(1)%fluid%omLvl( targetLevel ), &
               'om fine ', fieldProp(1)%fluid%omLvl( sourceLevel )

    allocate( glob%levelDesc( targetLevel )%depFromFiner(1))
    glob%levelDesc( targetLevel   )%intpFromFiner_high%nVals = 1
    allocate( glob%levelDesc( targetLevel  )%intpFromFiner_high%val( 1))

    do iDimn = 2, 3
      write(logUnit(1),*) 'Dimension ', iDimn
      intp = refIntp
      layout = refLayout
      call mus_init_stencil( stencil = layout%stencil, nBC = 0 )
      select case( iDimn)
      case( 3 )
        nSourceElems = 8
        nDOfs = 4 ! 3d polynomial with 10 coefficients
        intp%nMaxSources = 8
        tree%nElems = nSourceElems+1
        ! Initialize the stencil
        call mus_define_d3q19( layout = layout, tree = tree ) 
      case(2)
        nSourceElems = 4
        nDOfs = 3 ! 2d polynomial with 6 coefficients
        intp%nMaxSources = 4
        tree%nElems = nSourceElems+1
        ! Initialize the stencil
        call mus_define_d2q9( layout = layout, tree = tree ) 
      end select
      ! Initialize the source from coarser for the dependencies
      call init( me = glob%levelDesc( targetLevel )%sourceFromFiner, unique = .true. ) 
      do iElem = 1, nSourceElems
        call append( me = glob%levelDesc( targetLevel )%sourceFromFiner, val = iElem, pos = dPos ) 
      end do
      glob%levelDesc( targetLevel   )%intpFromFiner_high%val( 1) = 1
      allocate(glob%levelDesc( targetLevel )%depFromFiner(1)%stencil(3, nSourceElems))
      do iElem = 1, nSourceElems
        glob%levelDesc( targetLevel )%depFromFiner(1)%stencil(:, iElem) = &
        & layout%fStencil%cxDir(:, iElem )
      end do

      call init( me = glob%levelDesc( targetlevel )%depFromFiner(1)%elem )
      call init( me = glob%levelDesc( targetlevel )%depFromFiner(1)%elemBuffer )
      do iElem = 1, nSourceElems
        call append( me = glob%levelDesc( targetlevel )%depFromFiner(1)%elem, &
          &          val = iElem )
        call append( me = glob%levelDesc( targetlevel )%depFromFiner(1)%elemBuffer, &
          &          val = iElem )
      enddo
      ! Indexing array for pointing to a position where the target ghost is stored
      call init( me = ind )
      call append( me = ind, val = 1 ) 
      ! initialize the transformation from and into moment space
      call mus_init_moments( me = layout%moment, layout = layout )
      ! Initialize the polynomial coefficients for each moment
      call set_coefficients( me = momentCoefficients, nDOFs = nDOFs, &
      &  QQ = layout%fStencil%QQ, valMin = 0.3_rk, valmax = 1.5_rk )
      allocate( srcMom(layout%fStencil%QQ ,layout%fStencil%QQ ))
      allocate( srcPdf(layout%fStencil%QQ ))
      allocate( refMom(layout%fStencil%QQ ))
      allocate( tgtMom(layout%fStencil%QQ ))
      allocate( tgtPdf(layout%fStencil%QQ ))
      allocate( targetEq(layout%fStencil%QQ ))
      ! Initialize the matrix
      allocate( intp%fromFiner(   minLevel:maxLevel ))
      allocate( intp%fromCoarser( minLevel:maxLevel ))

      ! Assign initial values
      do iVal = minLevel, maxLevel
        if( iVal == targetLevel ) then ! assign QQ source elements
          nElems = 1
        else ! .. and 1 target element
          nElems = 2**layout%fStencil%nDims
        endif
        ! Allocate and reset the state array
        if( allocated( pdf( iVal )%state )) deallocate( pdf( iVal )%state )
        allocate( pdf( iVal )%state( nElems*layout%fStencil%QQ, 1 ))
        pdf( iVal )%state(:,:) = -1._rk
        ! Allocate the neighbor array. Is it needed?
        if( allocated( pdf( iVal )%neigh )) deallocate( pdf( iVal )%neigh )
        allocate( pdf( iVal )%neigh( nElems*layout%fStencil%QQ ))
        ! Allocate moment buffer with size ( QQ:nSourceELems )
        if( allocated( glob%levelDesc( iVal )%intpBufFromFiner )) &
          deallocate( glob%levelDesc( iVal )%intpBufFromFiner )
        allocate( glob%levelDesc( iVal )%intpBufFromFiner(   &
          &          layout%fStencil%QQ , nSourceElems ))
        ! Offset for all element types ( actually only needed for gfC here)
        glob%levelDesc( iVal )%offset(:,:) = 0
      enddo

      ! Set states for the source elements
      do iElem = 1, nSourceElems
        sourceCoord = 0.25_rk*real( childPosition( iElem, :), kind=rk) + 0.5_rk
        write(logUnit(1),*) 'sourceElem', iElem, 'crd', sourceCoord
        ! Set each moment
        do iVal = 1, layout%fStencil%QQ
          ! Set the source moments one after another.
          ! The order is the same as in Z-curve (sources are the children)
          srcMom( iVal, iElem ) = polynomial( a = momentCoefficients(:,iVal), &
          & coord = sourceCoord, order = 1, nDims = iDimn )
        end do
        ! shear stress:
        ! set the shear stress sigma = -(1-omega/2)*Pi1
        ! But the resulting moment Pi = Pi0 + Pi1 (+ Pi2...) 
        !                  Pi>1 is not available, so neglect.
        ! Pi = Pi0 + sigma/(-(1-omega/2))
        if( layout%fStencil%QQ .eq. 19 ) then
          srcMom(5:10, iElem) = srcMom(5:10, iElem)/(-(1._rk - fieldProp( 1)%fluid%omLvl( sourceLevel )*0.5_rk))
          ! Add the Pi(0) part to the shear stress
          ! pressure to the diagonal components ...
          srcMom( 5:7, iElem )  = srcMom(5:7, iElem ) + srcMom(1, iElem)/3._rk
          ! ... and momentum flux to all components
          srcMom( 5, iELem ) = srcMom( 5, iElem ) + srcMom(2, iElem)*srcMom(2, iElem)/srcMom(1, iElem)
          srcMom( 6, iELem ) = srcMom( 6, iElem ) + srcMom(3, iElem)*srcMom(3, iElem)/srcMom(1, iElem)
          srcMom( 7, iELem ) = srcMom( 7, iElem ) + srcMom(4, iElem)*srcMom(4, iElem)/srcMom(1, iElem)
          srcMom( 8, iELem ) = srcMom( 8, iElem ) + srcMom(2, iElem)*srcMom(3, iElem)/srcMom(1, iElem)
          srcMom( 9, iELem ) = srcMom( 9, iElem ) + srcMom(3, iElem)*srcMom(4, iElem)/srcMom(1, iElem)
          srcMom(10, iELem ) = srcMom(10, iElem ) + srcMom(2, iElem)*srcMom(4, iElem)/srcMom(1, iElem)
        else if( layout%fStencil%QQ .eq. 9 ) then
          srcMom(4:6, iElem) = srcMom(4:6, iElem)/(-(1._rk - fieldProp( 1)%fluid%omLvl( sourceLevel )*0.5_rk))
          ! Add the Pi(0) part to the shear stress
          ! pressure to the diagonal components ...
          srcMom( 4:5, iElem )  = srcMom(4:5, iElem ) + srcMom(1, iElem)/3._rk
          ! ... and momentum flux to all components
          srcMom( 4, iELem ) = srcMom( 4, iElem ) + srcMom(2, iElem)*srcMom(2, iElem)/srcMom(1, iElem)
          srcMom( 5, iELem ) = srcMom( 5, iElem ) + srcMom(3, iElem)*srcMom(3, iElem)/srcMom(1, iElem)
          srcMom( 6, iELem ) = srcMom( 6, iElem ) + srcMom(2, iElem)*srcMom(3, iElem)/srcMom(1, iElem)
        end if
        write(logUnit(1),'(i4,a,19f8.4)') iElem, 'gen srcMom', srcMom(:, iElem )
        ! Transform to pdf and store to the state vector
        srcPdf( : ) = matmul( layout%moment%toPdf%A, srcMom(:, iElem ))
        offset = layout%fStencil%QQ*( iElem -1)
        pdf( sourceLevel )%state( offset+1:offset+layout%fStencil%QQ, &
        pdf(sourceLevel)%nNext ) = srcPdf
      end do

      !------------append variables to variable system------------------------------
      call tem_var_sys_init( sys = varSys )
      call tem_var_append( sys = varSys, varname = 'pdf',     &
        &                  nComponents=layout%fStencil%QQ, pos = dPos )
      call tem_var_sys_mark( sys = varSys, systemName = 'lbm')
      !-----------------------------------------------------------------------------

      targetCoord(:) = [ 0.5_rk, 0.5_rk, 0.5_rk  ]
      !-----------------------------------------------------------------------------
      ! Do the interpolation
      if( layout%fStencil%nDims .eq. 2 ) then
        call fillMyGhostsFromFiner_linMoment2d(                         &
          &                          pdf       = pdf,                   &
          &                          scheme    = scheme,                &
          &                          glob      = glob,                  &
          &                          fieldProp = fieldProp,             &
          &                          level     = targetLevel,           &
          &                          layout    = layout,                &
          &                          intp      = intp,                  &
          &                          ind       = ind,                   &
          &                          varSys    = varSys,                &
          &                          derVarPos = derVarPos  )  
      else
        call fillMyGhostsFromFiner_linMoment3d(                         &
          &                          pdf       = pdf,                   &
          &                          scheme    = scheme,                &
          &                          glob      = glob,                  &
          &                          fieldProp = fieldProp,             &
          &                          level     = targetLevel,           &
          &                          layout    = layout,                &
          &                          intp      = intp,                  &
          &                          ind       = ind,                   &
          &                          varSys    = varSys,                &
          &                          derVarPos = derVarPos  )  
      end if
      !-----------------------------------------------------------------------------

      tgtPdf = pdf( targetLevel )%state( 1:layout%fStencil%QQ, pdf(targetLevel )%nNext )
      ! and now compare the results
      tgtMom( : ) = matmul( layout%moment%toMoments%A, tgtPdf(: ))
      write(logUnit(1),*)  ' Moment results'
      error = 0._rk
      if( layout%fStencil%QQ .eq. 19 ) then
        ! Loop over all the moments
        do iElem = 1, layout%fStencil%QQ
          ! Determine the analytical solution from the given polynomial
          ! Be aware that the higher moments have to be re-scaled!
          refMom(iElem) = polynomial( a = momentCoefficients(:,iElem), &
          & coord = targetCoord, order = 1, nDims = iDimn)
        end do

        ! Extract the shear stress
        tgtMom(5:10)  = get_shearstress3d( tgtMom, fieldProp( 1)%fluid%omLvl( targetLevel ))
        ! And calculate the equilibrium distribution from the macroscopic 
        ! quantities in tMom
        targetEq = getEquilibrium( tgtMom( 1),    & ! density
                        (/ tgtMom( 2)/tgtMom(1), & ! velX
                           tgtMom( 3)/tgtMom(1), & ! velY
                           tgtMom( 4)/tgtMom(1)/), & ! velZ
                           layout )

        ! Set the further higher moments to their equilibrium distribution
        refMom(11 ) = get_moment( targetEq, layout, (/2,1,0/) ) !mxxy
        refMom(12 ) = get_moment( targetEq, layout, (/2,0,1/) ) !mxxz
        refMom(13 ) = get_moment( targetEq, layout, (/1,2,0/) ) !myyx
        refMom(14 ) = get_moment( targetEq, layout, (/0,2,1/) ) !myyz
        refMom(15 ) = get_moment( targetEq, layout, (/1,0,2/) ) !mzzx
        refMom(16 ) = get_moment( targetEq, layout, (/0,1,2/) ) !mzzy
        refMom(17 ) = get_moment( targetEq, layout, (/2,2,0/) ) !mxxyy
        refMom(18 ) = get_moment( targetEq, layout, (/0,2,2/) ) !myyzz
        refMom(19 ) = get_moment( targetEq, layout, (/2,0,2/) ) !mzzxx      
      elseif( layout%fStencil%QQ .eq. 9 ) then
        ! Loop over all the moments
        do iElem = 1, layout%fStencil%QQ
          ! Determine the analytical solution from the given polynomial
          ! Be aware that the higher moments have to be re-scaled!
          refMom(iElem) = polynomial( a = momentCoefficients(:,iElem), &
          & coord = targetCoord, order = 1, nDims = iDimn)
        end do

        ! Extract the shear stress
        tgtMom(4:6)  = get_shearstress2d( tgtMom, fieldProp( 1)%fluid%omLvl( targetLevel ))
        ! And calculate the equilibrium distribution from the macroscopic 
        ! quantities in tMom
        targetEq = getEquilibrium( tgtMom( 1),    & ! density
                        (/ tgtMom( 2)/tgtMom(1), & ! velX
                           tgtMom( 3)/tgtMom(1), & ! velY
                           0._rk/), & ! velZ
                           layout )

        ! Set the further higher moments to their equilibrium distribution
        refMom(7) = get_moment( targetEq, layout, (/2,1,0/) ) !mxxy
        refMom(8) = get_moment( targetEq, layout, (/1,2,0/) ) !mxxz
        refMom(9) = get_moment( targetEq, layout, (/2,2,0/) ) !myyx
      else
        write(logUnit(1),*) 'Interpolation not completely implemented yet. Abort.'
        stop
      end if

      ! Evaluate the error for each moment
      do iElem = 1,layout%fStencil%QQ
        write(logUnit(1),*) 'ref moment', refMom(iElem), 'intp mom', tgtMom( iElem ), &
        &  'E', abs((refMom(iElem) - tgtMom( iElem ))/refMom(iElem))
        error = error + abs(refMom(iELem) - tgtMom( iElem ))
        write(logUnit(1),*) 'error', error
      end do
      if( error .gt. 1.E-12 ) then
        write(logUnit(1),*) 'ERROR!'
        write(logUnit(1),*) 'Error from the fillMyGhostFromFiner  routine:', error
        errCode = max( errCode, 10)
        return errCode
      else
        write(logUnit(1),*) 'finished with SUCCESS.'
        errCode = max( errCode, 0)
      end if
      deallocate( srcMom)
      deallocate( srcPdf)
      deallocate( refMom)
      deallocate( tgtMom)
      deallocate( tgtPdf)
      deallocate( targetEq)
      ! Initialize the matrix
      deallocate( intp%fromFiner)
      deallocate( intp%fromCoarser)
      deallocate(glob%levelDesc( targetLevel )%depFromFiner(1)%stencil )
      call tem_horizontalSpacer(fUnit=logUnit(1))
    end do ! iDim

    write(logUnit(1),*) 'returning error code', errCode
    return errCode

  end subroutine check_fromFiner
  !****************************************************************************

  !****************************************************************************
  ! Get the shear stress from the second order moment Pi
  ! Pi = Pi0 + Pi1 + Pi2 + ...
  ! approximate Pi1 with Pneq = Pi1 + Pi2 + ...
  ! sigma = -(1-omega/2)*(Pi - Pi0)
  function get_shearstress3d( mom, omega ) result( res )
    real(kind=rk), intent(in) :: mom(:)
    real(kind=rk), intent(in) :: omega
    real(kind=rk) :: res(6)

    res( 1) = mom(5) - mom(1)/3._rk - mom(2)*mom(2)/mom(1)
    res( 2) = mom(6) - mom(1)/3._rk - mom(3)*mom(3)/mom(1)
    res( 3) = mom(7) - mom(1)/3._rk - mom(4)*mom(4)/mom(1)
    res( 4) = mom(8)                - mom(2)*mom(3)/mom(1)
    res( 5) = mom(9)                - mom(3)*mom(4)/mom(1)
    res( 6) = mom(10)               - mom(2)*mom(4)/mom(1)
    res = -(1._rk - 0.5_rk*omega)*res
    res = convPrePost( omega )*res
  end function
  !****************************************************************************

  !****************************************************************************
  ! Get the shear stress from the second order moment Pi
  ! Pi = Pi0 + Pi1 + Pi2 + ...
  ! approximate Pi1 with Pneq = Pi1 + Pi2 + ...
  ! sigma = -(1-omega/2)*(Pi - Pi0)
  function get_shearstress2d( mom, omega ) result( res )
    real(kind=rk), intent(in) :: mom(:)
    real(kind=rk), intent(in) :: omega
    real(kind=rk) :: res(3)

    res( 1) = mom(4) - mom(1)/3._rk - mom(2)*mom(2)/mom(1)
    res( 2) = mom(5) - mom(1)/3._rk - mom(3)*mom(3)/mom(1)
    res( 3) = mom(6)                - mom(2)*mom(3)/mom(1)
    res = -(1._rk - 0.5_rk*omega)*res
  end function
  !****************************************************************************

end program mus_interpolation_test
!******************************************************************************
