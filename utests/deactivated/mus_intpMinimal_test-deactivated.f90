program utest_interpolation
  use env_module, only: rk, long_k, rk_mpi, PathLen, LabelLen, minLength
  use tem_stencil_module, only: d3q19_cxDir
  use treelmesh_module
  use mus_interpolate_header_module, only: interpolation_type,  &
    &                                      determine_timeLayer, &
    &                                      refinement_type
  use mus_field_prop_module, only: mus_field_prop_type 
  use mus_pdf_module,    only: pdf_data_type, pdf_global_type, pdf_globInfo_type
  use mus_scheme_layout_module, only: mus_scheme_layout_type, mus_define_d3q19
  use mus_interpolate_quadratic_module, only: mus_init_intpMatrices, &
    &    mus_interpolate_quadratic3d_leastSq
  use tem_polynomial_module, only: polynomial, set_divU, set_coefficients

  implicit none

  type( treelmesh_type ) :: tree
  real(kind=rk), allocatable :: a(:,:) !< coefficients
  real(kind=rk), allocatable :: sourceCoord(:,:)
  type( pdf_global_type ) :: glob           !< global parameters 
  type( mus_field_prop_type ), allocatable ::  fieldProp(:)
  type( mus_scheme_layout_type ) :: layout
  real(kind=rk), allocatable :: sourceMom(:)
  real(kind=rk) :: targetMom, error, limit, rand
  integer :: nSourceElems, iElem, iLevel
  integer :: sourceLevel, targetLevel, nDofs
  type(interpolation_type) :: intp

  nDofs = 10 ! 3d polynomial has 10 dofs 0, x, y,z, x^2, ...
  write(*,*) 'hi'
  limit = 1.E-12
  write(*,*) 'loading stencil d3q19'
  call mus_define_d3q19( layout, tree%nElems )
  nSourceElems = layout%fStencil%QQ
  write(*,*) 'setting inteprolation options'
  intp%compact = .false.
  iLevel = 1
  sourceLevel = iLevel
  targetLevel = iLevel + 1
  allocate( intp%fromFiner(   iLevel   ))
  allocate( intp%fromCoarser( targetLevel ))
  allocate( glob%levelDesc( iLevel:targetLevel))
  allocate( fieldProp(1))
  fieldProp(1)%fluid%omLvl = 1.8_rk
  glob%levelDesc( iLevel + 1 )%intpFromCoarser_high%nVals = 1
  allocate( glob%levelDesc( iLevel + 1 )%intpFromCoarser_high%val( 1))
  glob%levelDesc( iLevel + 1 )%intpFromCoarser_high%val( 1) = 1
  allocate( glob%levelDesc( iLevel + 1 )%depFromCoarser(1))
  glob%levelDesc( iLevel + 1 )%depFromCoarser(1)%elem%nVals = layout%fStencil%QQ
  allocate(glob%levelDesc( iLevel + 1 )%depFromCoarser(1)%stencil(3, nSourceElems))
  do iElem = 1, nSourceElems
    glob%levelDesc( iLevel + 1 )%depFromCoarser(1)%stencil(:, iElem) = &
    & layout%fStencil%cxDir(:, iElem )
  end do
  tree%nElems = 1

  write(*,*) 'testing the scalar least square interpolation in 3d'
  call mus_init_intpMatrices( me = intp%fromCoarser( targetLevel ),        &
    &                            intp = intp,                           &
    &                            layout = layout,                       &
    &                 dep  = glob%levelDesc( targetLevel )%depFromCoarser, &
    &                            glob  = glob,                          &
    &                            fieldProp = fieldProp(1),       &
    &                            compact   = intp%compact,              &
    &                            debug = .false.,&
    &                            targetLevel = targetLevel,                &
    &                            sourceLevel = iLevel )


  write(*,*) 'choosing polynomial'
  !!   f(x) = a0 + a1*x + a2*y + a3*z + a4*x^2 + a5*y^2 + a6*z^2 + a7*x*y 
  !!             + a8*y*z + a9*z*x
  write(*,*) 'setting random coefficients'
  !  a(:) = (/ 1._rk, 2._rk, 3._rk, 4._rk, 5._rk, 6._rk, &
  !  &         7._rk, 8._rk, 9._rk, 10._rk/)
  call set_Coefficients( me = a, nDOFs = nDOFs, &
  &  QQ = layout%fStencil%QQ, valMin = 0.5_rk, valMax = 1.6_rk )
  write(*,*) a(:,:)
  write(*,*) 'getting the coordinates of the source elements'

  allocate( sourceCoord(3, nSourceElems ) )
  allocate( sourceMom( nSourceElems ) )
  do iElem = 1, nSourceElems
    sourceCoord( :, iElem ) = d3q19_cxDir(:, iElem )
    sourceMom(iELem) = polynomial( a = a(:, iElem), coord = sourceCoord(:,iElem ), &
    &             order = 2, nDims = 3 )
  enddo
  error = 0._rk
  do iElem = 1, nSourceElems
    ! Do the interpolation
    targetMom = mus_interpolate_quadratic3d_leastSq(            &
    &     srcMom =    sourceMom(1:nSourceElems),                 &
    &     targetCoord =  sourceCoord(:, iElem), &
    &     nSourceElems =  nSourceElems,                             &
    &     matrix =    intp%fromCoarser( targetLevel )%matrixInv(1)%A )
    write(*,*) 'quantity which must be recovered:', sourceMom(iElem), &
    &   'result from interpolation', targetMom
    write(*,*) 'diff:', targetMom - sourceMom(iElem)
    error = error + abs(targetMom - sourceMom(iElem))
  end do
  write(*,*) 'accumulated error', error
  if( error .gt. limit ) then
    write(*,*) 'ERROR! Error above ', limit
    write(*,*) 'FAILED'
    stop -1
  else
    write(*,*) 'SUCCESS! Error below ', limit
    write(*,*) 'PASSED'
    stop 0
  end if

  contains 

end program utest_interpolation
