program mus_interpolation_test
  use aotus_module, only: flu_State, aot_get_val
  use aot_table_module, only: aot_table_open, aot_table_close, aot_table_length

  use env_module, only: rk, labelLen
  use treelmesh_module, only: treelmesh_type
  use tem_aux_module, only: tem_abort
  use tem_element_module, only: eT_ghostFromCoarserElem,                 &
    &                                eT_ghostFromFinerElem
  use tem_geometry_module, only: tem_baryOfId, tem_intp_bilinear
  use tem_grow_array_module, only: grw_intArray_type, init, append
  use tem_dyn_array_module, only: dyn_intArray_type, init, append
  use tem_param_module, only: cs2inv, cs2
  use tem_tools_module, only: tem_write
  use tem_var_system_module, only: tem_var_system_type, tem_var_append, tem_var_sys_init
  use tem_variable_module, only: tem_variable_type, tem_variable_out,  &
    &                            tem_variable_from, tem_variable_dump, &
    &                            assignment(=), dyn_varArray_type, init, &
    &                            append, PositionOfVal

  use mus_field_prop_module, only: mus_field_prop_type 
  use mus_param_module,  only: mus_param_type
  use mus_interpolate_header_module, only: interpolation_type,                 &
                                           determine_timeLayer
  use mus_interpolate_quadratic_module, only:                                  &
          fillFinerGhostsFromMe_quadraticCompact2d,                            &
          fillMyGhostsFromFiner_quadraticCompact2d,                            &
          fillFinerGhostsFromMe_quadraticCompact3d,                            &
          fillMyGhostsFromFiner_quadraticCompact3d,                            &
          fillFinerGhostsFromMe_quadratic2d,                           &
          fillFinerGhostsFromMe_quadratic3d,                           &
          fillFinerGhostsFromMe_quadraticReduced3d,                    &
          mus_init_intpMatrices
  use mus_moments_module, only: mus_init_moments
  use mus_pdf_module,    only: pdf_data_type, pdf_global_type, pdf_globInfo_type
  use mus_scheme_layout_module, only: mus_scheme_layout_type
  use mus_moments_module, only: init_transformation_matrix
  use mus_variable_module, only: mus_derVarPos_type
  use mus_scheme_type_module, only: mus_scheme_type
  use mus_scheme_layout_module, only: mus_scheme_layout_type, mus_load_newLayout,&
    &                          mus_weights_out, mus_define_d3q19

  implicit none

  write(*,*) 'Running mus_interpolation_test...'
  call check_intp_quadratic

contains

  subroutine check_intp_quadratic()

    integer :: depNeigh(16)
    integer :: nDepNeighs
    integer :: minLevel, maxLevel
    logical :: failure
    type( interpolation_type ) :: intp
    type( mus_field_prop_type ), allocatable :: fieldProp(:)
    type( mus_scheme_layout_type ) :: layout
    type( mus_scheme_type ) :: scheme
    type( pdf_data_type ), allocatable :: pdf(:)
    type( pdf_global_type ) :: glob
    type( tem_var_system_type ) :: varSys
    integer :: dPos, iVal, nElems, nSourceElems, iElem
    type( grw_intArray_type )   :: ind  
    type( mus_derVarPos_type ), allocatable :: derVarPos(:)
    type( treelmesh_type )    :: tree !< treelm mesh
    type( tem_variable_type ) :: tVar
    real(kind=rk), allocatable :: refVal(:), val(:), error
    real(kind=rk), allocatable :: mult_refVal(:), mult_val(:), error2
    integer :: errCode

    minLevel = 3
    maxLevel = 4
    nSourceElems = 19

    intp%timeLayer = 1
    allocate( pdf( 1:maxLevel ))
    allocate( glob%levelDesc( 1:maxLevel ))
    allocate( fieldProp(1))
    allocate( derVarPos(1))
    pdf(:)%nNext = 1
    pdf(:)%nNow  = 1
    tree%nElems = 1
    ! Initialize the source from coarser for the dependencies
    call init( me = glob%levelDesc( maxLevel )%sourceFromCoarser, unique = .true. ) 
    call append( me = glob%levelDesc( maxLevel )%sourceFromCoarser, val = 1, pos = dPos ) 
    ! Indexing array for pointing to a position where the ghost is stored
    call init( me = ind )
    call append( me = ind, val = 1, pos = dPos ) 
    ! Initialize the stencil
    call mus_define_d3q19( layout = layout, nElems = tree%nElems )

    ! Assign initial values
    do iVal = minLevel, maxLevel
      if( iVal == minLevel ) then ! assign QQ source elements
        nElems = layout%fStencil%QQ
      else ! .. and 1 target element
        nElems = 1
      endif
      ! Allocate and reset the state array
      allocate( pdf( iVal )%state( nElems*layout%fStencil%QQ, 1 ))
      pdf( iVal )%state(:,:) = -1._rk
      ! Allocate the neighbor array. Is it needed?
      allocate( pdf( iVal )%neigh( nElems*layout%fStencil%QQ ))
      ! Assign initial values to 
      pdf( iVal )%state(:,1) = (/  1.0_rk,       1.2_rk, &
    1.04_rk,       1.05_rk,       1.03_rk, &
    0.99_rk,       1.055_rk,        1.08_rk, &
    0.86_rk,        0.98_rk,      0.94_rk, &
    0.79_rk,        0.87_rk,      0.95_rk, &
    0.99_rk,        0.65_rk,      1.04_rk, &
    1.01_rk,       1.03_rk  /)
      ! Allocate moment buffer
      allocate( glob%levelDesc( iVal )%intpBufFromCoarser(   &
        &          layout%fStencil%QQ ,layout%fStencil%QQ ))
      glob%levelDesc( iVal )%offset = 0
      allocate( glob%levelDesc( iVal )%depFromCoarser(1))
      glob%levelDesc( iVal )%depFromCoarser(1)%coord = (/0.25_rk, 0.25_rk, 0.25_rk/)

      call init( me = glob%levelDesc( iVal )%depFromCoarser(1)%elem )
      do iElem = 1, nSourceElems
        call append( me = glob%levelDesc( iVal )%depFromCoarser(1)%elem, &
          &          val = iElem, pos = dPos )
      enddo
      fieldProp( 1)%fluid%omLvl( iVal ) = 1.8
    enddo

    call tem_var_sys_init( sys = varSys )
    call tem_var_append( sys = varSys, varname = 'pdf',     &
      &                  nComponents=19, pos = dPos )

    call mus_init_moments( me = layout%moment,  &
                           layout = layout      )

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

    ! Compare result values
    allocate( val( nElems*layout%fStencil%QQ ))
    allocate( refVal( nElems*layout%fStencil%QQ ))
    allocate( mult_val( nElems*layout%fStencil%QQ ))
    allocate( mult_refVal( nElems*layout%fStencil%QQ ))
    refVal = pdf( minLevel )%state( 1:layout%fStencil%QQ, 1)
    val = pdf( maxLevel )%state( :, 1)
    write(*,*) 'Fill Finer Ghosts from Me routine' 
    write(*,*) 'input' 
    write(*,*) refVal
    write(*,*) 'Finer ghost result' 
    write(*,*) val
    mult_val(1:19) = matmul( layout%moment%toMoments%A, pdf( maxLevel )%state( 1:19, 1 ))
    mult_refVal(1:19) = (/    18.555000000000003_rk,      0.39000000000000035_rk,&
    -0.83499999999999974_rk,      0.18500000000000016_rk,       9.2900000000000009_rk,&
    9.8950000000000014_rk,       9.5549999999999997_rk,      0.30999999999999994_rk,&
    9.4999999999999862E-002_rk, 0.22999999999999987_rk,     -0.36999999999999988_rk,&
    8.9999999999999969E-002_rk, 0.40999999999999992_rk,      0.14500000000000013_rk,&
    -7.0000000000000062E-002_rk,-0.29499999999999993_rk,       3.6900000000000004_rk, &
    3.9749999999999996_rk,       3.5499999999999998_rk  /)
    write(*,*) 'Matrix-Vector Multiplication ' 
    write(*,*) 'input' 
    write(*,*) mult_refVal
    write(*,*) 'Finer ghost result' 
    write(*,*) mult_val
    write(*,*) 'result'
!    write(*,'(19f8.4)') pdf( maxLevel )%state( 1:19, 1 )
    error = sum( (refVal(:)-val(:))**2) 
    error2 = sum( (mult_refVal(:)-mult_val(:))**2) 
    write(*,*) 'Error norm', error

    if( error .gt. 1.E-15 .or. error2 .gt. 1.E-15) then
      write(*,*) 'ERROR!'
      write(*,*) 'Error from the fillFinerGhostsFromMe routine:', error
      write(*,*) 'Error from the matrix-vector multiplication: ', error2
      errCode = 10
    else
      write(*,*) 'finished with SUCCESS.'
      errCode = 0
    end if
    call exit(errCode )

  end subroutine check_intp_quadratic
  
  ! Evaluate a 3d polynomial with given coefficients at a coordinate 
  function polynomial( a, coord ) result( res )
    real(kind=rk), intent(in) :: a(:) !< coefficients
    real(kind=rk), intent(in) :: coord(3)
    real(kind=rk) :: res
    real(kind=rk) :: x,y,z

    x = coord(1)
    y = coord(2)
    z = coord(3)

    res = a(1) + a(2)*x + a(3)*y + a(4)*z + a(5)*x**2 + a(6)*y**2 &
    &   + a(7)*z**2 + a(8)*x*y + a(9)*y*z + a(10)*z*x

  end function


end program mus_interpolation_test
