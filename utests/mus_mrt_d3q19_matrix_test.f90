! See copyright notice in the COPYRIGHT file.
! This utest program tests the compact quadratic velocity interpolation 2D
program mus_mrt_d3q19_matrix_test
  use env_module,         only: rk, eps, init_random_seed
  use tem_general_module, only: tem_start, tem_finalize, tem_general_type
  use tem_logging_module, only: logUnit

  use mus_mrtInit_module, only: MMtrD3Q19, MMivD3Q19

  implicit none

  logical :: error = .true.
  real(kind=rk) :: tolerance

  type( tem_general_type ) :: general
  integer :: clock
  integer :: iDir

  real(kind=rk) :: f(19)
  real(kind=rk) :: m_matrix(19)
  real(kind=rk) :: m_explicit(19)
  real(kind=rk) :: Imatrix(19,19)

  call tem_start( 'MRT D3Q19 Matrix test', 'utest', general )
  tolerance = eps * 2500._rk
  write(*,*) 'tolerance = ', tolerance

  ! initialize f
  CALL SYSTEM_CLOCK( COUNT = clock )
  call init_random_seed( clock )
  do iDir = 1, 19
  call random_number( f(iDir) )
  end do
  ! f(:) = [ (iDir, iDir = 1, 19) ]

  ! calculate moment by m_matrix   = f * matrix
  do iDir = 1, 19
    m_matrix(iDir) = sum( f(1:19) * MMtrD3Q19(iDir,1:19) )
  end do

  ! caluclate moment by m_explicit = subroutine(f)
  m_explicit(1:19) = get_moment_d3q19( f(1:19) )

  ! compare differece between these two
  call check_error( m_matrix(:), m_explicit(:), 19, tolerance, error )

  ! multiply transformation matrix by its inverse to check if gives identity.
  Imatrix = matmul( MMtrD3Q19, MMIvD3Q19 )
  do iDir = 1, 19
    write(logUnit(1), "(19F7.2)") Imatrix( iDir, : )
  end do

  call tem_finalize(general)

  if (.not. error) then
    write(*,*) 'PASSED'
  else
    write(*,*) 'FAILED'
  end if

contains

  subroutine check_error( a, b, nVals, tolerance, error )
    !---------------------------------------------------------------------------
    integer, intent(in) :: nVals
    real(kind=rk), intent(in) :: a(nVals)
    real(kind=rk), intent(in) :: b(nVals)
    real(kind=rk), intent(in) :: tolerance
    logical, intent(out) :: error
    !---------------------------------------------------------------------------
    real(kind=rk) :: diff(nVals), max_error
    integer :: iVal
    !---------------------------------------------------------------------------

    write( logUnit(1), *) 'Calculating errors between array A and B:'
    diff(:)   = abs(a(:) - b(:))
    max_error = maxval( diff(:) )
    if ( max_error > tolerance ) then
      error = .true.
      write(logUnit(1), *) "Max error: ", max_error
      write(logUnit(1), *) "tolerance: ", tolerance
      write(logUnit(1), *) "Max error exceeds tolerance!"
      write(logUnit(1), *) "  a                 b               diff"
      do iVal = 1, nVals
        if ( diff(iVal) > tolerance ) &
          & write(logUnit(1),*) iVal, a(iVal), b(iVal), diff(iVal)
      end do
    else
      error = .false.
      write(logUnit(1), *) "Max error: ", max_error
      write(logUnit(1), *) "tolerance: ", tolerance
      write(logUnit(1), *) "Max error within tolerance!"
      do iVal = 1, nVals
        write(logUnit(1),*) iVal, a(iVal), b(iVal), diff(iVal)
      end do
    end if

  end subroutine check_error

  function get_moment_d3q19( f ) result ( m )

    real(kind=rk), intent(in) :: f(19)
    real(kind=rk) :: m(19)

    m(1) =   f(18) + f( 5) + f(16) &
      &    + f( 1) + f(15) + f( 2) &
      &    + f(17) + f( 4) + f( 6) &
      &    + f(14) + f(10) + f(13) &
      &    + f( 8) + f( 3) + f(12) &
      &    + f( 9) + f(11) + f( 7) &
      &    + f(19)

    m(2) =   -f(19) + f(18)+ f(15)+ f(17) &
      &     + f(16) + f(14)+ f(11)+ f(12) &
      &     + f(13) + f(10)+ f( 7)+ f( 9) &
      &     + f( 8)

    !epsilon
    m(3) =  f(19)- 2._rk*(f( 4)+ f( 1)+ f( 5) &
      &     + f( 2)+ f( 6)+ f( 3))+ f(18)     &
      &     + f(15)+ f(17)                    &
      &     + f(16)+ f(14)+ f(11)+ f(12)      &
      &     + f(13)+ f(10)+ f( 7)+ f( 9)      &
      &     + f( 8)

    !jx
    m(4) =  f( 4)- f( 1)+ f(18)          &
      &     - f(15)+ f(17)               &
      &     - f(16)+ f(14)- f(11)+ f(12) &
      &     - f(13)

    !qx
    m(5) =  -2.0_rk*( f( 4)- f( 1))+f(18) &
      &     - f(15)+ f(17)                &
      &     - f(16)+ f(14)- f(11)+ f(12)  &
      &     - f(13)
    !jy
    m(6) =  f( 5)- f( 2)+ f(18)           &
      &     - f(15)- f(17)                &
      &     + f(16)+ f(10)- f( 7)+ f( 9)  &
      &     - f( 8)
    !qy
    m(7) =  -2.0_rk*(f( 5)- f( 2))+ f(18) &
      &     - f(15)- f(17)                &
      &     + f(16)+ f(10)- f( 7)+ f( 9)  &
      &     - f( 8)
    !jz
    m(8) =  f( 6)- f( 3)                   &
      &     + f(14)- f(11)- f(12)  + f(13) &
      &     + f(10)- f( 7)- f( 9)  + f( 8)
    !qz
    m(9) =  -2.0_rk*(f( 6)- f( 3))         &
      &     + f(14)- f(11)- f(12)  + f(13) &
      &     + f(10)- f( 7)- f( 9)  + f( 8)
    !3pxx
    m(10) = 2.0_rk*(f( 4)+ f( 1))-(f( 5)         &
      &      + f( 2)+ f( 6)+ f( 3))              &
      &      + f(18)+ f(15)+ f(17)+ f(16)        &
      &      + f(14)+ f(11)+ f(12)+ f(13)        &
      &      -2.0_rk* (f(10)+ f( 7)+ f( 9)+ f( 8))
    !3pixx
    m(11) = -2._rk*(f( 4)+ f( 1))+ f( 5)         &
      &      + f( 2)+ f( 6)+ f( 3)               &
      &      + f(18)+ f(15)+ f(17)+ f(16)        &
      &      + f(14)+ f(11)+ f(12)+ f(13)        &
      &      -2.0_rk*( f(10)+ f( 7)+ f( 9)+ f( 8))
    !pww
    m(12) = f( 5)                         &
      &     + f( 2)- f( 6)- f( 3)         &
      &     + f(18)+ f(15)+ f(17)+ f(16)  &
      &     -(f(14)+ f(11)+ f(12)+ f(13))
    !piww
    m(13) = - f( 5)                       &
      &     - f( 2)+ f( 6)+ f( 3)         &
      &     + f(18)+ f(15)+ f(17)+ f(16)  &
      &     -(f(14)+ f(11)+ f(12)+ f(13))
    !pxy
    m(14) = f(18)+ f(15)- f(17)- f(16)
    !pyz
    m(15) = f(10)+ f( 7)- f( 9)- f( 8)
    !pxz
    m(16) = f(14)+ f(11)- f(12)- f(13)
    !mx
    m(17) = f(18)- f(15)+ f(17)- f(16)  &
      &      - f(14)+ f(11)- f(12)+ f(13)
    !my
    m(18) = -f(18)+ f(15)+ f(17)- f(16) &
      &      + f(10)- f( 7)+ f( 9)- f( 8)
    !mz
    m(19) = +f(14)- f(11)- f(12)+ f(13) &
      &      - f(10)+ f( 7)+ f( 9)- f( 8)

  end function get_moment_d3q19

end program mus_mrt_d3q19_matrix_test
!******************************************************************************
