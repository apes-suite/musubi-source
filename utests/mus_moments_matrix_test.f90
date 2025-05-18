! See copyright notice in the COPYRIGHT file.
! This utest program tests the moment conversion matrix
! A vector of pdf is converted to moments, then density and the first order
! moments are checked. And then the moments are converted back to pdf and are
! checked with original pdf vector.
program mus_moments_matrix_test
  use env_module,         only: rk, eps
  use tem_general_module, only: tem_start, tem_finalize, tem_general_type

  use mus_moments_module,       only: mus_init_moments
  use mus_moments_type_module,  only: mus_moment_type
  use mus_scheme_layout_module, only: mus_scheme_layout_type, &
    &                                 mus_define_d3q19
  use mus_scheme_header_module, only: mus_scheme_header_type  

  implicit none

  logical :: error
  real(kind=rk) :: tolerance, max_error

  type( tem_general_type ) :: general
  type( mus_scheme_layout_type ) :: layout
  type( mus_scheme_header_type) :: schemeHeader
  type( mus_moment_type ) :: moment
  real(kind=rk), allocatable :: valOldPDF(:)
  real(kind=rk), allocatable :: valNewPDF(:)
  real(kind=rk), allocatable :: valMoments(:)
  real(kind=rk) :: mx, my, mz, rho ! , Sxx
  real(kind=rk) :: omega, rho0
  integer :: ii, QQ

  call tem_start( 'Moments matrix utest', general )
  error = .false.
  tolerance = eps * 2500._rk
  write(*,*) 'tolerance = ', tolerance

  ! Initialize stencil
  write(*,*) "Initialize stencil D3Q19"
  call mus_define_d3q19( layout = layout, nElems = 1 )
  QQ = layout%fStencil%QQ

  omega = 1.8_rk
  rho0  = 1.0_rk

  schemeHeader%kind = 'fluid'
  schemeHeader%relaxation = 'bgk'
  schemeHeader%layout = 'd3q19'

  ! initialize moment type
  write(*,*) "Initialize moment matrix"
  call mus_init_moments( moment, QQ, layout%fStencil%cxDir, &
    &                    layout%fStencil%label, schemeHeader )

  ! Output matrix ------------------------------------------------------
    write(*, "(A)") ' toMoment matrix:--------------------------'
  do ii = 1, QQ
    write(*, "(19F5.1)") moment%toMoments%A(ii,1:QQ)
  end do
    write(*, "(A)") ''

    write(*, "(A)") ' toPDF    matrix:--------------------------'
  do ii = 1, QQ
    write(*, "(19F6.2)") moment%toPDF%A(ii,1:QQ)
  end do
  ! Output matrix ------------------------------------------------------

  allocate( valOldPdf(  QQ ) )
  allocate( valNewPdf(  QQ ) )
  allocate( valMoments( QQ ) )

  ! Initialize pdf
  write(*,*) "Initialize pdf (oldPDF)"
  do ii = 1, QQ
    valOldPDF( ii ) = real(ii,rk) * 0.1_rk
  enddo

  ! Convert to moments
  write(*,*) "Convert to moments"
  valMoments = matmul( moment%toMoments%A, valOldPDF )

  ! Check zero and first moments
  write(*,*) "Calculating rho, momX, momY, momZ from oldPDF"
  rho = sum( valOldPDF )
  mx  = sum( valOldPDF(:) * layout%fStencil%cxDirRK(1,:) )
  my  = sum( valOldPDF(:) * layout%fStencil%cxDirRK(2,:) )
  mz  = sum( valOldPDF(:) * layout%fStencil%cxDirRK(3,:) )
  ! Sxx = getShearRateTensor_acoustic( valOldPDF, omega, layout, rho0 )

  write(*,*) "Calculating the error of rho, momX, momY, momZ"
  if ( abs( rho - valMoments(1)) > tolerance ) then
    error = .true.
    write(*,*) " Discrepency in density exceeds tolerance!"
    write(*,*) "   Discrepency = ", abs( rho - valMoments(1))
    write(*,*) "   tolerance   = ", tolerance
  else if ( abs( mx - valMoments(moment%first_moments(1)) ) > tolerance ) then
    error = .true.
    write(*,*) " Discrepency in momX exceeds tolerance!"
    write(*,*) "   Discrepency = ", abs( mx - valMoments(moment%first_moments(1)))
    write(*,*) "   tolerance   = ", tolerance
  else if ( abs( my - valMoments(moment%first_moments(2)) ) > tolerance ) then
    error = .true.
    write(*,*) " Discrepency in momY exceeds tolerance!"
    write(*,*) "   Discrepency = ", abs( mx - valMoments(moment%first_moments(2)))
    write(*,*) "   tolerance   = ", tolerance
  else if ( abs( mz - valMoments(moment%first_moments(3)) ) > tolerance ) then
    error = .true.
    write(*,*) " Discrepency in momY exceeds tolerance!"
    write(*,*) "   Discrepency = ", abs( mx - valMoments(moment%first_moments(3)))
    write(*,*) "   tolerance   = ", tolerance
  end if
  ! todo: high order moments can also be checked

  ! convert back to pdf
  write(*,*) "Convert moments back to pdf (newPDF)"
  valNewPDF = matmul( moment%toPdf%A, valMoments)

  ! check originalPdf and newPdf
  write(*,*) "Calculate difference between newPDF and oldPDF"
  max_error = maxval( abs(valNewPDF - valOldPDF) )
  write(*, "(A, ES20.10)") 'max error = ', max_error
  if ( max_error > tolerance ) then
    error = .true.
    do ii = 1, layout%fStencil%QQ
      write(*, "(A5, 2A20)") 'iDir', 'oldPDF', 'newPDF'
      write(*, "(I5, 2ES20.10)") ii, valOldPDF(ii), valNewPDF(ii)
    end do
  end if

  call tem_finalize(general)
  if (.not. error) then
    write(*,'(A)') 'PASSED'
  else
    write(*,'(A)') 'FAILED'
  end if

end program mus_moments_matrix_test
!******************************************************************************
