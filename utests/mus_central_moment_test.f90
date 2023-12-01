! See copyright notice in the COPYRIGHT file.
program mus_central_moment_test
  use env_module,          only: rk, eps
  use tem_general_module,  only: tem_start, tem_finalize, tem_general_type
  use mus_compute_cumulant_module, only: central_moment_split, central_moment

  implicit none

  logical :: error
  real(kind=rk) :: tolerance

  type( tem_general_type ) :: general
  integer :: ii, jj, kk
  real(kind=rk) :: f(-1:1,-1:1,-1:1)
  real(kind=rk) :: ux, uy, uz, rho, inv_rho
  real(kind=rk) :: res_cm, res_split, err_sum

  call tem_start( 'Central Moment Test', 'utest', general )
  error = .false.
  tolerance = eps * 2500._rk
  write(*,"(A,F20.14)") 'tolerance = ', tolerance

  write(*,"(A)") 'Initialize PDF'
  do kk = -1, 1
    do jj = -1, 1
      do ii = -1, 1
        f( ii, jj, kk ) = abs(0.1d0 * ii + 0.2d0 * jj + 0.3d0 * kk * kk)
        write(*,"(A,3I3,A,F6.3)") 'f(',ii,jj,kk,') = ', f( ii, jj, kk )
      end do
    end do
  end do
  f = f / sum(f)

  ! calc velocity
  rho = sum( f )
  inv_rho = 1.0_rk / rho
  ux  = ( sum(f(1,:,:)) - sum(f(-1, :, :)) ) * inv_rho
  uy  = ( sum(f(:,1,:)) - sum(f( :,-1, :)) ) * inv_rho
  uz  = ( sum(f(:,:,1)) - sum(f( :, :,-1)) ) * inv_rho
  write(*,"(A,F6.3)") 'rho = ', rho
  write(*,"(A,F6.3)") 'ux = ', ux
  write(*,"(A,F6.3)") 'uy = ', uy
  write(*,"(A,F6.3)") 'uz = ', uz

  err_sum = 0.0d0
  ! compare two ways of central moments calculation
  do kk = 0,2
    do jj = 0,2
      do ii = 0,2
        res_cm    = central_moment(f, ii,jj,kk, ux, uy, uz)
        res_split = central_moment_split(f, ii,jj,kk, ux, uy, uz)
        write(*,"(A,3I2,2(A,F12.7))") 'moment ', ii, jj, kk, &
          &                          ', normal: ', res_cm, &
          &                          ', split: ', res_split
        err_sum = err_sum + abs( res_cm - res_split )
      end do
    end do
  end do

  if ( err_sum < tolerance ) then
    error = .false.
  else
    error = .true.
  end if

  call tem_finalize(general)
  if (error) then
    write(*,"(A)") 'FAILED'
  else
    write(*,'(A)') 'PASSED'
  end if

contains

end program mus_central_moment_test
