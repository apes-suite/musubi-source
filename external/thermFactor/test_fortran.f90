module eNRTL
  use, intrinsic :: iso_c_binding, only: c_ptr, c_int, c_double, c_float,      &
    &                                    c_int_least8_t, c_int_least64_t,      &
    &                                    c_char, c_loc, c_bool
  
  implicit none

  interface
    function init_enrtl_f90( filename ) bind(c, name='init_enrtl')
    !function init_enrtl_f90( filename ) bind(c, name='_Z10init_enrtlPKc')
      use, intrinsic :: iso_c_binding
      character(kind=c_char), dimension(*) :: filename
      !> Result, indicating the status of encode
      integer(kind=c_int) :: success
    end function init_enrtl_f90
  end interface  

  interface
    subroutine calc_therm_factor( nSpc, Temp, Press, Mole_frac, Therm_factors )&
      & bind(c, name='calc_therm_factor_C')
      use, intrinsic :: iso_c_binding
      integer(kind=c_int), value, intent(in) :: nSpc
      real(kind=c_double), value, intent(in) :: Temp
      real(kind=c_double), value, intent(in) :: Press
      real(kind=c_double), dimension(*), intent(in) :: Mole_frac
      real(kind=c_double), dimension(*), intent(out) :: Therm_factors
    end subroutine calc_therm_factor
  end interface  

  interface
    subroutine calc_ms_diff_matrix(nSpc, Temp, Press, Mole_frac, D_ij_out)     &
      & bind(c, name='calc_ms_diff_matrix_C')
      use, intrinsic :: iso_c_binding
      integer(kind=c_int), value, intent(in) :: nSpc
      real(kind=c_double), value, intent(in) :: Temp
      real(kind=c_double), value, intent(in) :: Press
      real(kind=c_double), dimension(*), intent(in) :: Mole_frac
      real(kind=c_double), dimension(*), intent(out) :: D_ij_out
    end subroutine calc_ms_diff_matrix  
  end interface

end module eNRTL

program test
  use eNRTL
  use, intrinsic :: iso_c_binding, only: c_ptr, c_int, c_double, c_float,      &
    &                                    c_int_least8_t, c_int_least64_t,      &
    &                                    c_char, c_loc, c_null_char

  implicit none

  character(kind=c_char, len=100) :: filename
  integer(kind=c_int) :: success
  integer(kind=c_int) :: nSpc
  real(kind=c_double) :: temp, press
  real(kind=c_double), dimension(:), allocatable :: mole_frac
  real(kind=c_double), dimension(:), allocatable :: therm_factors
  real(kind=c_double), dimension(:), allocatable :: D_ij_out
 
  nSpc = 3
  temp = 298.15
  press = 1.01325e5
  allocate(mole_frac(nSpc))
  allocate(therm_factors(nSpc*nSpc))
  allocate(D_ij_out(nSpc*nSpc))
  mole_frac(1) = 0.9
  mole_frac(2) = 0.05
  mole_frac(3) = 0.05

  filename = 'H2O_NaCl.dat' // C_NULL_CHAR
  write(*,*) 'initializing eNrtl parameter from file: ', trim(filename)
  write(*,*) 'success ', success
  success = init_enrtl_f90(filename)
  write(*,*) 'success ', success
  
  write(*,*) 'calculating therm_factors for nspc ', nSpc
  call calc_therm_factor( nSpc, temp, press, mole_frac, therm_factors )
  write(*,*) 'therm_factors ', therm_factors

  write(*,*) 'Calculating MS diff coeff'
  call calc_ms_diff_matrix( nSpc, temp, press, mole_frac, D_ij_out)
  write(*,*) 'Diff_Coeff', D_ij_out

end program test
