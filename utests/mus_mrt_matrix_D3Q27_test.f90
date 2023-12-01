! See copyright notice in the COPYRIGHT file.
! This utest program tests the correctnes of M and Minv for WMRT D3Q27
program mus_mrt_matrix_D3Q27_test
  use env_module,         only: rk
  use tem_general_module, only: tem_start, tem_finalize, tem_general_type

  use mus_mrtInit_module, only: check_mrt_matrix_d3q27

  implicit none

  logical :: error = .true.

  type( tem_general_type ) :: general
  integer :: clock


  call tem_start( 'MRT D3Q27 Matrix test', 'utest', general )

  ! initialize f
  CALL SYSTEM_CLOCK( COUNT = clock )

  error = check_mrt_matrix_d3q27()

  call tem_finalize(general)

  if (.not. error) then
    write(*,*) 'PASSED'
  else
    write(*,*) 'FAILED'
  end if

contains

end program mus_mrt_matrix_D3Q27_test
!******************************************************************************
