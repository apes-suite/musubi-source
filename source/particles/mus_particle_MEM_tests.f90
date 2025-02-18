module mus_particle_MEM_tests_module
  ! Use statements go here
use env_module,                        only : rk, long_k, stdOutUnit, newUnit
use tem_geometry_module,               only : tem_CoordOfReal
use mus_geom_module,                   only : mus_geom_type
use mus_scheme_type_module,            only : mus_scheme_type
use mus_param_module,                  only : mus_param_type
use mus_particle_type_module,          only : mus_particle_MEM_type, &
                                            & mus_particle_group_type
use mus_particle_MEM_module,           only : updateParticleOwner
use mus_particle_comm_type_module,     only : mus_particles_communication_type 

implicit none

contains

! All tests in this module should be run with the test case in the lua file 
! at /utests_particles/env/musubi.lua
subroutine runtests_MEM(particleGroup, scheme, geometry, params)
  !> particle group
  type(mus_particle_group_type), intent(inout) :: particleGroup
  !> Scheme for access to fluid data
  type(mus_scheme_type), intent(inout) :: scheme
  !> Geometry for access to tree
  type(mus_geom_type), intent(in) :: geometry
  !> Params for access to dt, dx, etc.
  type(mus_param_type), intent(in) :: params
  !> Communication data type for particles
  type(mus_particles_communication_type) :: comm
  !-----------------------------------------!
  logical :: test_updateParticleOwner_success = .FALSE.
  integer :: myRank
  !-----------------------------------------!
  myRank = params%general%proc%rank

  call test_updateParticleOwner( scheme   = scheme,                          &
                               & geometry = geometry,                        &
                               & myRank   = myRank,                          &
                               & comm     = comm,                            &
                               & passed   = test_updateParticleOwner_success )

  if( test_updateParticleOwner_success ) then
    write(stdOutUnit, *) "test_updateParticleOwner SUCCEEDED"
  else
    write(stdOutUnit, *) "test_updateParticleOwner FAILED"
  end if


end subroutine runtests_MEM

! Routine to test the updating of particle owner.
subroutine test_updateParticleOwner(scheme, geometry, myRank, comm, passed)
  !> Scheme for access to fluid data
  type(mus_scheme_type), intent(inout) :: scheme
  !> Geometry for access to tree
  type(mus_geom_type), intent(in) :: geometry
  !> This process' rank
  integer, intent(in) :: myRank
  !> Communication data type for particles
  type(mus_particles_communication_type) :: comm
  !> Logical to indicate whether test passed
  logical, intent(out) :: passed
  ! -------------------------------------- !
  ! Mock particle to update owner of
  type(mus_particle_MEM_type) :: particle
  integer :: coordOfOrigin(4)
  integer :: lev
  ! -------------------------------------- !
  ! Set passed to true initially
  ! Will be set to false as soon as one of the 
  ! tests fails
  passed = .TRUE.
  lev = geometry%tree%global%maxLevel
  ! Initialize particle position on process 0
  particle%pos(1:3) = (/ 32.0, 16.0, 16.0 /) 

  ! Set particle owner to this rank for now
  particle%owner = myRank

  particle%coordOfOrigin = tem_CoordOfReal( mesh  = geometry%tree, &
                                          & point = particle%pos(1:3), &
                                          & level = lev            )

  call updateParticleOwner( this = particle, &
                          & scheme = scheme, &
                          & geometry = geometry, &
                          & myRank = myRank, &
                          & procs = comm%proc, &
                          & Nprocs = comm%nProcs )

  if( particle%owner /= 0 ) then
      passed = .FALSE.
  end if

  ! Set particle position on process 1
  particle%pos(1:3) = (/ 32.0, 48.0, 16.0 /) 
  particle%coordOfOrigin = tem_CoordOfReal( mesh  = geometry%tree, &
                                          & point = particle%pos(1:3), &
                                          & level = lev            )

  call updateParticleOwner( this = particle, &
                          & scheme = scheme, &
                          & geometry = geometry, &
                          & myRank = myRank, &
                          & procs = comm%proc, &
                          & Nprocs = comm%nProcs )

  if( particle%owner /= 1 ) then
      passed = .FALSE.
  end if

  ! Set particle position on process 2
  particle%pos(1:3) = (/ 32.0, 16.0, 48.0 /) 
  particle%coordOfOrigin = tem_CoordOfReal( mesh  = geometry%tree, &
                                          & point = particle%pos(1:3), &
                                          & level = lev            )

  call updateParticleOwner( this = particle, &
                          & scheme = scheme, &
                          & geometry = geometry, &
                          & myRank = myRank, &
                          & procs = comm%proc, &
                          & Nprocs = comm%nProcs )

  if( particle%owner /= 2 ) then
      passed = .FALSE.
  end if

  ! Set particle position on process 3
  particle%pos(1:3) = (/ 32.0, 48.0, 48.0 /) 
  particle%coordOfOrigin = tem_CoordOfReal( mesh  = geometry%tree, &
                                          & point = particle%pos(1:3), &
                                          & level = lev            )

  call updateParticleOwner( this = particle, &
                          & scheme = scheme, &
                          & geometry = geometry, &
                          & myRank = myRank, &
                          & procs = comm%proc, &
                          & Nprocs = comm%nProcs )

  if( particle%owner /= 3 ) then
      passed = .FALSE.
  end if

end subroutine test_updateParticleOwner

end module mus_particle_MEM_tests_module 