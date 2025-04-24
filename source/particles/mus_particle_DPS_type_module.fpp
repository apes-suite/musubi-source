! Copyright (c) 2025 Tristan Vlogman <t.g.vlogman@utwente.nl>
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are met:
!
! 1. Redistributions of source code must retain the above copyright notice,
! this list of conditions and the following disclaimer.
!
! 2. Redistributions in binary form must reproduce the above copyright notice,
! this list of conditions and the following disclaimer in the documentation
! and/or other materials provided with the distribution.
!
! THIS SOFTWARE IS PROVIDED BY THE UNIVERSITY OF SIEGEN “AS IS” AND ANY EXPRESS
! OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
! OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
! IN NO EVENT SHALL UNIVERSITY OF SIEGEN OR CONTRIBUTORS BE LIABLE FOR ANY
! DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
! (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
! LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
! ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
! (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
! SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
! **************************************************************************** !
?? include 'particleArrayMacros.inc'
!> Provides the data type for the discrete particles
module mus_particle_DPS_type_module
  use env_module, only: rk, zeroLength, minLength

  use tem_dyn_array_module, only: dyn_intArray_type
  use tem_grow_array_module, only: grw_intarray_type
  use tem_logging_module, only: logUnit

  use mus_particle_array_module, only: maxContainerSize

  implicit none

  private

  !> Basic particle type for Discrete Particle Simulations
  !! In DPS particles are unresolved i.e. they are not represented
  !! on the grid. Instead, their effect on the flow is considered
  !! using the volume fraction in the fluid_GNS scheme's equilibrium
  !! distribution functions and forcing terms. Hence the scheme kind
  !! should be set to 'fluid_GNS' when using these kinds of particles.
  type mus_particle_DPS_type
    !- Unique particleID for identification across all processes
    integer :: particleID
  
    !- Owner process of this particle
    !  Process is owner if coordOfOrigin is local to this process
    integer :: owner
  
    !> Previous owner of this particle. Used to determine when to 
    !! Send data like particle%momInc to new owner
    integer :: previousOwner
    
    !- Array matching the shape of the send%proc array
    !  the ith entry indicates whether this particle exists on
    !  send%proc(i) and hence we need to send velocity updates there
    logical, allocatable :: existsOnProc(:)
    
    !- Array matching the shape of the send%proc array
    !  the ith entry indicates whether this particle is new on
    !  send%proc(i) in which case we need to send over data
    !  to that proc so it can be added to its particleGroup
    logical, allocatable :: addToProc(:)
    
    !- Array matching the shape of the send%proc array
    !  the ith entry indicates whether this particle is new on
    !  send%proc(i) in which case we need to send over data
    !  to that proc so it can be added to its particleGroup
    logical, allocatable :: removeFromProc(:)
    
    !> Logical which tells us whether to initialize this particle or not
    !! is set to true only immediately after receiving this particle
    !! from a neighboring process
    logical :: newForMe = .FALSE.
  
    !> hasCollided tells us whether particle has just had its velocity modified
    !! in a collision and that this information needs to be sent to other processes
    logical :: hasCollided = .FALSE.
  
    !> removeParticle_global indicates that this particle needs to be removed after e.g.
    !! hitting an open boundary. This information is first sent to all other procs
    !! that know about this particle, then the particle is actually removed from 
    !! the particleGroup.
    logical :: removeParticle_global = .FALSE.
  
    !> removeParticle_local indicates that this particle needs to be removed from 
    !! this process's particleGroup because its coordOfOrigin is not within one 
    !! lattice site of a local fluid cell
    logical :: removeParticle_local = .FALSE.
  
    !- Radius of (spherical) particle
    real(kind=rk) :: radius
  
    ! radius in lattice units, rounded up
    integer :: Rn                   
  
    !- mass and rotational inertia
    real(kind=rk) :: mass
    real(kind=rk) :: rotInertia
  
    !- Particle origin position and translational + angular velocity (x,y,z,rx,ry,rz)
    real(kind=rk) :: vel(6)
    real(kind=rk) :: pos(6)
    real(kind=rk) :: oldPos(6)
    
    ! integer coordinate of desired point
    integer :: coordOfOrigin(4)     
    integer :: oldCoordOfOrigin(4)     
  
    ! position of origin in levelDesc%total list
    ! integer :: posOfOrigin
    integer :: posOfOrigin
    
    ! Hydrodynamic force and torque acting on the particle (fx,fy,fz,mx,my,mz) 
    real(kind=rk) :: F(6) = 0.0_rk
  
    ! Momentum to be transferred to fluid by modification of pdfs over 1 LBM time step 
    ! = the sum of all momentum transferred in the DEM subcycles.
    real(kind=rk) :: Favg(3) = 0.0_rk
  
    ! Fluid velocity interpolated to particle location
    ! For computation of hydrodynamic forces 
    real(kind=rk) :: u_fluid(3) = 0.0_rk
  
    ! Fluid density interpolated to particle location
    ! For computation of hydrodynamic forces 
    real(kind=rk) :: rho_fluid = 0.0_rk
  
    ! Fluid volume fraction interpolated to particle location
    ! For computation of hydrodynamic forces 
    real(kind=rk) :: eps_f_fluid = 0.0_rk
  
    ! Pressure gradient interpolated to particle location
    real(kind=rk) :: grad_p_fluid(3) = 0.0_rk
  
    ! Velocity curl interpolated to particle location
    real(kind=rk) :: curl_u_fluid(3) = 0.0_rk
  
    ! External force vector which is loaded from the lua script.
    ! Can only be constant as of now (for example to simulate gravity)
    real(kind=rk) :: Fext(6) = 0.0_rk 
    
    ! Buffer for the DEM collision force and torque
    ! We need to store forces at last two times for velocity verlet integration
    real(kind=rk) :: F_DEM(2,6) = 0.0_rk
  
    ! Index pointing to row in Fcoll for current DEM time step (t)
    integer :: F_DEM_now = 1
    ! Index pointing to row in Fcoll at next DEM time step (t + dt_dem)
    integer :: F_DEM_next = 2
  
    !> List containing indices of nearby particles to check for 
    !! collisions during the DEM subcycles
    type(dyn_intarray_type) :: DEM_neighborlist
  
    ! For DEM treatment of particle-wall collisions we store the distance
    ! to the nearest wall and the normal vector (pointing towards the wall)
    ! of that wall
    integer :: nWallPos = 0 ! number of elements in wallPosSum
    real(kind=rk) :: avgWallPos(3) = 0.0_rk
    real(kind=rk) :: rwall(3) = 0.0_rk
  
    ! Logical that indicates whether this particle is close enough to a wall
    ! that we need to compute wall interactions during the subcycling loop
    logical :: interactWithWall = .FALSE. 
  
  end type mus_particle_DPS_type

?? copy :: DPA_decltxt(particle_DPS)

  interface allocateProcessMasks
    module procedure allocateProcessMasks_DPS
  end interface

  public :: mus_particle_DPS_type
  public :: allocateProcessMasks
  public :: dyn_particle_DPS_array_type
  public :: init_da_particle_DPS
  public :: destroy_da_particle_DPS
  public :: append_da_particle_DPS
  public :: expand_da_particle_DPS
  public :: truncate_da_particle_DPS
  public :: swap_da_particle_DPS
  public :: remove_particle_from_da_particle_DPS
  public :: sortposofval_particle_DPS


contains


! ---- Dynamic particle array methods for Discrete Particle Simulations ---- !
?? copy :: DPA_inittxt(particle_DPS)
?? copy :: DPA_destroytxt(particle_DPS)
?? copy :: DPA_appendtxt(particle_DPS)
?? copy :: DPA_expandtxt(particle_DPS)
?? copy :: DPA_swaptxt(particle_DPS) 
?? copy :: DPA_removetxt(particle_DPS)
?? copy :: DPA_truncatetxt(particle_DPS)
?? copy :: DPA_sortposofval_txt(particle_DPS)

  subroutine allocateProcessMasks_DPS( particle, nProcs )
    !> Particle to initialize
    type(mus_particle_DPS_type), intent(inout) :: particle
    !> Number of processes to communicate particle data with
    integer :: nProcs
    ! -----------------------------------------------!
    ! Allocate space for the existsOnProc mask which tells us on which other procs
    ! this particle lives at the current time step
    if( allocated(particle%existsOnProc) ) deallocate( particle%existsOnProc )
    allocate( particle%existsOnProc( nProcs ) )
    particle%existsOnProc( 1:nProcs ) = .FALSE.
    
    ! addToProc is used to determine whether to send over data needed to add this
    ! particle to the receiving process's particle group
    if( allocated(particle%addToProc) ) deallocate( particle%addToProc )
    allocate( particle%addToProc( nProcs ) )
    particle%addToProc( 1:nProcs ) = .FALSE.
    
    ! removeFromProc is used to determine whether to send over the signal that
    ! a particle needs to be removed from the receiving proc's particle group
    if( allocated(particle%removeFromProc) ) deallocate( particle%removeFromProc )
    allocate( particle%removeFromProc( nProcs ) )
    particle%removeFromProc( 1:nProcs ) = .FALSE.
  
  end subroutine allocateProcessMasks_DPS 

end module mus_particle_DPS_type_module
