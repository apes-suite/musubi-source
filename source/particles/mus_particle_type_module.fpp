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
!> mus_particle_type_module contains the main particle data types for LBM-DEM
!! simulations of solid particles in the flow

?? include 'particleArrayMacros.inc'
module mus_particle_type_module

  use env_module, only: rk, long_k, zeroLength, &
    &   stdOutUnit, minLength, labelLen

  use tem_stencil_module, only: tem_stencilHeader_type, &
    &                           tem_stencil_findIndexOfDir
  use tem_topology_module, only: tem_idOfCoord
  use tem_dyn_array_module, only:  &
    &   init, append, destroy,     &
    &   empty, dyn_intArray_type,  &
    &   dyn_longArray_type,        &
    &   PosOfVal_long,             &
    &   SortPosOfVal_long
  use tem_grow_array_module, only:      &
    &   init, append, destroy, empty,   &
    &   grw_int2darray_type,            &
    &   grw_logical2darray_type,        &
    &   grw_intarray_type,              &
    &   grw_longarray_type,             &
    &   grw_real2darray_type
  
  use mus_geom_module, only: mus_geom_type
  use mus_scheme_type_module, only: mus_scheme_type
  use mus_param_module, only: mus_param_type
  use mus_particle_comm_type_module, only: mus_particles_communication_type
  use mus_particle_interpolation_module, only: mus_particle_interpolator_type

  implicit none

  ! --------- INTERFACE BLOCKS FOR GENERIC PROCEDURES -------- !
  ! This is used so that names for similar procedures for different kinds of 
  ! particles can be referenced using the same name 
  
  interface allocateProcessMasks
    module procedure allocateProcessMasks_MEM
    module procedure allocateProcessMasks_DPS
  end interface
  ! ----------------------------------------------------------!

  ! use this logUnit for writing to file containing particle velocities for
  ! debug
  integer, parameter, public :: particleLogUnit = 55  
  ! max container size for particle dynamic array
  integer, public :: maxContainerSize = 10000


  ! Basic particle type
  type mus_particle_MEM_type
    !- Unique particleID for identification across all processes
    integer :: particleID
  
    !- Owner process of this particle
    !  Process is owner if coordOfOrigin is local to this process
    integer :: owner
    
    !> Process who was owner in last time step. We need this for the 
    !! averaging of forces over two time steps
    integer :: previousOwner = -1
    
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
  
    !> removeParticle indicates that this particle needs to be removed after e.g.
    !! hitting an open boundary. This information is first sent to all other procs
    !! that know about this particle, then the particle is actually removed from 
    !! the particleGroup.
    logical :: removeParticle_global = .FALSE.
  
    !- Radius of (spherical) particle
    real(kind=rk) :: radius
  
    !- mass and rotational inertia
    real(kind=rk) :: mass
    real(kind=rk) :: rotInertia
  
    ! radius in lattice units, rounded up
    integer :: Rn                   
    
    !- Particle origin position and translational + angular velocity (x,y,z,rx,ry,rz)
    real(kind=rk) :: vel(6)
    real(kind=rk) :: pos(6)
    real(kind=rk) :: oldPos(6)
    
    ! integer coordinate of desired point
    integer :: coordOfOrigin(4)     
    integer :: oldCoordOfOrigin(4)     
  
  
    ! Hydrodynamic force and torque acting on the particle (fx,fy,fz,mx,my,mz) 
    ! F is the average of hydrodynamic force computation in current 
    ! and last time step F = 0.5 * ( Fbuff(1,:) + Fbuff(2,:) )
    real(kind=rk) :: F(6) = 0.0_rk
    real(kind=rk) :: Fbuff(2,6) = 0.0_rk         
  
    ! External force vector which is loaded from the lua script.
    ! Can only be constant as of now (for example to simulate gravity)
    real(kind=rk) :: Fext(6) = 0.0_rk 
    
    ! Index pointing to which row of Fbuff to fill at this time step 
    integer :: Fnow = 1                 
    ! Index pointing to which row of Fbuff contains force data of last time step
    integer :: Flast = 2                 
  
    ! Buffer for the DEM collision force and torque
    ! We need to store forces at last two times for velocity verlet integration
    real(kind=rk) :: F_DEM(2,6) = 0.0_rk
  
    ! Index pointing to row in Fcoll for current DEM time step (t)
    integer :: F_DEM_now = 1
    ! Index pointing to row in Fcoll at next DEM time step (t + dt_dem)
    integer :: F_DEM_next = 2
  
    ! For DEM treatment of particle-wall collisions we store the distance
    ! to the nearest wall and the normal vector (pointing towards the wall)
    ! of that wall
    integer :: nWallPos = 0 ! number of elements in wallPosSum
    real(kind=rk) :: avgWallPos(3) = 0.0_rk
    real(kind=rk) :: rwall(3) = 0.0_rk
  
    ! Logical that indicates whether this particle is close enough to a wall
    ! that we need to compute wall interactions during the subcycling loop
    logical :: interactWithWall = .FALSE. 
  
    !- Dynamic array containing indices of the currently covered fluid elements
    !- These need to be excluded in the loop over elements in the kernel
    !- Pertains to levelDesc total list (from which kernel lists are generated) 
    type(dyn_intArray_type) :: exclusionList 
    
    !> Buffer for exclusion list used in moveParticle routine 
    !  Used to determine newly uncovered fluid neighbors
    type(dyn_intArray_type) :: exclusionListBuffer 
  
    !> Number of fluid neighbors for this particle
    integer :: NfluidNeighbors
  
    !> Indices in levelDesc total list of elements that 
    !  need to be turned to fluid after moving particle
    type(grw_intArray_type) :: makeFluidList 
  
  end type mus_particle_MEM_type

?? copy :: DPA_decltxt(particle_MEM)

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

  !> Data type representing a collection of particles, typically all 
  !  particles on a single process/rank
  type mus_particle_group_type
      ! holds the particle objects
      ! particles(1:nLocal) holds local particles
      ! particles(nLocal+1:nParticles) are remote particles
      
      ! --- FOR MOMENTUM-EXCHANGE PARTICLES --- !
      ! Dynamic array holding particle objects in this group
      type(dyn_particle_mem_array_type) :: particles_MEM
  
      ! --------------------------------------- !
  
      ! --- FOR DPS PARTICLES --- !
      ! Dynamic array holding particle DPS objects in this group
      type(dyn_particle_DPS_array_type) :: particles_DPS
  
      !> Object containing information on interpolation boundaries 
      !! and weight functions to interpolate and distribute fluid 
      !! properties from lattice to particle and vice-versa.
      type(mus_particle_interpolator_type) :: interpolator
      ! ------------------------- !
  
      !> Kind of particles to simulate. At present this can be
      !! * Momentum-exchange Method (particle_kind = 'MEM')
      character(len=labelLen) :: particle_kind
  
      ! Number of particles in this group
      integer :: nParticles
    
      ! ---- DATA FOR PARTICLE COLLISIONS ---- !
      !> Logical set to TRUE if particle collision are enabled, false otherwise
      logical :: enableCollisions = .FALSE.
  
      !> Collision time
      real(kind=rk) :: collision_time
  
      !> Threshold gap at which to consider two particles colliding
      real(kind=rk) :: collision_tol
  
      !> Density to use for hydrodynamic force computation when local density is not available
      !! This is in lattice units!
      real(kind=rk) :: rho0_lat
  
      !> Number of DEM subcycles to execute per LBM time step
      integer :: Nsubcycles
  
      ! ---- Do we need this or is this just old stuff? --- !
      ! Integer array that tells us whether particles should interact with 
      ! a boundary given with boundaryID = bcID as a
      ! * BC_interaction(bcID) = 0 : in this case particles bounce off this boundary
      ! * BC_interaction(bcID) = 1 : particles disappear from domain when hitting this BC
      !                              use this for open boundaries
      integer, allocatable :: BC_interaction(:) 
  
  
      !> Log particle data every this many time steps, i.e. for particleLogInterval = 1
      !! the particle data is logged every time step.
      integer :: particleLogInterval
  
      !> Buffer size for the force, velocity and other particle data buffers.
      !! Should be set in musubi.lua particles table.
      integer :: particleBufferSize
  
      !> Particles are communicated to processes once they come within halo_distance 
      !! of the boundary between two processes. For MEM (fully resolved) particles we 
      !! usually choose this to be one particle diameter. For DPS (unresolved) particles 
      !! we choose halo_distance = the mesh size.
      real(kind=rk) :: halo_distance
  
      ! -- FOR PARALLELIZATION -- !
      ! Communication type for sending operations
      type(mus_particles_communication_type) :: send
      ! Communication type for receiving operations
      type(mus_particles_communication_type) :: recv
    
      ! Pointers to the main particle routines needed in each time step
      procedure(pg_applyHydroForceFunc), pointer :: applyHydrodynamicForces => null()
      procedure(pg_moveFunc), pointer :: moveParticles => null()
      procedure(pg_mapFunc), pointer :: mapParticles => null()
      procedure(pg_momentumTransferFunc), pointer :: transferMomentumToFluid => null()
      procedure(pg_modifyAuxFieldFunc), pointer :: addParticleSourcesToAuxField => null()
  
      ! Pointers to the specific routines used to transfer momentum to the fluid 
      ! and interpolate flow field. These get assigned to different routines 
      ! depending on whether we use DPS, DPS_twoway, DPS_oneway particle kind etc.
      procedure(momTransFunc), pointer, nopass :: transfer_momentum => null()
      procedure(modAuxFieldFunc), pointer, nopass :: modify_auxfield => null()
      procedure(interpolateFluidPropFunc), pointer, nopass :: intp => null()
      procedure(calcVelAndPGradFunc), pointer, nopass :: calc_vel_and_p_grad => null()
  
      ! Pointers to the procedures for calculating hydrodynamic forces
      procedure(pg_dragForceFunc), pointer, nopass :: calcDragForce => null()
      procedure(pg_liftForceFunc), pointer, nopass :: calcLiftForce => null()
      procedure(pg_pressureForceFunc), pointer, nopass :: calcPressureForce => null()
  
      ! Procedure pointer to force contribution exchange routine
      procedure(exchangepIDvectorFunc),  pointer :: exchangeForces => null()
  
      ! Procedure pointer to particle state exchange routine
      procedure(exchangepIDvectorFunc),  pointer :: exchangeParticleStates => null()
  
      ! Procedure pointer to routine which communicates new particles to/from neighboring procs
      procedure(exchangepIDvectorFunc),  pointer :: exchangeNewParticles => null()
  
  end type mus_particle_group_type


  ! --- INTERFACE DEFINITION FOR THE PROCEDURE POINTERS --- !
  abstract interface
    subroutine connectivityFunc(this, scheme, stencil)
      import :: mus_particle_MEM_type
      import :: rk
      import :: mus_scheme_type
      import :: tem_stencilHeader_type
  
      class(mus_particle_MEM_type), intent(inout) :: this
      type(mus_scheme_type), intent(inout) :: scheme
      type(tem_stencilHeader_type), intent(in) :: stencil
    end subroutine connectivityFunc
  
    subroutine hydroForceFunc(this, scheme, stencil, params)
      import :: mus_particle_MEM_type
      import :: rk
      import :: mus_scheme_type
      import :: tem_stencilHeader_type
      import :: mus_param_type
  
      class(mus_particle_MEM_type), intent(inout) :: this
      type(mus_scheme_type), intent(inout) :: scheme
      type(tem_stencilHeader_type), intent(in) :: stencil
      type(mus_param_type), intent(in) :: params
    end subroutine hydroForceFunc
  
    subroutine wallForceFunc(this, BCinteraction, scheme, stencil, geometry, params, rmflag)
      import :: mus_particle_MEM_type
      import :: rk
      import :: mus_scheme_type
      import :: tem_stencilHeader_type
      import :: mus_geom_type
      import :: mus_param_type
  
      class(mus_particle_MEM_type), intent(inout) :: this
      integer, intent(in) :: BCinteraction(:)
      type(mus_scheme_type), intent(inout) :: scheme
      type(tem_stencilHeader_type), intent(in) :: stencil
      type(mus_geom_type), intent(in) :: geometry
      type(mus_param_type), intent(in) :: params
      logical, intent(out) :: rmflag
    end subroutine wallForceFunc
  
    subroutine interactForceFunc(this, particleGroup, params)
      import :: mus_particle_MEM_type
      import :: mus_particle_group_type
      import :: mus_param_type
  
      class(mus_particle_MEM_type), intent(inout) :: this
      type(mus_particle_group_type), intent(inout) :: particleGroup
      type(mus_param_type), intent(in) :: params
  
    end subroutine interactForceFunc
  
    subroutine applyBCFunc(this, scheme, stencil, params)
      import :: mus_particle_MEM_type
      import :: rk
      import :: mus_scheme_type
      import :: tem_stencilHeader_type
      import :: mus_param_type
  
      class(mus_particle_MEM_type), intent(inout) :: this
      type(mus_scheme_type), intent(inout) :: scheme
      type(tem_stencilHeader_type), intent(in) :: stencil
      type(mus_param_type), intent(in) :: params
    end subroutine applyBCFunc
  
    subroutine mapFunc(this, particleGroup, scheme, stencil, geometry, params, rmflag)
      import :: mus_particle_MEM_type
      import :: mus_particle_group_type
      import :: rk
      import :: mus_scheme_type
      import :: mus_geom_type
      import :: tem_stencilHeader_type
      import :: mus_param_type
  
      class(mus_particle_MEM_type), intent(inout) :: this
      type(mus_particle_group_type), intent(inout) :: particleGroup
      type(mus_scheme_type), intent(inout) :: scheme
      type(tem_stencilHeader_type), intent(in) :: stencil
      type(mus_geom_type), intent(in) :: geometry
      type(mus_param_type), intent(in) :: params
      logical, intent(out) :: rmflag
    end subroutine mapFunc
  
    subroutine updateVelocityFunc(this, dt)
      import :: mus_particle_MEM_type
      import :: rk
  
      class(mus_particle_MEM_type), intent(inout) :: this
      real(kind=rk), intent(in) :: dt
    end subroutine updateVelocityFunc
  
    subroutine updatePositionFunc(this, geometry, dt)
      import :: mus_particle_MEM_type
      import :: rk
      import :: mus_geom_type
  
      class(mus_particle_MEM_type), intent(inout) :: this
      type(mus_geom_type), intent(in) :: geometry
      real(kind=rk), intent(in) :: dt
    end subroutine updatePositionFunc
  
    subroutine exclusionListFunc(this, scheme, geometry, myRank, procs, nProcs, dx, rmflag)
      import :: mus_particle_MEM_type
      import :: rk
      import :: mus_scheme_type
      import :: mus_geom_type
  
      class(mus_particle_MEM_type), intent(inout) :: this
      type(mus_scheme_type), intent(inout) :: scheme
      type(mus_geom_type), intent(in) :: geometry
      integer, intent(in) :: myRank
      integer, intent(in) :: procs(:)
      integer, intent(in) :: nProcs
      real(kind=rk), intent(in) :: dx
      logical, intent(out) :: rmflag
    end subroutine exclusionListFunc
  
    subroutine swapForceBufferFunc(this)
      import :: mus_particle_MEM_type
  
      class(mus_particle_MEM_type), intent(inout) :: this
  
    end subroutine swapForceBufferFunc
  
    ! ----- INTERFACES FOR PARTICLE GROUP PROCEDURE POINTERS ---- !
    subroutine pg_applyHydroForceFunc(particleGroup, scheme, geometry, params)
      import :: mus_particle_group_type
      import :: mus_scheme_type
      import :: mus_geom_type
      import :: mus_param_type
  
      class(mus_particle_group_type), intent(inout) :: particleGroup
      type(mus_scheme_type), intent(inout) :: scheme
      type(mus_geom_type), intent(in) :: geometry
      type(mus_param_type), intent(in) :: params
    end subroutine pg_applyHydroForceFunc
  
    subroutine pg_moveFunc(particleGroup, scheme, geometry, params)
      import :: mus_particle_group_type
      import :: mus_scheme_type
      import :: mus_geom_type
      import :: mus_param_type
  
      class(mus_particle_group_type), intent(inout) :: particleGroup
      type(mus_scheme_type), intent(inout) :: scheme
      type(mus_geom_type), intent(in) :: geometry
      type(mus_param_type), intent(in) :: params
    end subroutine pg_moveFunc
  
    subroutine pg_mapFunc(particleGroup, scheme, geometry, params)
      import :: mus_particle_group_type
      import :: mus_scheme_type
      import :: mus_geom_type
      import :: mus_param_type
  
      class(mus_particle_group_type), intent(inout) :: particleGroup
      type(mus_scheme_type), intent(inout) :: scheme
      type(mus_geom_type), intent(in) :: geometry
      type(mus_param_type), intent(in) :: params
    end subroutine pg_mapFunc
  
    subroutine pg_momentumTransferFunc(particleGroup, scheme, geometry, params)
      import :: mus_particle_group_type
      import :: mus_scheme_type
      import :: mus_geom_type
      import :: mus_param_type
  
      class(mus_particle_group_type), intent(inout) :: particleGroup
      type(mus_scheme_type), intent(inout) :: scheme
      type(mus_geom_type), intent(in) :: geometry
      type(mus_param_type), intent(in) :: params
    end subroutine pg_momentumTransferFunc
  
    subroutine pg_modifyAuxFieldFunc( particleGroup, scheme, geometry, params )
      import :: mus_particle_group_type
      import :: mus_scheme_type
      import :: mus_geom_type
      import :: mus_param_type
  
      class(mus_particle_group_type), intent(inout) :: particleGroup
      type(mus_scheme_type), intent(inout) :: scheme
      type(mus_geom_type), intent(in) :: geometry
      type(mus_param_type), intent(in) :: params
    end subroutine pg_modifyAuxFieldFunc
  
    subroutine momTransFunc( particle, interpolator, scheme, &
                                          & geometry, params, Ftot          )
      import :: mus_particle_DPS_type
      import :: rk
      import :: mus_particle_interpolator_type
      import :: mus_scheme_type
      import :: mus_geom_type
      import :: mus_param_type
    
      type(mus_particle_DPS_type), intent(inout) :: particle
      type(mus_particle_interpolator_type), intent(in) :: interpolator
      type(mus_scheme_type), intent(inout) :: scheme
      type(mus_geom_type), intent(in) :: geometry
      type(mus_param_type), intent(in) :: params
      real(kind=rk), intent(out) :: Ftot(3)
    end subroutine
  
    subroutine modAuxFieldFunc( particle, interpolator, scheme, &
      & geometry, params, Ftot )
  
      import :: mus_particle_DPS_type
      import :: rk
      import :: mus_particle_interpolator_type
      import :: mus_scheme_type
      import :: mus_geom_type
      import :: mus_param_type
  
      type(mus_particle_DPS_type), intent(inout) :: particle
      type(mus_particle_interpolator_type), intent(in) :: interpolator
      type(mus_scheme_type), intent(inout) :: scheme
      type(mus_geom_type), intent(in) :: geometry
      type(mus_param_type), intent(in) :: params
      real(kind=rk), intent(out) :: Ftot(3)
    end subroutine modAuxFieldFunc
  
    subroutine interpolateFluidPropFunc( xp, coord_xp, scheme, geom_origin, dx, &
                        & interpolator, vel_xp, rho_xp, eps_f_xp, posOfCoord )
      import :: rk
      import :: mus_scheme_type
      import :: mus_particle_interpolator_type
      real(kind=rk), intent(in) :: xp(3)
      integer, intent(in) :: coord_xp(4)
      type(mus_scheme_type), intent(in) :: scheme
      real(kind=rk), intent(in) :: geom_origin(3)
      real(kind=rk), intent(in) :: dx
      type(mus_particle_interpolator_type), intent(in) :: interpolator
      real(kind=rk), intent(out) :: vel_xp(3)
      real(kind=rk), intent(out) :: rho_xp
      real(kind=rk), intent(out) :: eps_f_xp
      integer, intent(in), optional :: posOfCoord
    end subroutine
  
    subroutine calcVelAndPGradFunc(coord, scheme, grad_p, curl_u, err, posOfCoord )
      import :: rk
      import :: mus_scheme_type
      integer, intent(in) :: coord(4)
      type(mus_scheme_type), intent(in) :: scheme
      real(kind=rk), intent(out) :: grad_p(3)
      real(kind=rk), intent(out) :: curl_u(3)
      logical, intent(out) :: err
      integer, intent(in), optional :: posOfCoord
    end subroutine
  
    subroutine exchangepIDvectorFunc(this, send, recv, comm, myRank, message_flag)
      import :: mus_particle_group_type
      import :: mus_particles_communication_type
      
      class(mus_particle_group_type), intent(inout) :: this
      type(mus_particles_communication_type), intent(inout) :: send
      type(mus_particles_communication_type), intent(inout) :: recv
      integer, intent(in) :: comm
      integer, intent(in) :: myRank
      integer, intent(in) :: message_flag
    end subroutine exchangepIDvectorFunc
  
    ! Routines for calculating hydrodynamic forces on particles
    subroutine pg_dragForceFunc( particle, eps_p, nu, Fd )
      import :: mus_particle_DPS_type
      import :: rk
      !> Particle to apply force to
      type(mus_particle_DPS_type), intent(inout) :: particle
      !> Solid volume fraction interpolated to location of the particle
      real(kind=rk), intent(in) :: eps_p
      !> Fluid kinematic viscosity (phy)
      real(kind=rk), intent(in) :: nu
      !> Output: drag force on particle
      real(kind=rk), intent(out) :: Fd(3)
    end subroutine pg_dragForceFunc
  
    subroutine pg_liftForceFunc( particle, nu, Flift )
      import :: mus_particle_DPS_type
      import :: rk
      !> Particle to apply force to
      type(mus_particle_DPS_type), intent(inout) :: particle
      !> Fluid kinematic viscosity (phy)
      real(kind=rk), intent(in) :: nu
      !> Output: drag force on particle
      real(kind=rk), intent(out) :: Flift(3)
    end subroutine pg_liftForceFunc
  
    subroutine pg_pressureForceFunc( particle, Fp )
      import :: mus_particle_DPS_type
      import :: rk
      !> Particle to apply force to
      type(mus_particle_DPS_type), intent(inout) :: particle
      !> Fluid kinematic viscosity (phy)
      real(kind=rk), intent(out) :: Fp(3)
    end subroutine pg_pressureForceFunc
  end interface


contains
  !--- ROUTINES FOR OPERATING ON PARTICLE DYNARRAYS ---!
  
! ---- Dynamic particle array methods for momentum exchange method ---- !
?? copy :: DPA_inittxt(particle_MEM)
?? copy :: DPA_destroytxt(particle_MEM)
?? copy :: DPA_appendtxt(particle_MEM)
?? copy :: DPA_expandtxt(particle_MEM)
?? copy :: DPA_swaptxt(particle_MEM) 
?? copy :: DPA_removetxt(particle_MEM)
?? copy :: DPA_truncatetxt(particle_MEM)
?? copy :: DPA_sortposofval_txt(particle_MEM)

! ---- Dynamic particle array methods for Discrete Particle Simulations ---- !
?? copy :: DPA_inittxt(particle_DPS)
?? copy :: DPA_destroytxt(particle_DPS)
?? copy :: DPA_appendtxt(particle_DPS)
?? copy :: DPA_expandtxt(particle_DPS)
?? copy :: DPA_swaptxt(particle_DPS) 
?? copy :: DPA_removetxt(particle_DPS)
?? copy :: DPA_truncatetxt(particle_DPS)
?? copy :: DPA_sortposofval_txt(particle_DPS)

  ! ************************************************************************ !
  !> Routines for allocating the existsOnProc, addToProc and removeFromProc masks
  !! used to determine when particles should be sent over to new processes or which 
  !! processes need to receive position, velocity updates etc.
  
  subroutine allocateProcessMasks_MEM( particle, scheme, geometry, nProcs, lev )
    !> Particle to initialize
    type(mus_particle_MEM_type), intent(inout) :: particle
    !> Scheme for access to level descriptor
    type(mus_scheme_type), intent(inout) :: scheme
    !> Geometry to determine TreeIDs of elements 'covered' by particle
    type(mus_geom_type), intent(in) :: geometry
    !> Number of processes to communicate particle data with
    integer :: nProcs
    !> Level of this particle
    integer :: lev
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
  
  end subroutine allocateProcessMasks_MEM 
  
  
  subroutine allocateProcessMasks_DPS( particle, scheme, geometry, nProcs, lev )
    !> Particle to initialize
    type(mus_particle_DPS_type), intent(inout) :: particle
    !> Scheme for access to level descriptor
    type(mus_scheme_type), intent(inout) :: scheme
    !> Geometry to determine TreeIDs of elements 'covered' by particle
    type(mus_geom_type), intent(in) :: geometry
    !> Number of processes to communicate particle data with
    integer :: nProcs
    !> Level of this particle
    integer :: lev
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
  
  ! ************************************************************************ !
  ! Routines for printing particle group data, to be used for debugging
  subroutine printParticleGroup(particleGroup, logUnit )
    ! Particle group to print
    type(mus_particle_group_type), intent(in) :: particleGroup
    integer, intent(in) :: logUnit
    !------------------------------------------------------------------------
    integer :: iP, i
    character(len=100) :: formatString

    !------------------------------------------------------------------------
    write(formatString,'(A,I0,A)') '(', particleGroup%send%nProcs, 'L2)'

    if(particleGroup%particles_MEM%nvals < 1) then
      write(logUnit,*) 'No particles in this group!'
    end if

    do iP = 1, particleGroup%particles_MEM%nvals
      write(logUnit, *) '-- Particle ID ', particleGroup%particles_MEM &
        &                                              %val(iP)       &
        &                                              %ParticleID, '--'
      
      write(logUnit, *) 'Particle owner', particleGroup%particles_MEM &
        &                                              %val(iP)       &
        &                                              %owner
      ! coordOfOrigin
      write(logUnit, '(A)', advance='no') 'coordOfOrigin = [ '
      write(logUnit, '(4I3)', advance = 'no') ( particleGroup%particles_MEM &
        &                                          %val(iP)       &
        &                                          %coordOfOrigin(i), i = 1,4) 
      write(logUnit, '(A)') ' ]'
      ! position
      write(logUnit, '(A)', advance='no') 'pos = [ '
      write(logUnit, '(6E10.3)', advance = 'no') ( particleGroup%particles_MEM &
        &                                                          %val(iP)       &
        &                                                          %pos(i), i = 1,6) 
      write(logUnit, '(A)') ' ]'
      
      ! velocity
      write(logUnit, '(A)', advance='no') 'vel = [ '
      write(logUnit, '(6E10.3)', advance = 'no') ( particleGroup%particles_MEM &
        &                                                          %val(iP)       &
        &                                                          %vel(i), i = 1,6) 
      write(logUnit, '(A)') ' ]'
      ! force
      write(logUnit, '(A)', advance='no') 'F = [ '
      write(logUnit, '(6E10.3)', advance = 'no') ( particleGroup%particles_MEM &
        &                                                          %val(iP)       &
        &                                                          %F(i), i = 1,6) 
      write(logUnit, '(A)') ' ]'
      ! radius
      write(logUnit, '(A)', advance = 'no') 'radius = '
      write(logUnit, '(E10.3)') particleGroup%particles_MEM &
        &                                        %val(iP)       &
        &                                        %radius
      ! Rn
      write(logUnit, '(A)', advance = 'no') 'Rn = '
      write(logUnit, '(I0.3)') particleGroup%particles_MEM &
        &                                        %val(iP)       &
        &                                        %Rn
      ! mass
      write(logUnit, '(A)', advance = 'no') 'mass = '
      write(logUnit, '(E10.3)') particleGroup%particles_MEM &
        &                                        %val(iP)       &
        &                                        %mass

      if( particleGroup%send%nProcs > 0 ) then
        write(logUnit, '(A)', advance = 'no') 'send%proc = ['
        write(logUnit, '(4I6)', advance = 'no') ( particleGroup%send%proc(i), i = 1,particleGroup%send%nProcs) 
        write(logUnit, '(A)') ']'
        ! mask to show which procs this particle also exists on
        write(logUnit, '(A)', advance = 'no') 'existsOnProc = ['
        write(logUnit, trim(formatString), advance = 'no') (particleGroup%particles_MEM &
          &                                        %val(iP)       &
          &                                        %existsOnProc(i),i=1,particleGroup%send%nProcs)
        write(logUnit, '(A)') ']'

        write(logUnit, '(A)', advance = 'no') 'addToProc = ['
        write(logUnit, trim(formatString), advance = 'no') (particleGroup%particles_MEM &
          &                                        %val(iP)       &
          &                                        %addToProc(i),i=1,particleGroup%send%nProcs)
        write(logUnit, '(A)') ']'

        write(logUnit, '(A)', advance = 'no') 'removeFromProc = ['
        write(logUnit, trim(formatString), advance = 'no') (particleGroup%particles_MEM &
          &                                        %val(iP)       &
          &                                        %removeFromProc(i),i=1,particleGroup%send%nProcs)
        write(logUnit, '(A)') ']'
      end if

      
      

    end do

  end subroutine printParticleGroup


  subroutine printParticleGroup2_MEM(particleGroup, logUnit, myRank, iter )
    ! Particle group to print
    type(mus_particle_group_type), intent(in) :: particleGroup
    integer, intent(in) :: logUnit
    integer, intent(in) :: myRank
    integer, intent(in) :: iter    ! current solver iteration
    !------------------------------------------------------------------------
    integer :: iP, i
    character(len=100) :: formatString

    !------------------------------------------------------------------------
    write(formatString,'(A,I0,A)') '(', particleGroup%send%nProcs, 'L2)'
      
    write(logUnit,*) '----- PARTICLEGROUP RANK', myRank, ' -----'
    write(logUnit,*) 'iter = ', iter
    write(logUnit, '(A)', advance = 'no') 'send%proc = ['
    write(logUnit, '(4I6)', advance = 'no') ( particleGroup%send%proc(i), i = 1,particleGroup%send%nProcs) 
    write(logUnit, '(A)') ']'
    write(logUnit, '(A)', advance = 'no') 'recv%proc = ['
    write(logUnit, '(4I6)', advance = 'no') ( particleGroup%recv%proc(i), i = 1,particleGroup%recv%nProcs) 
    write(logUnit, '(A)') ']'

    if(particleGroup%particles_MEM%nvals < 1) then
      write(logUnit,*) 'No particles in this group!'
    end if

    do iP = 1, particleGroup%particles_MEM%nvals
      write(logUnit, *) 'Particle ID', particleGroup%particles_MEM &
        &                                              %val(iP)       &
        &                                              %ParticleID
      
      write(logUnit, *) 'Particle owner', particleGroup%particles_MEM &
        &                                              %val(iP)       &
        &                                              %owner
      ! coordOfOrigin
      write(logUnit, '(A)', advance='no') 'coordOfOrigin = [ '
      write(logUnit, '(4I3)', advance = 'no') ( particleGroup%particles_MEM &
        &                                          %val(iP)       &
        &                                          %coordOfOrigin(i), i = 1,4) 
      write(logUnit, '(A)') ' ]'
      ! position
      write(logUnit, '(A)', advance='no') 'pos = [ '
      write(logUnit, '(6E15.8)', advance = 'no') ( particleGroup%particles_MEM &
        &                                                          %val(iP)       &
        &                                                          %pos(i), i = 1,3) 
      write(logUnit, '(A)') ' ]'

      if( particleGroup%send%nProcs > 0 ) then
        ! mask to show which procs this particle also exists on
        write(logUnit, '(A)', advance = 'no') 'existsOnProc = ['
        write(logUnit, trim(formatString), advance = 'no') (particleGroup%particles_MEM &
          &                                        %val(iP)       &
          &                                        %existsOnProc(i),i=1,particleGroup%send%nProcs)
        write(logUnit, '(A)') ']'

        write(logUnit, '(A)', advance = 'no') 'addToProc = ['
        write(logUnit, trim(formatString), advance = 'no') (particleGroup%particles_MEM &
          &                                        %val(iP)       &
          &                                        %addToProc(i),i=1,particleGroup%send%nProcs)
        write(logUnit, '(A)') ']'

        write(logUnit, '(A)', advance = 'no') 'removeFromProc = ['
        write(logUnit, trim(formatString), advance = 'no') (particleGroup%particles_MEM &
          &                                        %val(iP)       &
          &                                        %removeFromProc(i),i=1,particleGroup%send%nProcs)
        write(logUnit, '(A)') ']'
      end if
    end do
    
    write(logUnit,*) '----------------------------------'

  end subroutine printParticleGroup2_MEM


  subroutine printParticleGroupData(particleGroup, logUnit)
    !> Particle group to print
    type(mus_particle_group_type), intent(in) :: particleGroup
    !> Unit to write output to
    integer, intent(in) :: logUnit
    !------------------------------------------------------------------------
    !------------------------------------------------------------------------
    write(logUnit,*) '-------- PARTICLEGROUP DATA ---------'
    write(logUnit,*) 'particle_kind = ', trim(particleGroup%particle_kind)
    write(logUnit,*) 'nParticles = ', particleGroup%nParticles
    write(logUnit,*) 'enableCollisions = ', particleGroup%enableCollisions
    if( particleGroup%enableCollisions ) then
      write(logUnit,*) 'collision_time = ', particleGroup%collision_time
      write(logUnit,*) 'collision_tol = ', particleGroup%collision_tol
    end if
    write(logUnit,*) 'particleLogInterval = ', particleGroup%particleLogInterval
    write(logUnit,*) 'particleBufferSize = ', particleGroup%particleBufferSize
    write(logUnit,*) 'Number of DEM subcycles = ', particleGroup%Nsubcycles

    write(logUnit,*) '-------- GLOBAL PARTICLE SETTINGS ---------'
    write(logUnit,*) 'Max particle dyn array size = ', maxContainerSize
    
  end subroutine printParticleGroupData

  subroutine printParticleGroup2_DPS(particleGroup, logUnit, myRank, iter )
    ! Particle group to print
    type(mus_particle_group_type), intent(in) :: particleGroup
    integer, intent(in) :: logUnit
    integer, intent(in) :: myRank
    integer, intent(in) :: iter    ! current solver iteration
    !------------------------------------------------------------------------
    integer :: iP, i
    integer(kind=long_k) :: TreeID

    !------------------------------------------------------------------------
    write(logUnit,*) 'DPS particles in group: ', particleGroup%particles_DPS%nvals
    write(logUnit,*) '----- PARTICLEGROUP RANK', myRank, ' -----'
    write(logUnit,*) 'iter = ', iter
    write(logUnit, '(A)', advance = 'no') 'send%proc = ['
    write(logUnit, '(4I6)', advance = 'no') ( particleGroup%send%proc(i), i = 1,particleGroup%send%nProcs) 
    write(logUnit, '(A)') ']'
    write(logUnit, '(A)', advance = 'no') 'recv%proc = ['
    write(logUnit, '(4I6)', advance = 'no') ( particleGroup%recv%proc(i), i = 1,particleGroup%recv%nProcs) 
    write(logUnit, '(A)') ']'

    if(particleGroup%particles_DPS%nvals < 1) then
      write(logUnit,*) 'No particles in this group!'
    end if

    do iP = 1, particleGroup%particles_DPS%nvals
      write(logUnit, *) 'Particle ID', particleGroup%particles_DPS &
        &                                              %val(iP)       &
        &                                              %ParticleID
      
      write(logUnit, *) 'Particle owner', particleGroup%particles_DPS &
        &                                              %val(iP)       &
        &                                              %owner
      ! coordOfOrigin
      write(logUnit, '(A)', advance='no') 'coordOfOrigin = [ '
      write(logUnit, '(4I3)', advance = 'no') ( particleGroup%particles_DPS &
        &                                          %val(iP)       &
        &                                          %coordOfOrigin(i), i = 1,4) 
      write(logUnit, '(A)') ' ]'
      ! TreeID of origin
      TreeID = tem_IdOfCoord( particleGroup%particles_DPS      &
        &                                  %val(iP)            &
        &                                  %coordOfOrigin(1:4) )
      write(logUnit, '(A)', advance='no') 'TreeIDOfOrigin = '
      write(logUnit, '(I10)', advance='no') TreeID
      write(logUnit, '(A)') 

      ! position
      write(logUnit, '(A)', advance='no') 'pos = [ '
      write(logUnit, '(6E15.8)', advance = 'no') ( particleGroup%particles_DPS &
        &                                                          %val(iP)       &
        &                                                          %pos(i), i = 1,3) 
      write(logUnit, '(A)') ' ]'
    end do
    write(logUnit,*) '----------------------------------'

  end subroutine printParticleGroup2_DPS

  subroutine printpIDlist(particleGroup)
    ! Particle group to print
    type(mus_particle_group_type), intent(in) :: particleGroup
    !------------------------------------------------------------------------
    integer :: iP, nvals
    !------------------------------------------------------------------------
    nvals = particleGroup%particles_MEM%nvals

    ! Unsorted pIDlist
    write(stdOutUnit, '(A)', advance='no') 'pIDlist = [ '
    write(stdOutUnit, '(6I3)', advance='no') ( particleGroup%particles_MEM  &
      &                                      %pIDlist(iP), &
      &                                      iP = 1,nvals )
    write(stdOutUnit, '(A)') ' ]'
    
    ! Sorted pIDlist
    write(stdOutUnit, '(A)', advance='no') 'pIDsort = [ '
    write(stdOutUnit, '(6I3)', advance='no') ( particleGroup%particles_MEM  &
      &                                      %pIDsort(iP), &
      &                                       iP = 1,nvals )
    write(stdOutUnit, '(A)') ' ]'

  end subroutine printpIDlist

  subroutine test_append_da_particle(particleGroup)
    type(mus_particle_group_type), intent(inout) :: particleGroup
    !-----------------------------------------------------------------------!
    ! Array of five particles to test whether sorting by pID's works
    type(mus_particle_MEM_type) :: testParticles(5)

    integer :: particleIDs(5)
    integer :: iParticle
    logical :: wasAdded

    !-----------------------------------------------------------------------!
    wasAdded = .false.
    call init_da_particle_MEM(particleGroup%particles_MEM, 5)

    ! Set test particle properties
    particleIDs = (/ 12, 8, 3, 10, 1 /)

    
    do iParticle = 1,5
      testParticles(iParticle)%pos = 0.0_rk  
      testParticles(iParticle)%vel = 0.0_rk  
      testParticles(iParticle)%radius = 0.05_rk  
      testParticles(iParticle)%mass = 1.0_rk  
      testParticles(iParticle)%particleID = particleIDs(iParticle)

      call append_da_particle_MEM( me       = particleGroup%particles_MEM,  &
                             & particle = testParticles(iParticle), &
                             & wasAdded = wasAdded                  )
    end do

    call printParticleGroup(particleGroup, stdOutUnit)
    call printpIDlist(particleGroup)
    call remove_particle_from_da_particle_mem(particleGroup%particles_MEM, 2)
    write(stdOutUnit,*) 'After removing 5th particle'
    call printpIDlist(particleGroup)

  end subroutine test_append_da_particle

end module mus_particle_type_module
