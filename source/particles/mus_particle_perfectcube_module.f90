module mus_particle_perfectcube_module
use env_module,               only : rk, long_k, stdOutUnit
use tem_topology_module,      only : tem_IdOfCoord, tem_FirstIdAtLevel, &
                                   & tem_coordOfId
use mus_scheme_type_module,   only : mus_scheme_type
use mus_geom_module,          only : mus_geom_type
use mus_param_module,         only : mus_param_type
use mus_particle_aux_module,  only : treeIDlocalOnMyRank, has_remote_neighbors, &
                                   & getBaryOfCoord
use mus_particle_type_module, only : mus_particle_DPS_type, &
                                    & mus_particle_group_type, &
                                    & init_da_particle_DPS, &
                                    & append_da_particle_DPS
implicit none

type perfectcube_params_type
  !> Particle radius 
  real(kind=rk) :: radius
  !> Particle mass
  real(kind=rk) :: mass
  !> Particle rotational inertia
  real(kind=rk) :: rotInertia
  !> Particles per side of the perfect cube
  integer :: particles_per_side
  !> Total number of particles per process
  integer :: particles_per_process
end type perfectcube_params_type

type(perfectcube_params_type), save :: perfectcube_params

contains

subroutine createParticles_perfectcube( particleGroup, origin, L, particles_per_side, &
                            & Rp, mass, Nparticles, myRank )
  !> Array of particles
  type(mus_particle_group_type), intent(inout) :: particleGroup
  !> Cartesian coordinate origin of cubical domain to spawn particles in
  real(kind=rk), intent(in) :: origin(3)
  !> Side length of the perfect cube
  real(kind=rk) :: L
  !> Number of particles per one side of the cube. In total there will be 
  !! particles_per_side^3 particles in the cubical domain
  integer, intent(in) :: particles_per_side
  !> particle radius
  real(kind=rk), intent(in) :: Rp
  real(kind=rk), intent(in) :: mass
  !> Output: number of particles
  integer, intent(out) :: Nparticles
  !> This process rank
  integer, intent(in) :: myRank
  ! ---------------------------------------------- !
  logical :: wasAdded
  real(kind=rk) :: spacing, x, y, z
  integer :: ix, iy, iz
  integer :: particleID
  type(mus_particle_DPS_type) :: particle
  ! ---------------------------------------------- !
  Nparticles = particles_per_side**3
  call init_da_particle_DPS( particleGroup%particles_DPS, Nparticles )
  spacing = L / particles_per_side

  particleID = myRank*Nparticles + 1 
  do ix = 1, particles_per_side
    x = origin(1) + 0.5*spacing + (ix - 1)*spacing
    do iy = 1, particles_per_side
      y = origin(2) + 0.5*spacing + (iy - 1)*spacing
      do iz = 1, particles_per_side
        z = origin(3) + 0.5*spacing + (iz - 1)*spacing
        
        particle%pos(1) = x
        particle%pos(2) = y
        particle%pos(3) = z
        particle%pos(4:6) = 0.0_rk

        particle%vel = 0.0_rk
        particle%Fext = 0.0_rk

        particle%radius = Rp
        particle%mass   = mass
        particle%rotInertia  = 0.4*particle%mass*Rp**2
        particle%particleID = particleID 

        call append_da_particle_DPS( me = particleGroup%particles_DPS, &
                                   & particle = particle, &
                                   & length = 1, &
                                   & wasAdded = wasAdded )

        particleID = particleID + 1

        end do
    end do
  end do
  particleGroup%nParticles = particleGroup%particles_DPS%nvals

end subroutine createParticles_perfectcube

!> This routine initializes a specified number of particles (nParticles)
!! on the first Nparticles fluid elements which satisfy
!! * Is a local element
!! * All its neighbor elements are also local elements
subroutine createParticles_onTreeIDs( particleGroup, scheme, geometry, params, &
                                    & radius, mass, Nparticles, myRank )
  type(mus_particle_group_type), intent(inout) :: particleGroup
  type(mus_scheme_type), intent(in) :: scheme
  type(mus_geom_type), intent(in) :: geometry
  type(mus_param_type), intent(in) :: params
  real(kind=rk), intent(in) :: radius
  real(kind=rk), intent(in) :: mass
  integer, intent(in) :: Nparticles
  integer, intent(in) :: myRank
  ! ---------------------------------------- !
  integer :: iDir, iElem, particleCounter, lev
  integer :: particleID
  logical :: hasRemoteNeighbors, wasAdded
  type(mus_particle_DPS_type) :: particle
  ! ---------------------------------------- !
  lev = geometry%tree%global%maxLevel
  particleCounter = 0
  particleID = myRank*Nparticles + 1 

  do iElem = 1, scheme%pdf(lev)%nElems_fluid
    ! For each element check whether it has remote neighbors
    hasRemoteNeighbors = has_remote_neighbors( iElem  = iElem,  &
                                             & scheme = scheme, &
                                             & lev    = lev     )
    ! We only place particles at cells which do not have
    ! any remote neighbors
    if( hasRemoteNeighbors ) then
      continue
    else
      ! Place a particle at this location

      ! Set the position to the baryCenter of this element
      call placeParticleAtTreeID(                                &
                & particle = particle,                           &
                & TreeID   = scheme%levelDesc(lev)%total(iElem), &
                & geometry = geometry,                           &
                & params   = params                              )

      ! Set initial forces and velocity
      particle%vel = 0.0_rk
      particle%Fext = 0.0_rk

      ! Set other particle properties
      particle%radius = radius
      particle%mass   = mass
      particle%rotInertia  = 0.4*particle%mass*radius**2
      particle%particleID = particleID 

      ! Increment the particle counter
      particleCounter = particleCounter + 1

      ! Append the particle to the particleGroup array
      call append_da_particle_DPS( me = particleGroup%particles_DPS, &
                                 & particle = particle, &
                                 & length = 1, &
                                 & wasAdded = wasAdded )

      particleID = particleID + 1

      ! If we have initialized the required number of particles, exit routine.
      if(particleCounter >= Nparticles) return
                              
    end if
  end do ! iElem
end subroutine createParticles_onTreeIDs

!> Changes the position of a particle to the barycenter of the given TreeID
subroutine placeParticleAtTreeID(particle, TreeID, geometry, params)
  type(mus_particle_DPS_type), intent(inout) :: particle
  integer(kind=long_k), intent(in) :: TreeID
  type(mus_geom_type), intent(in) :: geometry
  !> Params for access to dt, dx, etc.
  type(mus_param_type), intent(in) :: params
  ! ---------------------------------------- !
  integer :: coord(4), lev
  integer(kind=long_k) :: TIDoffset
  real(kind=rk) :: bary(3), geom_origin(3), dx
  
  ! ---------------------------------------- !
  lev = geometry%tree%global%maxLevel
  ! First get the barycenter of the element with this TreeID
  TIDoffset = tem_firstIdAtLevel(lev)
  geom_origin  = geometry%tree%global%origin
  dx = params%physics%dxLvl(lev)

  coord = tem_coordOfID( TreeID = TreeID, &
                       & offset = TIDoffset)

  bary = getBaryOfCoord( coord  = coord(1:3), &
                       & origin = geom_origin,             &
                       & dx     = dx                       )

  ! Change particle position to this location
  particle%pos(1:3) = bary
  particle%pos(4:6) = 0.0_rk
end subroutine placeParticleAtTreeID

end module mus_particle_perfectcube_module
