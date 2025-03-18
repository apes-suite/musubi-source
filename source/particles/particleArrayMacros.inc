?? text :: DPA_decltxt(tname, tstring)
! \brief Smart Dynamic Array (DA) for ?tstring?
!
! This datatype implements a dynamic array for particles, 
! which is capable of growing and adding of unique elements. 
! It is available for various types of particles. 
! here we deal with $tstring$.
!> dynamic array (DA) type for mus_?tname?_type
type dyn_?tname?_array_type
  ! Current number of particles in the array
  integer :: nvals = 0

  ! Current container size
  integer :: containersize = 0

  ! The actual array holding the particles
  type(mus_?tname?_type), allocatable :: val(:)

  ! Array of particle ID's which mirrors the actual particle array val(:)
  ! So order of the ID's in this list is identical to order of val(:)
  integer, allocatable :: pIDlist(:) 

  ! Sorted positions of the particle ID list
  ! example, just like the tem dyn array types
  ! pIDlist:  8, 6, 7, 4, 5
  ! pIDsort:  4, 5, 2, 3, 1
  ! So we can for example traverse the pIDlist in sorted order using
  ! do i = 1, nvals
  !   val = pIDlist(pIDsort(i))
  ! end do

  integer, allocatable :: pIDsort(:) 

end type
?? end text DPA_decltxt


?? text :: DPA_inittxt(tname)

  subroutine init_da_?tname?(me, length)
    type(dyn_?tname?_array_type), intent(inout) :: me !< dynamic array to init
    integer, intent(in), optional :: length !< initial length of the container

    !----------------------------------------------!
    if (present(length)) then
      me%containersize = length
    else
      me%containersize = zeroLength
    end if

    ! deallocate ...
    if( allocated( me%val ) ) deallocate(me%val)
    if( allocated( me%pIDlist ) ) deallocate(me%pIDlist)
    if( allocated( me%pIDsort ) ) deallocate(me%pIDsort)

    ! ... and reallocate
    allocate(me%val(me%containersize))
    allocate(me%pIDlist(me%containersize))
    allocate(me%pIDsort(me%containersize))
    me%nvals = 0

  end subroutine init_da_?tname?

?? end text DPA_inittxt

?? text :: DPA_destroytxt(tname)
  ! Destroy particle dynamic array
  subroutine destroy_da_?tname?(me)
    type(dyn_?tname?_array_type), intent(inout) :: me !< dynamic array to destroy
    !----------------------------------------------!

    me%containersize = 0
    me%nvals         = 0

    if( allocated( me%val ) ) deallocate(me%val)
    if( allocated( me%pIDlist ) ) deallocate(me%pIDlist)
    if( allocated( me%pIDsort ) ) deallocate(me%pIDsort)
  end subroutine destroy_da_?tname?
?? end text DPA_destroytxt

?? text :: DPA_appendtxt(tname)

  subroutine append_da_?tname?(me, particle, length,  wasAdded)
    ! Dynamic array to append to
    type(dyn_?tname?_array_type), intent(inout) :: me 
    ! Particles to append
    type(mus_?tname?_type), intent(inout) :: particle
    ! Number of elements to increase containersize with
    integer, intent(in), optional :: length
    ! Logical will be set to 1 is appending was succesful, 0 if not
    logical :: wasAdded !< 
    !----------------------------------------------!
    integer :: foundPos
    integer :: i

    ! do a binary search on existing entries (returns closest entry next to
    ! it if not found).
    foundpos = sortposofval_?tname?(me, particle%particleID, .true.)
    wasadded = .false.

    ! if it found the value, the position is smaller than nvals
    if (foundpos <= me%nvals) then

      ! Check if particle ID is already in particle array 
      ! If so, do nothing
      if ( me%pIDlist(me%pIDsort(foundpos)) == particle%particleID ) then
        ! write(stdOutUnit,*) "Warning append_da_particle: particle already in array "
        return
      else
        ! need to append a new value!

        if (me%nvals == huge(me%nvals)) then
           write(stdOutUnit,*) "reached end of integer range for dynamic particle array!"
           write(stdOutUnit,*) "aborting!!"
           stop
        end if

        wasadded = .true.
        if (me%nvals == me%containersize) then
          ! container is full, need to expand it
          call expand_da_?tname?(me = me, length = length)
        end if
        me%nvals = me%nvals + 1

        ! put the new value into the last position in the
        ! array.
        me%val(me%nvals) = particle
        me%pIDlist(me%nvals) = particle%particleID
        do while( foundpos < me%nvals )
          if(me%pIDlist(me%pIDsort(foundpos)) /= particle%particleID) then
            exit
          end if
          ! in case of multiple entries with the same value
          ! move on to the first differing entry.
          foundpos = foundpos + 1
        end do
        ! shift the sorted list of indices, to create a
        ! whole for the value to be inserted, at position
        ! foundpos.
        do i=me%nvals-1,foundpos,-1
          me%pIDsort(i+1) = me%pIDsort(i)
        end do
        ! put the index of the new value into the
        ! sorted list at the now freed position.
        me%pIDsort(foundpos) = me%nvals

      end if

    else

      ! value to append is larger than all existing ones,
      ! just put it to the end of the list, this captures
      ! also the case of empty lists.
      ! in this case foundpos = me%nvals + 1 holds.
      wasadded = .true.
      if (foundpos > me%containersize) then
        ! expand the array, if its boundary is reached
        call expand_da_?tname?(me = me, length = length)
      end if
      me%nvals = foundpos
      me%val(foundpos) = particle
      me%pIDlist(foundpos) = particle%particleID
      me%pIDsort(foundpos) = foundpos

    end if
  end subroutine append_da_?tname?

?? end text DPA_appendtxt


?? text :: DPA_expandtxt(tname)

  subroutine expand_da_?tname?(me, length)
    !------------------------------------------------------------------------
    type(dyn_?tname?_array_type), intent(inout) :: me !< array to resize
    !> optional length to expand the array with
    integer, intent(in), optional :: length
    !------------------------------------------------------------------------
    type(mus_?tname?_type), allocatable :: swpval(:)
    integer, allocatable :: swpIDlist(:)
    integer, allocatable :: swpIDsort(:)
    integer :: explen
    !------------------------------------------------------------------------

    ! if length is present, use that, otherwise double the size
    if( present( length ) ) then
      explen = length
    else
      ! set the global minimum length, if doubling would be smaller than that
      explen = max(me%containersize, minLength)
    end if


    ! check whether the new size will exceed the max container size.
    if( (me%containersize + explen) >= maxContainerSize ) then
      ! if so, expand to the maximum size
      me%containersize = maxContainerSize
    else
      ! if not, expand to the calculated size
      me%containersize = me%containersize + explen
    end if

    ! Now make a larger array and copy all the current values into it.
    ! only need to copy values, if there are actually values to append.
    if (me%nvals > 0) then
      allocate(swpval(me%containersize))
      swpval(1:me%nvals) = me%val(1:me%nvals)
      call move_alloc( swpval, me%val )
      
      allocate(swpIDlist(me%containersize))
      swpIDlist(1:me%nvals) = me%pIDlist(1:me%nvals)
      call move_alloc( swpIDlist, me%pIDlist )
      
      allocate(swpIDsort(me%containersize))
      swpIDsort(1:me%nvals) = me%pIDsort(1:me%nvals)
      call move_alloc( swpIDsort, me%pIDsort )
    
    else ! me%nvals == 0
      if( allocated(me%val) ) then
        deallocate(me%val)
        allocate(me%val(me%containersize))
      end if
      
      if( allocated(me%pIDlist) ) then
        deallocate(me%pIDlist)
        allocate(me%pIDlist(me%containersize))
      end if
      
      if( allocated(me%pIDsort) ) then
        deallocate(me%pIDsort)
        allocate(me%pIDsort(me%containersize))
      end if
    end if

  end subroutine expand_da_?tname?


?? end text DPA_expandtxt

?? text :: DPA_swaptxt(tname) 
  !> Swaps the position of two particles in particle dynamic array
  !! New position of iElem1 = old position of iElem2 and vice-versa
  !! Also updates the pIDlist and PIDsort arrays
  subroutine swap_da_?tname?( me, iElem1, iElem2 )
    !> Particle array to operate on
    type(dyn_?tname?_array_type), intent(inout) :: me 
    !> Current index of one element
    integer, intent(in) :: iElem1
    !> Current index of other element
    integer, intent(in) :: iElem2
    !------------------------------------------------------------------------
    type(mus_?tname?_type) :: tmp
    integer :: tmpPID
    integer :: k, k1, k2
    !------------------------------------------------------------------------

    k1 = -1
    k2 = -1

    if(iElem1 == iElem2) then
      return
    end if
    ! 1. Update particle array and pIDlist
    ! Copy first element to tmp
    tmp = me%val(iElem1)
    tmpPID = me%pIDlist(iElem1)

    ! Put second element at position of first
    me%val(iElem1) = me%val(iElem2)
    me%pIDlist(iElem1) = me%pIDlist(iElem2)

    ! Put first element at old position of second
    me%val(iElem2) = tmp
    me%pIDlist(iElem2) = tmpPID

    ! 2. Update the sorted list 
    ! Find the indices of PIDsort pointing to iElem1 and iElem2
    do k = 1, me%nvals
      if (me%pIDsort(k) == iElem1) then
        k1 = k
        exit
      end if
    end do
    
    do k = 1, me%nvals
      if (me%pIDsort(k) == iElem2) then
        k2 = k
        exit
      end if
    end do

    ! Swap the values
    me%pIDsort(k1) = iElem2
    me%pIDsort(k2) = iElem1

  end subroutine swap_da_?tname?
?? end text DPA_swaptxt


?? text :: DPA_truncatetxt(tname) 
  !> Truncate the dynamic particle array to only fit the actual entries
  subroutine truncate_da_?tname?( me )
    !> Particle array to operate on
    type(dyn_?tname?_array_type), intent(inout) :: me
    !------------------------------------------------------------------------
    type(mus_?tname?_type), allocatable :: swpval(:)
    !------------------------------------------------------------------------

    if (me%nvals < me%containersize) then
      allocate(swpval(me%nvals))

      swpval = me%val(:me%nvals)

      call move_alloc(swpval, me%val)

      me%containersize = me%nvals
    end if

  end subroutine truncate_da_?tname?
?? end text DPA_truncatetxt

?? text :: DPA_removetxt(tname)
  ! Destroy particle dynamic array
  subroutine remove_particle_from_da_?tname?(particles, iElem)
     !> Particle group to operate on
     type(dyn_?tname?_array_type), intent(inout) :: particles
     !> Current index of element to remove
     integer, intent(in) :: iElem
     !------------------------------------------------------------------------
     integer :: ks, k
     !------------------------------------------------------------------------

     ! First check if this is a positive particle ID, so indicates an actual particle
     if( particles%val(iElem)%particleID < 0 ) return

     ! Place element to remove at the end of the array
     call swap_da_?tname?( me = particles,          &
                         & iElem1 = iElem,          &
                         & iElem2 = particles%nvals )

     ! Set particle ID to negative to indicate this particle no longer belongs
     ! to dyn_array
     particles%val( particles%nvals )%particleID = &
       &  -1*particles%val( particles%nvals )%particleID

     ! Update PIDsort: first do linear search to find 
     ! the element of pIDsort corresponding to iElem
     searchLoop: do ks = 1, particles%nvals
       if( particles%pIDsort(ks) == particles%nvals ) then
         ! pIDsort(ks) points to the element to remove
         do k = ks, particles%nvals - 1
           ! Shift all elements to the right of ks one to the left.
           particles%pIDsort(k) &
             & = particles%pIDsort(k + 1)
         end do
         exit searchLoop
       end if
     end do searchLoop

     ! After successful removal of particle, decrease nvals 
     particles%nvals = particles%nvals - 1

  end subroutine remove_particle_from_da_?tname?
?? end text DPA_removetxt


?? text :: DPA_sortposofval_txt(tname)

  !> return the sorted position of a value in the given dynamic array
  !!
  !! if the value was not found,
  !!  - return 0 if nextifnotfound = .false.
  !!  - return position at the end if nextifnotfound = .true.
  function sortposofval_?tname?(me, pID, nextifnotfound, lower, upper) result(pos)
    !------------------------------------------------------------------------
    type(dyn_?tname?_array_type), intent(in) :: me !< dynamic array
    integer, intent(in) :: pID !< particle ID to look for
    !> flag to indicate, if the next entry in the list should be returned,
    !! if the searched one is not found.
    logical, intent(in), optional :: nextifnotfound
    integer, intent(in), optional :: lower !< lower search limit
    integer, intent(in), optional :: upper !< upper search limit
    integer :: pos !< position of val in the sorted particle ID list, 0 if not found
    !------------------------------------------------------------------------
    logical :: retnext
    integer :: lb, ub
    integer :: mid
    integer :: lb_val, ub_val
    integer :: mid_val
    !------------------------------------------------------------------------

    retnext = .false.
    if (present(nextifnotfound)) retnext = nextifnotfound

    lb = 1
    ub = me%nvals

    if( present( lower ) ) lb = lower
    if( present( upper ) ) ub = upper

    pos = 0
    if (retnext) pos = lb

    !> binary search on sorted list
    do while(ub >= lb)
      lb_val = me%pIDlist(me%pIDsort(lb))

      ! If pID is smaller than smallest val in pIDlist, exit
      if (pID < lb_val) then
        if (retnext) pos = lb
        exit
      end if

      ub_val = me%pIDlist(me%pIDsort(ub))

      ! Also if pID is greater than greatest value in list, exit
      if (pID > ub_val) then
        if (retnext) pos = ub+1
        exit
      end if

      ! safe guard against integer limit overflow
      mid = lb + (ub-lb) / 2
      mid_val = me%pIDlist(me%pIDsort(mid))
      if (pID == mid_val) then
        pos = mid
        exit
      end if
      if (pID > mid_val) then
        lb = mid + 1
      else
        ub = mid - 1
      end if
    end do
  end function sortposofval_?tname?
?? end text DPA_sortposofval_txt

! --- Macros for loading particle data from lua script --- !
?? text :: loadParticleCreatorData_txt(tname)

  subroutine load_?tname?_creator_data( conf, parent, particle_creator,     &
                                   &  nParticles, chunkSize, scheme,        &
                                   & geometry, myRank  )
    type(flu_State), intent(in) :: conf

    !> Handle to parent table if position, velocity tables etc are inside
    !! another table
    integer, intent(in), optional :: parent

    !> Dynamic particle array to append particles read from file to
    type(mus_particle_creator_type), intent(inout) :: particle_creator

    !> Total number of particles to read
    integer, intent(in) :: nParticles

    ! Number of particles to read as one "chunk"
    integer, intent(in) :: chunkSize

    !> Scheme
    type(mus_scheme_type), intent(in) :: scheme

    !> Geometry 
    type(mus_geom_type), intent(in) :: geometry

    !> This MPI process rank
    integer, intent(in) :: myRank
    !--------------------------------------------!
    type(mus_?tname?_type), allocatable :: particleBuffer(:)
    integer :: iStart, iChunk, iParticle, kparticle
    integer :: nChunks, nChunkVals
    logical :: wasAdded, isLocal
    !--------------------------------------------!
    nChunks = nParticles/chunkSize + 1
    iStart = 1
    kparticle = 1  
    allocate( particleBuffer(chunkSize) )
    particle_creator%global_Nparticles = 0
    particle_creator%N_times_called = 0

    chunkLoop: do iChunk = 1, nChunks
      ! Load one chunk into the particleBuffer
      write(stdOutUnit,*) 'loading chunk ', iChunk
      ! Note: loading position values also sets the particle ID's!
      call mus_load_?tname?_data_chunk( conf           = conf,           &
                                       & parent         = parent,      &
                                       & particleBuffer = particleBuffer, &
                                       & nChunkVals     = nChunkVals,     &
                                       & key            = 'position',     &
                                       & iStart         = iStart          )

      call mus_load_?tname?_data_chunk( conf           = conf,           &
                                       & parent         = parent,      &
                                       & particleBuffer = particleBuffer, &
                                       & nChunkVals     = nChunkVals,     &
                                       & key            = 'velocity',     &
                                       & iStart         = iStart          )

      call mus_load_?tname?_data_chunk( conf           = conf,           &
                                       & parent         = parent,      &
                                       & particleBuffer = particleBuffer, &
                                       & nChunkVals     = nChunkVals,     &
                                       & key            = 'force',        &
                                       & iStart         = iStart          )

      call mus_load_?tname?_data_chunk( conf           = conf,           &
                                       & parent         = parent,      &
                                       & particleBuffer = particleBuffer, &
                                       & nChunkVals     = nChunkVals,     &
                                       & key            = 'radius',       &
                                       & iStart         = iStart          )

      call mus_load_?tname?_data_chunk( conf           = conf,           &
                                       & parent         = parent,      &
                                       & particleBuffer = particleBuffer, &
                                       & nChunkVals     = nChunkVals,     &
                                       & key            = 'mass',         &
                                       & iStart         = iStart          )


      ! nChunkVals is the number of particles actually read from the lua file
      ! Can be less than chunkSize if there are fewer than chunkSize particles
      ! left in lua table
      ! write(stdOutUnit,*) 'nChunkVals = ', nChunkVals

      !-- SET PARTICLE ROTATIONAL INERTIA --!
      do iParticle = 1, nChunkVals
        particleBuffer( iParticle )%rotInertia &
          & = 0.4 * particleBuffer( iParticle )%mass &
          &       * particleBuffer( iParticle )%radius**2

      end do

      ! Append particles in particleBuffer to particleGroup
      write(stdOutUnit,*) 'Loading particles: nChunkVals = ', nChunkVals
      do iParticle = 1, nChunkVals
        wasAdded = .FALSE.
        ! write(stdOutUnit,*) 'Appending particle'

        ! Check if this particle position is local to my process
        isLocal = positionLocalOnMyRank( pos = particlebuffer(iparticle)%pos(1:3), &
                                       & geometry = geometry,                      &
                                       & scheme = scheme,                          & 
                                       & myRank = myRank                           )

        if(isLocal) then
          call append( me = particle_creator%position, val = particlebuffer(iparticle)%pos(1:6) )
          call append( me = particle_creator%velocity, val = particlebuffer(iparticle)%vel(1:6) )
          call append( me = particle_creator%force, val = particlebuffer(iparticle)%Fext(1:6) )
          call append( me = particle_creator%radius, val = particlebuffer(iparticle)%radius )
          call append( me = particle_creator%mass, val = particlebuffer(iparticle)%mass )
          call append( me = particle_creator%IDoffset, val = kparticle )
        end if

        ! Increment global number of particles
        particle_creator%global_Nparticles = particle_creator%global_Nparticles + 1 

        ! kparticle keeps track of the number of particles that have been added to the 
        ! particle creators on all processes. So after we've loaded all the particles 
        ! it will be equal to the number of particles specified in the lua script.
        kparticle = kparticle + 1

      end do

      iStart = iStart + chunkSize
    end do chunkLoop

    particle_creator%nParticles = particle_creator%radius%nvals
    deallocate( particleBuffer )

  end subroutine load_?tname?_creator_data

?? end text loadParticleCreatorData_txt

?? text :: loadParticleData_txt(tname)

  subroutine load_?tname?_data( conf, parent, particles,     &
                                   &  nParticles, chunkSize  )
    type(flu_State), intent(in) :: conf

    !> Handle to parent table if position, velocity tables etc are inside
    !! another table
    integer, intent(in), optional :: parent

    !> Dynamic particle array to append particles read from file to
    type(dyn_?tname?_array_type), intent(inout) :: particles

    !> Total number of particles to read
    integer, intent(in) :: nParticles

    ! Number of particles to read as one "chunk"
    integer, intent(in) :: chunkSize
    !--------------------------------------------!
    type(mus_?tname?_type), allocatable :: particleBuffer(:)
    integer :: iStart, iChunk, iParticle
    integer :: nChunks, nChunkVals
    logical :: wasAdded
    !--------------------------------------------!
    nChunks = nParticles/chunkSize + 1
    iStart = 1
    allocate( particleBuffer(chunkSize) )

    chunkLoop: do iChunk = 1, nChunks
      ! Load one chunk into the particleBuffer
      ! write(stdOutUnit,*) 'loading chunk ', iChunk
      ! Note: loading position values also sets the particle ID's!
      call mus_load_?tname?_data_chunk( conf           = conf,           &
                                       & parent         = parent,      &
                                       & particleBuffer = particleBuffer, &
                                       & nChunkVals     = nChunkVals,     &
                                       & key            = 'position',     &
                                       & iStart         = iStart          )

      call mus_load_?tname?_data_chunk( conf           = conf,           &
                                       & parent         = parent,      &
                                       & particleBuffer = particleBuffer, &
                                       & nChunkVals     = nChunkVals,     &
                                       & key            = 'velocity',     &
                                       & iStart         = iStart          )

      call mus_load_?tname?_data_chunk( conf           = conf,           &
                                       & parent         = parent,      &
                                       & particleBuffer = particleBuffer, &
                                       & nChunkVals     = nChunkVals,     &
                                       & key            = 'force',        &
                                       & iStart         = iStart          )

      call mus_load_?tname?_data_chunk( conf           = conf,           &
                                       & parent         = parent,      &
                                       & particleBuffer = particleBuffer, &
                                       & nChunkVals     = nChunkVals,     &
                                       & key            = 'radius',       &
                                       & iStart         = iStart          )

      call mus_load_?tname?_data_chunk( conf           = conf,           &
                                       & parent         = parent,      &
                                       & particleBuffer = particleBuffer, &
                                       & nChunkVals     = nChunkVals,     &
                                       & key            = 'mass',         &
                                       & iStart         = iStart          )


      ! nChunkVals is the number of particles actually read from the lua file
      ! Can be less than chunkSize if there are fewer than chunkSize particles
      ! left in lua table
      ! write(stdOutUnit,*) 'nChunkVals = ', nChunkVals

      !-- SET PARTICLE ROTATIONAL INERTIA --!
      do iParticle = 1, nChunkVals
        particleBuffer( iParticle )%rotInertia &
          & = 0.4 * particleBuffer( iParticle )%mass &
          &       * particleBuffer( iParticle )%radius**2

      end do

      ! Append particles in particleBuffer to particleGroup
      ! write(stdOutUnit,*) 'Loading particles: nChunkVals = ', nChunkVals
      do iParticle = 1, nChunkVals
        wasAdded = .FALSE.
        ! write(stdOutUnit,*) 'Appending particle'

        call append_DA_?tname?( me        = particles,                       &
                              & particle  = particleBuffer(iParticle),   &
                              & length    = 1,                               &
                              & wasAdded  = wasAdded                         )

      end do

      iStart = iStart + chunkSize
    end do chunkLoop

    deallocate( particleBuffer )

  end subroutine load_?tname?_data


?? end text loadParticleData_txt

?? text :: loadParticleDataChunk_txt(tname)

  subroutine mus_load_?tname?_data_chunk( conf, parent, particleBuffer, &
                                         & nChunkVals, key, iStart       )
    
    type(flu_State), intent(in) :: conf

    !> Handle to parent table if position, velocity tables etc are inside
    !! another table
    integer, intent(in), optional :: parent

    ! Buffer to hold the data read from chunk
    type(mus_?tname?_type), allocatable, intent(inout) :: particleBuffer(:)
    
    ! Number of values actually read from table
    ! If entire buffer is filled, nChunkVals = size(particleBuffer)
    integer, intent(out) :: nChunkVals

    ! key = 'position', 'velocity', 'radius' or 'mass'
    character(len=*), intent(in) :: key
    ! Index to start reading the lua array at
    integer, intent(in) :: iStart

    !--------------------------------------------!
    integer :: iError
    integer :: vErr(6)
    real(kind=rk) :: vecBuffer(6)
    real(kind=rk) :: realBuffer

    integer :: nChunk

    integer :: nVals
    integer :: thandle
    integer :: iParticle, iChunk
    !--------------------------------------------!
    iError = 0
    nChunkVals = 0
    
    if( .NOT.( allocated( particleBuffer ) ) ) then
      write(stdOutUnit,*) 'ERROR mus_load_particle_data_chunk:'
      call tem_abort( 'particleBuffer not allocated')
    end if

    ! Number of elements to read = size of particleBuffer
    nChunk = size(particleBuffer)

    ! Determine whether we're loading a vector or scalar quantity
    if( key == 'position' .OR. key == 'velocity' .OR. key == 'force' ) then
      ! initialize vector particle quantity
      ! call aot_table_open(L = conf, thandle = thandle, key = key)
      call aot_table_open(L       = conf,    &
        &                 thandle = thandle, &
        &                 parent  = parent,  &
        &                 key     = key      )
      if (thandle /= 0) then
        ! get the number of position vectors, should be equal to Nparticles
        nVals = aot_table_length(L=conf, thandle=thandle)
        
        ! Loop over all the position vals
        iChunk = 1
        do iParticle=iStart, iStart + nChunk - 1

          ! Check to make sure we haven't reached the end of the lua table
          if( iParticle > nVals ) exit 

          ! Get the entire position vector
          call aot_get_val( L = conf, thandle = thandle, &
                          & pos = iParticle, val = vecBuffer, &
                          & ErrCode = vErr, &
                          & default = [0.0_double_k, 0.0_double_k, &
                          &            0.0_double_k, 0.0_double_k, &
                          &            0.0_double_k, 0.0_double_k ])
          
          if (btest(iError, aoterr_Fatal)) then
            write(*,*) 'FATAL Error occured, while retrieving particle pos'
            if (btest(iError, aoterr_NonExistent)) write(*,*) &
              &  'Variable not existent!'
            if (btest(iError, aoterr_WrongType)) write(*,*) &
              &  'Variable has wrong type!'
          else
            if (btest(iError, aoterr_NonExistent)) write(*,*) &
              &  'Variable not set in' &
              &  // ' config, Using default value!'
          end if

          ! Copy property to the particleBuffer
          if( key == 'position' ) then
            particleBuffer(iChunk)%pos(1:6) = vecBuffer(1:6)
            particleBuffer(iChunk)%particleID = iParticle
          else if( key == 'velocity' ) then
            particleBuffer(iChunk)%vel(1:6) = vecBuffer(1:6)
          else if( key == 'force' ) then
            particleBuffer(iChunk)%Fext(1:6) = vecBuffer(1:6)
          end if

          iChunk = iChunk + 1

          ! Increment nChunkVals to signify we've read another chunk
          nChunkVals = nChunkVals + 1
        end do
        
      end if
      call aot_table_close(L = conf, thandle = thandle)


    else if ( key == 'radius' .OR. key == 'mass' ) then
      ! initialize scalar particle quantity
      ! call aot_table_open(L = conf, thandle = thandle, key = key)
      call aot_table_open(L       = conf,    &
        &                 thandle = thandle, &
        &                 parent  = parent,  &
        &                 key     = key      )
      if (thandle /= 0) then
        ! get the number of position vectors, should be equal to Nparticles
        nVals = aot_table_length(L=conf, thandle=thandle)
        
        iChunk = 1
        do iParticle = iStart, iStart + nChunk - 1
          
          ! Check to make sure we haven't reached the end of the lua table
          if( iParticle > nVals ) exit 
          
          call aot_get_val(L = conf, thandle = thandle, &
            &              val = realBuffer, ErrCode = iError, &
            &              pos = iParticle)
          if (btest(iError, aoterr_Fatal)) then
            write(*,*) 'FATAL Error occured, while retrieving radius'
            if (btest(iError, aoterr_NonExistent)) write(*,*) &
              &  'Variable not existent!'
            if (btest(iError, aoterr_WrongType)) write(*,*) &
              &  'Variable has wrong type!'
          else
            if (btest(iError, aoterr_NonExistent)) write(*,*) &
              &  'Variable not set in' &
              &  // ' config, Using default value!'
          end if

          ! Copy property to the particleBuffer
          if( key == 'radius' ) then
            particleBuffer(iChunk)%radius = realBuffer
          else if ( key == 'mass' ) then
            particleBuffer(iChunk)%mass = realBuffer
          end if
          
          iChunk = iChunk + 1
          nChunkVals = nChunkVals + 1
        
        end do
      end if
      call aot_table_close(L = conf, thandle = thandle)

    else
      write(stdOutUnit,'(A)') &
        & 'Error MUS_LOAD_PARTICLE_DATA: unrecognized property'
    end if ! particle property


  end subroutine mus_load_?tname?_data_chunk


?? end text loadParticleDataChunk_txt

?? text :: createParticlesInBlob_txt(tname)
   subroutine mus_create_predefined_?tname?( particles, particleblob, Nparticles )
    !> Dynamic particle array to append particles read from file to
    type(dyn_?tname?_array_type), intent(inout) :: particles
    !> Particleblob object which represents the shape and location at which the 
    !! blob of particles is positioned
    type(mus_particle_blob_cylinder_type), intent(in) :: particleblob 
    !> Desired number of particles
    integer, intent(in) :: Nparticles
    ! ------------------------------------------------ !
    ! ------------------------------------------------ !
    



   end subroutine mus_create_predefined_?tname?
?? end text createParticlesInBlob_txt