! See copyright notice in the COPYRIGHT file.
program memAlloc
  use env_module, only: rk, long_k
  use tem_debug_module, only: tem_reportStatus, tem_debug_type
  use tem_logging_module, only: tem_logging_init
  implicit none
  real(kind=rk), allocatable :: rArray(:)
  integer(kind=long_k), allocatable :: iArray(:)
  real(kind=rk) :: memSize, nBytes
  type(tem_debug_type) :: debug
  integer :: iSize, iType
  integer(kind=long_k) :: nEntries
  logical :: error
  character(len=64) :: buffer, bufType
  real(kind=rk), allocatable :: rTest(:)
  
  allocate( rTest( 4))
  rTest = [ 4._rk, 16._rk, 24._rk, 32._rk ]
  error = .false.

  ! initialize the debug type
  call tem_logging_init( me    = debug%logger,                                 &
    &                    level = 10,                                           &
    &                    rank  = 0,                                            &
    &                    filename = 'memTrace.res' )

  do iType = 2,1, -1
    if( iType .eq. 1 ) then
      bufType = 'rk'
      nBytes = 8._rk
    else
      bufType = 'longk'
      nBytes = 8._rk
    end if
!    do iSize =14, 32, 1
    do iSize = 1, size( rTest )
!      nEntries = int(2,kind=long_k)**iSize
      nEntries = int( rTest(iSize)*1.E9_rk/nBytes, long_k )
      memSize = rTest(iSize )*1000._rk !real(nEntries,kind=rk)*8._rk/1000000._rk
      write(*,*) 'trying to allocate nWords', nEntries, ' size MB', memSize
      write(buffer , '(2a,f20.3)') trim(bufType), ' before allocating MB ', memSize
      call tem_reportStatus( debug = debug, level = 1,                         &
                             text = buffer, string = 'VmRSS')
      flush( debug%logger%funit )
   
      if ( iType .eq. 1 ) then
        allocate( rarray( nEntries ))
        rarray( 1:nEntries )  = 0._rk
        if( .not. allocated( rarray )) then
          error = .true. 
          write(*,*) 'Error: not allocated! '
        else
        end if
      else
        allocate( iarray( nEntries ))
        iarray( 1:nEntries ) = 0_long_k
        if( .not. allocated( iarray )) then
          error = .true. 
          write(*,*) 'Error: not allocated! '
        else
        end if
      end if
      write(buffer , '(2a,f20.3)') trim( bufType ), ' after  allocating MB ', memSize
      call tem_reportStatus( debug = debug, level = 1,                         &
                             text = buffer, string = 'VmRSS')
      flush( debug%logger%funit )
   
      if ( iType .eq. 1 ) then
        deallocate( rarray )
      else
        deallocate( iarray )
      end if
      write(*,*) 'done.'
    end do
  end do
  close( debug%logger%funit )
  if( error ) then
    write(*,*) 'FAILED'
    stop -1
  else
    write(*,*) 'PASSED'
    stop 0 
  end if

end program memAlloc
