! See copyright notice in the COPYRIGHT file.
!> This module contains macro function to test macro access
!! definition defined in lbm_macros.inc
?? include 'header/lbm_macros.inc'
module mus_macroAccess_module

  implicit none

  contains

  ! MACRO functions
  pure function NGPOS(iDir, node, nElems) 
    integer, intent(in) :: iDir, node, nElems
    integer :: NGPOS
    NGPOS = ?NGPOS?(iDir, node, nElems)
  end function

  pure function IDX(varPosiDir, node, nScalars, nElems)
    integer, intent(in) :: varPosiDir, node, nScalars, nElems
    integer :: IDX
    IDX = 0
    if (nElems < 0) RETURN
    IDX = ?IDX?(varPosiDir, node, nScalars, nElems)
  end function

  pure function FETCH(iDir, iField, node, QQ, nScalars, nElems, neigh)
    integer, intent(in) :: iDir, iField, node, QQ, nScalars, neigh(:), nElems
    integer :: FETCH
!    FETCH = neigh( NGPOS(iDir, node, QQ) ) + (iField-1)*QQ
    FETCH = ?FETCH?(iDir, iField, node, QQ, nScalars, nElems, neigh)
  end function

  pure function SAVE(iDir, iField, node, QQ, nScalars, nElems, neigh)
    integer, intent(in) :: iDir, iField, node, QQ, nScalars, neigh(:), nElems
    integer :: SAVE
    SAVE = ?SAVE?(iDir, iField, node, QQ, nScalars, nElems, neigh)
  end function

!   pure function ELEMIDX(index, QQ, nElems)
!     integer, intent(in) :: index, QQ, nElems
!     integer :: ELEMIDX
!     ELEMIDX = ?ELEMIDX?(index, QQ, nElems)
!   end function
!
!   pure function DIRIDX(index, QQ, nElems)
!     integer, intent(in) :: index, QQ, nElems
!     integer :: DIRIDX
!     DIRIDX = ?DIRIDX?(index, QQ, nElems )
!   end function
  ! END MACRO functions

end module mus_macroAccess_module
