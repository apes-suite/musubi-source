! Copyright (c) 2020 Harald Klimach <harald.klimach@uni-siegen.de>
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are met:
!
! 1. Redistributions of source code must retain the above copyright notice, this
! list of conditions and the following disclaimer.
!
! 2. Redistributions in binary form must reproduce the above copyright notice,
! this list of conditions and the following disclaimer in the documentation
! and/or other materials provided with the distribution.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
! AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
! DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
! FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
! DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
! SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
! OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
! OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
! ---------------------------------------------------------------------------- !
?? include 'header/lbm_macros.inc'

!> Module to provide information about pre-processing settings of the
!! current executable.
module mus_ppInfo_module
  use tem_logging_module, only: logUnit
  use tem_tools_module, only: tem_horizontalSpacer

  implicit none

  private

  public :: mus_print_ppInfo


contains


  !> Print information on the pre-processor options of the executable.
  subroutine mus_print_ppInfo()
    call tem_horizontalSpacer(fUnit = logUnit(1))
    write(logUnit(1),*) '| PRE-PROCESSING SETTINGS:'
??IF (SOA) THEN
    write(logUnit(1),*) '| Using SOA (structure of arrays data layout)'
??ELSE
    write(logUnit(1),*) '| Using AOS (array of structures data layout)'
??ENDIF

??IF (PUSH) THEN
    write(logUnit(1),*) '| Using PUSH for streaming'
??ELSE
    write(logUnit(1),*) '| Using PULL for streaming'
??ENDIF
    call tem_horizontalSpacer(fUnit = logUnit(1))
  end subroutine mus_print_ppInfo

end module mus_ppInfo_module
