! Copyright (c) 2013, 2015-2016, 2019-2020 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2013-2015 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
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
! ****************************************************************************** !
!> This module contains routines which initiliaze advection relaxation and
!! flow field for lbm incompressible model.
module mus_initFluidIncomp_module
  ! include treelm modules
  use env_module,         only: labelLen
  use tem_aux_module,     only: tem_abort
  use tem_logging_module, only: logUnit

  ! include musubi modules
  use mus_bgk_module,   only: bgk_advRel_generic
  use mus_d3q19_module, only: bgk_advRel_d3q19_incomp, &
    &                         trt_advRel_d3q19_incomp
  use mus_d3q27_module, only: bgk_advRel_d3q27
  use mus_d2q9_module,  only: bgk_advRel_d2q9_incomp, mrt_advRel_d2q9_incomp
  use mus_mrt_d3q19_module,   only: mrt_advRel_d3q19_incomp,         &
    &                         mrt_advRel_d3q19_incomp_generic, &
    &                         mrt_advRel_generic
  use mus_mrt_d3q27_module, only: weighted_mrt_advRel_d3q27_incomp, &
    &                             weighted_mrt_advRel_d3q27_generic
  use mus_scheme_type_module, only: kernel

  implicit none

  private

  public :: mus_init_advRel_fluidIncomp

contains

! **************************************************************************** !
  !> Initialize the relaxation model for lbm incompressible model
  subroutine mus_init_advRel_fluidIncomp( relaxation, layout, compute )
    ! --------------------------------------------------------------------------
    character(len=labelLen), intent(inout) :: relaxation
    character(len=labelLen), intent(in) :: layout
    procedure( kernel ), pointer, intent(out) :: compute
    ! --------------------------------------------------------------------------
    write(logUnit(1),*) 'Choosing fluid_incompressible relaxation model: ' &
      &                 // trim(relaxation) // ' for layout ' // trim(layout)

    select case (trim(relaxation))
    case ('bgk_generic')
      compute => bgk_advRel_generic

    case ('mrt_bgk')
      compute => mrt_advRel_generic

    case ('bgk')
      select case (trim(layout))
      case ('d3q27')
        compute => bgk_advRel_d3q27
      case ('d3q19')
        compute => bgk_advRel_d3q19_incomp
      case ('d2q9')
        compute => bgk_advRel_d2q9_incomp
      case default
        compute => bgk_advRel_generic
      end select

    case ('mrt_generic')
      select case (trim(layout))
      case ('d3q27')
        compute => weighted_mrt_advRel_d3q27_generic
      case ('d3q19')
        compute => mrt_advRel_d3q19_incomp_generic
      case default
        compute => mrt_advRel_generic
      end select

    case ('mrt_weighted')
      if ( trim(layout) == 'd3q27' ) then
        compute => weighted_mrt_advRel_d3q27_incomp
      else
        write(logUnit(1),*) 'Stencil '//trim(layout)//' is not supported yet!'
        call tem_abort()
      end if

    case ('mrt')
      select case (trim(layout))
      case ('d3q27')
        write(logUnit(1),*) 'Stencil '//trim(layout)//' is not supported yet!'
        write(logUnit(1),*) 'But you can use "mrt_weighted" !'
        call tem_abort()
      case ('d3q19')
        compute => mrt_advRel_d3q19_incomp
      case ('d2q9')
        compute => mrt_advRel_d2q9_incomp
      case default
        write(logUnit(1),*) 'Stencil '//trim(layout)//' is not supported yet!'
        call tem_abort()
      end select

    case ('trt')
      select case (trim(layout))
      case ('d3q19')
        compute => trt_advRel_d3q19_incomp
      case default
        write(logUnit(1),*) 'Stencil '//trim(layout)//' is not supported yet!'
        call tem_abort()
      end select
    case default
      write(logUnit(1),*) 'Relaxation '//trim(relaxation)//' is not supported!'
      call tem_abort()
    end select

  end subroutine mus_init_advRel_fluidIncomp
! **************************************************************************** !

end module mus_initFluidIncomp_module
