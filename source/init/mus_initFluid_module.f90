! Copyright (c) 2013, 2017, 2019-2020 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2014-2016 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2014 Kartik Jain <kartik.jain@uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2017 Raphael Haupt <raphael.haupt@uni-siegen.de>
! Copyright (c) 2021 Gregorio Gerardo Spinelli <gregoriogerardo.spinelli@dlr.de>
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
!> This module contains routines which initiliaze advection relaxation and
!! flow field for lbm model.
module mus_initFluid_module
  use env_module,         only: labelLen
  use tem_aux_module,     only: tem_abort
  use tem_logging_module, only: logUnit

  use mus_bgk_module,       only: bgk_advRel_generic
  use mus_compute_cumulant_module,  only: cumulant_d3q27, cascaded_d3q27,  &
    &                                     cumulant_d3q27_extended_generic, &
    &                                     cumulant_d3q27_extended_fast
  use mus_d3q27_module,     only: bgk_advRel_d3q27, &
    &                             trt_advRel_d3q27, &
    &                             bgk_improved_advRel_d3q27, &
    &                             bgk_Regularized_d3q27, &
    &                             bgk_RecursiveRegularized_d3q27, &
    &                             bgk_ProjectedRecursiveRegularized_d3q27, &
    &                             bgk_HybridRecursiveRegularized_d3q27, &
    &                             bgk_DualRelaxationTime_RR_d3q27, &
    &                             bgk_HybridRecursiveRegularizedCorr_d3q27
  use mus_d3q19_module,     only: bgk_advRel_d3q19,       &
    &                             bgk_advRel_d3q19_block, &
    &                             trt_advRel_d3q19,  &
    &                             bgk_Regularized_d3q19, &
    &                             bgk_RecursiveRegularized_d3q19, &
    &                             bgk_ProjectedRecursiveRegularized_d3q19, &
    &                             bgk_HybridRecursiveRegularized_d3q19, &
    &                             bgk_DualRelaxationTime_RR_d3q19, &
    &                             bgk_HybridRecursiveRegularizedCorr_d3q19
  use mus_d2q9_module,      only: mrt_advRel_d2q9, bgk_advRel_d2q9, &
    &                             bgk_improved_advRel_d2q9, &
    &                             bgk_Regularized_d2q9, &
    &                             bgk_RecursiveRegularized_d2q9, &
    &                             bgk_ProjectedRecursiveRegularized_d2q9, &
    &                             bgk_HybridRecursiveRegularized_d2q9, &
    &                             bgk_DualRelaxationTime_RR_d2q9, &
    &                             bgk_HybridRecursiveRegularizedCorr_d2q9
  use mus_mrt_d3q19_module,       only: mrt_advRel_d3q19,         &
    &                             mrt_advRel_d3q19_generic, &
    &                             mrt_advRel_generic
  use mus_mrt_d3q27_module, only: weighted_mrt_advRel_d3q27, &
    &                             weighted_mrt_advRel_d3q27_generic
  use mus_test_module,      only: vec_fma
  use mus_scheme_type_module, only: kernel

  implicit none

  private

  public :: mus_init_advRel_fluid

contains

  ! ************************************************************************** !
  !> Assigning compute kernel routine by scheme relaxation type for fluid kind.
  !!
  subroutine mus_init_advRel_fluid( relaxation, layout, compute )
    ! --------------------------------------------------------------------------
    character(len=labelLen), intent(inout) :: relaxation
    character(len=labelLen), intent(in) :: layout
    procedure( kernel ), pointer, intent(out) :: compute
    ! --------------------------------------------------------------------------

    write(logUnit(1),*) 'Choosing fluid relaxation model: '                 &
      &                 // trim(relaxation) // ' for layout ' // trim(layout)

    select case (trim(relaxation))
    case ('bgk_generic')
      compute => bgk_advRel_generic

    case ('bgk_improved')
      select case (trim(layout))
      case ('d3q27')
        compute => bgk_improved_advRel_d3q27
      case ('d2q9')
        compute => bgk_improved_advRel_d2q9
      case default
        write(logUnit(1),*) 'Stencil '//trim(layout)//' is not supported yet!'
        call tem_abort()
      end select

    case ('bgk_block')
      if ( trim(layout) == 'd3q19' ) then
        compute => bgk_advRel_d3q19_block
      else
        write(logUnit(1),*) 'Stencil '//trim(layout)//' is not supported yet!'
        call tem_abort()
      end if

    case ('mrt_bgk')
      compute => mrt_advRel_generic

    case ('bgk')
      select case (trim(layout))
      case ('d3q27')
        compute => bgk_advRel_d3q27
      case ('d3q19')
        compute => bgk_advRel_d3q19
      case ('d2q9')
        compute => bgk_advRel_d2q9
      case default
        compute => bgk_advRel_generic
      end select

    case ('mrt_generic')
      select case (trim(layout))
      case ('d3q27')
        compute => weighted_mrt_advRel_d3q27_generic
      case ('d3q19')
        compute => mrt_advRel_d3q19_generic
      case default
        compute => mrt_advRel_generic
      end select

    case ('mrt')
      select case (trim(layout))
      case ('d3q27')
        compute => weighted_mrt_advRel_d3q27
      case ('d3q19')
        compute => mrt_advRel_d3q19
      case ('d2q9')
        compute => mrt_advRel_d2q9
      case default
        write(logUnit(1),*) 'Stencil '//trim(layout)//' is not supported yet!'
        call tem_abort()
      end select

    case ('trt')
      select case (trim(layout))
      case ('d3q19')
        compute => trt_advRel_d3q19
      case ('d3q27')
        compute => trt_advRel_d3q27
      case default
        write(logUnit(1),*) 'Stencil '//trim(layout)//' is not supported yet!'
        call tem_abort()
      end select

    case ('cumulant')
      select case (trim(layout))
      case ('d3q27')
        compute => cumulant_d3q27
      case default
        write(logUnit(1),*) 'Stencil '//trim(layout)//' is not supported yet!'
        call tem_abort()
      end select

    case ('cumulant_extended')
      select case (trim(layout))
      case ('d3q27')
        compute => cumulant_d3q27_extended_fast
      case default
        write(logUnit(1),*) 'Stencil '//trim(layout)//' is not supported yet!'
        call tem_abort()
      end select
    
    case ('cumulant_extended_generic')
      select case (trim(layout))
      case ('d3q27')
        compute => cumulant_d3q27_extended_generic
      case default
        write(logUnit(1),*) 'Stencil '//trim(layout)//' is not supported yet!'
        call tem_abort()
      end select
      
    case ('cascaded')
      select case (trim(layout))
      case ('d3q27')
        compute => cascaded_d3q27
      case default
        write(logUnit(1),*) 'Stencil '//trim(layout)//' is not supported yet!'
        call tem_abort()
      end select

    case ('vec_fma', 'test')
      compute => vec_fma

    case ('hrr_bgk')
      select case (trim(layout))
      case ('d3q27')
        compute => bgk_HybridRecursiveRegularized_d3q27
      case ('d3q19')
        compute => bgk_HybridRecursiveRegularized_d3q19
      case ('d2q9')
        compute => bgk_HybridRecursiveRegularized_d2q9
      case default
        write(logUnit(1),*) 'Stencil '//trim(layout)//' is not supported yet!'
        call tem_abort()
      end select

    case ('hrr_bgk_corrected', 'prr_bgk_corrected', 'rr_bgk_corrected')
      select case (trim(layout))
      case ('d3q27')
        compute => bgk_HybridRecursiveRegularizedCorr_d3q27
      case ('d3q19')
        compute => bgk_HybridRecursiveRegularizedCorr_d3q19
      case ('d2q9')
        compute => bgk_HybridRecursiveRegularizedCorr_d2q9
      case default
        write(logUnit(1),*) 'Stencil '//trim(layout)//' is not supported yet!'
        call tem_abort()
      end select

    case ('drt_bgk')
      select case (trim(layout))
      case ('d3q27')
        compute => bgk_DualRelaxationTime_RR_d3q27
      case ('d3q19')
        compute => bgk_DualRelaxationTime_RR_d3q19
      case ('d2q9')
        compute => bgk_DualRelaxationTime_RR_d2q9
      case default
        write(logUnit(1),*) 'Stencil '//trim(layout)//' is not supported yet!'
        call tem_abort()
      end select

    case ('rr_bgk')
      select case (trim(layout))
      case ('d3q27')
        compute => bgk_RecursiveRegularized_d3q27
      case ('d3q19')
        compute => bgk_RecursiveRegularized_d3q19
      case ('d2q9')
        compute => bgk_RecursiveRegularized_d2q9
      case default
        write(logUnit(1),*) 'Stencil '//trim(layout)//' is not supported yet!'
        call tem_abort()
      end select

    case ('prr_bgk')
      select case (trim(layout))
      case ('d3q27')
        compute => bgk_ProjectedRecursiveRegularized_d3q27
      case ('d3q19')
        compute => bgk_ProjectedRecursiveRegularized_d3q19
      case ('d2q9')
        compute => bgk_ProjectedRecursiveRegularized_d2q9
      case default
        write(logUnit(1),*) 'Stencil '//trim(layout)//' is not supported yet!'
        call tem_abort()
      end select

    case ('r_bgk')
      select case (trim(layout))
      case ('d3q27')
        compute => bgk_Regularized_d3q27
      case ('d3q19')
        compute => bgk_Regularized_d3q19
      case ('d2q9')
        compute => bgk_Regularized_d2q9
      case default
        write(logUnit(1),*) 'Stencil '//trim(layout)//' is not supported yet!'
        call tem_abort()
      end select

    case default
      write(logUnit(1),*) 'Relaxation '//trim(relaxation)//' is not supported!'
      call tem_abort()
    end select

  end subroutine mus_init_advRel_fluid
  ! **************************************************************************** !

end module mus_initFluid_module
! ****************************************************************************** !
