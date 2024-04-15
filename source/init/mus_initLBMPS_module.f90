! Copyright (c) 2013 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2015 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
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
!> This module contains routines which assigns compute kernel for passive scalar
!! model
module mus_initLBMPS_module

  ! include treelm modules
  use env_module,         only: labelLen
  use tem_aux_module,     only: tem_abort
  use tem_logging_module, only: logUnit

  ! include musubi modules
  use mus_bgk_module,         only: bgk_advrel_flekkoy,                        &
    &                               bgk_advrel_flekkoy_noFluid,                &
    &                               bgk_advRel_flekkoy_2nd,                    &
    &                               trt_advRel_flekkoy_general
  use mus_scheme_type_module, only: kernel

  implicit none

  private

  public :: mus_init_advRel_lbm_ps

contains

! ****************************************************************************** !
  !> Initialize the relaxation model for lbm passive scalar scheme kind
  subroutine mus_init_advRel_lbm_ps( relaxation, layout, order, compute )
    ! ---------------------------------------------------------------------------
    character(len=labelLen), intent(in) :: relaxation
    character(len=labelLen), intent(in) :: layout
    character(len=labelLen), intent(in) :: order
    procedure( kernel ), pointer, intent(out) :: compute
    ! ---------------------------------------------------------------------------

    write(logUnit(1),*) 'Choosing LBM Passive Scalar relaxation model: '//     &
      &                 trim(relaxation)

    select case( trim(relaxation) )
    case( 'bgk' )
      select case( trim(order) )
      case( 'first' )
        compute => bgk_advRel_flekkoy
      case( 'second' )
        compute => bgk_advRel_flekkoy_2nd
      case default
        write(logUnit(1),*) 'Order '//trim(order)//' is not supported yet!'
        call tem_abort()
      end select
    case( 'trt' )
      write(*, *) 'using trt_advRel scheme.'
      compute => trt_advRel_flekkoy_general
    case( 'bgk_noFluid' )
      select case( trim(layout) )
      case( 'flekkoy' )
        compute => bgk_advRel_flekkoy_noFluid
      case default
        write(logUnit(1),*) 'Stencil '//trim(layout)//' is not supported yet!'
        call tem_abort()
      end select
    case default
      write(logUnit(1),*) 'The selected relaxation model is not supported: '// &
        &                  trim(relaxation)
      call tem_abort()
    end select

  end subroutine mus_init_advRel_lbm_ps
! ****************************************************************************** !

end module mus_initLBMPS_module
