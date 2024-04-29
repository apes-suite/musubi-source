! Copyright (c) 2012-2013, 2015-2017, 2020 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2012-2014 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2012, 2014 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2012 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2012-2013 Manuel Hasert <m.hasert@grs-sim.de>
! Copyright (c) 2014 Kartik Jain <kartik.jain@uni-siegen.de>
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
!> This module contains scheme property type and module related to scheme prop
module mus_scheme_header_module

  ! include treelm modules
  use env_module,         only: LabelLen
  use tem_tools_module,   only: tem_horizontalSpacer
  use tem_aux_module,     only: tem_abort
  use tem_logging_module, only: logUnit

  ! include aotus modules
  use aotus_module,     only: flu_State, aot_get_val, aoterr_NonExistent
  use aot_table_module, only: aot_table_open, aot_table_close
  use aot_out_module,   only: aot_out_type, aot_out_val, aot_out_open_table,   &
    &                         aot_out_close_table

  implicit none

  private

  public :: mus_scheme_header_type, mus_scheme_header_out
  public :: mus_load_scheme_header

  !> Datatype containing information to identify the scheme
  !!
  !! Combination of scheme kind, relaxation and layout%stencilKind
  !! are used to choose the correct compute kernel for the
  !! scheme
  !! 
  !!> | type | options |
  !!> |:-----------------|:--------------|
  !!> | kind | **fluid** (default)             |
  !!> |      | **fluid_incompressible**        |
  !!> |      | **isotherm_acEq**               |
  !!> |      | **multispecies_gas**            |
  !!> |      | **multispecies_liquid**         |
  !!> |      | **nernst_planck**               |
  !!> |      | **passive_scalar**              |
  !!> |      | **poisson**                     |
  !!> |      | **poisson_boltzmann_linear**    |
  !!> |      | **poisson_boltzmann_nonlinear** |
  !!> |------|---------------------------------|
  !!> | layout | **D2Q9**           |
  !!> |        | **D3Q19** (default)|
  !!> |        | **D3Q27**          |
  !!> |        | _D1Q3_             |
  !!> |        | _D2Q5_             |
  !!> |        | _D3Q6_             |
  !!> |        | _D3Q7_             |
  !!> |        | _D3Q13_            |
  !!> |        | _D3Q15_            |
  !!> |        | _flekkoy_          |
  !!> |--------|--------------------|
  !!> | relaxation | **BGK**  (default)       |
  !!> |            | **MRT**                  |
  !!> |            | **TRT**                  |
  !!> |            | **bgk_withthermodynfac** |
  !!> |            | **mrt_withthermodynfac** |
  !!> |            | _mrt_bgk_                |
  !!> |            | _mrt_generic_            |
  !!> |            | _bgk_generic_            |
  !!> |            | _bgk_improved_           |
  !!> |            | _bgk_block_              |
  !!> |            | _cumulant_               |
  !!> |            | _cascaded_               |
  !!> |            | _vec_fma_                |
  !!> |            | _test_                   |
  !!> |            | _bgk_noFluid_            |
  !!> |------------|--------------------------|
  type mus_scheme_header_type
    !> scheme kind, Ex: fluid, fluid_incompressible, multispecies_gas, 
    !! multispecies_liquid, poisson, poisson_boltzmann_linear, 
    !! poisson_boltzmann_nonlinear, nernst_planck, isotherm_acEq
    character(len=labelLen) :: kind
    !> scheme layout, Ex: d3q19
    character(len=labelLen) :: layout
    !> scheme relaxation type Ex: BGK, MRT, bgk_pl, bgk_cy, bgk_cs...
    character(len=labelLen) :: relaxation
    !> equilibrium order Ex: first, second
    character(len=labelLen) :: order
  end type mus_scheme_header_type

contains

! ****************************************************************************** !
  !> load scheme header info from lua file identify table or from scheme table
  !! or from config
  !!
  !! Load scheme label, kind, layoutKind and relaxation
  !!```lua
  !! identify = { kind = 'simType',
  !!              layout = 'stencilLayout',
  !!              relaxation = 'relaxationType' }
  !!```
  !! For a possible selection of the other parameters
  !! - simType: fluid, fluid_incompressible, multispecies_liquid
  !! - [[mus_scheme_layout_module]]: d2q9, d3q19, ...
  !! - relaxationType: bgk, mrt, ...
  !!
  subroutine mus_load_scheme_header( me, conf, parent, scaling )
  ! -----------------------------------------------------------------------------
    !> returns scheme identify information
    type( mus_scheme_header_type ), intent(out) :: me
    type( flu_State ) :: conf !< flu state
    !> parent handle if scheme table is defined
    integer, intent(in), optional :: parent
    !> scaling: diffusive or acoustic?
    character(len=*), intent(in) :: scaling
    ! ---------------------------------------------------------------------------
    integer :: thandle !< handle for scheme identify table
    integer :: iError
    ! ---------------------------------------------------------------------------

    call tem_horizontalSpacer(fUnit = logUnit(1))
    write(logUnit(1),*) 'Loading Scheme identify table: '
    !default values
    me%kind = 'fluid'
    me%layout = 'd3q19'
    me%relaxation = 'bgk'

    call aot_table_open( L       = conf,      &
      &                  parent  = parent,    &
      &                  thandle = thandle,   &
      &                  key     = 'identify' )

    if( thandle .gt. 0 ) then
      ! get schemekind
      call aot_get_val( L       = conf,    &
        &               thandle = thandle, &
        &               key     = 'kind',  &
        &               val     = me%kind, &
        &               default = 'fluid', &
        &               ErrCode = iError   )

      ! get layoutkind
      call aot_get_val( L       = conf,                                        &
        &               thandle = thandle,                                     &
        &               key     = 'layout',                                    &
        &               val     = me%layout,                                   &
        &               default = 'd3q19',                                     &
        &               ErrCode = iError )

      ! get relaxation
      call aot_get_val( L       = conf,                                        &
        &               thandle = thandle,                                     &
        &               key     = 'relaxation',                                &
        &               val     = me%relaxation,                               &
        &               default = 'bgk',                                       &
        &               ErrCode = iError )
      
      ! get order of equilibrium
      call aot_get_val( L       = conf,                                        &
        &               thandle = thandle,                                     &
        &               key     = 'order',                                     &
        &               val     = me%order,                                    &
        &               default = 'second',                                    &
        &               ErrCode = iError ) 

    else
      write(logUnit(1),*) 'Scheme Identify table not defined.'
      write(logUnit(1),*) 'Setting default values for scheme..'
    endif

    call aot_table_close( L=conf, thandle=thandle )

    write(logUnit(1),*) 'kind: '// trim( me%kind )
    write(logUnit(1),*) 'Layout: '// trim( me%layout )
    write(logUnit(1),*) 'relaxation: '// trim( me%relaxation )
    write(logUnit(1),*) 'order: '// trim( me%order )
    call tem_horizontalSpacer(fUnit = logUnit(1))

    ! Both multispeciees and poisson equation must have diffusive scaling
    ! since diffusive scaling is used to recover macroscopic equations
    ! from asymptotic analysis
    select case(trim(me%kind))
    case ('fluid', 'fluid_incompressible', 'isotherm_acEq' )
      if(trim(scaling) /= 'acoustic') then
         write(logUnit(1),*)'ERROR: Choose scaling = "acoustic" for ' &
           &                // trim(me%kind)
         write(logUnit(1),*)'Aborting'
         call tem_abort()
      end if
    case ( 'multispecies_gas', 'multispecies_liquid', 'nernst_planck', &
      &    'passive_scalar, ''poisson', 'poisson_boltzmann_linear',    &
      &    'poisson_boltzmann_nonlinear'                               )
      if(trim(scaling) /= 'diffusive') then
         write(logUnit(1),*)'ERROR: Choose scaling = "diffusive" for ' &
           &                // trim(me%kind)
         write(logUnit(1),*)'Aborting'
         call tem_abort()
      end if
    end select

  end subroutine mus_load_scheme_header
! ****************************************************************************** !

! ****************************************************************************** !
  !> Dumps scheme header
  subroutine mus_scheme_header_out( me, conf )
  ! -----------------------------------------------------------------------------
    !> returns scheme identify information
    type( mus_scheme_header_type ), intent(in) :: me
    type(aot_out_type) :: conf
    ! ---------------------------------------------------------------------------

    ! the label does not have to be outputted.
    ! because it is the name of the outputted varSys
    call aot_out_open_table( put_conf = conf, tname = 'identify' )
    call aot_out_val( put_conf = conf,                                         &
      &               vname    = 'kind',                                       &
      &               val      = trim( me%kind ))
    call aot_out_val( put_conf = conf,                                         &
      &               vname    = 'relaxation',                                 &
      &               val      = trim( me%relaxation ))
    call aot_out_val( put_conf = conf,                                         &
      &               vname    = 'layout',                                     &
      &               val      = trim( me%layout ))
    call aot_out_close_table( put_conf = conf )
  end subroutine mus_scheme_header_out
! ****************************************************************************** !

end module mus_scheme_header_module
! ****************************************************************************** !
