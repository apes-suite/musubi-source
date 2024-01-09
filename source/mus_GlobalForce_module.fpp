! Copyright (c) 2023 Gregorio Gerardo Spinelli <gregoriogerardo.spinelli@dlr.de>
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
!> This module keeps all information about the Global force
!! It applies the glocalForce to the outstate.
!!
!! author: Gregorio Gerardo Spinelli
?? include 'header/lbm_macros.inc'
module mus_GlobalForce_module
    use iso_c_binding, only: c_f_pointer
    
    use mpi
    ! include treelm modules
    use env_module,           only: rk
    use tem_param_module,     only: cs2, cs2inv, rho0, cs4inv, div1_2
    use tem_varSys_module,    only: tem_varSys_type
  
    ! include musubi modules
    use mus_scheme_header_module,   only: mus_scheme_header_type
    use mus_derVarPos_module,       only: mus_derVarPos_type
    use mus_physics_module,         only: mus_convertFac_type
    use mus_fluid_module,           only: proc_applyGlobalForce
    use mus_varSys_module,          only: mus_varSys_data_type
    use mus_scheme_layout_module,   only: mus_scheme_layout_type
  
    implicit none
  
    private
  
    public :: mus_assign_globalForce_ptr
    public :: applyGlobalForce
    public :: applyGlobalForce_MRT
    public :: applyGlobalForce_MRT_d2q9
    public :: applyGlobalForce_MRT_d3q19
    public :: applyGlobalForce_MRT_d3q27
    public :: applyGlobalForce_1stOrd
    public :: force_discretization
  
  contains
  
    ! ************************************************************************** !
    !> This routines assigns the pointer applyGlobalForce_ptr to the proper
    !! applyGlobalForce function according to the kind of relaxation and layout
    subroutine mus_assign_globalForce_ptr( schemeHeader, applyGlobalForce_ptr )
      ! ---------------------------------------------------------------------------
      !> scheme header
      type( mus_scheme_header_type ), intent(in) :: schemeHeader
  
      !> pointer to applyGlobalForce
      procedure(proc_applyGlobalForce), pointer, intent(out) :: applyGlobalForce_ptr
      ! ---------------------------------------------------------------------------
  
      if (trim(schemeHeader%kind) == 'fluid' .or. &
        & trim(schemeHeader%kind) == 'fluid_incompressible') then

        select case (trim(schemeHeader%relaxation))

        case ('mrt', 'mrt_bgk','mrt_generic')
          select case (trim(schemeHeader%layout))
          case ('d2q9')
            applyGlobalForce_ptr => applyGlobalForce_MRT_d2q9
          case ('d3q19')
            applyGlobalForce_ptr => applyGlobalForce_MRT_d3q19
          case ('d3q27')
            applyGlobalForce_ptr => applyGlobalForce_MRT_d3q27
          case default
            applyGlobalForce_ptr => applyGlobalForce_MRT
          end select

        case ('r_bgk', 'rr_bgk', 'prr_bgk', 'hrr_bgk', 'rr_bgk_corrected', &
          &   'prr_bgk_corrected', 'hrr_bgk_corrected', 'drt_bgk')
          applyGlobalForce_ptr => applyGlobalForce_HRRAlike

        case default
          applyGlobalForce_ptr => applyGlobalForce
        end select

      else
        applyGlobalForce_ptr => applyGlobalForce_dummy
      end if
  
    end subroutine mus_assign_globalForce_ptr
    ! ************************************************************************** !

! ****************************************************************************** !
  !> Update state with global force variable "force_phy".
  !!
  !! Force term used here is from:
  !! "Discrete lattice effects on the forcing term in the lattice Boltzmann
  !!  method", Zhaoli Guo, Chugung Zheng and Baochang Shi.
  !! In the paper, use force term is referred as Method 2 as:
  !! \[ F_i = w_i( (\vec{e}_i-\vec{u}*)/cs2 +
  !!       (\vec{e}_i \cdot \vec{u}*)\vec{e}_i/cs4) \cdot \vec{F} \]
  !! Force must be defined as body force per unit volume
  !! KM: If this force formula is used then velocity needs to be
  !! computed as u = \sum c_i f_i + \vec{F}/2
  !!
  !! Simuilar to derive routine but it updates the state whereas derive
  !! is used for tracking.
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[proc_applyGlobalForce]].
  subroutine applyGlobalForce_dummy( force_phy, inState, outState, neigh,      &
    &                          auxField, nPdfSize, iLevel, varSys, phyConvFac, &
    &                          derVarPos                                       )
    ! -------------------------------------------------------------------- !
    !> Global Force to be applied on state
    real(kind=rk), intent(in) :: force_phy(3)

    !> input  pdf vector
    real(kind=rk), intent(in) :: inState(:)

    !> output pdf vector
    real(kind=rk), intent(inout) :: outState(:)

    !> connectivity Array corresponding to state vector
    integer,intent(in) :: neigh(:)

    !> auxField array
    real(kind=rk), intent(in) :: auxField(:)

    !> number of elements in state Array
    integer, intent(in) :: nPdfSize

    !> current level
    integer, intent(in) :: iLevel

    !> variable system
    type(tem_varSys_type), intent(in) :: varSys

    !> Physics conversion factor for current level
    type(mus_convertFac_type), intent(in) :: phyConvFac

    !> position of derived quantities in varsys
    type(mus_derVarPos_type), intent(in) :: derVarPos
    ! -------------------------------------------------------------------- !

  end subroutine applyGlobalForce_dummy
! ****************************************************************************** !


! ****************************************************************************** !
  !> Discretize force in velocity space
  !!
  !! Force term used here is from:
  !! "Discrete lattice effects on the forcing term in the lattice Boltzmann
  !!  method", Zhaoli Guo, Chugung Zheng and Baochang Shi.
  !! In the paper, use force term is referred as Method 2 as:
  !! \[ F_i = w_i( (\vec{e}_i-\vec{u}*)/cs2 +
  !!       (\vec{e}_i \cdot \vec{u}*)\vec{e}_i/cs4) \cdot \vec{F} \]
  !! Force must be defined as body force per unit volume
  !! KM: If this force formula is used then velocity needs to be
  !! computed as u = \sum c_i f_i + \vec{F}/2
  !!
  pure function force_discretization( force_phy, velocity, omega, QQ, layout ) &
    &  result ( forceTerm )
    ! -------------------------------------------------------------------- !
    !> Global Force in Lattice units
    real(kind=rk), intent(in) :: force_phy(3)

    !> Velocity
    real(kind=rk), intent(in) :: velocity(3)

    !> Collision frequency
    real(kind=rk), intent(in) :: omega

    !> Stencil size
    integer, intent(in) :: QQ

    !> the layout used
    type( mus_scheme_layout_type ), intent(in) :: layout

    !> moment vector
    real(kind=rk) :: forceTerm(QQ)
    ! -------------------------------------------------------------------- !
    real(kind=rk) :: ucx, uMinusCX(3), omega_fac
    integer :: iDir
    ! ---------------------------------------------------------------------------

    ! adjust omega factor for force
    omega_fac = -div1_2 !1.0_rk - omega * 0.5_rk

    ! force term:
    ! F_i = w_i( (\vec{e}_i-\vec{u}*)/cs2 +
    !       (\vec{e}_i \cdot \vec{u}*)\vec{e}_i/cs4) \cdot \vec{F}
    do iDir = 1, QQ
      ucx = dot_product( layout%fStencil%cxDirRK(:, iDir), &
        &                velocity )
      uMinusCx = layout%fStencil%cxDirRK(:, iDir) - velocity

      forceTerm(iDir) = omega_fac * layout%weight( iDir ) &
        &       * dot_product( uMinusCx * cs2inv          &
        &       + ucx * layout%fStencil%cxDirRK(:,iDir)   &
        &       * cs4inv, force_phy )

    end do

  end function force_discretization
! ****************************************************************************** !


! ****************************************************************************** !
  !> Update state with global force variable "force_phy".
  !!
  !! Force term used here is from:
  !! "Discrete lattice effects on the forcing term in the lattice Boltzmann
  !!  method", Zhaoli Guo, Chugung Zheng and Baochang Shi.
  !! In the paper, use force term is referred as Method 2 as:
  !! \[ F_i = w_i( (\vec{e}_i-\vec{u}*)/cs2 +
  !!       (\vec{e}_i \cdot \vec{u}*)\vec{e}_i/cs4) \cdot \vec{F} \]
  !! Force must be defined as body force per unit volume
  !! KM: If this force formula is used then velocity needs to be
  !! computed as u = \sum c_i f_i + \vec{F}/2
  !!
  !! Simuilar to derive routine but it updates the state whereas derive
  !! is used for tracking.
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[proc_applyGlobalForce]].
  subroutine applyGlobalForce( force_phy, inState, outState, neigh, auxField, &
    &                          nPdfSize, iLevel, varSys, phyConvFac,          &
    &                          derVarPos                                      )
    ! -------------------------------------------------------------------- !
    !> Global Force to be applied on state
    real(kind=rk), intent(in) :: force_phy(3)

    !> input  pdf vector
    real(kind=rk), intent(in) :: inState(:)

    !> output pdf vector
    real(kind=rk), intent(inout) :: outState(:)

    !> connectivity Array corresponding to state vector
    integer,intent(in) :: neigh(:)

    !> auxField array
    real(kind=rk), intent(in) :: auxField(:)

    !> number of elements in state Array
    integer, intent(in) :: nPdfSize

    !> current level
    integer, intent(in) :: iLevel

    !> variable system
    type(tem_varSys_type), intent(in) :: varSys

    !> Physics conversion factor for current level
    type(mus_convertFac_type), intent(in) :: phyConvFac

    !> position of derived quantities in varsys
    type(mus_derVarPos_type), intent(in) :: derVarPos
    ! -------------------------------------------------------------------- !
    type(mus_varSys_data_type), pointer :: fPtr
    real(kind=rk) :: G(3), velocity(3), ucx, uMinusCX(3), forceTerm
    integer :: iElem, iDir, QQ, nScalars, statePos, elemOff
    integer :: vel_pos(3)
    real(kind=rk) :: omega, omega_fac
    ! ---------------------------------------------------------------------------
    ! convert c pointer to solver type fortran pointer
    call c_f_pointer( varSys%method%val( derVarPos%pdf )%method_data, fPtr )

    ! convert physical to lattice
    ! For incompressible model: this forceField should be divided by rho0.
    ! Since rho0 =1, this term is also valid for incompressible model
    G = force_phy / phyConvFac%body_force
    
    associate (layout => fPtr%solverData%scheme%layout,                      &
      &        omega_vec => fPtr%solverData%scheme%field(1)%fieldProp%fluid  &
      &                              %viscKine%omLvl(iLevel)%val,            &
      &        nSolve => fPtr%solverData%scheme%pdf(iLevel)%nElems_local     )

      ! constant parameter
      QQ = layout%fStencil%QQ
  
      nScalars = varSys%nScalars
      ! Position of velocity variable in auxField
      vel_pos = varSys%method%val(derVarPos%velocity)%auxField_varPos(1:3)
  
      do iElem = 1, nSolve
  
        ! element offset
        elemoff = (iElem-1)*varSys%nAuxScalars
        ! obtain velocity from auxField
        velocity(1) = auxField(elemOff + vel_pos(1))
        velocity(2) = auxField(elemOff + vel_pos(2))
        velocity(3) = auxField(elemOff + vel_pos(3))
  
        ! get the correct omega value
        omega = omega_vec(iElem)
        omega_fac = 1.0_rk - omega * 0.5_rk
  
        ! force term:
        ! F_i = w_i( (\vec{e}_i-\vec{u}*)/cs2 +
        !       (\vec{e}_i \cdot \vec{u}*)\vec{e}_i/cs4) \cdot \vec{F}
        do iDir = 1, QQ
          ucx = dot_product( layout%fStencil%cxDirRK(:, iDir), &
            &                velocity )
          uMinusCx = layout%fStencil%cxDirRK(:, iDir) - velocity
  
          forceTerm = dot_product( uMinusCx * cs2inv               &
            &       + ucx * layout%fStencil%cxDirRK(:,iDir) &
            &       * cs4inv, G )
  
          ! position in state array
          statePos = ?SAVE?(iDir, 1, iElem, QQ, nScalars, nPdfSize, neigh)
          ! update outstate
          outState(statePos) = outState(statePos)                         &
            & + omega_fac * layout%weight( iDir ) * forceTerm
  
        end do
  
      end do !iElem

    end associate

  end subroutine applyGlobalForce
! ****************************************************************************** !

! ****************************************************************************** !
  !> Update state with global force variable "force_phy".
  !!
  !! Force term used here is from:
  !! "Discrete lattice effects on the forcing term in the lattice Boltzmann
  !!  method", Zhaoli Guo, Chugung Zheng and Baochang Shi.
  !! In the paper, use force term is referred as Method 2 as:
  !! \[ F_i = w_i( (\vec{e}_i-\vec{u}*)/cs2 +
  !!       (\vec{e}_i \cdot \vec{u}*)\vec{e}_i/cs4) \cdot \vec{F} \]
  !! Force must be defined as body force per unit volume
  !! KM: If this force formula is used then velocity needs to be
  !! computed as u = \sum c_i f_i + \vec{F}/2
  !!
  !! Simuilar to derive routine but it updates the state whereas derive
  !! is used for tracking.
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[proc_applyGlobalForce]].
  subroutine applyGlobalForce_HRRAlike( force_phy, inState, outState, neigh,   &
    &                          auxField, nPdfSize, iLevel, varSys, phyConvFac, &
    &                          derVarPos                                       )
    ! -------------------------------------------------------------------- !
    !> Global Force to be applied on state
    real(kind=rk), intent(in) :: force_phy(3)

    !> input  pdf vector
    real(kind=rk), intent(in) :: inState(:)

    !> output pdf vector
    real(kind=rk), intent(inout) :: outState(:)

    !> connectivity Array corresponding to state vector
    integer,intent(in) :: neigh(:)

    !> auxField array
    real(kind=rk), intent(in) :: auxField(:)

    !> number of elements in state Array
    integer, intent(in) :: nPdfSize

    !> current level
    integer, intent(in) :: iLevel

    !> variable system
    type(tem_varSys_type), intent(in) :: varSys

    !> Physics conversion factor for current level
    type(mus_convertFac_type), intent(in) :: phyConvFac

    !> position of derived quantities in varsys
    type(mus_derVarPos_type), intent(in) :: derVarPos
    ! -------------------------------------------------------------------- !
    type(mus_varSys_data_type), pointer :: fPtr
    real(kind=rk) :: G(3), velocity(3), ucx, uMinusCX(3), forceTerm
    integer :: iElem, iDir, QQ, nScalars, statePos, elemOff
    integer :: vel_pos(3)
    real(kind=rk) :: omega_fac
    ! ---------------------------------------------------------------------------
    ! convert c pointer to solver type fortran pointer
    call c_f_pointer( varSys%method%val( derVarPos%pdf )%method_data, fPtr )

    ! convert physical to lattice
    ! For incompressible model: this forceField should be divided by rho0.
    ! Since rho0 =1, this term is also valid for incompressible model
    G = force_phy / phyConvFac%body_force
    
    associate (layout => fPtr%solverData%scheme%layout,                      &
      &        nSolve => fPtr%solverData%scheme%pdf(iLevel)%nElems_local     )

      ! constant parameter
      QQ = layout%fStencil%QQ
  
      nScalars = varSys%nScalars
      ! Position of velocity variable in auxField
      vel_pos = varSys%method%val(derVarPos%velocity)%auxField_varPos(1:3)
  
      do iElem = 1, nSolve
  
        ! element offset
        elemoff = (iElem-1)*varSys%nAuxScalars
        ! obtain velocity from auxField
        velocity(1) = auxField(elemOff + vel_pos(1))
        velocity(2) = auxField(elemOff + vel_pos(2))
        velocity(3) = auxField(elemOff + vel_pos(3))
  
        ! get the correct omega value
        ! HRR alike collision schemes already take into account 
        ! the part (div1_2  - omega * div1_2) with f_neq
        omega_fac = div1_2 !1.0_rk - omega * div1_2
  
        ! force term:
        ! F_i = w_i( (\vec{e}_i-\vec{u}*)/cs2 +
        !       (\vec{e}_i \cdot \vec{u}*)\vec{e}_i/cs4) \cdot \vec{F}
        do iDir = 1, QQ
          ucx = dot_product( layout%fStencil%cxDirRK(:, iDir), &
            &                velocity )
          uMinusCx = layout%fStencil%cxDirRK(:, iDir) - velocity
  
          forceTerm = dot_product( uMinusCx * cs2inv               &
            &       + ucx * layout%fStencil%cxDirRK(:,iDir) &
            &       * cs4inv, G )
  
          ! position in state array
          statePos = ?SAVE?(iDir, 1, iElem, QQ, nScalars, nPdfSize, neigh)
          ! update outstate
          outState(statePos) = outState(statePos)                         &
            & + omega_fac * layout%weight( iDir ) * forceTerm
  
        end do
  
      end do !iElem

    end associate

  end subroutine applyGlobalForce_HRRAlike
! ****************************************************************************** !

! ****************************************************************************** !
  !> Update state with source variable "force" for MRT collision model.
  !!
  !! Force term used here is from:
  !! Chai, Z., & Zhao, T. (2012). Effect of the forcing term in the
  !! multiple-relaxation-time lattice Boltzmann equation on the shear stress
  !! or the strain rate tensor. Physical Review E, 86(1), 1–11.
  !! Force term for MRT is
  !! \[ \bar{F} = M^-1 (I-0.5 S) M F' \] and
  !! \[ F'_i = w_i( (\vec{e}_i-\vec{u}*)/cs2 +
  !!       (\vec{e}_i \cdot \vec{u}*)\vec{e}_i/cs4) \cdot \vec{F} \]
  !!
  !! \vec{F} is the force that must be defined as body force per unit volume
  !!
  !! Simuilar to derive routine but it updates the state whereas derive
  !! is used for tracking.
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[proc_applyGlobalForce]].
  subroutine applyGlobalForce_MRT( force_phy, inState, outState, neigh,        &
    &                          auxField, nPdfSize, iLevel, varSys, phyConvFac, &
    &                          derVarPos                                       )
    ! -------------------------------------------------------------------- !
    !> Global Force to be applied on state
    real(kind=rk), intent(in) :: force_phy(3)

    !> input  pdf vector
    real(kind=rk), intent(in) :: inState(:)

    !> output pdf vector
    real(kind=rk), intent(inout) :: outState(:)

    !> connectivity Array corresponding to state vector
    integer,intent(in) :: neigh(:)

    !> auxField array
    real(kind=rk), intent(in) :: auxField(:)

    !> number of elements in state Array
    integer, intent(in) :: nPdfSize

    !> current level
    integer, intent(in) :: iLevel

    !> variable system
    type(tem_varSys_type), intent(in) :: varSys

    !> Physics conversion factor for current level
    type(mus_convertFac_type), intent(in) :: phyConvFac

    !> position of derived quantities in varsys
    type(mus_derVarPos_type), intent(in) :: derVarPos
    ! -------------------------------------------------------------------- !
    type(mus_varSys_data_type), pointer :: fPtr
    real(kind=rk) :: G(3), velocity(3), ucx, uMinusCX(3)
    integer :: iElem, iDir, QQ, nScalars, statePos, elemOff
    integer :: vel_pos(3)
    real(kind=rk) :: omegaKine, discForce
    real(kind=rk) :: forceTerm(27), momForce(27), s_mrt(27)
    real(kind=rk) :: mInvXOmega(27,27)
    ! ---------------------------------------------------------------------------
    ! convert c pointer to solver type fortran pointer
    call c_f_pointer( varSys%method%val( derVarPos%pdf )%method_data, fPtr )

    ! convert physical to lattice
    ! For incompressible model: this forceField should be divided by rho0.
    ! Since rho0 =1, this term is also valid for incompressible model
    G = force_phy / phyConvFac%body_force
    
    associate (layout => fPtr%solverData%scheme%layout,                      &
      &        omega_vec => fPtr%solverData%scheme%field(1)%fieldProp%fluid  &
      &                              %viscKine%omLvl(iLevel)%val,            &
      &        nSolve => fPtr%solverData%scheme%pdf(iLevel)%nElems_local,    &
      &        fluid =>  fPtr%solverData%scheme%field(1)%fieldProp%fluid     )

      ! constant parameter
      QQ = layout%fStencil%QQ
  
      nScalars = varSys%nScalars
      ! Position of velocity variable in auxField
      vel_pos = varSys%method%val(derVarPos%velocity)%auxField_varPos(1:3)
  
  
      do iElem = 1, nSolve
        ! element offset
        elemoff = (iElem-1)*varSys%nAuxScalars
        ! obtain velocity from auxField
        velocity(1) = auxField(elemOff + vel_pos(1))
        velocity(2) = auxField(elemOff + vel_pos(2))
        velocity(3) = auxField(elemOff + vel_pos(3))
  
        ! get the correct omega value
        omegaKine = omega_vec(iElem)
        ! MRT omegas
        ! overwrite omegaKine term in the element loop
        s_mrt(1:QQ) = fluid%mrtPtr( omegaKine=omegaKine,                  &
          &                         omegaBulk=fluid%omegaBulkLvl(iLevel), &
          &                         QQ=QQ                                 )
  
        ! M^-1 * (I-0.5 S)
        s_mrt(1:QQ) = 1.0_rk - 0.5_rk * s_mrt(1:QQ)
        do iDir = 1, QQ
          mInvXOmega(1:QQ,iDir) = layout%moment%toPDF%A(1:QQ,iDir) &
            &                   * s_mrt(iDir)
        end do
  
        ! force term:
        ! F_i = w_i( (\vec{e}_i-\vec{u}*)/cs2 +
        !       (\vec{e}_i \cdot \vec{u}*)\vec{e}_i/cs4) \cdot \vec{F}
        do iDir = 1, QQ
          ucx = dot_product( layout%fStencil%cxDirRK(:, iDir), &
            &                velocity )
          uMinusCx = layout%fStencil%cxDirRK(:, iDir) - velocity
  
          forceTerm(iDir) = layout%weight(iDir)             &
            &       * dot_product( uMinusCx * cs2inv               &
            &       + ucx * layout%fStencil%cxDirRK(:,iDir) &
            &       * cs4inv, G )
        end do
  
        ! Force moments: M * F
        !do iDir = 1, QQ
        !  momForce(iDir) = sum(layout%moment%toMoments%A(iDir,:) * forceTerm)
        !end do
        momForce(1:QQ) = matmul( layout%moment%toMoments%A(1:QQ, 1:QQ), &
          &                forceTerm(1:QQ) )
  
        !discForce = matmul( omegaTerm, forceTerm )
        do iDir = 1, QQ
          ! discrete force
          ! \bar{F} =  M^-1 (I-S/2) M F
          discForce = dot_product(mInvXOmega(iDir,1:QQ),  momForce(1:QQ))
          ! position in state array
          statePos = ?SAVE?(iDir, 1, iElem, QQ, nScalars, nPdfSize, neigh)
          ! update outstate
          outState(statePos) = outState(statePos) + discForce
        end do
  
      end do !iElem

    end associate

  end subroutine applyGlobalForce_MRT
! ****************************************************************************** !

! ****************************************************************************** !
  !> Update state with source variable "force" for MRT collision model.
  !!
  !! Force term used here is from:
  !! Chai, Z., & Zhao, T. (2012). Effect of the forcing term in the
  !! multiple-relaxation-time lattice Boltzmann equation on the shear stress
  !! or the strain rate tensor. Physical Review E, 86(1), 1–11.
  !! Force term for MRT is
  !! \[ \bar{F} = M^-1 (I-0.5 S) M F' \] and
  !! \[ F'_i = w_i( (\vec{e}_i-\vec{u}*)/cs2 +
  !!       (\vec{e}_i \cdot \vec{u}*)\vec{e}_i/cs4) \cdot \vec{F} \]
  !!
  !! \vec{F} is the force that must be defined as body force per unit volume
  !!
  !! Simuilar to derive routine but it updates the state whereas derive
  !! is used for tracking.
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[proc_applyGlobalForce]].
  subroutine applyGlobalForce_MRT_d3q27( force_phy, inState, outState, neigh,  &
    &                          auxField, nPdfSize, iLevel, varSys, phyConvFac, &
    &                          derVarPos                                       )
    ! -------------------------------------------------------------------- !
    !> Global Force to be applied on state
    real(kind=rk), intent(in) :: force_phy(3)

    !> input  pdf vector
    real(kind=rk), intent(in) :: inState(:)

    !> output pdf vector
    real(kind=rk), intent(inout) :: outState(:)

    !> connectivity Array corresponding to state vector
    integer,intent(in) :: neigh(:)

    !> auxField array
    real(kind=rk), intent(in) :: auxField(:)

    !> number of elements in state Array
    integer, intent(in) :: nPdfSize

    !> current level
    integer, intent(in) :: iLevel

    !> variable system
    type(tem_varSys_type), intent(in) :: varSys

    !> Physics conversion factor for current level
    type(mus_convertFac_type), intent(in) :: phyConvFac

    !> position of derived quantities in varsys
    type(mus_derVarPos_type), intent(in) :: derVarPos
    ! -------------------------------------------------------------------- !
    type(mus_varSys_data_type), pointer :: fPtr
    real(kind=rk) :: F(3), velocity(3)
    integer :: iElem, iDir, QQ, nScalars, statePos, elemOff
    integer :: vel_pos(3)
    real(kind=rk) :: discForce
    real(kind=rk) :: momForce(27), s_mrt(27)
    real(kind=rk) :: mInvXOmega(27,27)
    ! ---------------------------------------------------------------------------
    ! convert c pointer to solver type fortran pointer
    call c_f_pointer( varSys%method%val( derVarPos%pdf )%method_data, fPtr )

    ! convert physical to lattice
    ! For incompressible model: this forceField should be divided by rho0.
    ! Since rho0 =1, this term is also valid for incompressible model
    F = force_phy / phyConvFac%body_force
    
    associate (layout => fPtr%solverData%scheme%layout,                      &
      &        omega_vec => fPtr%solverData%scheme%field(1)%fieldProp%fluid  &
      &                              %viscKine%omLvl(iLevel)%val,            &
      &        nSolve => fPtr%solverData%scheme%pdf(iLevel)%nElems_local,    &
      &        fluid =>  fPtr%solverData%scheme%field(1)%fieldProp%fluid     )

      ! constant parameter
      QQ = layout%fStencil%QQ
  
      nScalars = varSys%nScalars
      ! Position of velocity variable in auxField
      vel_pos = varSys%method%val(derVarPos%velocity)%auxField_varPos(1:3)
  
      ! MRT omegas
      s_mrt = fluid%mrtPtr( omegaKine=1._rk,                      &
        &                   omegaBulk=fluid%omegaBulkLvl(iLevel), &
        &                   QQ=QQ                                 )
      ! M^-1 * (I-0.5 S)
      s_mrt(2:4) = 1.0_rk - 0.5_rk * s_mrt(2:4)
      s_mrt(10) = 1.0_rk - 0.5_rk * s_mrt(10)
  
      do iElem = 1, nSolve
        ! element offset
        elemoff = (iElem-1)*varSys%nAuxScalars
        ! obtain velocity from auxField
        velocity(1) = auxField(elemOff + vel_pos(1))
        velocity(2) = auxField(elemOff + vel_pos(2))
        velocity(3) = auxField(elemOff + vel_pos(3))
  
        ! MRT omegas
        ! overwrite omegaKine term in the element loop
        ! get the correct omega value
        s_mrt(5:9) = omega_vec(iElem)
  
        ! M^-1 * (I-0.5 S)
        s_mrt(5:9) = 1.0_rk - 0.5_rk * s_mrt(5:9)
        do iDir = 2, 10
          mInvXOmega(1:QQ,iDir) = layout%moment%toPDF%A(1:QQ,iDir) &
            &                   * s_mrt(iDir)
        end do
  
        ! force term:
        ! F_i = w_i( (\vec{e}_i-\vec{u}*)/cs2 +
        !       (\vec{e}_i \cdot \vec{u}*)\vec{e}_i/cs4) \cdot \vec{F}
        ! Force moments: M * F
        momForce(1:QQ) = 0._rk
        momForce(2:4) = F(1:3)
        momForce(5) = F(1) * velocity(2) + F(2) * velocity(1)
        momForce(6) = F(2) * velocity(3) + F(3) * velocity(2)
        momForce(7) = F(1) * velocity(3) + F(3) * velocity(1)
        momForce(8) = -2._rk * ( F(2) * velocity(2) - 2._rk * F(1) * velocity(1) &
          &                      + F(3) * velocity(3) )
        momForce(9) = 2._rk * ( F(2) * velocity(2) - F(3) * velocity(3) )
        momForce(10) = 2._rk * ( F(1) * velocity(1) + F(2) * velocity(2) &
          &                      + F(3) * velocity(3) )
  
        do iDir = 1, QQ
          ! discrete force
          ! \bar{F} =  M^-1 (I-S/2) M F
          discForce = dot_product(mInvXOmega(iDir,2:10),  momForce(2:10))
          ! position in state array
          statePos = ?SAVE?(iDir, 1, iElem, QQ, nScalars, nPdfSize, neigh)
          ! update outstate
          outState(statePos) = outState(statePos) + discForce
        end do
  
      end do !iElem

    end associate

  end subroutine applyGlobalForce_MRT_d3q27
! ****************************************************************************** !


! ****************************************************************************** !
  !> Update state with source variable "force" for d3q19 MRT collision model.
  !!
  !! Force term used here is from:
  !! Chai, Z., & Zhao, T. (2012). Effect of the forcing term in the
  !! multiple-relaxation-time lattice Boltzmann equation on the shear stress
  !! or the strain rate tensor. Physical Review E, 86(1), 1–11.
  !! Force term for MRT is
  !! \[ \bar{F} = M^-1 (I-0.5 S) M F' \] and
  !! \[ F'_i = w_i( (\vec{e}_i-\vec{u}*)/cs2 +
  !!       (\vec{e}_i \cdot \vec{u}*)\vec{e}_i/cs4) \cdot \vec{F} \]
  !!
  !! Force must be defined as body force per unit volume
  !!
  !! Simuilar to derive routine but it updates the state whereas derive
  !! is used for tracking.
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[proc_applyGlobalForce]].
  subroutine applyGlobalForce_MRT_d3q19( force_phy, inState, outState, neigh,  &
    &                          auxField, nPdfSize, iLevel, varSys, phyConvFac, &
    &                          derVarPos                                       )
    ! -------------------------------------------------------------------- !
    !> Global Force to be applied on state
    real(kind=rk), intent(in) :: force_phy(3)

    !> input  pdf vector
    real(kind=rk), intent(in) :: inState(:)

    !> output pdf vector
    real(kind=rk), intent(inout) :: outState(:)

    !> connectivity Array corresponding to state vector
    integer,intent(in) :: neigh(:)

    !> auxField array
    real(kind=rk), intent(in) :: auxField(:)

    !> number of elements in state Array
    integer, intent(in) :: nPdfSize

    !> current level
    integer, intent(in) :: iLevel

    !> variable system
    type(tem_varSys_type), intent(in) :: varSys

    !> Physics conversion factor for current level
    type(mus_convertFac_type), intent(in) :: phyConvFac

    !> position of derived quantities in varsys
    type(mus_derVarPos_type), intent(in) :: derVarPos
    ! -------------------------------------------------------------------- !
    type(mus_varSys_data_type), pointer :: fPtr
    real(kind=rk) :: F(3), velocity(3)
    integer :: iElem, iDir, QQ, nScalars, statePos, elemOff
    integer :: vel_pos(3)
    real(kind=rk) :: omegaKine, discForce
    real(kind=rk) :: momForce(19), s_mrt(19)
    real(kind=rk) :: mInvXOmega(19,19)
    ! ---------------------------------------------------------------------------
    ! convert c pointer to solver type fortran pointer
    call c_f_pointer( varSys%method%val( derVarPos%pdf )%method_data, fPtr )

    ! convert physical to lattice
    ! For incompressible model: this forceField should be divided by rho0.
    ! Since rho0 =1, this term is also valid for incompressible model
    F = force_phy / phyConvFac%body_force
    
    associate (layout => fPtr%solverData%scheme%layout,                      &
      &        omega_vec => fPtr%solverData%scheme%field(1)%fieldProp%fluid  &
      &                              %viscKine%omLvl(iLevel)%val,            &
      &        nSolve => fPtr%solverData%scheme%pdf(iLevel)%nElems_local,    &
      &        fluid =>  fPtr%solverData%scheme%field(1)%fieldProp%fluid     )

      ! constant parameter
      QQ = 19
  
      nScalars = varSys%nScalars
      ! Position of velocity variable in auxField
      vel_pos = varSys%method%val(derVarPos%velocity)%auxField_varPos(1:3)
  
      ! MRT omegas
      ! overwrite omegaKine term in the element loop
      ! KM: For incompressible model: omegaBulk is unused in mrtPtr
      s_mrt = fluid%mrtPtr( omegaKine=1._rk,                  &
        &                   omegaBulk=fluid%omegaBulkLvl(iLevel), &
        &                   QQ=QQ                                 )
  
      ! F = M^-1 (I-0.5 S) M F
      ! (I-0.5 S) - omega for force term
      s_mrt = 1.0_rk - 0.5_rk * s_mrt
  
      do iElem = 1, nSolve
        ! element offset
        elemoff = (iElem-1)*varSys%nAuxScalars
        ! obtain velocity from auxField
        velocity(1) = auxField(elemOff + vel_pos(1))
        velocity(2) = auxField(elemOff + vel_pos(2))
        velocity(3) = auxField(elemOff + vel_pos(3))
  
        ! get the correct omega value
        omegaKine = omega_vec(iElem)
        ! MRT omegas
        ! overwrite omegaKine term in the element loop
        s_mrt(10) = 1.0_rk - 0.5_rk * omegaKine
        s_mrt(12) = s_mrt(10)
        s_mrt(14) = s_mrt(10)
        s_mrt(15) = s_mrt(10)
        s_mrt(16) = s_mrt(10)
  
        ! M^-1 (1-0.5 S)
        do iDir = 1, QQ
          mInvXOmega(:,iDir) = layout%moment%toPDF%A(:,iDir) * s_mrt(iDir)
        end do
  
        ! force term:
        ! F_i = w_i( (\vec{e}_i-\vec{u}*)/cs2 +
        !       (\vec{e}_i \cdot \vec{u}*)\vec{e}_i/cs4) \cdot \vec{F}
        ! Force moments: M * F
        momForce(1:QQ) = 0._rk
        momForce(2) = 2._rk * ( F(1) * velocity(1) + F(2) * velocity(2) &
          &                     + F(3) * velocity(3) )
        momForce(4) = F(1)
        momForce(6) = F(2)
        momForce(8) = F(3)
        momForce(10) = -2._rk * ( F(2) * velocity(2) - 2._rk * F(1) * velocity(1) &
          &                       + F(3) * velocity(3) )
        momForce(12) = 2._rk * ( F(2) * velocity(2) - F(3) * velocity(3) )
        momForce(14) = F(1) * velocity(2) + F(2) * velocity(1)
        momForce(15) = F(2) * velocity(3) + F(3) * velocity(2)
        momForce(16) = F(1) * velocity(3) + F(3) * velocity(1)
  
        do iDir = 1, QQ
          ! discrete force
          ! \bar{F} =  M^-1 (I-S/2) M F
          discForce = sum(mInvXOmega(iDir,2:16) * momForce(2:16))
          ! position in state array
          statePos = ?SAVE?(iDir, 1, iElem, QQ, nScalars, nPdfSize, neigh)
          ! update outstate
          outState(statePos) = outState(statePos) + discForce
        end do
  
      end do !iElem

    end associate

  end subroutine applyGlobalForce_MRT_d3q19
! ****************************************************************************** !


! ****************************************************************************** !
  !> Update state with source variable "force" for d3q19 MRT collision model.
  !!
  !! Force term used here is from:
  !! Chai, Z., & Zhao, T. (2012). Effect of the forcing term in the
  !! multiple-relaxation-time lattice Boltzmann equation on the shear stress
  !! or the strain rate tensor. Physical Review E, 86(1), 1–11.
  !! Force term for MRT is
  !! \[ \bar{F} = M^-1 (I-0.5 S) M F' \] and
  !! \[ F'_i = w_i( (\vec{e}_i-\vec{u}*)/cs2 +
  !!       (\vec{e}_i \cdot \vec{u}*)\vec{e}_i/cs4) \cdot \vec{F} \]
  !!
  !! Force must be defined as body force per unit volume
  !!
  !! Simuilar to derive routine but it updates the state whereas derive
  !! is used for tracking.
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[proc_applyGlobalForce]].
  subroutine applyGlobalForce_MRT_d2q9( force_phy, inState, outState, neigh,   &
    &                          auxField, nPdfSize, iLevel, varSys, phyConvFac, &
    &                          derVarPos                                       )
    ! -------------------------------------------------------------------- !
    !> Global Force to be applied on state
    real(kind=rk), intent(in) :: force_phy(3)

    !> input  pdf vector
    real(kind=rk), intent(in) :: inState(:)

    !> output pdf vector
    real(kind=rk), intent(inout) :: outState(:)

    !> connectivity Array corresponding to state vector
    integer,intent(in) :: neigh(:)

    !> auxField array
    real(kind=rk), intent(in) :: auxField(:)

    !> number of elements in state Array
    integer, intent(in) :: nPdfSize

    !> current level
    integer, intent(in) :: iLevel

    !> variable system
    type(tem_varSys_type), intent(in) :: varSys

    !> Physics conversion factor for current level
    type(mus_convertFac_type), intent(in) :: phyConvFac

    !> position of derived quantities in varsys
    type(mus_derVarPos_type), intent(in) :: derVarPos
    ! -------------------------------------------------------------------- !
    type(mus_varSys_data_type), pointer :: fPtr
    real(kind=rk) :: F(3), velocity(3)
    integer :: iElem, iDir, QQ, nScalars, statePos, elemOff
    integer :: vel_pos(3)
    real(kind=rk) :: omegaKine, discForce
    real(kind=rk) :: momForce(9), s_mrt(9)
    real(kind=rk) :: mInvXOmega(9,9)
    ! ---------------------------------------------------------------------------
    ! convert c pointer to solver type fortran pointer
    call c_f_pointer( varSys%method%val( derVarPos%pdf )%method_data, fPtr )

    ! convert physical to lattice
    ! For incompressible model: this forceField should be divided by rho0.
    ! Since rho0 =1, this term is also valid for incompressible model
    F = force_phy / phyConvFac%body_force
    
    associate (layout => fPtr%solverData%scheme%layout,                      &
      &        omega_vec => fPtr%solverData%scheme%field(1)%fieldProp%fluid  &
      &                              %viscKine%omLvl(iLevel)%val,            &
      &        nSolve => fPtr%solverData%scheme%pdf(iLevel)%nElems_local,    &
      &        fluid =>  fPtr%solverData%scheme%field(1)%fieldProp%fluid     )

      ! constant parameter
      QQ = 9
  
      nScalars = varSys%nScalars
      ! Position of velocity variable in auxField
      vel_pos = varSys%method%val(derVarPos%velocity)%auxField_varPos(1:3)
  
      ! MRT omegas
      ! overwrite omegaKine term in the element loop
      ! KM: For incompressible model: omegaBulk is unused in mrtPtr
      s_mrt = fluid%mrtPtr( omegaKine=1._rk,                      &
        &                   omegaBulk=fluid%omegaBulkLvl(iLevel), &
        &                   QQ=QQ                                 )
  
      ! F = M^-1 (I-0.5 S) M F
      ! (I-0.5 S) - omega for force term
      s_mrt = 1.0_rk - 0.5_rk * s_mrt
  
      do iElem = 1, nSolve
        ! element offset
        elemoff = (iElem-1)*varSys%nAuxScalars
        ! obtain velocity from auxField
        velocity(1) = auxField(elemOff + vel_pos(1))
        velocity(2) = auxField(elemOff + vel_pos(2))
  
        ! get the correct omega value
        omegaKine = omega_vec(iElem)
        ! MRT omegas
        ! overwrite omegaKine term in the element loop
        s_mrt(8:9) = 1.0_rk - 0.5_rk * omegaKine
  
        ! M^-1 (1-0.5 S)
        do iDir = 1, QQ
          mInvXOmega(:,iDir) = layout%moment%toPDF%A(:,iDir) * s_mrt(iDir)
        end do
  
        ! force term:
        ! F_i = w_i( (\vec{e}_i-\vec{u}*)/cs2 +
        !       (\vec{e}_i \cdot \vec{u}*)\vec{e}_i/cs4) \cdot \vec{F}
        ! Force moments: M * F
        momForce(1) = 0._rk
        momForce(2) = 6._rk * ( F(1) * velocity(1) + F(2) * velocity(2) )
        momForce(3) = -momForce(2)
        momForce(4) = F(1)
        momForce(5) = -F(1)
        momForce(6) = F(2)
        momForce(7) = -F(2)
        momForce(8) = 2._rk * ( F(1) * velocity(1) - F(2) * velocity(2) )
        momForce(9) = F(1) * velocity(2) + F(2) * velocity(1)
  
        do iDir = 1, QQ
          ! discrete force
          ! \bar{F} =  M^-1 (I-S/2) M F
          discForce = sum(mInvXOmega(iDir,2:QQ) * momForce(2:QQ))
          ! position in state array
          statePos = ?SAVE?(iDir, 1, iElem, QQ, nScalars, nPdfSize, neigh)
          ! update outstate
          outState(statePos) = outState(statePos) + discForce
        end do
  
      end do !iElem

    end associate

  end subroutine applyGlobalForce_MRT_d2q9
! ****************************************************************************** !

! ************************************************************************** !
  !> Update state with source variable "force_1stOrd"
  !!
  !! Force term used here is from:
  !! "A D3Q27 multiple-relaxation-time lattice Boltzmann method for
  !! turbulent flows", K. Suga, Y. Kuwata, K. Takashima, R. Chikasue
  !!
  !! \[ F_i = w_i/c_s^2 ( \vec{e}_i \cdot \vec{F} ) \]
  !! Force must be defined as body force per unit volume
  !! This force term can be applied for both compressible and incompressible
  !! LBM models
  !!
  !! Similar to derive routine but it updates the state whereas derive
  !! is used for tracking
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[proc_applyGlobalForce]].
  subroutine applyGlobalForce_1stOrd( force_phy, inState, outState, neigh,     &
    &                          auxField, nPdfSize, iLevel, varSys, phyConvFac, &
    &                          derVarPos                                       )
    ! -------------------------------------------------------------------- !
    !> Global Force to be applied on state
    real(kind=rk), intent(in) :: force_phy(3)

    !> input  pdf vector
    real(kind=rk), intent(in) :: inState(:)

    !> output pdf vector
    real(kind=rk), intent(inout) :: outState(:)

    !> connectivity Array corresponding to state vector
    integer,intent(in) :: neigh(:)

    !> auxField array
    real(kind=rk), intent(in) :: auxField(:)

    !> number of elements in state Array
    integer, intent(in) :: nPdfSize

    !> current level
    integer, intent(in) :: iLevel

    !> variable system
    type(tem_varSys_type), intent(in) :: varSys

    !> Physics conversion factor for current level
    type(mus_convertFac_type), intent(in) :: phyConvFac

    !> position of derived quantities in varsys
    type(mus_derVarPos_type), intent(in) :: derVarPos
    ! -------------------------------------------------------------------- !
    type(mus_varSys_data_type), pointer :: fPtr
    real(kind=rk) :: F(3)
    real(kind=rk) :: forceTerm
    integer :: iElem, iDir, QQ, nScalars
    ! ---------------------------------------------------------------------- !
    ! convert c pointer to solver type fortran pointer
    call c_f_pointer( varSys%method%val( derVarPos%pdf )%method_data, fPtr )

    ! convert physical to lattice
    ! For incompressible model: this forceField should be divided by rho0.
    ! Since rho0 =1, this term is also valid for incompressible model
    F = force_phy / phyConvFac%body_force
    
    associate (layout => fPtr%solverData%scheme%layout,                      &
      &        nSolve => fPtr%solverData%scheme%pdf(iLevel)%nElems_local     )

      ! constant parameter
      QQ = layout%fStencil%QQ
      nScalars = varSys%nScalars
  
      do iElem = 1, nSolve
        ! force term:
        ! F_i = w_i/cs2 (\vec{e}_i \cdot \vec{F}
        do iDir = 1, QQ
  
          forceTerm = dot_product( layout%fStencil    &
            &                      %cxDirRK(:,iDir), F )
  
          outState( ?SAVE?(iDir,1,iElem,QQ,nScalars,nPdfSize,neigh) ) &
            & = outState(                                             &
            & ?SAVE?(iDir,1,iElem,QQ,nScalars,nPdfSize,neigh) )       &
            & + layout%weight( iDir ) * cs2inv * forceTerm
  
        end do
  
      end do !iElem

    end associate

  end subroutine applyGlobalForce_1stOrd
! ************************************************************************** !
  
  end module mus_GlobalForce_module
  ! ****************************************************************************** !
  
