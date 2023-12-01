! See copyright notice in the COPYRIGHT file.
!> Module to create Treelm environment for unit testcases.
module mus_utestEnv_module
  use, intrinsic :: iso_c_binding, only: C_NEW_LINE, c_loc
  use env_module,           only: rk, labelLen
  use tem_param_module,     only: rho0, cs2
  use tem_bc_prop_module,   only: tem_bc_prop_type, load_tem_bc_prop,         &
    &                             tem_empty_bc_prop
  use tem_property_module,  only: prp_hasbnd
  use flu_binding,          only: flu_state
  use treelmesh_module,     only: treelmesh_type
  use aotus_module,         only: open_config_chunk
  use treelmesh_module,     only: load_tem
  use tem_grow_array_module, only: grw_labelarray_type, init
  use tem_general_module,   only: tem_general_type, tem_load_general, tem_start
  use tem_logging_module,   only: tem_logging_load_primary
  use tem_varSys_module,    only: tem_varSys_type, tem_varSys_init, &
    &                             tem_varSys_append_stateVar,       &
    &                             tem_varSys_proc_point,            &
    &                             tem_varSys_proc_element,          &
    &                             tem_varSys_proc_setParams,        &
    &                             tem_varSys_proc_getParams,        &
    &                             tem_varSys_proc_setupIndices,     &
    &                             tem_varSys_proc_getValOfIndex

  use mus_scheme_header_module,  only: mus_scheme_header_type
  use mus_mrtRelaxation_module,  only: mus_assign_mrt_ptr
  use mus_fluid_module,  only: mus_fluid_type
  use mus_varSys_module, only: mus_varSys_solverData_type, &
    &                          mus_get_new_solver_ptr
  use mus_derVarPos_module, only: mus_derVarPos_type
  use mus_variable_module, only: mus_append_auxField, mus_store_derVarPos

  implicit none

  private

  public :: init_fluid
  public :: init_varSys
  public :: load_env

  character, parameter :: nl = C_NEW_LINE

  character(len=200), parameter :: cubeconf =                                  &
    &   'mesh = {' // nl                                                       &
    & //'         predefined = "cube",' // nl                                  &
    & //'         origin = {0.0, 0.0, 0.0},' // nl                             &
    & //'         length=1.0,' // nl                                           &
    & //'         refinementLevel = 1' // nl                                   &
    & //'       }' // nl

contains

  !> Create the treelm environment for unit tests.
  subroutine load_env(tree, boundary, general )
    !---------------------------------------------------------------------------
    type(treelmesh_type), intent(out) :: tree
    type(tem_bc_prop_type), intent(out) :: boundary 
    type(tem_general_type), intent(out) :: general
    !---------------------------------------------------------------------------
    type(flu_state) :: conf
    integer :: iProp
    !---------------------------------------------------------------------------

    ! Init the Treelm environment
    call tem_start('MUSUBI unit test', 'utest', general)

    print *, "Hello from loadenv" 
    ! Open the configuration file 
    call open_config_chunk(L = conf, chunk = trim(cubeconf))

    ! load and initialize logUnit
    call tem_logging_load_primary(conf = conf,              &
      &                           rank = general%proc%rank  )

    call tem_load_general( me = general, conf = conf )

    ! Load the mesh first.
    call load_tem(me = tree, conf = conf, myPart = general%proc%rank, &
      &           nParts = general%proc%comm_size, comm = general%proc%comm )

    ! Load the properties now.
    do iprop = 1, size(tree%Property)
      ! Is the current property setting the prp_hasBnd bit?
      if (tree%global%property(iprop)%bitpos == prp_hasBnd) then
        ! Found the property setting the prp_hasBnd, load the
        ! Boundaries accordingly and leave the loop.
        call load_tem_bc_prop(me = boundary,                &
          &    offset = tree%Property(iprop)%Offset,        &
          &    nElems = tree%Property(iprop)%nElems,        &
          &    basename = trim(tree%global%dirname)//'bnd', &
          &    comm = general%proc%comm,                    & 
          &    mypart = general%proc%rank )
        exit
      end if
    end do

  end subroutine load_env

  ! fill fluid properties with values of omega and rho0 for utest
  subroutine init_fluid( fluid, omega, nSolve, schemeHeader )
    type(mus_fluid_type), intent(out) :: fluid
    real(kind=rk), intent(in) :: omega
    integer, intent(in) :: nSolve
    type(mus_scheme_header_type), intent(in) :: schemeHeader

    allocate( fluid%viscKine%dataOnLvl(1) )
    call init(fluid%viscKine%dataOnLvl(1), nSolve)

    allocate( fluid%viscKine%omLvl(1) )
    ! Assume constant omega
    allocate( fluid%viscKine%omLvl(1)%val(nSolve) ) 
    allocate( fluid%omegaBulkLvl(1) )

    fluid%viscKine%omLvl(1)%val(:) = omega
    fluid%viscKine%dataOnLvl(1)%val(:) = cs2 * (1.0_rk/omega - 0.5_rk) 
    fluid%omegaBulkLvl = omega

    ! assign function pointer
    call mus_assign_mrt_ptr(fluid%mrtPtr, schemeHeader)

  end subroutine init_fluid

  ! init varSys and append pdf variable to varSys
  ! this routine is used by all compute kernel utests 
  subroutine init_varSys( varSys, sysName, QQ, solverData, schemeHeader, &
    &                     derVarPos )
    type( tem_varSys_type ) :: varSys
    character(len=labelLen), intent(in) :: sysName
    integer, intent(in) :: QQ
    type(mus_varSys_solverData_type), target :: solverData
    type(mus_scheme_header_type), intent(in) :: schemeHeader
    type(mus_derVarPos_type), allocatable, intent(out)  :: derVarPos(:)

    integer :: pdfPos
    logical :: wasAdded
    procedure(tem_varSys_proc_point), pointer :: get_point => NULL()
    procedure(tem_varSys_proc_element), pointer :: get_element => NULL()
    procedure(tem_varSys_proc_setParams), pointer :: set_params => null()
    procedure(tem_varSys_proc_getParams), pointer :: get_params => null()
    procedure(tem_varSys_proc_setupIndices), pointer :: &
      &                                      setup_indices => null()
    procedure(tem_varSys_proc_getValOfIndex), pointer :: &
      &                                       get_valOfIndex => null()
    character(len=labelLen) :: fldLabel(1)
    ! array of derive variable names depends on scheme kind
    type(grw_labelarray_type) ::  derVarName

    call tem_varSys_init( varSys, sysName )

    call tem_varSys_append_stateVar(                         &
      & me             = varSys,                             &
      & varName        = 'pdf',                              &
      & nComponents    = QQ,                                 &
      & method_data    = mus_get_new_solver_ptr(solverData), &
      & get_point      = get_point,                          &
      & get_element    = get_element,                        &
      & set_params     = set_params,                         &
      & get_params     = get_params,                         &
      & setup_indices  = setup_indices,                      &
      & get_valOfIndex = get_valOfIndex,                     &
      & pos            = pdfPos,                             &
      & wasAdded       = wasAdded                            )

    fldLabel = ''
    ! initialize derVarname list
    call init(me = derVarName, length=1)
    ! append auxField variable depending on scheme kinds
    call mus_append_auxField( varSys       = varSys,       &
      &                       solverData   = solverData,   &
      &                       schemeHeader = schemeHeader, &
      &                       nFields      = 1,            & 
      &                       fldLabel     = fldLabel,     &
      &                       derVarName   = derVarName    )


    allocate( derVarPos(1) ) ! allocate for nFields = 1
    ! store derVarPos
    call mus_store_derVarPos( derVarPos  = derVarPos,  &
      &                       derVarname = derVarname, &
      &                       varSys     = varSys,     &
      &                       nFields    = 1,          &
      &                       fldLabel   = fldLabel    )
    write(*,*) 'scheme ', derVarPos(1)%density
    write(*,*) 'dens_pos ', varSys%method%val(derVarPos(1)%density)%auxField_varPos(1)

  end subroutine init_varsys

end module mus_utestEnv_module
