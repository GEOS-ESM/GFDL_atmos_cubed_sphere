module fv_arrays_mod
  
  use fv_mp_mod, only: halo_size => ng
  
  implicit none
  
  public fv_grid_type, fv_grid_bounds_type
  
  private

  integer, parameter :: R_GRID = 8

  type fv_grid_type
     logical, pointer :: nested
     integer, pointer :: grid_type
     logical:: sw_corner, se_corner, ne_corner, nw_corner
     real, allocatable, dimension(:,:) :: dxa, dya
     real, allocatable, dimension(:,:) :: area
     ! rarea, del6_u/v are used in subroutine deln_flux which is called
     ! if some of the optional arguments mfx, mfy, damp_c and mass are present
     real, allocatable, dimension(:,:) :: rarea ! NOT INITIALIZED AT THE MOMENT
     real, allocatable :: del6_u(:,:), del6_v(:,:) ! NOT INITIALIZED AT THE MOMENT
     ! da_min is used in fv_tp_2d if damp_c and mass are present
     real(kind=R_GRID) :: da_min ! NOT INITIALIZED AT THE MOMENT
  end type fv_grid_type

  interface fv_grid_type
     procedure fv_grid_type_initialize ! constructor
  end interface fv_grid_type

  type fv_grid_bounds_type
     integer :: is,  ie,  js,  je
     integer :: isd, ied, jsd, jed
  end type fv_grid_bounds_type

  interface fv_grid_bounds_type
     procedure fv_grid_bounds_type_initialize ! constructor
  end interface fv_grid_bounds_type
  
contains

  function fv_grid_bounds_type_initialize(n) result(bd)

    ! Arguments
    integer, intent(in) :: n
    type(fv_grid_bounds_type) :: bd ! output

    ! Start
    bd%is = 1
    bd%ie = n
    bd%isd = bd%is - halo_size
    bd%ied = bd%ie + halo_size

    bd%js = 1
    bd%je = n
    bd%jsd = bd%js - halo_size
    bd%jed = bd%je + halo_size

  end function fv_grid_bounds_type_initialize

  function fv_grid_type_initialize(bd, npx, npy, nested, grid_type) result(gridstruct)

    ! Arguments
    type(fv_grid_bounds_type), intent(in) :: bd
    integer, intent(in) :: npx, npy
    logical, intent(in), target :: nested
    integer, intent(in), target :: grid_type
    type(fv_grid_type) :: gridstruct ! output

    ! NOTE: We are NOT initializing rarea, del6_u/v, da_min
    ! These variables are used if some of the optional arguments (mfx, mfy,
    ! damp_c, mass) to fv_tp_2d are present. At the moment, for the purpose of
    ! top_core driver, we are not passing the optional arguments

    ! Start
    gridstruct%grid_type => grid_type
    gridstruct%nested => nested
    
    call initialize_corners_(bd, npx, npy, gridstruct)
    allocate(gridstruct%dxa(bd%isd:bd%ied, bd%jsd:bd%jed), source = 1.0)
    allocate(gridstruct%dya(bd%isd:bd%ied, bd%jsd:bd%jed), source = 1.0)
    allocate(gridstruct%area(bd%isd:bd%ied, bd%jsd:bd%jed), source = 1.0)

  end function fv_grid_type_initialize

  subroutine initialize_corners_(bd, npx, npy, gridstruct)

    ! Arguments
    type(fv_grid_bounds_type), intent(in) :: bd
    integer, intent(in) :: npx, npy
    type(fv_grid_type), intent(inout) :: gridstruct

    ! Start
    if ( (gridstruct%grid_type < 3) .and. (.not. gridstruct%nested) ) then
       if ( bd%is==1 .and. bd%js==1 ) gridstruct%sw_corner = .true.
       if ( (bd%ie+1)==npx .and. bd%js==1 ) gridstruct%se_corner = .true.
       if ( (bd%ie+1)==npx .and. (bd%je+1)==npy ) gridstruct%ne_corner = .true.
       if (     bd%is==1   .and. (bd%je+1)==npy ) gridstruct%nw_corner = .true.
    endif

  end subroutine initialize_corners_

end module fv_arrays_mod
