module input_arrays_mod

  use fv_arrays_mod, only: fv_grid_bounds_type, fv_grid_type

  implicit none

  private

  public InputArrays_T

  type InputArrays_T
     real, allocatable :: q(:, :)
     real, allocatable :: crx(:, :), cry(:, :)
     real, allocatable :: xfx(:, :), yfx(:, :)
     real, allocatable :: ra_x(:, :), ra_y(:, :)
  end type InputArrays_T

  interface InputArrays_T
     procedure :: input_arrays_initialize_
  end interface InputArrays_T

contains

  function input_arrays_initialize_(bd, npx, npy, gridstruct) result(arrays)

    ! Arguments
    type(fv_grid_bounds_type), intent(in) :: bd
    integer, intent(in) :: npx, npy
    type(fv_grid_type), intent(in) :: gridstruct
    type(InputArrays_T) :: arrays ! output
    ! Locals
    integer :: i, j

    call initialize_q_(bd, npx, npy, arrays%q)
    allocate(arrays%crx(bd%is:bd%ie+1, bd%jsd:bd%jed), source = 0.5)
    allocate(arrays%cry(bd%isd:bd%ied, bd%js:bd%je+1), source = 0.5)
    allocate(arrays%xfx(bd%is:bd%ie+1, bd%jsd:bd%jed), source = 0.5)
    allocate(arrays%yfx(bd%isd:bd%ied, bd%js:bd%je+1), source = 0.5)
    call initialize_ra_x_(bd, gridstruct%area, arrays%xfx, arrays%ra_x)
    call initialize_ra_y_(bd, gridstruct%area, arrays%yfx, arrays%ra_y)

  end function input_arrays_initialize_

  subroutine initialize_q_(bd, npx, npy, q)

    ! Arguments
    type(fv_grid_bounds_type), intent(in) :: bd
    integer, intent(in) :: npx, npy
    real, allocatable, intent(out) :: q(:,:)
    ! Locals
    integer :: j, i

    allocate(q(bd%isd:bd%ied, bd%jsd:bd%jed))
    do j = bd%jsd, bd%jed
       do i = bd%isd, bd%ied
          q(i,j) = sin(3.1415927 * float(i*j) / float( (npx-1) * (npy-1) ) )
       end do
    end do

  end subroutine initialize_q_

  subroutine initialize_ra_x_(bd, area, xfx, ra_x)

    ! Arguments
    type(fv_grid_bounds_type), intent(in) :: bd
    real, intent(in) :: area(:,:)
    real, intent(in) :: xfx(:,:)
    real, allocatable, intent(out) :: ra_x(:,:)
    ! Locals
    integer :: j, i

    allocate(ra_x(bd%is:bd%ie, bd%jsd:bd%jed))
    ! do j = bd%jsd, bd%jed
    !    do i = bd%is, bd%ie
    !       ra_x(i,j) = area(i,j) + (xfx(i,j) - xfx(i+1,j))
    !    end do
    ! end do
    ra_x = 1.0

  end subroutine initialize_ra_x_

  subroutine initialize_ra_y_(bd, area, yfx, ra_y)

    ! Arguments
    type(fv_grid_bounds_type), intent(in) :: bd
    real, intent(in) :: area(:,:)
    real, intent(in) :: yfx(:,:)
    real, allocatable, intent(out) :: ra_y(:,:)
    ! Locals
    integer :: j, i

    allocate(ra_y(bd%isd:bd%ied, bd%js:bd%je))
    ! do j = bd%js, bd%je
    !    do i = bd%isd, bd%ied
    !       ra_y(i,j) = area(i,j) + (yfx(i,j) - yfx(i,j+1))
    !    end do
    ! end do
    ra_y = 1.0

  end subroutine initialize_ra_y_

end module input_arrays_mod
