module output_arrays_mod

  use fv_arrays_mod, only: fv_grid_bounds_type

  implicit none
  
  private

  public OutputArrays_T

  type OutputArrays_T
     real, allocatable :: fx(:, :) ! flux in x (E)
     real, allocatable :: fy(:, :) ! flux in y (N)
  end type OutputArrays_T

  interface OutputArrays_T
     procedure :: output_arrays_allocate_mem_
  end interface OutputArrays_T

contains
  
  function output_arrays_allocate_mem_(bd) result(arrays)

    ! Arguments
    type(fv_grid_bounds_type), intent(in) :: bd
    type(OutputArrays_T) :: arrays ! output

    ! Start
    allocate(arrays%fx(bd%is:bd%ie+1, bd%js:bd%je))
    allocate(arrays%fy(bd%is:bd%ie, bd%js:bd%je+1))

  end function output_arrays_allocate_mem_
  
end module output_arrays_mod
