module driver_cpu_mod

  use tp_core_mod, only: fv_tp_2d
  use fv_arrays_mod, only: fv_grid_bounds_type, fv_grid_type
  use input_arrays_mod, only: InputArrays_T
  use output_arrays_mod, only: OutputArrays_T

  implicit none

  private

  public run_driver

  integer, parameter :: hord = 8
  real, parameter :: lim_fac = 1.0

contains

  subroutine run_driver(resolution, n_iterations) bind(C, name='run_driver')

    ! Arguments
    integer, intent(in) :: resolution
    integer, intent(in) :: n_iterations
    ! Locals
    type(fv_grid_bounds_type) :: bd
    type(fv_grid_type) :: gridstruct
    type(InputArrays_T) :: in_arrays
    type(OutputArrays_T) :: out_arrays
    integer :: npx, npy
    integer :: iter
    real :: start, finish

    print *, 'resolution: ', resolution
    print *, 'number of times fv_tp_2d is called: ', n_iterations
    
    ! Initialize
    bd = fv_grid_bounds_type(resolution)
    npx = resolution + 1; npy = resolution + 1
    gridstruct = fv_grid_type(bd, npx, npy, .false., 0)
    in_arrays = InputArrays_T(bd, npx, npy, gridstruct)
    out_arrays = OutputArrays_T(bd)

    ! Run fv_tp_2d
    call cpu_time(start)
    do iter = 1, n_iterations
       call fv_tp_2d( &
            in_arrays%q, in_arrays%crx, in_arrays%cry, &
            npx, npy, hord, &
            out_arrays%fx, out_arrays%fy, &
            in_arrays%xfx, in_arrays%yfx, &
            gridstruct, bd, &
            in_arrays%ra_x, in_arrays%ra_y, &
            lim_fac)
    end do
    call cpu_time(finish)
    
    print *, 'time taken: ', finish - start, 's'
    print *, 'sum(fx): ', sum(out_arrays%fx), ', sum(fy): ', sum(out_arrays%fy)

  end subroutine run_driver

end module driver_cpu_mod
