module output_mod

  implicit none
  
  private

  public Output_T, get_output

  type Output_T
     real :: sum_pkz
     real :: sum_pt
  end type Output_T

contains

  subroutine get_output(filename, output)

    ! Arguments
    character(len=*), intent(in) :: filename
    type(Output_T), intent(out) :: output
    ! Locals
    integer :: file_handle

    ! Start
    print *, 'Reading output data from file: ', trim(filename)

    open(newunit = file_handle, file = filename, form = 'unformatted', status = 'old')

    read(file_handle) output%sum_pkz
    read(file_handle) output%sum_pt
    
    close(file_handle)

  end subroutine get_output

end module output_mod
