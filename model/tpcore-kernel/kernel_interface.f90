module tpcore_kernel_interface

  use iso_c_binding, only: c_int, c_float

  private
  public :: kernel1_orig, kernel1
  
  interface
     subroutine kernel1_orig(nrows, ncols, yfx, fy2, fyy) bind(C, name="kernel1_orig")
       import c_int, c_float
       integer(kind=c_int), value, intent(in) :: nrows, ncols
       real(kind=c_float), dimension(*), intent(in) :: yfx, fy2
       real(kind=c_float), dimension(*), intent(out) :: fyy
     end subroutine kernel1_orig
  end interface

  interface
     subroutine kernel1(nrows, ncols, yfx, fy2, fyy) bind(C, name="kernel1")
       import c_int, c_float
       integer(kind=c_int), value, intent(in) :: nrows, ncols
       real(kind=c_float), dimension(*), intent(in) :: yfx, fy2
       real(kind=c_float), dimension(*), intent(out) :: fyy
     end subroutine kernel1
  end interface

end module tpcore_kernel_interface
