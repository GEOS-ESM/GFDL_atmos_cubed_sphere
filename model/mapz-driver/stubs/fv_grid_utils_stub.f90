!***********************************************************************
!*                   GNU Lesser General Public License                 
!*
!* This file is part of the FV3 dynamical core.
!*
!* The FV3 dynamical core is free software: you can redistribute it 
!* and/or modify it under the terms of the
!* GNU Lesser General Public License as published by the
!* Free Software Foundation, either version 3 of the License, or 
!* (at your option) any later version.
!*
!* The FV3 dynamical core is distributed in the hope that it will be 
!* useful, but WITHOUT ANYWARRANTY; without even the implied warranty 
!* of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
!* See the GNU General Public License for more details.
!*
!* You should have received a copy of the GNU Lesser General Public
!* License along with the FV3 dynamical core.  
!* If not, see <http://www.gnu.org/licenses/>.
!***********************************************************************

!>@brief The module 'fv_grid_utils' contains routines for setting up and
!! computing grid-related quantities. 
!>@details Many of these are useful for computing diagnostics or setting up initial conditions.

module fv_grid_utils_mod

  use iso_fortran_env, only: real64
  use mpp_domains_mod, only: domain2d, BITWISE_EXACT_SUM, BITWISE_EFP_SUM
  use mpp_domains_mod, only: mpp_global_sum
  use fv_mp_mod, only: is_master, mp_reduce_sum
  
  implicit none
  
  private
  
  real, parameter:: ptop_min=1.d-8
  integer, public, parameter :: R_GRID = real64  
  
  public ptop_min, g_sum
  
contains
  
  !>@brief The function 'g_sum' is the fast version of 'globalsum'. 
  real function g_sum(domain, p, ifirst, ilast, jfirst, jlast, ngc, area, mode, reproduce)
    integer, intent(IN) :: ifirst, ilast
    integer, intent(IN) :: jfirst, jlast, ngc
    integer, intent(IN) :: mode  ! if ==1 divided by area
    logical, intent(in), optional :: reproduce
    real, intent(IN) :: p(ifirst:ilast,jfirst:jlast)      ! field to be summed
    real(kind=R_GRID), intent(IN) :: area(ifirst-ngc:ilast+ngc,jfirst-ngc:jlast+ngc)
    type(domain2d), intent(IN) :: domain
    integer :: i,j
    real gsum
    logical, SAVE :: g_sum_initialized = .false.
    real(kind=R_GRID), SAVE :: global_area
    real :: tmp(ifirst:ilast,jfirst:jlast) 
    integer :: err
    
    if ( .not. g_sum_initialized ) then
       global_area = mpp_global_sum(domain, area, flags=BITWISE_EFP_SUM)
       if ( is_master() ) write(*,*) 'Global Area=',global_area
       g_sum_initialized = .true.
    end if
    
    !-------------------------
    ! FMS global sum algorithm:
    !-------------------------
    if ( present(reproduce) ) then
       if (reproduce) then
          gsum = mpp_global_sum(domain, p(:,:)*area(ifirst:ilast,jfirst:jlast), &
               flags=BITWISE_EFP_SUM)
       else
          gsum = mpp_global_sum(domain, p(:,:)*area(ifirst:ilast,jfirst:jlast))
       endif
    else
       !-------------------------
       ! Quick local sum algorithm
       !-------------------------
       gsum = 0.
       do j=jfirst,jlast
          do i=ifirst,ilast
             gsum = gsum + p(i,j)*area(i,j)
          enddo
       enddo
       call mp_reduce_sum(gsum)
    endif
    
    if ( mode==1 ) then
       g_sum = gsum / global_area
    else
       g_sum = gsum
    endif
  end function g_sum

end module fv_grid_utils_mod
