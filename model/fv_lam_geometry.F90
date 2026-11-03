! Isolated NASA port: full halo geometry, one bounded spherical tile.
module fv_lam_geometry_mod
 use fv_arrays_mod, only: fv_atmos_type, R_GRID
 use netcdf
 use constants_mod, only: pi
 use mpp_mod, only: mpp_error,FATAL
 use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
 implicit none
 private
 public :: lam_read_geometry, lam_nc_check
contains
 subroutine lam_nc_check(status)
 integer,intent(in)::status
 if(status/=nf90_noerr) call mpp_error(FATAL,'LAM NetCDF: '//trim(nf90_strerror(status)))
 end subroutine
 subroutine lam_read_geometry(a)
 type(fv_atmos_type),intent(inout)::a
 integer::nc,v,d,nx,ny,i,j,h,ng
 real(R_GRID),allocatable::x(:,:),y(:,:)
 character(64)::units
 ng=a%ng;h=3
 call lam_nc_check(nf90_open(trim(a%neststruct%lam_directory)//'/grid.tile7.halo3.nc',nf90_nowrite,nc))
 call lam_nc_check(nf90_inq_dimid(nc,'nxp',d));call lam_nc_check(nf90_inquire_dimension(nc,d,len=nx))
 call lam_nc_check(nf90_inq_dimid(nc,'nyp',d));call lam_nc_check(nf90_inquire_dimension(nc,d,len=ny))
 if(ng/=h.or.nx/=2*(a%npx-1+2*h)+1.or.ny/=2*(a%npy-1+2*h)+1) call mpp_error(FATAL,'LAM halo3 supergrid shape mismatch')
 allocate(x(nx,ny),y(nx,ny))
 call lam_nc_check(nf90_inq_varid(nc,'x',v));call lam_nc_check(nf90_get_att(nc,v,'units',units))
 if(trim(units)/='degrees_east') call mpp_error(FATAL,'LAM x must be degrees_east')
 call lam_nc_check(nf90_get_var(nc,v,x))
 call lam_nc_check(nf90_inq_varid(nc,'y',v));call lam_nc_check(nf90_get_att(nc,v,'units',units))
 if(trim(units)/='degrees_north') call mpp_error(FATAL,'LAM y must be degrees_north')
 call lam_nc_check(nf90_get_var(nc,v,y));call lam_nc_check(nf90_close(nc))
 if(.not.all(ieee_is_finite(x)).or..not.all(ieee_is_finite(y)).or.any(abs(y)>90.)) call mpp_error(FATAL,'LAM invalid coordinates')
 do j=1-ng,a%npy+ng
 do i=1-ng,a%npx+ng
 a%grid_global(i,j,1,1)=x(2*(i+h)-1,2*(j+h)-1)*pi/180.
 a%grid_global(i,j,2,1)=y(2*(i+h)-1,2*(j+h)-1)*pi/180.
 enddo
 enddo
 a%gridstruct%grid_64=a%grid_global(a%bd%isd:a%bd%ied+1,a%bd%jsd:a%bd%jed+1,:,1)
 end subroutine
end module
