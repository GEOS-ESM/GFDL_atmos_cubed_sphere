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

!>@brief The module 'fv_tracer2d.F90' performs sub-cycled tracer advection.
!>@see \cite lin2004vertically

! Modules Included:
! <table>
! <tr>
!     <th>Module Name</th>
!     <th>Functions Included</th>
!   </tr>
! <table>
!   <tr>
!     <td>boundary_mod</td>
!     <td>nested_grid_BC_apply_intT</td>
!   </tr>
!   <tr>
!     <td>fv_arrays_mod</td>
!     <td>fv_grid_type, fv_nest_type, fv_atmos_type, fv_grid_bounds_type</td>
!   </tr>
!   <tr>
!   <tr>
!     <td>fv_mp_mod</td>
!     <td>mp_reduce_max, ng, mp_gather, is_master, group_halo_update_type, 
!         start_group_halo_update, complete_group_halo_update</td>
!   </tr>
!    <tr>
!     <td>fv_timing_mod</td>
!     <td>timing_on, timing_off</td>
!   </tr>
!  <tr>
!     <td>mpp_mod</td>
!     <td>mpp_error, FATAL, mpp_broadcast, mpp_send, mpp_recv, mpp_sum, mpp_max</td>
!   </tr>
!   <tr>
!     <td>mpp_domains_mod</td>
!     <td>mpp_update_domains, CGRID_NE, domain2d</td>
!   </tr>
!   <tr>
!     <td>tp_core_mod</td>
!     <td>fv_tp_2d, copy_corners</td>
!   </tr>
! </table>

module fv_tracer2d_mod
!============ GEOS-LAM Update Start =============
! Record dry-pilot tracer-stage diagnostics; remove or isolate before production use.
 use lam_probe_mod,only:probe_mark,probe_array,probe_context,probe_level
!============ GEOS-LAM Update End =============
   use tp_core_mod,       only: fv_tp_2d, copy_corners
   use fv_mp_mod,         only: mp_reduce_max
   use fv_mp_mod,         only: ng, mp_gather, is_master
   use fv_mp_mod,         only: group_halo_update_type
   use fv_mp_mod,         only: start_group_halo_update, complete_group_halo_update
   use mpp_domains_mod,   only: mpp_update_domains, CGRID_NE, domain2d, mpp_get_boundary
   use mpp_domains_mod,   only: NON_BITWISE_EXACT_SUM, BITWISE_EXACT_SUM, BITWISE_EFP_SUM
   use fv_timing_mod,     only: timing_on, timing_off
   use boundary_mod,      only: nested_grid_BC_apply_intT
   use fv_arrays_mod,     only: fv_grid_type, fv_flags_type, fv_nest_type, fv_atmos_type, fv_grid_bounds_type, REAL8
!============ GEOS-LAM Update Start =============
! Support bounded tracer transport and dry-pilot diagnostic contracts.
   use mpp_mod,           only: mpp_error, FATAL, mpp_broadcast, mpp_send, mpp_recv, mpp_sum, mpp_max, mpp_pe
!============ GEOS-LAM Update End =============
   use fv_grid_utils_mod, only: g_sum_r8

implicit none
private

public :: tracer_2d, tracer_2d_nested, tracer_2d_1L, offline_tracer_advection
!============ GEOS-LAM Update Start =============
! Exercise the historical GCHP buffer probe; this is not the selected persistent-cell implementation.
public :: set_buffer_probe

real, allocatable, dimension(:,:,:) :: nest_fx_west_accum, nest_fx_east_accum, nest_fx_south_accum, nest_fx_north_accum
! Isolated harness snapshots. Never used unless GBUF_REFRESH=1 is set.
real, allocatable :: gbuf_low(:,:,:,:), gbuf_high(:,:,:,:)
real :: gbuf_t0, gbuf_t1, gbuf_tcap

contains

subroutine set_buffer_probe(low, high, t0, t1, tcap)
   real, intent(in) :: low(:,:,:,:), high(:,:,:,:), t0, t1, tcap
   if (any(shape(low)/=shape(high))) error stop 'GBUF snapshot shapes differ'
   if (t1<=t0.or.tcap<t0.or.tcap>t1) error stop 'GBUF invalid time bracket'
   if (allocated(gbuf_low)) deallocate(gbuf_low,gbuf_high)
   allocate(gbuf_low,source=low)
   allocate(gbuf_high,source=high)
   gbuf_t0=t0
   gbuf_t1=t1
   gbuf_tcap=tcap
end subroutine set_buffer_probe
!============ GEOS-LAM Update End =============

!>@brief The subroutine 'tracer_2d_1L' performs 2-D horizontal-to-lagrangian transport.
!>@details This subroutine is called if 'z_tracer = .true.'
!! It modifies 'tracer_2d' so that each layer uses a different diagnosed number 
!! of split tracer timesteps. This potentially accelerates tracer advection when there
!! is a large difference in layer-maximum wind speeds (cf. polar night jet).
subroutine tracer_2d_1L(q, dp1, mfx, mfy, cx, cy, gridstruct, bd, domain, npx, npy, npz,   &
                        nq,  hord, dt, id_divg, q_pack, nord_tr, trdm, lim_fac, dpA)

      type(fv_grid_bounds_type), intent(IN) :: bd
      integer, intent(IN) :: npx
      integer, intent(IN) :: npy
      integer, intent(IN) :: npz
      integer, intent(IN) :: nq    !< number of tracers to be advected
      integer, intent(IN) :: hord, nord_tr
      integer, intent(IN) :: id_divg
      real   , intent(IN) :: dt, trdm
      real   , intent(IN) :: lim_fac
      type(group_halo_update_type), intent(inout) :: q_pack
      real   , intent(INOUT) :: q(bd%isd:bd%ied,bd%jsd:bd%jed,npz,nq)   !< Tracers
      real   , intent(INOUT) :: dp1(bd%isd:bd%ied,bd%jsd:bd%jed,npz)    !< DELP before dyn_core
      real   , intent(IN   ) :: mfx(bd%is:bd%ie+1,bd%js:bd%je,  npz)    !< Mass Flux X-Dir
      real   , intent(IN   ) :: mfy(bd%is:bd%ie  ,bd%js:bd%je+1,npz)    !< Mass Flux Y-Dir
      real   , intent(IN   ) ::  cx(bd%is:bd%ie+1,bd%jsd:bd%jed  ,npz)  !< Courant Number X-Dir
      real   , intent(IN   ) ::  cy(bd%isd:bd%ied,bd%js :bd%je +1,npz)  !< Courant Number Y-Dir
      real   , optional, intent(OUT) :: dpA(bd%is:bd%ie,bd%js:bd%je,npz)    ! DELP after advection
      type(fv_grid_type), intent(IN), target :: gridstruct
      type(domain2d), intent(INOUT) :: domain

! Local Arrays
      real :: qn2(bd%isd:bd%ied,bd%jsd:bd%jed,nq)   !< 3D tracers
      real :: dp2(bd%is:bd%ie,bd%js:bd%je)
      real :: fx(bd%is:bd%ie+1,bd%js:bd%je )
      real :: fy(bd%is:bd%ie , bd%js:bd%je+1)
      real :: ra_x(bd%is:bd%ie,bd%jsd:bd%jed)
      real :: ra_y(bd%isd:bd%ied,bd%js:bd%je)
      real :: xfx(bd%is:bd%ie+1,bd%jsd:bd%jed  ,npz)
      real :: yfx(bd%isd:bd%ied,bd%js: bd%je+1, npz)
      real :: mfx2(bd%is:bd%ie+1,bd%js:bd%je,npz)
      real :: mfy2(bd%is:bd%ie,bd%js:bd%je+1,npz)
      real ::  cx2(bd%is:bd%ie+1,bd%jsd:bd%jed, npz)
      real ::  cy2(bd%isd:bd%ied,bd%js :bd%je +1, npz)
      real :: c2d(bd%is:bd%ie,bd%js:bd%je)
      real :: cmax(npz), cmin(npz)
!============ GEOS-LAM Update Start =============
! Support bounded tracer transport and dry-pilot diagnostic contracts.
      real(kind=8) :: lam_face_budget(nq)
!============ GEOS-LAM Update End =============
      real :: qmax(npz*(nq+1)), qmin(npz*(nq+1))
      integer :: icount(npz,nq)
      real :: frac
      integer :: nsplt
      integer :: i,j,k,n,it,iq
!============ GEOS-LAM Update Start =============
! Exercise the historical GCHP buffer probe; this is not the selected persistent-cell implementation.
      real :: gbuf_halo_factor
      character(len=32) :: gbuf_halo_string
      integer :: gbuf_halo_status
      character(len=32) :: gbuf_refresh_string
      character(len=32) :: gbuf_stage_label
      character(len=32) :: gbuf_force_string
      integer :: gbuf_refresh_status
      integer :: gbuf_force_status
      logical :: refresh_buffer
      real(kind=8) :: gbuf_time, gbuf_alpha, gbuf_value, gbuf_mass(nq), innerq
      integer :: is0, ie0, js0, je0
!============ GEOS-LAM Update End =============

      real, pointer, dimension(:,:) :: area, rarea
      real, pointer, dimension(:,:,:) :: sin_sg
      real, pointer, dimension(:,:) :: dxa, dya, dx, dy

      integer :: is,  ie,  js,  je
      integer :: isd, ied, jsd, jed

      is  = bd%is
      ie  = bd%ie
      js  = bd%js
      je  = bd%je
      isd = bd%isd
      ied = bd%ied
      jsd = bd%jsd
      jed = bd%jed
!============ GEOS-LAM Update Start =============
! Support bounded tracer transport and dry-pilot diagnostic contracts.
      is0=max(is,4)
      ie0=min(ie,npx-4)
      js0=max(js,4)
      je0=min(je,npy-4)
!============ GEOS-LAM Update End =============

       area => gridstruct%area
      rarea => gridstruct%rarea

      sin_sg => gridstruct%sin_sg
      dxa    => gridstruct%dxa 
      dya    => gridstruct%dya 
      dx     => gridstruct%dx  
      dy     => gridstruct%dy  
!============ GEOS-LAM Update Start =============
! Record dry-pilot tracer-stage diagnostics; remove or isolate before production use.
      lam_face_budget=0.d0
      print *, 'LAM_TRACER_1L_ENTER',is,ie,js,je,isd,ied,jsd,jed
      print *, 'LAM_TRACER_1L_BOUNDS_DXA',lbound(dxa),ubound(dxa)
      print *, 'LAM_TRACER_1L_BOUNDS_DYA',lbound(dya),ubound(dya)
      print *, 'LAM_TRACER_1L_BOUNDS_DX',lbound(dx),ubound(dx)
      print *, 'LAM_TRACER_1L_BOUNDS_DY',lbound(dy),ubound(dy)
      print *, 'LAM_TRACER_1L_BOUNDS_SIN',lbound(sin_sg),ubound(sin_sg)
      flush(6)
!============ GEOS-LAM Update End =============

!$OMP parallel do default(none) shared(is,ie,js,je,isd,ied,jsd,jed,npz,cx,xfx,dxa,dy, &
!$OMP                                  sin_sg,cy,yfx,dya,dx) private(i,j,k)
  do k=1,npz
     do j=jsd,jed
        do i=is,ie+1
           if (cx(i,j,k) > 0.) then
              xfx(i,j,k) = cx(i,j,k)*dxa(i-1,j)*dy(i,j)*sin_sg(i-1,j,3)
           else
              xfx(i,j,k) = cx(i,j,k)*dxa(i,  j)*dy(i,j)*sin_sg(i,  j,1)
           endif
        enddo
     enddo
     do j=js,je+1
        do i=isd,ied
           if (cy(i,j,k) > 0.) then
              yfx(i,j,k) = cy(i,j,k)*dya(i,j-1)*dx(i,j)*sin_sg(i,j-1,4)
           else
              yfx(i,j,k) = cy(i,j,k)*dya(i,j  )*dx(i,j)*sin_sg(i,j,  2)
           endif
        enddo
     enddo
  enddo  ! k-loop
!============ GEOS-LAM Update Start =============
! Record dry-pilot tracer-stage diagnostics; remove or isolate before production use.
  print *, 'LAM_TRACER_1L_FLUX_GEOMETRY_DONE'
  flush(6)
!============ GEOS-LAM Update End =============

                               call timing_on('COMM_TOTAL')
                         call timing_on('COMM_TRACER')
  call complete_group_halo_update(q_pack, domain)
!============ GEOS-LAM Update Start =============
! Exercise the historical GCHP buffer probe; this is not the selected persistent-cell implementation.
  gbuf_halo_factor=1.0
  call get_environment_variable('GBUF_HALO_FACTOR',gbuf_halo_string,status=gbuf_halo_status)
  if (gbuf_halo_status==0) read(gbuf_halo_string,*,iostat=gbuf_halo_status) gbuf_halo_factor
  if (gbuf_halo_factor/=1.0) then
     q(isd:is-1,:,:,:)=q(isd:is-1,:,:,:)*gbuf_halo_factor
     q(ie+1:ied,:,:,:)=q(ie+1:ied,:,:,:)*gbuf_halo_factor
     q(:,jsd:js-1,:,:)=q(:,jsd:js-1,:,:)*gbuf_halo_factor
     q(:,je+1:jed,:,:)=q(:,je+1:jed,:,:)*gbuf_halo_factor
  endif
  print *, 'LAM_TRACER_1L_HALO_DONE'
  print *, 'LAM_TRACER_1L_SEAM',mpp_pe(),q(is-1,js+7,1,1),q(is,js+7,1,1), &
       q(ie,js+7,1,1),q(ie+1,js+7,1,1)
  flush(6)
!============ GEOS-LAM Update End =============
                        call timing_off('COMM_TRACER')
                              call timing_off('COMM_TOTAL')

! Check for levels where Q does not need to be advected
!$OMP parallel do default(none) shared(nq,npz,is,ie,js,je,q,qmax) private(iq,i,j,k,n)
   do iq=1,nq
      do k=1,npz
         n=(iq-1)*npz + k
         qmax(n) = 0.0
         do j=js,je
            do i=is,ie
               qmax(n) = max(qmax(n),q(i,j,k,iq))
            enddo
         enddo
      enddo
   enddo
!$OMP parallel do default(none) shared(is,ie,js,je,npz,nq,&
!$OMP                                  sin_sg,cx,cy,cmax,qmax) private(i,j,k,n)
   do k=1,npz
      cmax(k) = 0.
      if ( k < npz/6 ) then
           do j=js,je
              do i=is,ie
                cmax(k) = max( cmax(k), abs(cx(i,j,k)), abs(cy(i,j,k)) )
              enddo
           enddo
      else
           do j=js,je
              do i=is,ie
                cmax(k) = max( cmax(k), max(abs(cx(i,j,k)),abs(cy(i,j,k)))+1.-sin_sg(i,j,5) )
              enddo
           enddo
      endif
 !!!  if ( is_master() )  write(*,*) 'tracer_2d_1L: k, nsplt, cmax =', k, int(1. + cmax(k)), cmax(k)
    ! add to qmax for allreduce
     n=nq*npz + k
     qmax(n) = cmax(k)
   enddo  ! k-loop

                        call timing_on('COMM_TOTAL')
                            call timing_on('COMM_TRACER_MAX')
  call mp_reduce_max(qmax,(nq+1)*npz)
!============ GEOS-LAM Update Start =============
! Record dry-pilot tracer-stage diagnostics; remove or isolate before production use.
  print *, 'LAM_TRACER_1L_REDUCE_DONE', maxval(qmax)
  flush(6)
!============ GEOS-LAM Update End =============
                           call timing_off('COMM_TRACER_MAX')
                       call timing_off('COMM_TOTAL')

  ! get cmax from allreduce array
  do k=1,npz
     n=nq*npz + k
     cmax(k) = qmax(n)
  enddo  ! k-loop
!============ GEOS-LAM Update Start =============
! Exercise the historical GCHP buffer probe; this is not the selected persistent-cell implementation.
  gbuf_force_string=''
  call get_environment_variable('GBUF_FORCE_SPLIT',gbuf_force_string,status=gbuf_force_status)
  if (gbuf_force_status==0.and.trim(gbuf_force_string)=='1') then
     ! Diagnostic schedule test only: retain the actual transport fluxes and
     ! Courant arrays, but force selected levels through nsplt=2 and 3.
     if(npz<3) error stop 'GBUF forced split requires three levels'
     if(mpp_pe()==0) write(*,'(A,2(1X,ES24.16))') 'GBUF_FORCE_ACTUAL',cmax(1),cmax(2)
     cmax(1)=max(cmax(1),1.25)
     cmax(2)=max(cmax(2),2.25)
  endif
  if (mpp_pe()==0) write(*,'(A,1X,ES24.16,1X,I0,1X,ES24.16)') &
       'GBUF_1L_COURANT',maxval(cmax),maxval(int(1.+cmax)),gbuf_halo_factor
  gbuf_refresh_string=''
  call get_environment_variable('GBUF_REFRESH',gbuf_refresh_string,status=gbuf_refresh_status)
  refresh_buffer=gbuf_refresh_status==0.and.trim(gbuf_refresh_string)=='1'
  if (refresh_buffer) then
     if (.not.allocated(gbuf_low)) error stop 'GBUF 1L refresh snapshots missing'
     if (any(shape(gbuf_low)/=[ie-is+1,je-js+1,npz,nq])) &
          error stop 'GBUF 1L snapshot extent mismatch'
  endif
!============ GEOS-LAM Update End =============

!$OMP parallel do default(none) shared(is,ie,js,je,isd,ied,jsd,jed,npz,cx,xfx, &
!$OMP                                  cy,yfx,mfx,mfy,cmax,mfx2,mfy2,cx2,cy2)   &
!$OMP                          private(nsplt, frac)
  do k=1,npz

     mfx2(:,:,k)=mfx(:,:,k)
     mfy2(:,:,k)=mfy(:,:,k)
     cx2(:,:,k)=cx(:,:,k)
     cy2(:,:,k)=cy(:,:,k)
     nsplt = int(1. + cmax(k))
     if ( nsplt > 1 ) then
        frac  = 1. / real(nsplt)
        do j=jsd,jed
           do i=is,ie+1
               cx2(i,j,k) =  cx(i,j,k) * frac
              xfx(i,j,k) = xfx(i,j,k) * frac
           enddo
        enddo
        do j=js,je
           do i=is,ie+1
              mfx2(i,j,k) = mfx(i,j,k) * frac
           enddo
        enddo
        do j=js,je+1
           do i=isd,ied
              cy2(i,j,k) =  cy(i,j,k) * frac
             yfx(i,j,k) = yfx(i,j,k) * frac
           enddo
        enddo
        do j=js,je+1
           do i=is,ie
              mfy2(i,j,k) = mfy(i,j,k) * frac
           enddo
        enddo
     endif

  enddo

! Begin k-independent tracer transport; can not be OpenMPed because the mpp_update call.
  do k=1,npz

!$OMP parallel do default(none) shared(k,is,ie,js,je,isd,ied,jsd,jed,xfx,area,yfx,ra_x,ra_y)
     do j=jsd,jed
        do i=is,ie
           ra_x(i,j) = area(i,j) + (xfx(i,j,k) - xfx(i+1,j,k))
        enddo
        if ( j>=js .and. j<=je ) then
           do i=isd,ied
              ra_y(i,j) = area(i,j) + (yfx(i,j,k) - yfx(i,j+1,k))
           enddo
        endif
     enddo

     nsplt = int(1. + cmax(k))
!============ GEOS-LAM Update Start =============
! Exercise the historical GCHP buffer probe; this is not the selected persistent-cell implementation.
     if (mpp_pe()==0) write(*,'(A,1X,I0,1X,I0,1X,ES24.16)') &
          'GBUF_1L_LEVEL',k,nsplt,cmax(k)
     if (refresh_buffer) qn2(:,:,:)=q(:,:,k,:)

     do it=1,nsplt
       if(k==1)then
         print *, 'LAM_TRACER_1L_SPLIT',k,it,nsplt
         flush(6)
       endif

       if (refresh_buffer) then
          gbuf_time=gbuf_tcap+real(it-1,kind=8)*dt/real(nsplt,kind=8)
          if (gbuf_time<gbuf_t0.or.gbuf_time>gbuf_t1) &
               error stop 'GBUF 1L substep crosses snapshot bracket'
          gbuf_alpha=(gbuf_time-gbuf_t0)/(gbuf_t1-gbuf_t0)
          gbuf_mass=0.
          do iq=1,nq
             do j=js,je
                do i=is,ie
                   if (i<=3.or.i>=npx-3.or.j<=3.or.j>=npy-3) then
                      gbuf_value=real(real((1.-gbuf_alpha)*gbuf_low(i-is+1,j-js+1,k,iq) &
                           +gbuf_alpha*gbuf_high(i-is+1,j-js+1,k,iq),kind=4),kind=8)
                      gbuf_mass(iq)=gbuf_mass(iq)+(gbuf_value-qn2(i,j,iq))*dp1(i,j,k)*area(i,j)/9.80665d0
                      qn2(i,j,iq)=gbuf_value
                   endif
                enddo
             enddo
          enddo
          call mpp_update_domains(qn2,domain)
          do iq=1,nq
             do j=jsd,jed
                do i=isd,ied
                   if (i<1.or.i>=npx.or.j<1.or.j>=npy) &
                        qn2(i,j,iq)=qn2(max(1,min(npx-1,i)),max(1,min(npy-1,j)),iq)
                enddo
             enddo
          enddo
          if (gbuf_halo_factor/=1..and.nq>1) then
             do iq=2,nq
                do j=jsd,jed
                   do i=isd,ied
                      if (i<1.or.i>=npx.or.j<1.or.j>=npy) &
                           qn2(i,j,iq)=gbuf_halo_factor*qn2(max(1,min(npx-1,i)), &
                                                               max(1,min(npy-1,j)),iq)
                   enddo
                enddo
             enddo
          endif
          if (mpp_pe()==0) then
             write(*,'(A,2(1X,I0),2(1X,ES24.16))') &
                  'GBUF_1L_REFRESH_TIME',k,it,gbuf_time,gbuf_alpha
             write(*,'(A,2(1X,I0),*(1X,ES24.16))') &
                  'GBUF_1L_REFRESH_MASS',k,it,gbuf_mass
          endif
       endif
!============ GEOS-LAM Update End =============

!$OMP parallel do default(none) shared(k,is,ie,js,je,rarea,mfx2,mfy2,dp1,dp2)
        do j=js,je
           do i=is,ie
              dp2(i,j) = dp1(i,j,k) + ((mfx2(i,j,k)-mfx2(i+1,j,k))+(mfy2(i,j,k)-mfy2(i,j+1,k)))*rarea(i,j)
           enddo
        enddo

!============ GEOS-LAM Update Start =============
! Record dry-pilot tracer-stage diagnostics; remove or isolate before production use.
!$OMP parallel do default(none) shared(k,nsplt,it,is,ie,js,je,isd,ied,jsd,jed,npx,npy,cx2,xfx,hord,trdm, &
!$OMP                                  nord_tr,nq,gridstruct,bd,cy2,yfx,mfx2,mfy2,qn2,q,ra_x,ra_y,dp1,dp2,rarea,lim_fac,lam_face_budget, &
!$OMP                                  refresh_buffer,is0,ie0,js0,je0,domain) &
!$OMP                          private(fx,fy,innerq)
        do iq=1,nq
        if(k==1)then
          print *, 'LAM_TRACER_1L_FVTP_BEGIN',iq
          flush(6)
        endif
        if ( nsplt /= 1 .or. refresh_buffer ) then
           if ( it==1 .and. .not.refresh_buffer ) then
!============ GEOS-LAM Update End =============
              do j=jsd,jed
                 do i=isd,ied
                    qn2(i,j,iq) = q(i,j,k,iq)
                 enddo
              enddo
           endif
           call fv_tp_2d(qn2(isd,jsd,iq), cx2(is,jsd,k), cy2(isd,js,k), &
                         npx, npy, hord, fx, fy, xfx(is,jsd,k), yfx(isd,js,k), &
                         gridstruct, bd, ra_x, ra_y, lim_fac, mfx=mfx2(is,js,k), mfy=mfy2(is,js,k))
           if ( it < nsplt ) then   ! not last call
              do j=js,je
              do i=is,ie
                 qn2(i,j,iq) = (qn2(i,j,iq)*dp1(i,j,k)+((fx(i,j)-fx(i+1,j))+(fy(i,j)-fy(i,j+1)))*rarea(i,j))/dp2(i,j)
              enddo
              enddo
           else
              do j=js,je
              do i=is,ie
                 q(i,j,k,iq) = (qn2(i,j,iq)*dp1(i,j,k)+((fx(i,j)-fx(i+1,j))+(fy(i,j)-fy(i,j+1)))*rarea(i,j))/dp2(i,j)
              enddo
              enddo
           endif
        else
           call fv_tp_2d(q(isd,jsd,k,iq), cx2(is,jsd,k), cy2(isd,js,k), &
                         npx, npy, hord, fx, fy, xfx(is,jsd,k), yfx(isd,js,k), &
                         gridstruct, bd, ra_x, ra_y, lim_fac, mfx=mfx2(is,js,k), mfy=mfy2(is,js,k))
           do j=js,je
              do i=is,ie
                 q(i,j,k,iq) = (q(i,j,k,iq)*dp1(i,j,k)+((fx(i,j)-fx(i+1,j))+(fy(i,j)-fy(i,j+1)))*rarea(i,j))/dp2(i,j)
              enddo
           enddo
        endif
!============ GEOS-LAM Update Start =============
! Exercise the historical GCHP buffer probe; this is not the selected persistent-cell implementation.
        innerq=0.
        if(js0<=je0) then
           if(is<=4.and.ie+1>=4) innerq=innerq+sum(dble(fx(4,js0:je0)))
           if(is<=14.and.ie+1>=14) innerq=innerq-sum(dble(fx(14,js0:je0)))
        endif
        if(is0<=ie0) then
           if(js<=4.and.je+1>=4) innerq=innerq+sum(dble(fy(is0:ie0,4)))
           if(js<=14.and.je+1>=14) innerq=innerq-sum(dble(fy(is0:ie0,14)))
        endif
        write(*,'(A,4(1X,I0),1X,ES24.16)')'GBUF_1L_Q_INNER',mpp_pe(),k,it,iq,innerq/9.80665d0
        lam_face_budget(iq)=lam_face_budget(iq)+sum(dble(fx(is,js:je)-fx(ie+1,js:je))) &
             +sum(dble(fy(is:ie,js)-fy(is:ie,je+1)))
        enddo   !  tracer-loop
        if(k==1)then
          print *, 'LAM_TRACER_1L_FVTP_DONE'
          flush(6)
        endif
        if (refresh_buffer .and. nq>1) then
           write(gbuf_stage_label,'("g1lk",I2.2,"s",I2.2)') k,it
           if (it<nsplt) then
              call probe_array(trim(gbuf_stage_label),qn2(is:ie,js:je,2),is,js)
           else
              call probe_array(trim(gbuf_stage_label),q(is:ie,js:je,k,2),is,js)
           endif
        endif
!============ GEOS-LAM Update End =============

        if ( it < nsplt ) then   ! not last call
             do j=js,je
                do i=is,ie
                   dp1(i,j,k) = dp2(i,j)
                enddo
             enddo
                               call timing_on('COMM_TOTAL')
                         call timing_on('COMM_TRACER')
             call mpp_update_domains(qn2, domain)
                        call timing_off('COMM_TRACER')
                              call timing_off('COMM_TOTAL')
        endif
     enddo  ! time-split loop

     if (present(dpA)) then
        dpA(:,:,k)=dp2
     endif

  enddo    ! k-loop
!============ GEOS-LAM Update Start =============
! Record dry-pilot tracer-stage diagnostics; remove or isolate before production use.
  do iq=1,nq
     print *, 'LAM_TRACER_FACE_FLUX',iq,lam_face_budget(iq)/9.80665d0
  enddo
!============ GEOS-LAM Update End =============

end subroutine tracer_2d_1L

!>@brief The subroutine 'tracer_2d' is the standard routine for sub-cycled tracer advection.
subroutine tracer_2d(q, dp1, mfx, mfy, cx, cy, gridstruct, bd, domain, npx, npy, npz,   &
                     nq,  hord, q_split, dt, id_divg, q_pack, nord_tr, trdm, lim_fac, dpA)

      type(fv_grid_bounds_type), intent(IN) :: bd
      integer, intent(IN) :: npx
      integer, intent(IN) :: npy
      integer, intent(IN) :: npz
      integer, intent(IN) :: nq    !< number of tracers to be advected
      integer, intent(IN) :: hord, nord_tr
      integer, intent(IN) :: q_split
      integer, intent(IN) :: id_divg
      real   , intent(IN) :: dt, trdm
      real   , intent(IN) :: lim_fac
      type(group_halo_update_type), intent(inout) :: q_pack
      real   , intent(INOUT) :: q(bd%isd:bd%ied,bd%jsd:bd%jed,npz,nq)   !< Tracers
      real   , intent(INOUT) :: dp1(bd%isd:bd%ied,bd%jsd:bd%jed,npz)    !< DELP before dyn_core
      real   , intent(IN   ) :: mfx(bd%is:bd%ie+1,bd%js:bd%je,  npz)    !< Mass Flux X-Dir
      real   , intent(IN   ) :: mfy(bd%is:bd%ie  ,bd%js:bd%je+1,npz)    !< Mass Flux Y-Dir
      real   , intent(IN   ) ::  cx(bd%is:bd%ie+1,bd%jsd:bd%jed  ,npz)  !< Courant Number X-Dir
      real   , intent(IN   ) ::  cy(bd%isd:bd%ied,bd%js :bd%je +1,npz)  !< Courant Number Y-Dir
      real   , optional, intent(OUT) :: dpA(bd%is:bd%ie,bd%js:bd%je,npz)! DELP after advection
      type(fv_grid_type), intent(IN), target :: gridstruct
      type(domain2d), intent(INOUT) :: domain

! Local Arrays
      real :: dp2(bd%is:bd%ie,bd%js:bd%je)
      real :: fx(bd%is:bd%ie+1,bd%js:bd%je )
      real :: fy(bd%is:bd%ie , bd%js:bd%je+1)
      real :: ra_x(bd%is:bd%ie,bd%jsd:bd%jed)
      real :: ra_y(bd%isd:bd%ied,bd%js:bd%je)
      real :: xfx(bd%is:bd%ie+1,bd%jsd:bd%jed  ,npz)
      real :: yfx(bd%isd:bd%ied,bd%js: bd%je+1, npz)
      real :: mfx2(bd%is:bd%ie+1,bd%js:bd%je,npz)
      real :: mfy2(bd%is:bd%ie,bd%js:bd%je+1,npz)
      real ::  cx2(bd%is:bd%ie+1,bd%jsd:bd%jed, npz)
      real ::  cy2(bd%isd:bd%ied,bd%js :bd%je +1, npz)
      real :: cmax(npz)
      real :: frac, rdt
      integer :: maxsplt, ksplt(npz)
!============ GEOS-LAM Update Start =============
! Exercise the historical GCHP buffer probe; this is not the selected persistent-cell implementation.
      integer :: i,j,k,n,it,iq,halo_status,refresh_status,is0,ie0,js0,je0
      character(32) :: halo_string, probe_label, refresh_string
      real :: halo_factor, faceq, innerq, gbuf_time, gbuf_alpha, gbuf_value, gbuf_mass(nq)
      logical :: halo_probe, refresh_buffer
!============ GEOS-LAM Update End =============

      real, pointer, dimension(:,:) :: area, rarea
      real, pointer, dimension(:,:,:) :: sin_sg
      real, pointer, dimension(:,:) :: dxa, dya, dx, dy

      integer :: is,  ie,  js,  je
      integer :: isd, ied, jsd, jed

      is  = bd%is
      ie  = bd%ie
      js  = bd%js
      je  = bd%je
      isd = bd%isd
      ied = bd%ied
      jsd = bd%jsd
      jed = bd%jed
!============ GEOS-LAM Update Start =============
! Support bounded tracer transport and dry-pilot diagnostic contracts.
      is0=max(is,4); ie0=min(ie,13)
      js0=max(js,4); je0=min(je,13)
!============ GEOS-LAM Update End =============

       area => gridstruct%area
      rarea => gridstruct%rarea

      sin_sg => gridstruct%sin_sg
      dxa    => gridstruct%dxa 
      dya    => gridstruct%dya 
      dx     => gridstruct%dx  
      dy     => gridstruct%dy  

!$OMP parallel do default(none) shared(is,ie,js,je,isd,ied,jsd,jed,npz,cx,xfx,dxa,dy, &  
!$OMP                                  sin_sg,cy,yfx,dya,dx,cmax,q_split)
    do k=1,npz
       do j=jsd,jed
          do i=is,ie+1
             if (cx(i,j,k) > 0.) then
                 xfx(i,j,k) = cx(i,j,k)*dxa(i-1,j)*dy(i,j)*sin_sg(i-1,j,3)
             else
                 xfx(i,j,k) = cx(i,j,k)*dxa(i,j)*dy(i,j)*sin_sg(i,j,1)
             endif
          enddo
       enddo
       do j=js,je+1
          do i=isd,ied
              if (cy(i,j,k) > 0.) then
                  yfx(i,j,k) = cy(i,j,k)*dya(i,j-1)*dx(i,j)*sin_sg(i,j-1,4)
              else
                  yfx(i,j,k) = cy(i,j,k)*dya(i,j)*dx(i,j)*sin_sg(i,j,2)
              endif
          enddo
       enddo
       if ( q_split == 0 ) then
         cmax(k) = 0.
         if ( k < npz/6 ) then
            do j=js,je
               do i=is,ie
                  cmax(k) = max( cmax(k), abs(cx(i,j,k)), abs(cy(i,j,k)) )
               enddo
            enddo
         else
            do j=js,je
               do i=is,ie
                  cmax(k) = max( cmax(k), max(abs(cx(i,j,k)),abs(cy(i,j,k)))+1.-sin_sg(i,j,5) )
               enddo
            enddo
         endif
       endif
    enddo
!============ GEOS-LAM Update Start =============
! Record dry-pilot tracer-stage diagnostics; remove or isolate before production use.
    print *, 'LAM_TRACER_2D_GEOMETRY_DONE',minval(xfx),maxval(xfx),minval(yfx),maxval(yfx)
    flush(6)
!============ GEOS-LAM Update End =============

!--------------------------------------------------------------------------------

! Determine global cmax on levels:
    if ( q_split == 0 ) then
       call mp_reduce_max(cmax,npz)
    else
       cmax = q_split-1.
    endif

! Determine maxsplt for outer iteration loop. 
! Advection will be done base on levels using ksplt
    maxsplt = 0
    do k=1,npz
       ksplt(k) = int(1. + cmax(k))
       maxsplt = max(maxsplt,ksplt(k))
    enddo
!============ GEOS-LAM Update Start =============
! Exercise the historical GCHP buffer probe; this is not the selected persistent-cell implementation.
    call get_environment_variable('GBUF_HALO_FACTOR',halo_string,status=halo_status)
    halo_probe = halo_status == 0
    halo_factor = 1.
    if (halo_probe) read(halo_string,*) halo_factor
    refresh_string=''
    call get_environment_variable('GBUF_REFRESH',refresh_string,status=refresh_status)
    refresh_buffer = refresh_status==0.and.trim(refresh_string)=='1'
    if (refresh_buffer) then
       if (.not.allocated(gbuf_low)) error stop 'GBUF refresh snapshots missing'
       if (any(shape(gbuf_low)/=[ie-is+1,je-js+1,npz,nq])) &
            error stop 'GBUF refresh snapshot extent mismatch'
    endif
    if (mpp_pe()==0) write(*,'(A,1X,I0,1X,ES16.8,1X,I0,1X,ES16.8)') &
         'GBUF_SPLIT',q_split,max(maxval(abs(cx)),maxval(abs(cy))),maxsplt,halo_factor
!============ GEOS-LAM Update End =============

!--------------------------------------------------------------------------------

!$OMP parallel do default(none) shared(is,ie,js,je,isd,ied,jsd,jed,npz,cx,cx2,xfx,mfx,mfx2,cy,cy2,yfx,mfy,mfy2,cmax,ksplt) &
!$OMP                          private( frac )
    do k=1,npz
        mfx2(:,:,k)=mfx(:,:,k)
        mfy2(:,:,k)=mfy(:,:,k)
        cx2(:,:,k)=cx(:,:,k)
        cy2(:,:,k)=cy(:,:,k)
        if( ksplt(k) /= 1 ) then
            frac  = 1. / real(ksplt(k))
            do j=jsd,jed
               do i=is,ie+1
                  cx2(i,j,k) =   cx(i,j,k) * frac
                  xfx(i,j,k) =  xfx(i,j,k) * frac
               enddo
            enddo
            do j=js,je
               do i=is,ie+1
                  mfx2(i,j,k) = mfx(i,j,k) * frac
               enddo
            enddo
            do j=js,je+1
               do i=isd,ied
                  cy2(i,j,k) =  cy(i,j,k) * frac
                  yfx(i,j,k) = yfx(i,j,k) * frac
               enddo
            enddo
            do j=js,je+1
               do i=is,ie
                  mfy2(i,j,k) = mfy(i,j,k) * frac
               enddo
            enddo
        endif
    enddo

    do it=1,maxsplt
                        call timing_on('COMM_TOTAL')
                            call timing_on('COMM_TRACER')
      call complete_group_halo_update(q_pack, domain)
!============ GEOS-LAM Update Start =============
! Exercise the historical GCHP buffer probe; this is not the selected persistent-cell implementation.
      if (refresh_buffer) then
         ! The first active substep starts at the AdvCore cap clock T.
         ! Level k's prescribed sample is the left edge of its own
         ! dt/ksplt(k) flux partition, not it/maxsplt for all levels.
         gbuf_mass=0.
         do k=1,npz
            if (it>ksplt(k)) cycle
            gbuf_time=gbuf_tcap+real(it-1)*dt/real(ksplt(k))
            if (gbuf_time<gbuf_t0.or.gbuf_time>gbuf_t1) &
                 error stop 'GBUF substep crosses snapshot bracket'
            gbuf_alpha=(gbuf_time-gbuf_t0)/(gbuf_t1-gbuf_t0)
            do iq=1,nq
            do j=js,je
            do i=is,ie
               if (i<=3.or.i>=npx-3.or.j<=3.or.j>=npy-3) then
                  gbuf_value=real(real((1.-gbuf_alpha)*gbuf_low(i-is+1,j-js+1,k,iq) &
                       +gbuf_alpha*gbuf_high(i-is+1,j-js+1,k,iq),kind=4))
                  gbuf_mass(iq)=gbuf_mass(iq)+(gbuf_value-q(i,j,k,iq))*dp1(i,j,k)*area(i,j)/9.80665d0
                  q(i,j,k,iq)=gbuf_value
               endif
            enddo
            enddo
            enddo
         enddo
         ! Refreshing after the previous asynchronous halo completion
         ! requires a new exchange before the next PPM stencil reads q.
         do iq=1,nq
            call mpp_update_domains(q(:,:,:,iq),domain)
         enddo
         do iq=1,nq
         do j=jsd,jed
         do i=isd,ied
            if (i<1.or.i>=npx.or.j<1.or.j>=npy) &
                 q(i,j,:,iq)=q(max(1,min(npx-1,i)),max(1,min(npy-1,j)),:,iq)
         enddo
         enddo
         enddo
         if (mpp_pe()==0) then
            write(*,'(A,1X,I0,1X,ES24.16,1X,ES24.16)') &
                 'GBUF_REFRESH_TIME',it,gbuf_tcap+real(it-1)*dt/real(ksplt(1)), &
                 (gbuf_tcap+real(it-1)*dt/real(ksplt(1))-gbuf_t0)/(gbuf_t1-gbuf_t0)
            write(*,'(A,1X,I0,*(1X,ES24.16))') 'GBUF_REFRESH_MASS',it,gbuf_mass
         endif
      endif
      if (halo_probe .and. nq>1) then
         do j=jsd,jed
         do i=isd,ied
            if(i<1.or.i>=npx.or.j<1.or.j>=npy) &
               q(i,j,:,2:nq)=halo_factor*q(max(1,min(npx-1,i)),max(1,min(npy-1,j)),:,2:nq)
         enddo
         enddo
      endif
      if (nq>1) then
         write(probe_label,'("halo",I2.2)')it
         call probe_array(trim(probe_label),q(:,:,:,2),isd,jsd)
      endif
      if (mpp_pe()==0) write(*,'(A,1X,I0,1X,I0,1X,ES16.8)') &
           'GBUF_SUBSTEP',it,maxsplt,maxval(abs(cx))/real(maxsplt)
      if (mpp_pe()==0) write(*,'(A,1X,I0,4(1X,ES16.8))') &
           'GBUF_DIRECTION',it,minval(cx2),maxval(cx2),minval(cy2),maxval(cy2)
      print *, 'LAM_TRACER_2D_HALO_DONE',minval(q),maxval(q),minval(dp1),maxval(dp1)
      flush(6)
!============ GEOS-LAM Update End =============
                           call timing_off('COMM_TRACER')
                       call timing_off('COMM_TOTAL')

!============ GEOS-LAM Update Start =============
! Handle regional mass-flux values used at dynamics and tracer stages.
!$OMP parallel do default(none) shared(is,ie,js,je,isd,ied,jsd,jed,npz,dp1,mfx2,mfy2,rarea,nq,&
!$OMP                                  area,xfx,yfx,q,cx2,cy2,npx,npy,hord,gridstruct,bd,it,ksplt,nord_tr,trdm,lim_fac,dpA, &
!$OMP                                  is0,ie0,js0,je0) &
!$OMP                          private(dp2, ra_x, ra_y, fx, fy, i, j, k, n, iq, faceq, innerq)
!============ GEOS-LAM Update End =============
     do k=1,npz

       if ( it <= ksplt(k) ) then

         do j=js,je
            do i=is,ie
               dp2(i,j) = dp1(i,j,k) + ((mfx2(i,j,k)-mfx2(i+1,j,k))+(mfy2(i,j,k)-mfy2(i,j+1,k)))*rarea(i,j)
            enddo
         enddo

         do j=jsd,jed
            do i=is,ie
               ra_x(i,j) = area(i,j) + (xfx(i,j,k) - xfx(i+1,j,k))
               if (cx2(i,j,k) > 1.0)  write(*,*) 'cx2(i,j,k) > 1.0 : ', cx2(i,j,k), i, j, k
            enddo
         enddo
         do j=js,je
            do i=isd,ied
               ra_y(i,j) = area(i,j) + (yfx(i,j,k) - yfx(i,j+1,k))
               if (cy2(i,j,k) > 1.0)  write(*,*) 'cy2(i,j,k) > 1.0 : ', cy2(i,j,k), i, j, k
            enddo
         enddo

         do iq=1,nq
           if ( it==1 .and. trdm>1.e-4 ) then
            call fv_tp_2d(q(isd,jsd,k,iq), cx2(is,jsd,k), cy2(isd,js,k), &
                          npx, npy, hord, fx, fy, xfx(is,jsd,k), yfx(isd,js,k), &
                          gridstruct, bd, ra_x, ra_y, lim_fac, mfx=mfx2(is,js,k), mfy=mfy2(is,js,k),   &
                          mass=dp1(isd,jsd,k), nord=nord_tr, damp_c=trdm)
           else
            call fv_tp_2d(q(isd,jsd,k,iq), cx2(is,jsd,k), cy2(isd,js,k), &
                          npx, npy, hord, fx, fy, xfx(is,jsd,k), yfx(isd,js,k), &
                          gridstruct, bd, ra_x, ra_y, lim_fac, mfx=mfx2(is,js,k), mfy=mfy2(is,js,k))
           endif
!============ GEOS-LAM Update Start =============
! Exercise the historical GCHP buffer probe; this is not the selected persistent-cell implementation.
           if (iq>1) then
              faceq=0.
              if(is==1) faceq=faceq+sum(fx(is,js:je))
              if(ie==npx-1) faceq=faceq-sum(fx(ie+1,js:je))
              if(js==1) faceq=faceq+sum(fy(is:ie,js))
              if(je==npy-1) faceq=faceq-sum(fy(is:ie,je+1))
              write(*,'(A,4(1X,I0),1X,ES24.16)')'GBUF_Q_FACE',mpp_pe(),it,k,iq,faceq/9.80665d0
           endif
           innerq=0.
           if (js0<=je0) then
              if (is<=4.and.ie+1>=4) innerq=innerq+sum(fx(4,js0:je0))
              if (is<=14.and.ie+1>=14) innerq=innerq-sum(fx(14,js0:je0))
           endif
           if (is0<=ie0) then
              if (js<=4.and.je+1>=4) innerq=innerq+sum(fy(is0:ie0,4))
              if (js<=14.and.je+1>=14) innerq=innerq-sum(fy(is0:ie0,14))
           endif
           write(*,'(A,4(1X,I0),1X,ES24.16)')'GBUF_Q_INNER',mpp_pe(),it,k,iq,innerq/9.80665d0
!============ GEOS-LAM Update End =============
           do j=js,je
              do i=is,ie
                 q(i,j,k,iq) = ( q(i,j,k,iq)*dp1(i,j,k) + &
                               ((fx(i,j)-fx(i+1,j))+(fy(i,j)-fy(i,j+1)))*rarea(i,j) )/dp2(i,j)
              enddo
           enddo
         enddo

         if ( it /= ksplt(k) ) then
            do j=js,je
               do i=is,ie
                  dp1(i,j,k) = dp2(i,j)
               enddo
            enddo
         elseif (it == ksplt(k)) then
            if (present(dpA)) then
               dpA(:,:,k)=dp2
            end if
         endif

       endif   ! ksplt(k)

     enddo ! npz
!============ GEOS-LAM Update Start =============
! Record dry-pilot tracer-stage diagnostics; remove or isolate before production use.

     if (nq>1) then
        write(probe_label,'("adv",I2.2)')it
        call probe_array(trim(probe_label),q(is:ie,js:je,:,2),is,js)
     endif
!============ GEOS-LAM Update End =============

     if ( it /= maxsplt ) then
                      call timing_on('COMM_TOTAL')
                          call timing_on('COMM_TRACER')
           call start_group_halo_update(q_pack, q, domain)
                          call timing_off('COMM_TRACER')
                      call timing_off('COMM_TOTAL')
     endif

   enddo  ! maxsplt


end subroutine tracer_2d


subroutine tracer_2d_nested(q, dp1, mfx, mfy, cx, cy, gridstruct, bd, domain, npx, npy, npz,   &
                     nq,  hord, q_split, dt, id_divg, q_pack, nord_tr, trdm, &
                     k_split, neststruct, parent_grid, lim_fac)

      type(fv_grid_bounds_type), intent(IN) :: bd
      integer, intent(IN) :: npx
      integer, intent(IN) :: npy
      integer, intent(IN) :: npz
      integer, intent(IN) :: nq    !< number of tracers to be advected
      integer, intent(IN) :: hord, nord_tr
      integer, intent(IN) :: q_split, k_split
      integer, intent(IN) :: id_divg
      real   , intent(IN) :: dt, trdm
      real   , intent(IN) :: lim_fac
      type(group_halo_update_type), intent(inout) :: q_pack
      real   , intent(INOUT) :: q(bd%isd:bd%ied,bd%jsd:bd%jed,npz,nq)   !< Tracers
      real   , intent(INOUT) :: dp1(bd%isd:bd%ied,bd%jsd:bd%jed,npz)    !< DELP before dyn_core
!============ GEOS-LAM Update Start =============
! Handle regional mass-flux values used at dynamics and tracer stages.
      real   , intent(IN) :: mfx(bd%is:bd%ie+1,bd%js:bd%je,  npz)    !< Full-step mass flux X
      real   , intent(IN) :: mfy(bd%is:bd%ie  ,bd%js:bd%je+1,npz)    !< Full-step mass flux Y
      real   , intent(IN) ::  cx(bd%is:bd%ie+1,bd%jsd:bd%jed  ,npz)  !< Full-step Courant X
      real   , intent(IN) ::  cy(bd%isd:bd%ied,bd%js :bd%je +1,npz)  !< Full-step Courant Y
      type(fv_grid_type), intent(IN), target :: gridstruct
      type(fv_nest_type), intent(INOUT) :: neststruct
      type(fv_atmos_type), pointer, intent(INOUT) :: parent_grid
!============ GEOS-LAM Update End =============
      type(domain2d), intent(INOUT) :: domain

! Local Arrays
      real :: dp2(bd%is:bd%ie,bd%js:bd%je)
      real :: fx(bd%is:bd%ie+1,bd%js:bd%je )
      real :: fy(bd%is:bd%ie , bd%js:bd%je+1)
      real :: ra_x(bd%is:bd%ie,bd%jsd:bd%jed)
      real :: ra_y(bd%isd:bd%ied,bd%js:bd%je)
      real :: xfx(bd%is:bd%ie+1,bd%jsd:bd%jed  ,npz)
      real :: yfx(bd%isd:bd%ied,bd%js: bd%je+1, npz)
!============ GEOS-LAM Update Start =============
! Keep internal tracer substep flux and Courant copies private to preserve full-step DYN values.
      real :: mfx_sub(bd%is:bd%ie+1,bd%js:bd%je,npz)
      real :: mfy_sub(bd%is:bd%ie,bd%js:bd%je+1,npz)
      real ::  cx_sub(bd%is:bd%ie+1,bd%jsd:bd%jed,npz)
      real ::  cy_sub(bd%isd:bd%ied,bd%js:bd%je+1,npz)
!============ GEOS-LAM Update End =============
      real :: cmax(npz)
      real :: cmax_t
      real :: c_global
      real :: frac, rdt
      integer :: nsplt, nsplt_parent, msg_split_steps = 1
      integer :: i,j,k,it,iq

      real, pointer, dimension(:,:) :: area, rarea
      real, pointer, dimension(:,:,:) :: sin_sg
      real, pointer, dimension(:,:) :: dxa, dya, dx, dy

      integer :: is,  ie,  js,  je
      integer :: isd, ied, jsd, jed

      is  = bd%is
      ie  = bd%ie
      js  = bd%js
      je  = bd%je
      isd = bd%isd
      ied = bd%ied
      jsd = bd%jsd
      jed = bd%jed

       area => gridstruct%area
      rarea => gridstruct%rarea

      sin_sg => gridstruct%sin_sg
      dxa    => gridstruct%dxa 
      dya    => gridstruct%dya 
      dx     => gridstruct%dx  
      dy     => gridstruct%dy  

!$OMP parallel do default(none) shared(is,ie,js,je,isd,ied,jsd,jed,npz,cx,xfx,dxa,dy, &
!$OMP                                  sin_sg,cy,yfx,dya,dx)
      do k=1,npz
         do j=jsd,jed
            do i=is,ie+1
               if (cx(i,j,k) > 0.) then
                  xfx(i,j,k) = cx(i,j,k)*dxa(i-1,j)*dy(i,j)*sin_sg(i-1,j,3)
               else
                  xfx(i,j,k) = cx(i,j,k)*dxa(i,j)*dy(i,j)*sin_sg(i,j,1)
               endif
            enddo
         enddo
         do j=js,je+1
            do i=isd,ied
               if (cy(i,j,k) > 0.) then
                  yfx(i,j,k) = cy(i,j,k)*dya(i,j-1)*dx(i,j)*sin_sg(i,j-1,4)
               else
                  yfx(i,j,k) = cy(i,j,k)*dya(i,j)*dx(i,j)*sin_sg(i,j,2)
               endif
            enddo
         enddo
      enddo

!--------------------------------------------------------------------------------
  if ( q_split == 0 ) then
! Determine nsplt

!$OMP parallel do default(none) shared(is,ie,js,je,npz,cmax,cx,cy,sin_sg) &
!$OMP                          private(cmax_t )
      do k=1,npz
         cmax(k) = 0.
         if ( k < 4 ) then
! Top layers: C < max( abs(c_x), abs(c_y) )
            do j=js,je
               do i=is,ie
                  cmax_t  = max( abs(cx(i,j,k)), abs(cy(i,j,k)) )
                  cmax(k) = max( cmax_t, cmax(k) )
               enddo
            enddo
         else
            do j=js,je
               do i=is,ie
                  cmax_t  = max(abs(cx(i,j,k)), abs(cy(i,j,k))) + 1.-sin_sg(i,j,5)
                  cmax(k) = max( cmax_t, cmax(k) )
               enddo
            enddo
         endif
      enddo
      call mp_reduce_max(cmax,npz)

!============ GEOS-LAM Update Start =============
! Support bounded tracer transport and dry-pilot diagnostic contracts.
! Find the global maximum Courant number for private substep copies.
!============ GEOS-LAM Update End =============
      c_global = cmax(1)
      if ( npz /= 1 ) then                ! if NOT shallow water test case
         do k=2,npz
            c_global = max(cmax(k), c_global)
         enddo
      endif
      nsplt = int(1. + c_global)
!!    if ( is_master() .and. nsplt > 3 )  write(*,*) 'Tracer_2d_split=', nsplt, c_global
   else
      nsplt = q_split
      if (gridstruct%nested .and. neststruct%nestbctype > 1) msg_split_steps = max(q_split/parent_grid%flagstruct%q_split,1)
   endif

!--------------------------------------------------------------------------------

   frac  = 1. / real(nsplt)

!============ GEOS-LAM Update Start =============
! Keep internal tracer substep flux and Courant copies private to preserve full-step DYN values.
! Keep the full-step capacitor arguments unchanged for fv_dynamics export.
! Only these private arrays carry one tracer substep's flux and Courant.
   mfx_sub = mfx
   mfy_sub = mfy
    cx_sub = cx
    cy_sub = cy

      if( nsplt /= 1 ) then
!$OMP parallel do default(none) shared(is,ie,js,je,isd,ied,jsd,jed,npz,cx_sub,frac,xfx,mfx_sub,cy_sub,yfx,mfy_sub)
          do k=1,npz
             do j=jsd,jed
                do i=is,ie+1
                   cx_sub(i,j,k) = cx_sub(i,j,k) * frac
!============ GEOS-LAM Update End =============
                   xfx(i,j,k) = xfx(i,j,k) * frac
                enddo
             enddo
             do j=js,je
                do i=is,ie+1
!============ GEOS-LAM Update Start =============
! Keep internal tracer substep flux and Courant copies private to preserve full-step DYN values.
                   mfx_sub(i,j,k) = mfx_sub(i,j,k) * frac
!============ GEOS-LAM Update End =============
                enddo
             enddo

             do j=js,je+1
                do i=isd,ied
!============ GEOS-LAM Update Start =============
! Support bounded tracer transport and dry-pilot diagnostic contracts.
                   cy_sub(i,j,k) = cy_sub(i,j,k) * frac
!============ GEOS-LAM Update End =============
                  yfx(i,j,k) = yfx(i,j,k) * frac
                enddo
             enddo

             do j=js,je+1
                do i=is,ie
!============ GEOS-LAM Update Start =============
! Keep internal tracer substep flux and Courant copies private to preserve full-step DYN values.
                  mfy_sub(i,j,k) = mfy_sub(i,j,k) * frac
!============ GEOS-LAM Update End =============
                enddo
             enddo
          enddo
      endif


    do it=1,nsplt
       if ( gridstruct%nested ) then
          neststruct%tracer_nest_timestep = neststruct%tracer_nest_timestep + 1
       end if
                        call timing_on('COMM_TOTAL')
                            call timing_on('COMM_TRACER')
      call complete_group_halo_update(q_pack, domain)
                           call timing_off('COMM_TRACER')
                       call timing_off('COMM_TOTAL')
	    
      if (gridstruct%nested) then
            do iq=1,nq
!============ GEOS-LAM Update Start =============
! Record dry-pilot tracer-stage diagnostics; remove or isolate before production use.
                 call probe_mark('BC_REQUEST_q')
                 call nested_grid_BC_apply_intT(q(isd:ied,jsd:jed,:,iq), &
                      0, 0, npx, npy, npz, bd, &
                      real(neststruct%tracer_nest_timestep)+real(merge(-1,nsplt*k_split,neststruct%regional)), real(nsplt*k_split), &
                 neststruct%q_BC(iq), bctype=neststruct%nestbctype  )
!============ GEOS-LAM Update End =============
           enddo
      endif


!============ GEOS-LAM Update Start =============
! Keep internal tracer substep flux and Courant copies private to preserve full-step DYN values.
!$OMP parallel do default(none) shared(is,ie,js,je,isd,ied,jsd,jed,npz,dp1,mfx_sub,mfy_sub,rarea,nq, &
!$OMP                                  area,xfx,yfx,q,cx_sub,cy_sub,npx,npy,hord,gridstruct,bd,it,nsplt,nord_tr,trdm,lim_fac) &
!$OMP                          private(dp2, ra_x, ra_y, fx, fy)
!============ GEOS-LAM Update End =============
      do k=1,npz

         do j=js,je
            do i=is,ie
!============ GEOS-LAM Update Start =============
! Keep internal tracer substep flux and Courant copies private to preserve full-step DYN values.
               dp2(i,j) = dp1(i,j,k) + ((mfx_sub(i,j,k)-mfx_sub(i+1,j,k))+(mfy_sub(i,j,k)-mfy_sub(i,j+1,k)))*rarea(i,j)
!============ GEOS-LAM Update End =============
            enddo
         enddo

         do j=jsd,jed
            do i=is,ie
               ra_x(i,j) = area(i,j) + (xfx(i,j,k) - xfx(i+1,j,k))
            enddo
         enddo
         do j=js,je
            do i=isd,ied
               ra_y(i,j) = area(i,j) + (yfx(i,j,k) - yfx(i,j+1,k))
            enddo
         enddo

         do iq=1,nq
         if ( it==1 .and. trdm>1.e-4 ) then
!============ GEOS-LAM Update Start =============
! Keep internal tracer substep flux and Courant copies private to preserve full-step DYN values.
            call fv_tp_2d(q(isd,jsd,k,iq), cx_sub(is,jsd,k), cy_sub(isd,js,k), &
                          npx, npy, hord, fx, fy, xfx(is,jsd,k), yfx(isd,js,k), &
                          gridstruct, bd, ra_x, ra_y, lim_fac, mfx=mfx_sub(is,js,k), mfy=mfy_sub(is,js,k),   &
                          mass=dp1(isd,jsd,k), nord=nord_tr, damp_c=trdm)
         else
            call fv_tp_2d(q(isd,jsd,k,iq), cx_sub(is,jsd,k), cy_sub(isd,js,k), &
                          npx, npy, hord, fx, fy, xfx(is,jsd,k), yfx(isd,js,k), &
                          gridstruct, bd, ra_x, ra_y, lim_fac, mfx=mfx_sub(is,js,k), mfy=mfy_sub(is,js,k))
!============ GEOS-LAM Update End =============
         endif
            do j=js,je
               do i=is,ie
                  q(i,j,k,iq) = ( q(i,j,k,iq)*dp1(i,j,k) + &
                                ((fx(i,j)-fx(i+1,j))+(fy(i,j)-fy(i,j+1)))*rarea(i,j) )/dp2(i,j)
               enddo
               enddo
          enddo
      enddo ! npz

      if ( it /= nsplt ) then
                      call timing_on('COMM_TOTAL')
                          call timing_on('COMM_TRACER')
           call start_group_halo_update(q_pack, q, domain)
                          call timing_off('COMM_TRACER')
                      call timing_off('COMM_TOTAL')
      endif
           !Apply nested-grid BCs
           if ( gridstruct%nested ) then
              do iq=1,nq


!============ GEOS-LAM Update Start =============
! Record dry-pilot tracer-stage diagnostics; remove or isolate before production use.
                 call probe_mark('BC_REQUEST_q')


!============ GEOS-LAM Update End =============
                 call nested_grid_BC_apply_intT(q(isd:ied,jsd:jed,:,iq), &
                      0, 0, npx, npy, npz, bd, &
                      real(neststruct%tracer_nest_timestep), real(nsplt*k_split), &
                 neststruct%q_BC(iq), bctype=neststruct%nestbctype  )
!============ GEOS-LAM Update Start =============
! Record dry-pilot tracer-stage diagnostics; remove or isolate before production use.
                 if (iq>=2.and.iq<=3) print *, 'LAM_AUX_TRANSPORT',mpp_pe(),iq, &
                      neststruct%q_BC(iq)%stage_seconds + &
                      real(neststruct%tracer_nest_timestep)/real(nsplt*k_split)* &
                      neststruct%q_BC(iq)%stage_dt, &
                      q(isd,jsd,1,iq),q(ied,jsd,1,iq), &
                      q(isd,jed,1,iq),q(ied,jed,1,iq)
!============ GEOS-LAM Update End =============

              end do
           end if


   enddo  ! nsplt

   if ( id_divg > 0 ) then
        rdt = 1./(frac*dt)

!$OMP parallel do default(none) shared(is,ie,js,je,npz,dp1,xfx,yfx,rarea,rdt)
        do k=1,npz
        do j=js,je
           do i=is,ie
              dp1(i,j,k) = ((xfx(i+1,j,k)-xfx(i,j,k)) + (yfx(i,j+1,k)-yfx(i,j,k)))*rarea(i,j)*rdt
           enddo
        enddo
        enddo
   endif

 end subroutine tracer_2d_nested

!============ GEOS-LAM Update Start =============
! Record the rejected offline exterior-face experiment for review, not production adoption.
subroutine offline_tracer_advection(q, pleB, pleA, mfx, mfy, cx, cy, &
                                    gridstruct, flagstruct, bd, domain, &
                                    npx, npy, npz,   &
                                    nq, dt, q_split, bounded_lam, native_mfx, native_mfy)
!============ GEOS-LAM Update End =============

      use fv_mapz_mod,        only: mapn_tracer
      use fv_fill_mod,        only: fillz

      integer, intent(IN) :: npx
      integer, intent(IN) :: npy
      integer, intent(IN) :: npz
      integer, intent(IN) :: nq    ! number of tracers to be advected
      integer, intent(IN) :: q_split
!============ GEOS-LAM Update Start =============
! Record the rejected offline exterior-face experiment for review, not production adoption.
      logical, intent(IN), optional :: bounded_lam
      type(fv_grid_bounds_type), intent(IN   ) :: bd
      real, intent(IN), optional :: native_mfx(bd%is:bd%ie+1,bd%js:bd%je,npz)
      real, intent(IN), optional :: native_mfy(bd%is:bd%ie,bd%js:bd%je+1,npz)
!============ GEOS-LAM Update End =============
      type(fv_flags_type), intent(INOUT) :: flagstruct
      type(fv_grid_type), intent(IN), target :: gridstruct
      type(domain2D), intent(INOUT) :: domain

      real, intent(IN   ) :: dt
      real, intent(IN   ) ::pleB(bd%is:bd%ie,bd%js:bd%je,npz+1)      ! DELP before dyn_core
      real, intent(INOUT) ::pleA(bd%is:bd%ie,bd%js:bd%je,npz+1)      ! DELP after dyn_core
      real, intent(IN   ) ::  cx(bd%is:bd%ie,bd%js:bd%je,npz)        ! Courant Number X-Dir
      real, intent(IN   ) ::  cy(bd%is:bd%ie,bd%js:bd%je,npz)        ! Courant Number Y-Dir
      real, intent(IN   ) :: mfx(bd%is:bd%ie,bd%js:bd%je,npz)        ! Mass Flux X-Dir
      real, intent(IN   ) :: mfy(bd%is:bd%ie,bd%js:bd%je,npz)        ! Mass Flux Y-Dir
      real, intent(INOUT) ::   q(bd%is:bd%ie,bd%js:bd%je,npz,nq)     ! Tracers
! Local Arrays
      real ::   xL(bd%isd:bd%ied+1,bd%jsd:bd%jed  ,npz)  ! X-Dir for MPP Updates
      real ::   yL(bd%isd:bd%ied  ,bd%jsd:bd%jed+1,npz)  ! Y-Dir for MPP Updates
      real ::  cxL(bd%is :bd%ie +1,bd%jsd:bd%jed  ,npz)  ! Courant Number X-Dir
      real ::  cyL(bd%isd:bd%ied  ,bd%js :bd%je +1,npz)  ! Courant Number Y-Dir
      real :: mfxL(bd%is :bd%ie +1,bd%js :bd%je   ,npz)  ! Mass Flux X-Dir
      real :: mfyL(bd%is :bd%ie   ,bd%js :bd%je +1,npz)  ! Mass Flux Y-Dir
!============ GEOS-LAM Update Start =============
! Handle regional mass-flux values used at dynamics and tracer stages.
      real :: cxa(bd%isd:bd%ied,bd%jsd:bd%jed,npz)
      real :: cya(bd%isd:bd%ied,bd%jsd:bd%jed,npz)
      real :: mfxa(bd%isd:bd%ied,bd%jsd:bd%jed,npz)
      real :: mfya(bd%isd:bd%ied,bd%jsd:bd%jed,npz)
      real ::  dpL(bd%isd:bd%ied  ,bd%jsd:bd%jed  ,npz)  ! Pressure Thickness
      real ::  dpA(bd%is :bd%ie   ,bd%js :bd%je   ,npz)  ! Pressure Thickness
      real :: dp_inner_max
!============ GEOS-LAM Update End =============
! Local Tracer Arrays
      real ::   q3(bd%isd:bd%ied,bd%jsd:bd%jed, npz,nq)! Ghosted 3D Tracers
      real ::   q2(bd%is :bd%ie ,               npz   )! 2D Tmp
      integer :: kord_tracers(nq)
! Local Buffer Arrarys
      real :: wbuffer(bd%js:bd%je,npz)
      real :: sbuffer(bd%is:bd%ie,npz)
      real :: ebuffer(bd%js:bd%je,npz)
      real :: nbuffer(bd%is:bd%ie,npz)
! Local Remap Arrays
      real  pe1(bd%is:bd%ie,npz+1)
      real  pe2(bd%is:bd%ie,npz+1)
      real  dp1(bd%is:bd%ie,npz)
      real  dp2(bd%is:bd%ie,npz)

! Local indices
      integer     :: i,j,k,n,iq

      real :: scalingFactor

      type(group_halo_update_type), save :: i_pack

      integer :: is,  ie,  js,  je
      integer :: isd, ied, jsd, jed

      is  = bd%is
      ie  = bd%ie
      js  = bd%js
      je  = bd%je
      isd = bd%isd
      ied = bd%ied
      jsd = bd%jsd
      jed = bd%jed
!============ GEOS-LAM Update Start =============
! Record dry-pilot tracer-stage diagnostics; remove or isolate before production use.
      print *, 'LAM_OFFLINE_TRACER_ENTER', nq, minval(q), maxval(q)
      flush(6)
!============ GEOS-LAM Update End =============

! Time-step
! Fill CX/CY C-Grid boundaries and update ghost regions
    xL(is:ie,js:je,:) = cx(:,:,:)
    yL(is:ie,js:je,:) = cy(:,:,:)
!============ GEOS-LAM Update Start =============
! Support bounded tracer transport and dry-pilot diagnostic contracts.
    do j=jsd,jed
       do i=isd,ied
          cxa(i,j,:)=cx(max(is,min(ie,i)),max(js,min(je,j)),:)
          cya(i,j,:)=cy(max(is,min(ie,i)),max(js,min(je,j)),:)
       enddo
    enddo
    call mpp_update_domains(cxa,domain)
    call mpp_update_domains(cya,domain)
    ebuffer = cx(ie,js:je,:)
    nbuffer = cy(is:ie,je,:)
!============ GEOS-LAM Update End =============
    call mpp_get_boundary(xL, yL, domain, &
                          wbufferx=wbuffer, ebufferx=ebuffer, &
                          sbuffery=sbuffer, nbuffery=nbuffer, &
                          gridtype=CGRID_NE )
!============ GEOS-LAM Update Start =============
! Record dry-pilot tracer-stage diagnostics; remove or isolate before production use.
    print *, 'LAM_CGRID_EDGE',mpp_pe(),cx(ie,js+7,1),ebuffer(js+7,1), &
         cy(is+min(7,ie-is),je,1),nbuffer(is+min(7,ie-is),1)
    xL(ie+1,js:je,:) = ebuffer
    yL(is:ie,je+1,:) = nbuffer
    ! Diagnostic bounded-grid repair: define physical exterior C-grid halos
    ! before the MPI update overwrites halos on inter-rank boundaries.
    do j=jsd,jed
       do i=isd,ied+1
          if (i < is .or. i > ie+1 .or. j < js .or. j > je) &
             xL(i,j,:) = xL(max(is,min(ie+1,i)),max(js,min(je,j)),:)
       enddo
    enddo
    do j=jsd,jed+1
       do i=isd,ied
          if (i < is .or. i > ie .or. j < js .or. j > je+1) &
             yL(i,j,:) = yL(max(is,min(ie,i)),max(js,min(je+1,j)),:)
       enddo
    enddo
    print *, 'LAM_CGRID_BEFORE_MPP',minval(xL),maxval(xL),minval(yL),maxval(yL)
    flush(6)
    call mpp_update_domains( xL, yL, domain, gridtype=CGRID_NE, complete=.true.)
    xL(ie+1,jsd:jed,:)=cxa(ie+1,jsd:jed,:)
    yL(isd:ied,je+1,:)=cya(isd:ied,je+1,:)
    cxL(is:ie+1,jsd:jed,:) = xL(is:ie+1,jsd:jed,:)
    cyL(isd:ied,js:je+1,:) = yL(isd:ied,js:je+1,:)
    print *, 'LAM_CGRID_AFTER_MPP',minval(cxL),maxval(cxL),minval(cyL),maxval(cyL)
    print *, 'LAM_CGRID_AFTER_MPP_EDGE',mpp_pe(),cxL(is,js+7,1),cxL(ie,js+7,1),cxL(ie+1,js+7,1)
    flush(6)
!============ GEOS-LAM Update End =============

! Fill MFX/MFY C-Grid boundaries
    xL(is:ie,js:je,:) = mfx(:,:,:)
    yL(is:ie,js:je,:) = mfy(:,:,:)
!============ GEOS-LAM Update Start =============
! Handle regional mass-flux values used at dynamics and tracer stages.
    do j=jsd,jed
       do i=isd,ied
          mfxa(i,j,:)=mfx(max(is,min(ie,i)),max(js,min(je,j)),:)
          mfya(i,j,:)=mfy(max(is,min(ie,i)),max(js,min(je,j)),:)
       enddo
    enddo
    call mpp_update_domains(mfxa,domain)
    call mpp_update_domains(mfya,domain)
    ebuffer = mfx(ie,js:je,:)
    nbuffer = mfy(is:ie,je,:)
!============ GEOS-LAM Update End =============
    call mpp_get_boundary(xL, yL, domain, &
                          wbufferx=wbuffer, ebufferx=ebuffer, &
                          sbuffery=sbuffer, nbuffery=nbuffer, &
                          gridtype=CGRID_NE )
    xL(ie+1,js:je,:) = ebuffer
    yL(is:ie,je+1,:) = nbuffer
    mfxL(is:ie+1,js:je,:) = xL(is:ie+1,js:je,:)
    mfyL(is:ie,js:je+1,:) = yL(is:ie,js:je+1,:)
!============ GEOS-LAM Update Start =============
! Exercise the historical GCHP buffer probe; this is not the selected persistent-cell implementation.
    mfxL(ie+1,js:je,:)=mfxa(ie+1,js:je,:)
    mfyL(is:ie,je+1,:)=mfya(is:ie,je+1,:)
    call probe_array('gchp_mfx',mfxL,is,js)
    call probe_array('gchp_mfy',mfyL,is,js)
    call probe_array('gchp_cx',cxL,is,jsd)
    call probe_array('gchp_cy',cyL,isd,js)
    if (mpp_pe()==0) write(*,'(A,1X,ES24.16)')'GBUF_DRY_FACE', &
         (sum(mfxL(is,js:je,:))-sum(mfxL(ie+1,js:je,:)) + &
          sum(mfyL(is:ie,js,:))-sum(mfyL(is:ie,je+1,:)))/9.80665d0
    if (present(native_mfx) .neqv. present(native_mfy)) &
         error stop 'LAM exterior mass-flux pair is incomplete'
    if (present(native_mfx) .and. present(native_mfy)) then
       ! The MAPL imports carry one value per active A-grid cell. On the
       ! physical east/north edges, their last value is an interior face;
       ! retain the FV native exterior C-grid face instead.
       if (ie==npx-1) then
          mfxL(ie+1,js:je,:)=native_mfx(ie+1,js:je,:)
       endif
       if (je==npy-1) then
          mfyL(is:ie,je+1,:)=native_mfy(is:ie,je+1,:)
       endif
    endif
!============ GEOS-LAM Update End =============

! Fill local tracers and pressure thickness
    dpL(bd%is:bd%ie,bd%js:bd%je,:) = pleB(:,:,2:npz+1) - pleB(:,:,1:npz)
    q3(is:ie,js:je,:,:) = q(is:ie,js:je,:,:)
!============ GEOS-LAM Update Start =============
! Record dry-pilot tracer-stage diagnostics; remove or isolate before production use.
    if (nq>1) call probe_array('gchp_q0',q(is:ie,js:je,:,2),is,js)
    call probe_array('gchp_dp0',dpL(is:ie,js:je,:),is,js)
    ! Diagnostic only: give unfilled physical-boundary halos defined values.
    do j=jsd,jed
       do i=isd,ied
          if (i < is .or. i > ie .or. j < js .or. j > je) then
             q3(i,j,:,:) = q3(max(is,min(ie,i)),max(js,min(je,j)),:,:)
             dpL(i,j,:) = dpL(max(is,min(ie,i)),max(js,min(je,j)),:)
          endif
       enddo
    enddo
    call mpp_update_domains(dpL,domain)
    call start_group_halo_update(i_pack, q3, domain)
    print *, 'LAM_OFFLINE_TRACER_HALO_STARTED'
    flush(6)
!============ GEOS-LAM Update End =============

    if ( flagstruct%z_tracer .and. q_split==0 ) then
         call tracer_2d_1L(q3, dpL, mfxL, mfyL, cxL, cyL, &
                         gridstruct, bd, domain, npx, npy, npz, nq,    &
                         flagstruct%hord_tr, dt, 0, i_pack, &
                         flagstruct%nord_tr, flagstruct%trdm2, flagstruct%lim_fac, dpA=dpA)
    else
         call tracer_2d(q3, dpL, mfxL, mfyL, cxL, cyL, gridstruct, bd, domain, npx, npy, npz, nq,    &
                        flagstruct%hord_tr, q_split, dt, 0, i_pack, &
                        flagstruct%nord_tr, flagstruct%trdm2, flagstruct%lim_fac, dpA=dpA)
    endif
!============ GEOS-LAM Update Start =============
! Exercise the historical GCHP buffer probe; this is not the selected persistent-cell implementation.
    print *, 'LAM_OFFLINE_TRACER_HORIZONTAL_DONE', minval(q3), maxval(q3)
    flush(6)
    if (nq>1) call probe_array('gchp_qh',q3(is:ie,js:je,:,2),is,js)
    call probe_array('gchp_dpa',dpA,is,js)
    call probe_array('gchp_dpd',pleA(:,:,2:npz+1)-pleA(:,:,1:npz),is,js)
    dp_inner_max=0.
    do k=1,npz
    do j=js,je
    do i=is,ie
       if (i>3.and.i<npx-3.and.j>3.and.j<npy-3) &
            dp_inner_max=max(dp_inner_max,abs(dpA(i,j,k)-(pleA(i,j,k+1)-pleA(i,j,k))))
    enddo
    enddo
    enddo
    if (mpp_pe()==0) write(*,'(A,2(1X,ES24.16))')'GBUF_PRESSURE', &
         dp_inner_max,maxval(abs(dpA-(pleA(:,:,2:npz+1)-pleA(:,:,1:npz))))
!============ GEOS-LAM Update End =============

!------------------------------------------------------------------
! Re-Map constituents
!------------------------------------------------------------------
       kord_tracers = flagstruct%kord_tr
!$OMP parallel do default(none) shared(is,ie,isd,ied,js,je,jsd,jed,npz,nq, &
!$OMP                                  pleB,dpA,pleA,q3,flagstruct,kord_tracers) &
!$OMP                          private(i,j,k,pe1,dp1,pe2,dp2)
       do j=js,je
        ! pressures mapping from (dpA is new delp after tracer_2d)
          pe1(:,1) = pleB(:,j,1)
          do k=2,npz+1
            pe1(:,k) = pe1(:,k-1) + dpA(:,j,k-1)
          enddo
          do k=1,npz
             dp1(:,k) = pe1(:,k+1) - pe1(:,k)
          enddo
        ! pressures mapping to
          pe2 = pleA(:,j,:)
          do k=1,npz
             dp2(:,k) = pe2(:,k+1) - pe2(:,k)
          enddo
          call mapn_tracer(nq, npz, pe1, pe2, q3, dp2, kord_tracers, j,     &
                           is, ie, isd, ied, jsd, jed, 0., flagstruct%fill)
       enddo
!============ GEOS-LAM Update Start =============
! Record dry-pilot tracer-stage diagnostics; remove or isolate before production use.
       print *, 'LAM_OFFLINE_TRACER_VERTICAL_DONE', minval(q3), maxval(q3)
       flush(6)
       if (nq>1) call probe_array('gchp_q1',q3(is:ie,js:je,:,2),is,js)
!============ GEOS-LAM Update End =============

       ! Rescale tracers based on pleA at destination timestep
       !------------------------------------------------------
       do iq=1,nq
            ! Scale tracers
            !---------------
            scalingFactor = calcScalingFactor(q(is:ie,js:je,1:npz,iq), q3(is:ie,js:je,1:npz,iq), pleB, pleA, &
                              npz, domain, gridstruct, flagstruct, bd)
!============ GEOS-LAM Update Start =============
! Record the rejected offline exterior-face experiment for review, not production adoption.
            print *, 'LAM_TRACER_GLOBAL_SCALE',iq,scalingFactor
            if(present(bounded_lam)) then
               if(bounded_lam) scalingFactor=1.0
            endif
!============ GEOS-LAM Update End =============
            q(is:ie,js:je,1:npz,iq) = q3(is:ie,js:je,1:npz,iq) * scalingFactor
       enddo

end subroutine offline_tracer_advection

!------------------------------------------------------------------------------------

         function calcScalingFactor(q1, q2, ple1, ple2, npz, domain, gridstruct, flagstruct, bd) result(scaling)
         integer, intent(in) :: npz
         type(fv_grid_bounds_type), intent(IN   ) :: bd
         real, intent(in) :: q1(bd%is:bd%ie,bd%js:bd%je,npz)
         real, intent(in) :: q2(bd%is:bd%ie,bd%js:bd%je,npz)
         real, intent(in) :: ple1(bd%is:bd%ie,bd%js:bd%je,npz+1)
         real, intent(in) :: ple2(bd%is:bd%ie,bd%js:bd%je,npz+1)
         type(domain2D), intent(INOUT) :: domain
         type(fv_grid_type), intent(IN   ) :: gridstruct
         type(fv_flags_type), intent(INOUT) :: flagstruct
         integer :: sflag
         real :: scaling

         integer :: k
         real(REAL8)     :: qsum1(bd%is:bd%ie,bd%js:bd%je)
         real(REAL8)     :: qsum2(bd%is:bd%ie,bd%js:bd%je)
         real(REAL8)     :: globalSums(2)
         real(REAL8), parameter :: TINY_DENOMINATOR = tiny(1.d0)
         real(REAL8)     :: scalingR8
         !-------
         ! Compute partial sum on local array first to minimize communication.
         ! This algorithm will not be strongly repdroducible under changes do domain
         ! decomposition, but uses far less communication bandwidth (and memory BW)
         ! then the preceding implementation.
         !-------
         qsum1(:,:) = 0.d0
         qsum2(:,:) = 0.d0
         do k=1,npz
            qsum1(:,:) = qsum1(:,:) + q1(:,:,k) * (ple1(:,:,k+1)-ple1(:,:,k))
            qsum2(:,:) = qsum2(:,:) + q2(:,:,k) * (ple2(:,:,k+1)-ple2(:,:,k))
         enddo

         ! numerator
         globalSums(1) = g_sum_r8(domain, qsum1, bd%is,bd%ie, bd%js,bd%je, 0, &
                                  gridstruct%area_64(bd%is:bd%ie,bd%js:bd%je), 1)
         ! denominator
         globalSums(2) = g_sum_r8(domain, qsum2, bd%is,bd%ie, bd%js,bd%je, 0, &
                                  gridstruct%area_64(bd%is:bd%ie,bd%js:bd%je), 1)
!============ GEOS-LAM Update Start =============
! Record dry-pilot tracer-stage diagnostics; remove or isolate before production use.
         print *, 'LAM_TRACER_SCALE_MASS',globalSums(1)/9.80665d0,globalSums(2)/9.80665d0
!============ GEOS-LAM Update End =============

         if (globalSums(2) > TINY_DENOMINATOR) then
            scalingR8 =  globalSums(1) / globalSums(2)
            !#################################################################
            ! This line was added to ensure strong reproducibility of the code
            !#################################################################
            scaling = REAL(scalingR8, KIND=kind(1.00))
         else
            scaling = 1.0
         end if

         end function calcScalingFactor

end module fv_tracer2d_mod
