! Isolated port. NASA SI v1 full halo snapshots, flat dry matched-eta scope.
! All clock/history belong to the grid's neststruct. No module SAVE buffers.
module fv_lam_provider_mod
 use fv_lam_aux_bridge_mod, only: aux_q, aux_lower, aux_upper, aux_ready, &
      met_ps,met_t,met_w,met_zh,met_uw,met_vw,met_us,met_vs, &
      met_lower,met_upper,met_ready,met_x0,met_y0
 use lam_probe_mod,only:probe_mark,probe_array,probe_context,probe_level
 use fv_arrays_mod, only: fv_atmos_type,fv_nest_BC_type_3D,allocate_fv_nest_BC_type,deallocate_fv_nest_BC_type
 use fv_grid_utils_mod, only: mid_pt_sphere,cell_center2,get_unit_vect2,get_latlon_vector,inner_prod, &
      latlon2xyz,vect_cross,normalize_vect
 use fv_lam_geometry_mod, only: lam_nc_check
 use constants_mod, only: rdgas,grav,kappa,cp_air
 use sw_core_mod, only: divergence_corner_nest
 use fv_restart_mod, only: d2c_setup
 use boundary_mod, only: nested_grid_BC_apply_intT
 use mpp_domains_mod, only: mpp_get_neighbor_pe,EAST,WEST,NORTH,SOUTH
 use netcdf
 use mpp_mod, only: mpp_error,FATAL,mpp_pe,mpp_transmit,NULL_PE
 use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
 implicit none
 private
 public::lam_prepare,lam_initial_state
contains
 subroutine read_snapshot(a,hour,pt,dp,dz,w,u,v,uc,vc,divg,initial,tracers,met_idx,tracer_only)
 type(fv_atmos_type),intent(inout)::a
 integer,intent(in)::hour
 logical,optional,intent(in)::initial
 integer,optional,intent(in)::met_idx
 logical,optional,intent(in)::tracer_only
 real,optional,intent(out)::tracers(:,:,:,:)
 real,intent(out)::pt(a%bd%isd:a%bd%ied+1,a%bd%jsd:a%bd%jed+1,a%npz)
 real,intent(out)::dp(a%bd%isd:a%bd%ied+1,a%bd%jsd:a%bd%jed+1,a%npz),dz(a%bd%isd:a%bd%ied+1,a%bd%jsd:a%bd%jed+1,a%npz)
 real,intent(out)::w(a%bd%isd:a%bd%ied+1,a%bd%jsd:a%bd%jed+1,a%npz)
 real,intent(out)::u(a%bd%isd:a%bd%ied,a%bd%jsd:a%bd%jed+1,a%npz),vc(a%bd%isd:a%bd%ied,a%bd%jsd:a%bd%jed+1,a%npz)
 real,intent(out)::v(a%bd%isd:a%bd%ied+1,a%bd%jsd:a%bd%jed,a%npz),uc(a%bd%isd:a%bd%ied+1,a%bd%jsd:a%bd%jed,a%npz)
 real,intent(out)::divg(a%bd%isd:a%bd%ied+1,a%bd%jsd:a%bd%jed+1,a%npz)
 real::ps(a%bd%isd:a%bd%ied+1,a%bd%jsd:a%bd%jed+1),zh(a%bd%isd:a%bd%ied+1,a%bd%jsd:a%bd%jed+1,a%npz+1)
 real::us(a%bd%isd:a%bd%ied+1,a%bd%jsd:a%bd%jed+1,a%npz),vs(a%bd%isd:a%bd%ied+1,a%bd%jsd:a%bd%jed+1,a%npz)
 real::uw(a%bd%isd:a%bd%ied+1,a%bd%jsd:a%bd%jed+1,a%npz),vw(a%bd%isd:a%bd%ied+1,a%bd%jsd:a%bd%jed+1,a%npz)
 real::ua(a%bd%isd:a%bd%ied,a%bd%jsd:a%bd%jed),va(a%bd%isd:a%bd%ied,a%bd%jsd:a%bd%jed)
 real::cuc(a%bd%isd:a%bd%ied+1,a%bd%jsd:a%bd%jed),cvc(a%bd%isd:a%bd%ied,a%bd%jsd:a%bd%jed+1)
 real::ak(a%npz+1),bk(a%npz+1),p0,p1,stamp
 real::moist(a%bd%isd:a%bd%ied+1,a%bd%jsd:a%bd%jed+1,a%npz)
 real(kind=8)::p(2),q(2),e1(3),e2(3),ex(3),ey(3),radial(3)
 real(kind=8) :: cvce(a%bd%isd:a%bd%ied,a%bd%jsd:a%bd%jed+1),cvcn(a%bd%isd:a%bd%ied,a%bd%jsd:a%bd%jed+1)
 real(kind=8) :: cuce(a%bd%isd:a%bd%ied+1,a%bd%jsd:a%bd%jed),cucn(a%bd%isd:a%bd%ied+1,a%bd%jsd:a%bd%jed)
 real::native_x(a%bd%isd:a%bd%ied+1,a%bd%jsd:a%bd%jed,3)
 real::native_y(a%bd%isd:a%bd%ied,a%bd%jsd:a%bd%jed+1,3)
 real::send_w(3,a%bd%jsd:a%bd%jed),send_e(3,a%bd%jsd:a%bd%jed)
 real::recv_w(3,a%bd%jsd:a%bd%jed),recv_e(3,a%bd%jsd:a%bd%jed)
 real::send_s(3,a%bd%isd:a%bd%ied),send_n(3,a%bd%isd:a%bd%ied)
 real::recv_s(3,a%bd%isd:a%bd%ied),recv_n(3,a%bd%isd:a%bd%ied)
 integer::nc,vid,i,j,k,n,isd,ied,jsd,jed,nd,dimids(3),len1,len2,len3
 integer::west_pe,east_pe,south_pe,north_pe,want_w,want_e,want_s,want_n
 character(256)::file
 character(64)::schema,contract
 character(32)::units
 character(4)::tracer_name
 logical::met_source_only
 isd=a%bd%isd;ied=a%bd%ied;jsd=a%bd%jsd;jed=a%bd%jed
 met_source_only=a%neststruct%lam_met_file_split.and.present(met_idx)
 if(met_source_only) then
 ak=a%ak;bk=a%bk
 else
 write(file,'(A,"/lam.",I3.3,".nc")')trim(a%neststruct%lam_directory),hour
 if(present(initial)) then
 if(initial) file=trim(a%neststruct%lam_directory)//'/lam.initial.nc'
 endif
 call lam_nc_check(nf90_open(trim(file),nf90_nowrite,nc))
 call lam_nc_check(nf90_get_att(nc,nf90_global,'schema',schema))
 if(trim(schema)/='NASA_LAM_SI_V1') call mpp_error(FATAL,'LAM unsupported snapshot schema')
 call lam_nc_check(nf90_get_att(nc,nf90_global,'halo',contract))
 if(trim(contract)/='4') call mpp_error(FATAL,'LAM requires full input halo 4')
 call lam_nc_check(nf90_get_att(nc,nf90_global,'vertical_order',contract))
 if(trim(contract)/='top_to_bottom') call mpp_error(FATAL,'LAM vertical order')
 call lam_nc_check(nf90_get_att(nc,nf90_global,'terrain',contract))
 if(trim(contract)/='flat') call mpp_error(FATAL,'LAM first unit requires flat terrain')
 call lam_nc_check(nf90_get_att(nc,nf90_global,'moisture',contract))
 if(trim(contract)/='dry') call mpp_error(FATAL,'LAM first unit requires dry input')
 call lam_nc_check(nf90_inq_varid(nc,'time_seconds',vid));call lam_nc_check(nf90_get_var(nc,vid,stamp))
 if(stamp/=hour*3600.) call mpp_error(FATAL,'LAM file timestamp does not match bracket')
 call lam_nc_check(nf90_inq_varid(nc,'ak',vid));call lam_nc_check(nf90_get_var(nc,vid,ak))
 call lam_nc_check(nf90_inq_varid(nc,'bk',vid));call lam_nc_check(nf90_get_var(nc,vid,bk))
 if(any(ak/=a%ak).or.any(bk/=a%bk)) call mpp_error(FATAL,'LAM matched-eta contract violated')
 if(present(tracer_only)) then
 if(tracer_only) then
 if(.not.present(tracers)) call mpp_error(FATAL,'LAM tracer-only read requires output')
 if(any(shape(tracers)/=[ied-isd+2,jed-jsd+2,a%npz,size(a%q,4)])) &
      call mpp_error(FATAL,'LAM tracer-only buffer shape')
 if(size(tracers,4)>3) call mpp_error(FATAL,'LAM dry named tracer schema')
 tracers=0.
 call get3('q',moist)
 if(any(moist/=0.)) call mpp_error(FATAL,'LAM dry unit rejects nonzero moisture')
 if(size(tracers,4)>0) tracers(:,:,:,1)=moist
 do n=2,size(tracers,4)
    write(tracer_name,'("q",I3.3)')n-1
    call get3(tracer_name,tracers(:,:,:,n))
    if(any(tracers(:,:,:,n)<0.)) call mpp_error(FATAL,'LAM negative tracer input '//tracer_name)
 enddo
 call lam_nc_check(nf90_close(nc))
 call probe_mark('TRACER_ONLY_SNAPSHOT',[real(hour),stamp])
 if(mpp_pe()==0) print *, 'LAM_TRACER_ONLY_SNAPSHOT',hour
 return
 endif
 endif
 call lam_nc_check(nf90_inq_varid(nc,'ps',vid))
 call shape_check('ps')
 call lam_nc_check(nf90_get_att(nc,vid,'units',units))
 if(trim(units)/='Pa') call mpp_error(FATAL,'LAM ps units')
 call lam_nc_check(nf90_get_var(nc,vid,ps,start=[isd+4,jsd+4]))
 call get3('t',pt);call get3('zh',zh);call get3('w',w);call get3('q',moist)
 if(present(tracers)) then
 if(any(shape(tracers)/=[ied-isd+2,jed-jsd+2,a%npz,size(a%q,4)])) &
      call mpp_error(FATAL,'LAM tracer snapshot buffer shape')
 if(size(tracers,4)>3) call mpp_error(FATAL,'LAM dry native tracer input supports Q,Q001,Q002 only')
 tracers=0.
 if(size(tracers,4)>0) tracers(:,:,:,1)=moist
 do n=2,size(tracers,4)
    write(tracer_name,'("q",I3.3)')n-1
    call get3(tracer_name,tracers(:,:,:,n))
    if(any(tracers(:,:,:,n)<0.)) call mpp_error(FATAL,'LAM negative tracer input '//tracer_name)
 enddo
 endif
 call probe_mark('SNAPSHOT_READ',[real(hour),stamp])
 call probe_array('SNAPSHOT_ZH',zh,isd,jsd)
 if(any(moist/=0.)) call mpp_error(FATAL,'LAM dry unit rejects nonzero moisture')
 if(any(abs(zh(:,:,a%npz+1))>1.e-10)) call mpp_error(FATAL,'LAM flat terrain surface height')
 call get3('u_s',us);call get3('v_s',vs);call get3('u_w',uw);call get3('v_w',vw)
 call lam_nc_check(nf90_close(nc))
 endif
 if(present(met_idx)) then
   if(.not.met_ready.or.met_idx<0.or.met_idx>1) call mpp_error(FATAL,'LAM MET bridge unavailable')
   n=met_idx+1
   if(isd+4<met_x0.or.ied+5>met_x0+min(size(met_ps,1),size(met_t,1), &
        size(met_w,1),size(met_zh,1),size(met_uw,1),size(met_vw,1))-1) &
        call mpp_error(FATAL,'LAM MET pressure local x window')
   if(isd+4<met_x0.or.ied+4>met_x0+min(size(met_us,1),size(met_vs,1))-1) &
        call mpp_error(FATAL,'LAM MET US_E local x window')
   if(jsd+4<met_y0.or.jed+5>met_y0+min(size(met_ps,2),size(met_t,2), &
        size(met_w,2),size(met_zh,2),size(met_us,2),size(met_vs,2))-1) &
        call mpp_error(FATAL,'LAM MET pressure local y window')
   if(jsd+4<met_y0.or.jed+4>met_y0+min(size(met_uw,2),size(met_vw,2))-1) &
        call mpp_error(FATAL,'LAM MET UW_E local y window')
   ps=met_ps(isd+5-met_x0:ied+6-met_x0,jsd+5-met_y0:jed+6-met_y0,1,n)
   pt=met_t(isd+5-met_x0:ied+6-met_x0,jsd+5-met_y0:jed+6-met_y0,:,n)
   zh=met_zh(isd+5-met_x0:ied+6-met_x0,jsd+5-met_y0:jed+6-met_y0,:,n)
   w=met_w(isd+5-met_x0:ied+6-met_x0,jsd+5-met_y0:jed+6-met_y0,:,n)
   us(isd:ied,jsd:jed+1,:)=met_us(isd+5-met_x0:ied+5-met_x0,jsd+5-met_y0:jed+6-met_y0,:,n)
   vs(isd:ied,jsd:jed+1,:)=met_vs(isd+5-met_x0:ied+5-met_x0,jsd+5-met_y0:jed+6-met_y0,:,n)
   uw(isd:ied+1,jsd:jed,:)=met_uw(isd+5-met_x0:ied+6-met_x0,jsd+5-met_y0:jed+5-met_y0,:,n)
   vw(isd:ied+1,jsd:jed,:)=met_vw(isd+5-met_x0:ied+6-met_x0,jsd+5-met_y0:jed+5-met_y0,:,n)
   call probe_mark('MET_BRIDGE_SNAPSHOT',[real(met_idx),real(met_lower),real(met_upper)])
 endif
 if(.not.all(ieee_is_finite(ps)).or.any(ps<=0.).or.any(pt<=0.).or.any(zh(:,:,1:a%npz)<=zh(:,:,2:a%npz+1))) call mpp_error(FATAL,'LAM nonpositive pressure/T/thickness')
 do k=1,a%npz
 do j=jsd,jed+1
 do i=isd,ied+1
 p0=ak(k)+bk(k)*ps(i,j);p1=ak(k+1)+bk(k+1)*ps(i,j)
 if(p0<=0..or.p1<=p0) call mpp_error(FATAL,'LAM invalid pressure interfaces')
 dp(i,j,k)=p1-p0;dz(i,j,k)=zh(i,j,k+1)-zh(i,j,k)
 enddo
 enddo
 enddo
 ! Request only the missing local endpoint normals from their owning MPI
 ! neighbor. The FMS C-grid halo update excludes these outermost endpoints.
 ! The exterior physical halo has no owner and retains the local fallback.
 if(present(met_idx)) then
   native_x=0.;native_y=0.
   send_w=0.;send_e=0.;send_s=0.;send_n=0.
   do j=jsd,jed
     do i=isd+1,ied
       native_x(i,j,:)=a%gridstruct%ew(:,i,j,1)
     enddo
   enddo
   do j=jsd+1,jed
     do i=isd,ied
       native_y(i,j,:)=a%gridstruct%es(:,i,j,2)
     enddo
   enddo
   call mpp_get_neighbor_pe(a%domain,WEST,west_pe)
   call mpp_get_neighbor_pe(a%domain,EAST,east_pe)
   call mpp_get_neighbor_pe(a%domain,SOUTH,south_pe)
   call mpp_get_neighbor_pe(a%domain,NORTH,north_pe)
   call mpp_transmit(ied+1,east_pe,want_w,west_pe)
   call mpp_transmit(isd,west_pe,want_e,east_pe)
   if(west_pe/=NULL_PE) then
     if(want_w<isd+1.or.want_w>ied) call mpp_error(FATAL,'LAM west metric owner bounds')
     do j=jsd,jed
       send_w(:,j)=native_x(want_w,j,:)
     enddo
   endif
   if(east_pe/=NULL_PE) then
     if(want_e<isd+1.or.want_e>ied) call mpp_error(FATAL,'LAM east metric owner bounds')
     do j=jsd,jed
       send_e(:,j)=native_x(want_e,j,:)
     enddo
   endif
   call mpp_transmit(send_w,size(send_w),west_pe,recv_e,size(recv_e),east_pe)
   call mpp_transmit(send_e,size(send_e),east_pe,recv_w,size(recv_w),west_pe)
   if(east_pe/=NULL_PE) then
     do j=jsd,jed
       native_x(ied+1,j,:)=recv_e(:,j)
     enddo
   endif
   if(west_pe/=NULL_PE) then
     do j=jsd,jed
       native_x(isd,j,:)=recv_w(:,j)
     enddo
   endif
   call mpp_transmit(jed+1,north_pe,want_s,south_pe)
   call mpp_transmit(jsd,south_pe,want_n,north_pe)
   if(south_pe/=NULL_PE) then
     if(want_s<jsd+1.or.want_s>jed) call mpp_error(FATAL,'LAM south metric owner bounds')
     do i=isd,ied
       send_s(:,i)=native_y(i,want_s,:)
     enddo
   endif
   if(north_pe/=NULL_PE) then
     if(want_n<jsd+1.or.want_n>jed) call mpp_error(FATAL,'LAM north metric owner bounds')
     do i=isd,ied
       send_n(:,i)=native_y(i,want_n,:)
     enddo
   endif
   call mpp_transmit(send_s,size(send_s),south_pe,recv_n,size(recv_n),north_pe)
   call mpp_transmit(send_n,size(send_n),north_pe,recv_s,size(recv_s),south_pe)
   if(north_pe/=NULL_PE) then
     do i=isd,ied
       native_y(i,jed+1,:)=recv_n(:,i)
     enddo
   endif
   if(south_pe/=NULL_PE) then
     do i=isd,ied
       native_y(i,jsd,:)=recv_s(:,i)
     enddo
   endif
 endif
 ! Covariant D and orthogonal C bases at the actual staggered points.
 do j=jsd,jed+1
 do i=isd,ied
 call mid_pt_sphere(a%gridstruct%grid_64(i,j,:),a%gridstruct%grid_64(i+1,j,:),p)
 call get_unit_vect2(a%gridstruct%grid_64(i,j,:),a%gridstruct%grid_64(i+1,j,:),e1)
 call get_latlon_vector(p,ex,ey)
 ! es(:,:,:,2) is the NASA C-grid normal at this edge.
 if(j>jsd.and.j<=jed) then
    e2=a%gridstruct%es(:,i,j,2)
 else if(present(met_idx).and.j>=1.and.j<=a%npy-1) then
    e2=native_y(i,j,:)
    if(abs(inner_prod(e2,e2)-1.)>1.e-10) call mpp_error(FATAL,'LAM owner y metric halo')
 else
    ! The outer local halo lacks a stored es vector. Define its C normal
    ! from the available edge tangent and midpoint radial direction.
    call latlon2xyz(p,radial)
    call vect_cross(e2,radial,e1)
    call normalize_vect(e2)
 endif
 cvce(i,j)=inner_prod(e2,ex);cvcn(i,j)=inner_prod(e2,ey)
 do k=1,a%npz
 u(i,j,k)=us(i,j,k)*inner_prod(e1,ex)+vs(i,j,k)*inner_prod(e1,ey)
 vc(i,j,k)=us(i,j,k)*inner_prod(e2,ex)+vs(i,j,k)*inner_prod(e2,ey)
 enddo
 enddo
 enddo
 do j=jsd,jed
 do i=isd,ied+1
 call mid_pt_sphere(a%gridstruct%grid_64(i,j,:),a%gridstruct%grid_64(i,j+1,:),p)
 call get_unit_vect2(a%gridstruct%grid_64(i,j,:),a%gridstruct%grid_64(i,j+1,:),e2)
 call get_latlon_vector(p,ex,ey)
 if(i>isd.and.i<=ied) then
    e1=a%gridstruct%ew(:,i,j,1)
 else if(present(met_idx).and.i>=1.and.i<=a%npx-1) then
    e1=native_x(i,j,:)
    if(abs(inner_prod(e1,e1)-1.)>1.e-10) call mpp_error(FATAL,'LAM owner x metric halo')
 else
    ! The outer local halo lacks a stored ew vector. The provider must
    ! supply a finite C normal there for boundary pressure transport.
    call latlon2xyz(p,radial)
    call vect_cross(e1,e2,radial)
    call normalize_vect(e1)
 endif
 cuce(i,j)=inner_prod(e1,ex);cucn(i,j)=inner_prod(e1,ey)
 do k=1,a%npz
 v(i,j,k)=uw(i,j,k)*inner_prod(e2,ex)+vw(i,j,k)*inner_prod(e2,ey)
 uc(i,j,k)=uw(i,j,k)*inner_prod(e1,ex)+vw(i,j,k)*inner_prod(e1,ey)
 enddo
 enddo
 enddo
 if (present(initial)) then
 if (initial) then
 block
 integer :: dump_unit
 character(64) :: dump_name
 write(dump_name,'("direct-vc-pre-d2c-r",I1,".bin")') mpp_pe()
 open(newunit=dump_unit,file=trim(dump_name),status='new',access='stream',form='unformatted')
 write(dump_unit) vc(a%bd%is:a%bd%ie,a%bd%js:a%bd%je+1,:)
 close(dump_unit)
 write(dump_name,'("direct-uc-pre-d2c-r",I1,".bin")') mpp_pe()
 open(newunit=dump_unit,file=trim(dump_name),status='new',access='stream',form='unformatted')
 write(dump_unit) uc(a%bd%is:a%bd%ie+1,a%bd%js:a%bd%je,:)
 close(dump_unit)
 write(dump_name,'("direct-vc-coeff-e-r",I1,".bin")') mpp_pe()
 open(newunit=dump_unit,file=trim(dump_name),status='new',access='stream',form='unformatted')
 write(dump_unit) cvce(a%bd%is:a%bd%ie,a%bd%js:a%bd%je+1)
 close(dump_unit)
 write(dump_name,'("direct-vc-coeff-n-r",I1,".bin")') mpp_pe()
 open(newunit=dump_unit,file=trim(dump_name),status='new',access='stream',form='unformatted')
 write(dump_unit) cvcn(a%bd%is:a%bd%ie,a%bd%js:a%bd%je+1)
 close(dump_unit)
 write(dump_name,'("direct-uc-coeff-e-r",I1,".bin")') mpp_pe()
 open(newunit=dump_unit,file=trim(dump_name),status='new',access='stream',form='unformatted')
 write(dump_unit) cuce(a%bd%is:a%bd%ie+1,a%bd%js:a%bd%je)
 close(dump_unit)
 write(dump_name,'("direct-uc-coeff-n-r",I1,".bin")') mpp_pe()
 open(newunit=dump_unit,file=trim(dump_name),status='new',access='stream',form='unformatted')
 write(dump_unit) cucn(a%bd%is:a%bd%ie+1,a%bd%js:a%bd%je)
 close(dump_unit)
 end block
 endif
 endif
 divg=0.
 do k=1,a%npz
 call d2c_setup(u(isd,jsd,k),v(isd,jsd,k),ua,va,cuc,cvc,a%flagstruct%nord>0, &
 isd,ied,jsd,jed,a%bd%is,a%bd%ie,a%bd%js,a%bd%je,a%npx,a%npy,a%gridstruct%grid_type,.true., &
 .false.,.false.,.false.,.false.,a%gridstruct%rsin_u,a%gridstruct%rsin_v,a%gridstruct%cosa_s,a%gridstruct%rsin2)
 uc(a%bd%is-1:a%bd%ie+2,a%bd%js-1:a%bd%je+1,k)=cuc(a%bd%is-1:a%bd%ie+2,a%bd%js-1:a%bd%je+1)
 vc(a%bd%is-1:a%bd%ie+1,a%bd%js-1:a%bd%je+2,k)=cvc(a%bd%is-1:a%bd%ie+1,a%bd%js-1:a%bd%je+2)
 call divergence_corner_nest(u(isd,jsd,k),v(isd,jsd,k),ua,va,divg(isd,jsd,k),a%gridstruct,a%flagstruct,a%bd)
 enddo
 ! Isolated one-time diagnostic of the direct-file provider's native
 ! projected D-grid components. Does not change the snapshot or dynamics.
 if (present(initial)) then
 if (initial) then
 block
 integer :: dump_unit
 character(64) :: dump_name
 write(dump_name,'("direct-u-yface-r",I1,".bin")') mpp_pe()
 open(newunit=dump_unit,file=trim(dump_name),status='new',access='stream',form='unformatted')
 write(dump_unit) u(a%bd%is:a%bd%ie,a%bd%js:a%bd%je+1,:)
 close(dump_unit)
 write(dump_name,'("direct-v-xface-r",I1,".bin")') mpp_pe()
 open(newunit=dump_unit,file=trim(dump_name),status='new',access='stream',form='unformatted')
 write(dump_unit) v(a%bd%is:a%bd%ie+1,a%bd%js:a%bd%je,:)
 close(dump_unit)
 print *, 'MET3_DIRECT_NATIVE',mpp_pe(),a%bd%is,a%bd%ie,a%bd%js,a%bd%je
 end block
 endif
 endif
contains
 subroutine get3(name,field)
 character(*),intent(in)::name
 real,intent(out)::field(:,:,:)
 character(32)::units
 call lam_nc_check(nf90_inq_varid(nc,name,vid))
 call shape_check(name)
 call lam_nc_check(nf90_get_att(nc,vid,'units',units))
 select case(name)
 case('t');if(trim(units)/='K') call mpp_error(FATAL,'LAM t units')
 case('q','q001','q002');if(trim(units)/='kg kg-1') call mpp_error(FATAL,'LAM q units')
 case('zh');if(trim(units)/='m') call mpp_error(FATAL,'LAM zh units')
 case default; if(trim(units)/='m s-1') call mpp_error(FATAL,'LAM wind units')
 end select
 call lam_nc_check(nf90_get_var(nc,vid,field,start=[isd+4,jsd+4,1]))
 if(.not.all(ieee_is_finite(field))) call mpp_error(FATAL,'LAM nonfinite file field '//name)
 end subroutine
 subroutine shape_check(name)
 character(*),intent(in)::name
 integer::rank,ids(3),lens(3),d,nx,ny,nz
 character(20)::dn,expected(3)
 call lam_nc_check(nf90_inquire_variable(nc,vid,ndims=rank,dimids=ids))
 nx=a%npx-1+8;ny=a%npy-1+8;nz=a%npz
 expected=['x                   ','y                   ','lev                 ']
 select case(name)
 case('u_s','v_s');ny=ny+1;expected(2)='ys'
 case('u_w','v_w');nx=nx+1;expected(1)='xs'
 case('zh');nz=nz+1;expected(3)='edge'
 end select
 if(rank/=merge(2,3,name=='ps'))call mpp_error(FATAL,'LAM snapshot rank '//name)
 do d=1,rank
 call lam_nc_check(nf90_inquire_dimension(nc,ids(d),name=dn,len=lens(d)))
 if(trim(dn)/=trim(expected(d)))call mpp_error(FATAL,'LAM snapshot dimension order '//name)
 enddo
 ! Explicit checks avoid relying on truncation of a larger malformed grid.
 if(lens(1)/=nx.or.lens(2)/=ny)call mpp_error(FATAL,'LAM snapshot horizontal shape '//name)
 if(rank==3)then
 if(lens(3)/=nz)call mpp_error(FATAL,'LAM snapshot vertical shape '//name)
 endif
 end subroutine
 end subroutine read_snapshot

 subroutine lam_prepare(a,time,dt)
 type(fv_atmos_type),intent(inout)::a
 real,intent(in)::time,dt
 integer::lower,upper,idx,i,j,k,n,isd,ied,jsd,jed
 real::pt(a%bd%isd:a%bd%ied+1,a%bd%jsd:a%bd%jed+1,a%npz),dp(a%bd%isd:a%bd%ied+1,a%bd%jsd:a%bd%jed+1,a%npz)
 real::dz(a%bd%isd:a%bd%ied+1,a%bd%jsd:a%bd%jed+1,a%npz),w(a%bd%isd:a%bd%ied+1,a%bd%jsd:a%bd%jed+1,a%npz)
 real::u(a%bd%isd:a%bd%ied,a%bd%jsd:a%bd%jed+1,a%npz),vc(a%bd%isd:a%bd%ied,a%bd%jsd:a%bd%jed+1,a%npz)
 real::v(a%bd%isd:a%bd%ied+1,a%bd%jsd:a%bd%jed,a%npz),uc(a%bd%isd:a%bd%ied+1,a%bd%jsd:a%bd%jed,a%npz)
 real::divg(a%bd%isd:a%bd%ied+1,a%bd%jsd:a%bd%jed+1,a%npz),pkz,alpha,p0,p1
 real::q_snapshot(a%bd%isd:a%bd%ied+1,a%bd%jsd:a%bd%jed+1,a%npz,size(a%q,4))
 if(allocated(a%neststruct%q_BC)) then
 if(size(a%neststruct%q_BC)/=size(a%q,4))then
 do n=1,size(a%neststruct%q_BC)
 call deallocate_fv_nest_BC_type(a%neststruct%q_BC(n))
 enddo
 deallocate(a%neststruct%q_BC)
 endif
 endif
 if(.not.allocated(a%neststruct%q_BC)) then
 allocate(a%neststruct%q_BC(size(a%q,4)))
 do n=1,size(a%q,4)
 call allocate_fv_nest_BC_type(a%neststruct%q_BC(n),a,1,0,0,.false.)
 enddo
 ! A cold host can establish the time bracket before it registers Q.
 ! Re-read the bracket when the tracer schema appears or changes.
 a%neststruct%lam_lower=-1.
 a%neststruct%lam_upper=-1.
 endif
 if(time<0..or.dt<0.) call mpp_error(FATAL,'LAM invalid elapsed time/dt')
 call probe_mark('PROVIDER_PREPARE',[time,dt,real(a%neststruct%lam_interval)])
 lower=int(floor(time/(3600.*a%neststruct%lam_interval)))*a%neststruct%lam_interval
 upper=lower+a%neststruct%lam_interval
 if(time+dt>upper*3600.+1.e-9) call mpp_error(FATAL,'LAM step crosses snapshot boundary; timestep must divide interval')
 isd=a%bd%isd;ied=a%bd%ied;jsd=a%bd%jsd;jed=a%bd%jed
 if(a%neststruct%lam_lower/=lower*3600..or.a%neststruct%lam_upper/=upper*3600.) then
 a%neststruct%lam_lower=lower*3600.;a%neststruct%lam_upper=upper*3600.
 do idx=0,1
 if(a%neststruct%lam_met_file_split) then
   call read_snapshot(a,merge(lower,upper,idx==0),pt,dp,dz,w,u,v,uc,vc,divg, &
        tracers=q_snapshot,tracer_only=.true.)
 else
   call read_snapshot(a,merge(lower,upper,idx==0),pt,dp,dz,w,u,v,uc,vc,divg,tracers=q_snapshot)
 endif
 if(.not.a%neststruct%lam_met_file_split) then
 do k=1,a%npz
 do j=jsd,jed+1
 do i=isd,ied+1
 if(a%flagstruct%hydrostatic) then
 p0=a%ptop+sum(dp(i,j,1:k-1));p1=p0+dp(i,j,k)
 pkz=(p1**kappa-p0**kappa)/(kappa*log(p1/p0))
 else
 pkz=exp(kappa*log((-rdgas/grav)*dp(i,j,k)*pt(i,j,k)/dz(i,j,k)))
 endif
 pt(i,j,k)=pt(i,j,k)/pkz
 enddo
 enddo
 enddo
 call store(a%neststruct%pt_BC,pt,0,0)
 call store(a%neststruct%delp_BC,dp,0,0)
 call store(a%neststruct%u_BC,u,0,1);call store(a%neststruct%vc_BC,vc,0,1)
 call store(a%neststruct%v_BC,v,1,0);call store(a%neststruct%uc_BC,uc,1,0)
 call store(a%neststruct%divg_BC,divg,1,1)
 if(.not.a%flagstruct%hydrostatic) then
 call store(a%neststruct%delz_BC,dz,0,0);call store(a%neststruct%w_BC,w,0,0)
 endif
 endif
 pt=0.;call store(a%neststruct%q_con_BC,pt,0,0)
 pt=kappa;call store(a%neststruct%cappa_BC,pt,0,0)
 do n=1,size(a%neststruct%q_BC)
    call store(a%neststruct%q_BC(n),q_snapshot(:,:,:,n),0,0)
 enddo
 enddo
 if(mpp_pe()==0) write(*,'(A,3F14.3)')'LAM_BRACKET ',time,a%neststruct%lam_lower,a%neststruct%lam_upper
 endif
 if(a%neststruct%lam_met_file_split.and..not.met_ready) then
   if(dt>0.) call mpp_error(FATAL,'LAM MET imports required before first DYN advance')
   call set_clock(a%neststruct%q_con_BC)
   call set_clock(a%neststruct%cappa_BC)
   do n=1,size(a%neststruct%q_BC)
      call set_clock(a%neststruct%q_BC(n))
      call apply(a%q(:,:,:,n),0,0,a%neststruct%q_BC(n))
   enddo
   a%neststruct%nest_timestep=0
   a%neststruct%tracer_nest_timestep=0
   if(mpp_pe()==0) print *, 'LAM_MET_DEFERRED_INIT',time
   return
 endif
 ! A host heartbeat supplies the two meteorological endpoints for this DYN
 ! call. Keep the file bracket above for named tracer boundaries only.
 if(met_ready) then
   if(time<met_lower-1.e-8.or.time+dt>met_upper+1.e-8) &
        call mpp_error(FATAL,'LAM MET stage outside host samples')
   do idx=0,1
     call read_snapshot(a,merge(lower,upper,idx==0),pt,dp,dz,w,u,v,uc,vc,divg,met_idx=idx)
     do k=1,a%npz
     do j=jsd,jed+1
     do i=isd,ied+1
       if(a%flagstruct%hydrostatic) then
         p0=a%ptop+sum(dp(i,j,1:k-1));p1=p0+dp(i,j,k)
         pkz=(p1**kappa-p0**kappa)/(kappa*log(p1/p0))
       else
         pkz=exp(kappa*log((-rdgas/grav)*dp(i,j,k)*pt(i,j,k)/dz(i,j,k)))
       endif
       pt(i,j,k)=pt(i,j,k)/pkz
     enddo
     enddo
     enddo
     call store(a%neststruct%pt_BC,pt,0,0)
     call store(a%neststruct%delp_BC,dp,0,0)
     call store(a%neststruct%u_BC,u,0,1);call store(a%neststruct%vc_BC,vc,0,1)
     call store(a%neststruct%v_BC,v,1,0);call store(a%neststruct%uc_BC,uc,1,0)
     call store(a%neststruct%divg_BC,divg,1,1)
     if(.not.a%flagstruct%hydrostatic) then
       call store(a%neststruct%delz_BC,dz,0,0)
       call store(a%neststruct%w_BC,w,0,0)
     endif
   enddo
   if(mpp_pe()==0) write(*,'(A,3F14.3)')'LAM_MET_BRACKET ',time,met_lower,met_upper
 endif
 ! Isolated MAPL two-grid probe: host-time brackets overwrite only tracer
 ! exterior side arrays. Physical tracer cells and meteorological files persist.
 if(aux_ready) then
   if(size(a%q,4)<3) error stop 'LAM AUX named tracer schema'
   if(time<aux_lower-1.e-8 .or. time+dt>aux_upper+1.e-8) &
        error stop 'LAM AUX stage outside host samples'
   do idx=0,1
     do n=2,3
       do k=1,a%npz
         q_snapshot(:,:,k,n)=aux_q(isd+4:ied+5,jsd+4:jed+5,n-1,idx+1)
       enddo
       call store(a%neststruct%q_BC(n),q_snapshot(:,:,:,n),0,0)
     enddo
   enddo
 endif
 ! Restore each invocation's exterior data in the representation expected before dynamics.
 call set_clock(a%neststruct%pt_BC);call set_clock(a%neststruct%delp_BC)
 call set_clock(a%neststruct%u_BC);call set_clock(a%neststruct%v_BC)
 call set_clock(a%neststruct%uc_BC);call set_clock(a%neststruct%vc_BC);call set_clock(a%neststruct%divg_BC)
 if(met_ready) then
   call set_met_clock(a%neststruct%pt_BC);call set_met_clock(a%neststruct%delp_BC)
   call set_met_clock(a%neststruct%u_BC);call set_met_clock(a%neststruct%v_BC)
   call set_met_clock(a%neststruct%uc_BC);call set_met_clock(a%neststruct%vc_BC)
   call set_met_clock(a%neststruct%divg_BC)
 endif
 call set_clock(a%neststruct%q_con_BC);call set_clock(a%neststruct%cappa_BC)
 do n=1,size(a%neststruct%q_BC)
   call set_clock(a%neststruct%q_BC(n))
   if(aux_ready.and.n>=2) then
      a%neststruct%q_BC(n)%lower_seconds=aux_lower
      a%neststruct%q_BC(n)%upper_seconds=aux_upper
   endif
 enddo
 if(.not.a%flagstruct%hydrostatic) then
 call set_clock(a%neststruct%w_BC);call set_clock(a%neststruct%delz_BC)
 if(met_ready) then
   call set_met_clock(a%neststruct%w_BC);call set_met_clock(a%neststruct%delz_BC)
 endif
 call probe_mark('PROVIDER_APPLY_w');call apply(a%w,0,0,a%neststruct%w_BC);call probe_mark('PROVIDER_APPLY_delz');call apply(a%delz,0,0,a%neststruct%delz_BC)
 endif
 call probe_mark('PROVIDER_APPLY_delp');call apply(a%delp,0,0,a%neststruct%delp_BC)
 call probe_mark('PROVIDER_APPLY_pt');call apply(a%pt,0,0,a%neststruct%pt_BC)
 call probe_mark('PROVIDER_APPLY_u');call apply(a%u,0,1,a%neststruct%u_BC);call probe_mark('PROVIDER_APPLY_v');call apply(a%v,1,0,a%neststruct%v_BC)
 do k=1,a%npz
 do j=jsd,jed
 do i=isd,ied
 if(i<1.or.i>=a%npx.or.j<1.or.j>=a%npy) then
 ! NH native theta -> Kelvin inversion, using the same NASA EOS.
 if(a%flagstruct%hydrostatic) then
 p0=a%ptop+sum(a%delp(i,j,1:k-1));p1=p0+a%delp(i,j,k)
 pkz=(p1**kappa-p0**kappa)/(kappa*log(p1/p0))
 a%pt(i,j,k)=a%pt(i,j,k)*pkz
 else
 a%pt(i,j,k)=exp((log(a%pt(i,j,k))+kappa*log((-rdgas/grav)*a%delp(i,j,k)/a%delz(i,j,k)))/(1.-kappa))
 endif
 endif
 enddo
 enddo
 enddo
 ! Time-matched tracer data occupy exterior halos only; physical cells persist.
 do n=1,size(a%neststruct%q_BC)
    call apply(a%q(:,:,:,n),0,0,a%neststruct%q_BC(n))
 enddo
 if(aux_ready.and.mpp_pe()==0) write(*,'(A,4F16.9)') 'LAM_AUX_STAGE ',time,dt, &
      a%q(isd,jsd,1,2),a%q(ied,jed,1,2)
 a%neststruct%nest_timestep=0;a%neststruct%tracer_nest_timestep=0
contains
 subroutine store(b,f,si,sj)
 type(fv_nest_BC_type_3D),intent(inout)::b
 integer,intent(in)::si,sj
 real,intent(in)::f(isd:,jsd:,:)
 if(idx==0) then
 if(a%bd%is==1) b%west_t0=f(lbound(b%west_t0,1):ubound(b%west_t0,1),lbound(b%west_t0,2):ubound(b%west_t0,2),:)
 if(a%bd%ie==a%npx-1) b%east_t0=f(lbound(b%east_t0,1):ubound(b%east_t0,1),lbound(b%east_t0,2):ubound(b%east_t0,2),:)
 if(a%bd%js==1) b%south_t0=f(lbound(b%south_t0,1):ubound(b%south_t0,1),lbound(b%south_t0,2):ubound(b%south_t0,2),:)
 if(a%bd%je==a%npy-1) b%north_t0=f(lbound(b%north_t0,1):ubound(b%north_t0,1),lbound(b%north_t0,2):ubound(b%north_t0,2),:)
 else
 if(a%bd%is==1) b%west_t1=f(lbound(b%west_t1,1):ubound(b%west_t1,1),lbound(b%west_t1,2):ubound(b%west_t1,2),:)
 if(a%bd%ie==a%npx-1) b%east_t1=f(lbound(b%east_t1,1):ubound(b%east_t1,1),lbound(b%east_t1,2):ubound(b%east_t1,2),:)
 if(a%bd%js==1) b%south_t1=f(lbound(b%south_t1,1):ubound(b%south_t1,1),lbound(b%south_t1,2):ubound(b%south_t1,2),:)
 if(a%bd%je==a%npy-1) b%north_t1=f(lbound(b%north_t1,1):ubound(b%north_t1,1),lbound(b%north_t1,2):ubound(b%north_t1,2),:)
 endif
 end subroutine
 subroutine set_clock(b)
 type(fv_nest_BC_type_3D),intent(inout)::b
 b%file_boundary=.true.;b%lower_seconds=a%neststruct%lam_lower;b%upper_seconds=a%neststruct%lam_upper
 b%stage_seconds=time;b%stage_dt=dt
 end subroutine
 subroutine set_met_clock(b)
 type(fv_nest_BC_type_3D),intent(inout)::b
 b%lower_seconds=met_lower;b%upper_seconds=met_upper
 end subroutine
 subroutine apply(f,si,sj,b)
 integer,intent(in)::si,sj
 real,intent(inout)::f(isd:,jsd:,:)
 type(fv_nest_BC_type_3D),intent(in)::b
 call nested_grid_BC_apply_intT(f,si,sj,a%npx,a%npy,a%npz,a%bd,0.,1.,b,1)
 end subroutine
 end subroutine lam_prepare

 subroutine lam_initial_state(a)
 type(fv_atmos_type),intent(inout)::a
 integer::nc,vid,isd,ied,jsd,jed,i,j,k
 real::pt(a%bd%isd:a%bd%ied+1,a%bd%jsd:a%bd%jed+1,a%npz),dp(a%bd%isd:a%bd%ied+1,a%bd%jsd:a%bd%jed+1,a%npz)
 real::dz(a%bd%isd:a%bd%ied+1,a%bd%jsd:a%bd%jed+1,a%npz),w(a%bd%isd:a%bd%ied+1,a%bd%jsd:a%bd%jed+1,a%npz)
 real::u(a%bd%isd:a%bd%ied,a%bd%jsd:a%bd%jed+1,a%npz),vc(a%bd%isd:a%bd%ied,a%bd%jsd:a%bd%jed+1,a%npz)
 real::v(a%bd%isd:a%bd%ied+1,a%bd%jsd:a%bd%jed,a%npz),uc(a%bd%isd:a%bd%ied+1,a%bd%jsd:a%bd%jed,a%npz)
 real::divg(a%bd%isd:a%bd%ied+1,a%bd%jsd:a%bd%jed+1,a%npz),p0,p1,pkz
 isd=a%bd%isd;ied=a%bd%ied;jsd=a%bd%jsd;jed=a%bd%jed
 call lam_nc_check(nf90_open(trim(a%neststruct%lam_directory)//'/lam.initial.nc',nf90_nowrite,nc))
 call lam_nc_check(nf90_inq_varid(nc,'ak',vid));call lam_nc_check(nf90_get_var(nc,vid,a%ak))
 call lam_nc_check(nf90_inq_varid(nc,'bk',vid));call lam_nc_check(nf90_get_var(nc,vid,a%bk))
 call lam_nc_check(nf90_close(nc))
 a%ptop=a%ak(1);a%ks=count(a%bk==0.)-1
 call read_snapshot(a,0,pt,dp,dz,w,u,v,uc,vc,divg,.true.)
 a%u=u;a%v=v;a%uc=uc;a%vc=vc;a%w=w(isd:ied,jsd:jed,:);a%phis=0.;a%q=0.;a%q_con=0.
 a%delp=dp(isd:ied,jsd:jed,:)
 a%pt=pt(isd:ied,jsd:jed,:)
 a%delz=dz(isd:ied,jsd:jed,:)
 do j=a%bd%js,a%bd%je
 do i=a%bd%is,a%bd%ie
 a%pe(i,1,j)=a%ptop
 do k=1,a%npz
 a%pe(i,k+1,j)=a%pe(i,k,j)+a%delp(i,j,k)
 enddo
 a%peln(i,:,j)=log(a%pe(i,:,j));a%pk(i,j,:)=a%pe(i,:,j)**kappa
 do k=1,a%npz
 if(a%flagstruct%hydrostatic) then
 a%pkz(i,j,k)=(a%pk(i,j,k+1)-a%pk(i,j,k))/(kappa*(a%peln(i,k+1,j)-a%peln(i,k,j)))
 else
 a%pkz(i,j,k)=((-rdgas/grav)*a%delp(i,j,k)*a%pt(i,j,k)/a%delz(i,j,k))**kappa
 endif
 enddo
 enddo
 enddo
 a%ps=a%ptop+sum(a%delp,dim=3)
 end subroutine
end module
