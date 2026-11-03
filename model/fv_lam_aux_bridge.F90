!============ GEOS-LAM Update Start =============
! Exchange dry-host auxiliary meteorology across tested MPI neighbors.
! Isolated host-to-FV3 tracer boundary probe; not a production input provider.
module fv_lam_aux_bridge_mod
  use mpi, only: MPI_REAL, MPI_DOUBLE_PRECISION, MPI_INTEGER, &
       MPI_SUM, MPI_Comm_rank, MPI_Comm_size, MPI_Allreduce, MPI_Sendrecv, &
       MPI_STATUS_SIZE, MPI_SUCCESS
  implicit none
  private
  real(kind=4), allocatable, public :: aux_q(:,:,:,:)
  real(kind=8), allocatable, public :: aux_lat_global(:,:)
  real(kind=8), public :: aux_lower=-1.d0, aux_upper=-1.d0
  logical, public :: aux_ready=.false.
  public :: set_aux_samples
  ! Isolated meteorological bridge: each rank keeps its own auxiliary slab
  ! and four exchanged neighbor columns. Slots 1/2 are host time and +dt.
  real(kind=4), allocatable, public :: met_ps(:,:,:,:),met_t(:,:,:,:),met_w(:,:,:,:),met_zh(:,:,:,:)
  real(kind=4), allocatable, public :: met_uw(:,:,:,:),met_vw(:,:,:,:),met_us(:,:,:,:),met_vs(:,:,:,:)
  ! Global auxiliary x index of local slab element 1.
  integer, public :: met_x0=1, met_y0=1
  integer :: bridge_comm=-1,bridge_rank=-1,bridge_npe=0,bridge_layout(2)=0
  integer :: bridge_i0=0,bridge_j0=0
  real(kind=8), public :: met_lower=-1.d0,met_upper=-1.d0
  logical, public :: met_ready=.false.
  integer :: met_fields=0
  public :: set_met_field,set_met_surface,set_bridge_layout
contains
  subroutine set_bridge_layout(comm,rank,npe,layout,i0,j0)
    integer, intent(in) :: comm,rank,npe,layout(2),i0,j0
    integer :: ierr,actual_rank,actual_npe
    call MPI_Comm_rank(comm,actual_rank,ierr)
    if(ierr/=MPI_SUCCESS.or.actual_rank/=rank) error stop 'LAM bridge communicator rank'
    call MPI_Comm_size(comm,actual_npe,ierr)
    if(ierr/=MPI_SUCCESS.or.actual_npe/=npe) error stop 'LAM bridge communicator size'
    if(product(layout)/=npe.or.i0<1.or.j0<1) error stop 'LAM bridge grid metadata'
    bridge_comm=comm;bridge_rank=rank;bridge_npe=npe
    bridge_layout=layout;bridge_i0=i0;bridge_j0=j0
  end subroutine

  subroutine exchange_met3(local,dest,nx,ny,nz,slot,first_x,first_y)
    real(kind=4), intent(in) :: local(:,:,:)
    real(kind=4), allocatable, intent(inout) :: dest(:,:,:,:)
    integer, intent(in) :: nx,ny,nz,slot
    integer, intent(out) :: first_x,first_y
    real(kind=4) :: send(4,ny,nz),recv(4,ny,nz)
    real(kind=4) :: send_y(nx,4,nz),recv_y(nx,4,nz)
    integer :: rank,npe,peer,ierr,status(MPI_STATUS_SIZE),width,local_width,height,local_height
    rank=bridge_rank;npe=bridge_npe
    if(bridge_comm<0.or.rank<0) error stop 'LAM bridge layout unconfigured'
    if(size(local,3)/=nz) error stop 'LAM MET local level extent'
    select case(npe)
    case(1)
      if(any(shape(local)/=[nx,ny,nz])) error stop 'LAM MET local one-rank extent'
      width=nx;height=ny;first_x=1;first_y=1
      if(.not.allocated(dest)) allocate(dest(width,height,nz,2))
      if(any(shape(dest)/=[width,height,nz,2])) error stop 'LAM MET local destination shape'
      dest(:,:,:,slot)=local
    case(2)
      peer=1-rank
      if(all(bridge_layout==[2,1])) then
        if(ny/=size(local,2)) error stop 'LAM MET x-split y extent'
        local_width=12+merge(1,0,nx==25.and.bridge_i0==13)
        if(size(local,1)/=local_width) error stop 'LAM MET x-split x extent'
        width=local_width+4;height=ny;first_x=bridge_i0-merge(4,0,rank==1);first_y=1
        if(.not.allocated(dest)) allocate(dest(width,height,nz,2))
        if(any(shape(dest)/=[width,height,nz,2])) error stop 'LAM MET local destination shape'
        if(rank==0) then
          send=local(9:12,:,:)
        else
          send=local(1:4,:,:)
        endif
        call MPI_Sendrecv(send,4*ny*nz,MPI_REAL,peer,2971, &
             recv,4*ny*nz,MPI_REAL,peer,2971,bridge_comm,status,ierr)
        if(ierr/=MPI_SUCCESS) error stop 'LAM MET x halo exchange'
        if(rank==0) then
          dest(1:local_width,:,:,slot)=local
          dest(local_width+1:width,:,:,slot)=recv
        else
          dest(1:4,:,:,slot)=recv
          dest(5:width,:,:,slot)=local
        endif
      else if(all(bridge_layout==[1,2])) then
        if(nx/=size(local,1)) error stop 'LAM MET y-split x extent'
        local_height=12+merge(1,0,ny==25.and.bridge_j0==13)
        if(size(local,2)/=local_height) error stop 'LAM MET y-split y extent'
        width=nx;height=local_height+4;first_x=1;first_y=bridge_j0-merge(4,0,rank==1)
        if(.not.allocated(dest)) allocate(dest(width,height,nz,2))
        if(any(shape(dest)/=[width,height,nz,2])) error stop 'LAM MET local destination shape'
        if(rank==0) then
          send_y=local(:,9:12,:)
        else
          send_y=local(:,1:4,:)
        endif
        call MPI_Sendrecv(send_y,4*nx*nz,MPI_REAL,peer,2971, &
             recv_y,4*nx*nz,MPI_REAL,peer,2971,bridge_comm,status,ierr)
        if(ierr/=MPI_SUCCESS) error stop 'LAM MET y halo exchange'
        if(rank==0) then
          dest(:,1:local_height,:,slot)=local
          dest(:,local_height+1:height,:,slot)=recv_y
        else
          dest(:,1:4,:,slot)=recv_y
          dest(:,5:height,:,slot)=local
        endif
      else
        error stop 'LAM MET unsupported two-rank layout'
      endif
    case default
      error stop 'LAM MET local transfer supports one or two ranks only'
    end select
  end subroutine

  subroutine set_met_field(name,local,seconds,step_seconds)
    character(*), intent(in) :: name
    real(kind=4), intent(in) :: local(:,:,:)
    real(kind=8), intent(in) :: seconds,step_seconds
    integer :: slot
    if(name=='MET3_UW_E') then
      met_ready=.false.;met_fields=0
      met_lower=seconds;met_upper=seconds+step_seconds
    endif
    if(step_seconds<=0.d0.or.seconds/=met_lower.or.seconds+step_seconds/=met_upper) &
         error stop 'LAM MET clock disagreement'
    slot=1
    if(index(name,'_NEXT')>0) slot=2
    select case(name)
    case('MET3_UW_E','MET3_UW_E_NEXT')
      call exchange_met3(local,met_uw,25,24,32,slot,met_x0,met_y0)
    case('MET3_VW_N','MET3_VW_N_NEXT')
      call exchange_met3(local,met_vw,25,24,32,slot,met_x0,met_y0)
    case('MET3_US_E','MET3_US_E_NEXT')
      call exchange_met3(local,met_us,24,25,32,slot,met_x0,met_y0)
    case('MET3_VS_N','MET3_VS_N_NEXT')
      call exchange_met3(local,met_vs,24,25,32,slot,met_x0,met_y0)
    case('MET3_T','MET3_T_NEXT')
      call exchange_met3(local,met_t,24,24,32,slot,met_x0,met_y0)
    case('MET3_W','MET3_W_NEXT')
      call exchange_met3(local,met_w,24,24,32,slot,met_x0,met_y0)
    case('MET3_ZH','MET3_ZH_NEXT')
      call exchange_met3(local,met_zh,24,24,33,slot,met_x0,met_y0)
    case default
      error stop 'LAM MET unexpected field'
    end select
    met_fields=met_fields+1
  end subroutine

  subroutine set_met_surface(now,next,seconds,step_seconds)
    real(kind=4), intent(in) :: now(:,:),next(:,:)
    real(kind=8), intent(in) :: seconds,step_seconds
    real(kind=4) :: tmp(size(now,1),size(now,2),1)
    if(met_fields/=14.or.seconds/=met_lower.or.seconds+step_seconds/=met_upper) &
         error stop 'LAM MET incomplete matched sample'
    if(any(shape(next)/=shape(now))) error stop 'LAM MET pressure pair shape'
    tmp(:,:,1)=now;call exchange_met3(tmp,met_ps,24,24,1,1,met_x0,met_y0)
    tmp(:,:,1)=next;call exchange_met3(tmp,met_ps,24,24,1,2,met_x0,met_y0)
    met_ready=.true.
    call dump_local_met(seconds)
  end subroutine

  subroutine dump_local_met(seconds)
    real(kind=8), intent(in) :: seconds
    character(len=1) :: enabled
    character(len=128) :: path
    integer :: env_status,rank,ierr,slot,unit
    call get_environment_variable('LAM_MET_LOCAL_DUMP',enabled,status=env_status)
    if(env_status/=0.or.enabled/='1') return
    if(seconds/=0.d0.and.seconds/=10800.d0.and.seconds/=10920.d0) return
    call MPI_Comm_rank(bridge_comm,rank,ierr)
    if(ierr/=MPI_SUCCESS) error stop 'LAM MET dump rank'
    do slot=1,2
      call dump_field('ps',met_ps(:,:,:,slot))
      call dump_field('t',met_t(:,:,:,slot))
      call dump_field('w',met_w(:,:,:,slot))
      call dump_field('zh',met_zh(:,:,:,slot))
      call dump_field('uw',met_uw(:,:,:,slot))
      call dump_field('vw',met_vw(:,:,:,slot))
      call dump_field('us',met_us(:,:,:,slot))
      call dump_field('vs',met_vs(:,:,:,slot))
    enddo
    print *, 'LAM_MET_LOCAL_DUMP',rank,int(seconds),met_x0,met_y0,shape(met_ps),shape(met_us)
  contains
    subroutine dump_field(name,field)
      character(*), intent(in) :: name
      real(kind=4), intent(in) :: field(:,:,:)
      write(path,'("met-local-",A,"-r",I0,"-t",I6.6,"-s",I0,".bin")') &
           name,rank,int(seconds),slot
      open(newunit=unit,file=trim(path),status='new',access='stream',form='unformatted')
      write(unit)field
      close(unit)
    end subroutine
  end subroutine

  subroutine set_aux_samples(q1_now,q1_next,q2_now,q2_next,lat,seconds,step_seconds)
    real(kind=4), intent(in) :: q1_now(:,:),q1_next(:,:),q2_now(:,:),q2_next(:,:)
    real(kind=8), intent(in) :: lat(:,:)
    real(kind=8), intent(in) :: seconds,step_seconds
    real(kind=4) :: local_q(24,24,2,2)
    real(kind=8) :: local_lat(24,24)
    integer :: local_mask(24,24),global_mask(24,24)
    integer :: rank,npe,i0,i1,j0,j1,ierr
    rank=bridge_rank;npe=bridge_npe
    if(bridge_comm<0.or.rank<0) error stop 'LAM AUX layout unconfigured'
    select case(npe)
    case(1)
      if(any(shape(q1_now)/=[24,24])) error stop 'LAM AUX one-rank shape'
      i0=1;i1=24;j0=1;j1=24
    case(2)
      i0=bridge_i0;j0=bridge_j0
      i1=i0+size(q1_now,1)-1;j1=j0+size(q1_now,2)-1
      if(all(bridge_layout==[2,1])) then
        if(any(shape(q1_now)/=[12,24])) error stop 'LAM AUX 2x1 local shape'
      else if(all(bridge_layout==[1,2])) then
        if(any(shape(q1_now)/=[24,12])) error stop 'LAM AUX 1x2 local shape'
      else
        error stop 'LAM AUX unsupported layout'
      endif
    case default
      error stop 'LAM AUX probe supports one rank or 2x1 only'
    end select
    if(any(shape(q1_next)/=shape(q1_now)).or.any(shape(q2_now)/=shape(q1_now)).or. &
       any(shape(q2_next)/=shape(q1_now)).or.any(shape(lat)/=shape(q1_now))) &
         error stop 'LAM AUX bridge field shape mismatch'
    if (step_seconds<=0.d0) error stop 'LAM AUX bridge timestep'
    if (.not.allocated(aux_q)) allocate(aux_q(24,24,2,2))
    if (.not.allocated(aux_lat_global)) allocate(aux_lat_global(24,24))
    local_q=0.;local_lat=0.;local_mask=0
    local_q(i0:i1,j0:j1,1,1)=q1_now;local_q(i0:i1,j0:j1,1,2)=q1_next
    local_q(i0:i1,j0:j1,2,1)=q2_now;local_q(i0:i1,j0:j1,2,2)=q2_next
    local_lat(i0:i1,j0:j1)=lat
    local_mask(i0:i1,j0:j1)=1
    call MPI_Allreduce(local_mask,global_mask,24*24,MPI_INTEGER,MPI_SUM,bridge_comm,ierr)
    if(ierr/=MPI_SUCCESS.or.any(global_mask/=1)) error stop 'LAM AUX MPI coverage/overlap'
    call MPI_Allreduce(local_q,aux_q,24*24*2*2,MPI_REAL,MPI_SUM,bridge_comm,ierr)
    if(ierr/=MPI_SUCCESS) error stop 'LAM AUX MPI tracer gather'
    call MPI_Allreduce(local_lat,aux_lat_global,24*24,MPI_DOUBLE_PRECISION,MPI_SUM,bridge_comm,ierr)
    if(ierr/=MPI_SUCCESS) error stop 'LAM AUX MPI latitude gather'
    print *, 'LAM_AUX_REDIST',rank,npe,i0,i1,j0,j1,shape(q1_now), &
         minval(global_mask),maxval(global_mask)
    aux_lower=seconds;aux_upper=seconds+step_seconds
    aux_ready=.true.
  end subroutine
end module
!============ GEOS-LAM Update End =============
