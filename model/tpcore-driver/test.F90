program main

#define VALIDATE

   use cudafor
   use tp_core_cuda, only : fv_tp_2d_cuda=>fv_tp_2d
   use tp_core_acc , only : fv_tp_2d_acc =>fv_tp_2d, alloc_memory, dealloc_memory
   use tp_core_cpu , only :                fv_tp_2d
   implicit none

   integer :: n        ! size of the vector

   real,dimension(:,:),allocatable        :: q
   real,dimension(:,:),allocatable        :: fxGPU
   real,dimension(:,:),allocatable        :: fyGPU

   real,dimension(:,:),allocatable        :: crx
   real,dimension(:,:),allocatable        :: cry
   real,dimension(:,:),allocatable        :: fx
   real,dimension(:,:),allocatable        :: fy
   real,dimension(:,:),allocatable        :: xfx
   real,dimension(:,:),allocatable        :: yfx
   real,dimension(:,:),allocatable        :: dxa
   real,dimension(:,:),allocatable        :: dya
   real,dimension(:,:),allocatable        :: area

   real,dimension(:,:),allocatable,device ::    q_device
   real,dimension(:,:),allocatable,device ::  crx_device
   real,dimension(:,:),allocatable,device ::  cry_device
   real,dimension(:,:),allocatable,device ::   fx_device
   real,dimension(:,:),allocatable,device ::   fy_device
   real,dimension(:,:),allocatable,device ::  xfx_device
   real,dimension(:,:),allocatable,device ::  yfx_device
   real,dimension(:,:),allocatable,device ::  dxa_device
   real,dimension(:,:),allocatable,device ::  dya_device
   real,dimension(:,:),allocatable,device :: area_device

   real,dimension(:,:),allocatable,device ::   qi_device
   real,dimension(:,:),allocatable,device ::   qj_device
   real,dimension(:,:),allocatable,device ::  fx2_device
   real,dimension(:,:),allocatable,device ::  fy2_device


   integer :: niters=1

   integer :: stream(2)
   integer :: npx, npy
   integer :: is, ie, isd, ied
   integer :: js, je, jsd, jed
   integer :: i, j, k, l, nt
   integer :: ng=3
   integer :: ord=12

   real*8 :: wtimer
   real*8 :: csync  , ccomp  , cgpu  , chost
   real*8 :: m1, m2, c0, c1, c2, c3, c4
   real   :: gpu_time
   type(cudaEvent) :: gpu_start, gpu_end
   real   :: comp_time
   type(cudaEvent) :: comp_start, comp_end
   real   :: sync_time
   type(cudaEvent) :: sync_start, sync_end

   integer :: ierr, istat

   integer, external :: iargc
   character(10) :: arg1

   if( iargc() .gt. 0 )then
     call getarg( 1, arg1 )
     read(arg1,'(i10)') n
     call getarg( 2, arg1 )
     read(arg1,'(i10)') niters
   else
     n = 180
   endif
   if( n .le. 0 ) n = 180
   npx = n+1   ; npy = n+1
   is  = 1     ; ie  = n
   js  = 1     ; je  = n
   isd = is-ng ; ied = ie+ng
   jsd = js-ng ; jed = je+ng

   allocate(   q(isd:ied,jsd:jed))
   allocate( crx(isd:ied,jsd:jed))
   allocate( cry(isd:ied,jsd:jed))
   allocate( fx (isd:ied,jsd:jed))
   allocate( xfx(isd:ied,jsd:jed))
   allocate( fy (isd:ied,jsd:jed))
   allocate( yfx(isd:ied,jsd:jed))
   allocate( dxa(isd:ied,jsd:jed))
   allocate( dya(isd:ied,jsd:jed))
   allocate(area(isd:ied,jsd:jed))

   allocate(fxGPU(isd:ied,jsd:jed))
   allocate(fyGPU(isd:ied,jsd:jed))

! Initial Data
   do j=jsd,jed
     do i=isd,ied
        q(i,j) = sin(3.1415927*float(i*j)/float((npx-1)*(npy-1)))
     enddo
   enddo
   xfx(:,:) = 0.5
   crx(:,:) = 0.5
   yfx(:,:) = 0.5
   cry(:,:) = 0.5
   dxa(:,:) = 1.0
   dya(:,:) = 1.0
  area(:,:) = 1.0

! CPU code
   c0 = wtimer()
!$omp region
   do nt=1,niters
      call fv_tp_2d(ord, q, crx, cry, xfx, yfx, dxa, dya, area, fx, fy, &
                    isd, ied, jsd, jed, is, ie, js, je, npx, npy, ng)
   enddo
!$omp end region
   c1 = wtimer()
   chost = c1-c0
    print *, chost, ' seconds on  HOST'

! Device Data
   allocate(   q_device(isd:ied,jsd:jed))
   allocate( crx_device(isd:ied,jsd:jed))
   allocate( cry_device(isd:ied,jsd:jed))
   allocate(  fx_device(isd:ied,jsd:jed))
   allocate(  fy_device(isd:ied,jsd:jed))
   allocate( xfx_device(isd:ied,jsd:jed))
   allocate( yfx_device(isd:ied,jsd:jed))
   allocate( dxa_device(isd:ied,jsd:jed))
   allocate( dya_device(isd:ied,jsd:jed))
   allocate(area_device(isd:ied,jsd:jed))

   allocate(  qi_device(isd:ied,jsd:jed))
   allocate(  qj_device(isd:ied,jsd:jed))
   allocate( fx2_device(isd:ied,jsd:jed))
   allocate( fy2_device(isd:ied,jsd:jed))

   q_device   = q
   crx_device = crx
   cry_device = cry
   xfx_device = xfx
   yfx_device = yfx
   dxa_device = dxa
   dya_device = dya
   area_device = area

   istat = cudaEventCreate(gpu_start)
   istat = cudaEventCreate(gpu_end)
   istat = cudaEventCreate(sync_start)
   istat = cudaEventCreate(sync_end)

! GPU Device code
   call alloc_memory(isd, ied, jsd, jed)
   q_device=q
   istat = cudaEventRecord(gpu_start, 0)
   do nt=1,niters
#ifdef CUDAFOR
      call fv_tp_2d_acc(ord, q_device, crx_device, cry_device, xfx_device, yfx_device, &
                    dxa_device, dya_device, area_device, fx_device, fy_device, &
                    isd, ied, jsd, jed, is, ie, js, je, npx, npy, ng)
#else
      call fv_tp_2d_acc(ord, q, crx, cry, xfx, yfx, dxa, dya, area, fxGPU, fyGPU, &
                    isd, ied, jsd, jed, is, ie, js, je, npx, npy, ng)
#endif
   enddo
   istat = cudaEventRecord(sync_start, 0)
   istat = cudaThreadSynchronize()
   istat = cudaEventRecord(sync_end, 0)
   istat = cudaEventRecord(gpu_end, 0)
   fxGPU=fx_device
   fyGPU=fy_device
   call dealloc_memory()

   istat = cudaEventElapsedTime(gpu_time, gpu_start, gpu_end)
   istat = cudaEventElapsedTime(sync_time, sync_start, sync_end)
   cgpu = gpu_time/1000.0
   csync = sync_time/1000.0
   ccomp = cgpu-csync

#ifdef VALIDATE
    ! check the results
    do j=js,je
      do i=is,ie+1
        if( abs(fx(i,j) - fxGPU(i,j)) .gt. 0.000001 )then
          print *, 'FX0', i, j, fx(i,j), fxGPU(i,j)
        endif
      enddo
    enddo
    do j=js,je+1
      do i=is,ie
        if( abs(fy(i,j) - fyGPU(i,j)) .gt. 0.000001 )then
          print *, 'FY0', i, j, fy(i,j), fyGPU(i,j)
        endif
      enddo
    enddo
#endif

    print *, chost, ' seconds on  HOST'
#ifdef CUDAFOR
    print *, cgpu , ' seconds on  ACC'
#else
    print *, cgpu , ' seconds on  ACC on HOST'
#endif
    print *, ccomp, ' seconds for COMPUTE'
    print *, csync, ' seconds for SYNC'
    print *, chost/cgpu, ' SPEEDUP '
    print *, chost/ccomp, ' SPEEDUP (compute only)'

#ifdef SKIP


! CUDA Fortran FV_TP_2D
#ifdef _CUDA
! Setup for multiple stream processing
   istat = cudaStreamCreate(stream(1))
   istat = cudaStreamCreate(stream(2))
   istat = cudaEventRecord(gpu_start, 0)
   q_device = q
   do nt=1,niters
   istat = cudaDeviceSynchronize()
#ifdef CUDAFOR
      call fv_tp_2d_cuda(ord, q_device, qi_device, qj_device, crx_device, cry_device, xfx_device, yfx_device, &
                    dxa_device, dya_device, area_device, fx_device, fy_device, fx2_device, fy2_device, &
                    isd, ied, jsd, jed, is, ie, js, je, npx, npy, ng, stream)
#else
      call fv_tp_2d(ord, q, crx, cry, xfx, yfx, dxa, dya, area, fxGPU, fyGPU, &
                    isd, ied, jsd, jed, is, ie, js, je, npx, npy, ng)
#endif
   enddo
#ifdef CUDAFOR
   fxGPU = fx_device
   fyGPU = fy_device
#endif
   istat = cudaEventRecord(sync_start, 0)
   istat = cudaThreadSynchronize()
   istat = cudaEventRecord(sync_end, 0)

   istat = cudaEventRecord(gpu_end, 0)
   istat = cudaEventElapsedTime(gpu_time, gpu_start, gpu_end)
   istat = cudaEventElapsedTime(comp_time, comp_start, comp_end)
   istat = cudaEventElapsedTime(sync_time, sync_start, sync_end)
   cgpu = gpu_time/1000.0
   csync = sync_time/1000.0
   ccomp = cgpu-csync

#ifdef VALIDATE
    ! check the results
    do j=js,je
      do i=is,ie+1
        if( abs(fx(i,j) - fxGPU(i,j)) .gt. 1.e-7)then
          print *, 'FX0', i, j, fx(i,j), fxGPU(i,j), abs(fx(i,j) - fxGPU(i,j))
        endif
      enddo
    enddo
    do j=js,je+1
      do i=is,ie
        if( abs(fy(i,j) - fyGPU(i,j)) .gt. 1.e-7)then
          print *, 'FY0', i, j, fy(i,j), fyGPU(i,j), abs(fy(i,j) - fyGPU(i,j))
        endif
      enddo
    enddo
#endif
#endif

    print *, '                        '
    print *, chost, ' seconds on  HOST'

#ifdef SKIP
#ifdef CUDAFOR
    print *, cgpu , ' seconds on  CUDA'
#else
    print *, cgpu , ' seconds on  CUDA on HOST'
#endif
    print *, ccomp, ' seconds for COMPUTE'
    print *, csync, ' seconds for SYNC'
    print *, chost/cgpu, ' SPEEDUP '
    print *, chost/ccomp, ' SPEEDUP (compute only)'
#endif
#endif

contains

  real*8 function wtimer()
  real*4 etime, tarray(2)
    wtimer = etime(tarray)
!   wtimer = MPI_Wtime()
!   wtiemr = omp_get_wtime()
  end function wtimer

end program

