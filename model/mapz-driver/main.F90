#define ASSERT_(a) if(.not. a) then; print *, __FILE__, __LINE__; call MPI_Abort(MPI_COMM_WORLD, 1, mpierr); endif

program Main

  use fms_mod, only: fms_init
  use fv_mapz_mod, only : Lagrangian_to_Eulerian
  use fv_arrays_mod, only: fv_flags_type, fv_grid_type
  use mpp_domains_mod, only: domain2d
  use external_ic_mod, only: mpp_domain_decomp
  
  use input_array_extents_mod, only: ArrayExtents_T, get_input_array_extents
  use input_scalars_mod, only: Scalars_T, get_input_scalars
  use input_arrays_mod, only: Arrays_T, get_input_arrays
  use input_flagstruct_mod, only: get_input_flagstruct
  use input_gridstruct_mod, only: get_input_gridstruct
  use input_domain_layout_mod, only: get_domain_layout, get_npx_npy
  use output_mod, only: Output_T, get_output
  
  implicit none

  include "mpif.h"
  
  integer, parameter :: halo_width = 3, nregions = 6, grid_type = 0

  type(Scalars_T) :: scalars
  type(ArrayExtents_T) :: arrext
  type(Arrays_T) :: arrays
  type(Output_T) :: output
  type(domain2d) :: mydomain
  type(fv_flags_type) :: flagstruct
  type(fv_grid_type) :: gridstruct

  integer :: nx, ny, npes, npx, npy, file_handle, mpierr, irank, nranks
  real :: start, finish
  character(len=256) :: scalars_file, arrext_file, arrays_file
  character(len=256) :: gridstruct_file, flagstruct_file, output_file
  
  call MPI_Init(mpierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, irank, mpierr)
  call MPI_Comm_size(MPI_COMM_WORLD, nranks, mpierr)

  ! Rank specific files
  write(scalars_file, '(a12, i2.2, a4)') 'mapz_scalars', irank, '.bin'
  write(arrext_file, '(a18, i2.2, a4)') 'mapz_array_extents', irank, '.bin'
  write(arrays_file, '(a11, i2.2, a4)') 'mapz_arrays', irank, '.bin'
  write(gridstruct_file, '(a15, i2.2, a4)') 'mapz_gridstruct', irank, '.bin'
  write(flagstruct_file, '(a15, i2.2, a4)') 'mapz_flagstruct', irank, '.bin'
  write(output_file, '(a11, i2.2, a4)') 'mapz_output', irank, '.bin'
  
  ! Domain layout
  call get_domain_layout('layout.nml', nx, ny)
  ASSERT_(nx*ny == nranks)

  ! Read input data
  call get_npx_npy('input.nml', npx, npy)
  call get_input_scalars(scalars_file, scalars)
  call get_input_array_extents(arrext_file, arrext)
  call get_input_arrays(arrays_file, arrext, arrays)
  call get_input_flagstruct(flagstruct_file, flagstruct)
  call get_input_gridstruct(gridstruct_file, arrext, gridstruct)
  
  ! Create domain2d
  call fms_init(MPI_COMM_WORLD)
  npes = nx*ny
  call mpp_domain_decomp( &
       mydomain, npx, npy, nregions, halo_width, grid_type, &
       npes, nx, ny/nregions)

  ! This is what we are testing
  call cpu_time(start)
  call Lagrangian_to_Eulerian( &
       scalars%last_step, scalars%consv, &
       arrays%ps, arrays%pe, arrays%delp, arrays%pkz, arrays%pk, &
       scalars%mdt, scalars%pdt, &
       arrext%km, &
       arrext%is, arrext%ie, arrext%js, arrext%je, &
       arrext%isd, arrext%ied, arrext%jsd, arrext%jed, &
       arrext%nq, &
       scalars%nwat, scalars%sphum, &
       arrays%q_con, &
       arrays%u, arrays%v, arrays%w, &
       arrays%delz, arrays%pt, arrays%q, arrays%hs, &
       scalars%r_vir, scalars%cp, scalars%akap, &
       arrays%cappa, &
       scalars%kord_mt, scalars%kord_wz, &
       arrays%kord_tr, &
       scalars%kord_tm, &
       arrays%peln, arrays%te0_2d, &
       scalars%ng, & ! halo width
       arrays%ua, arrays%va, arrays%omga, arrays%te, arrays%ws, &
       scalars%fill, scalars%reproduce_sum, scalars%out_dt, &
       arrays%dtdt, &
       scalars%ptop, &
       arrays%ak, arrays%bk, arrays%pfull, &
       flagstruct, &
       gridstruct, &
       mydomain, &
       scalars%do_sat_adj, scalars%hydrostatic, scalars%hybrid_z, &
       scalars%do_omega, scalars%adiabatic, scalars%do_adiabatic_init, &
       arrays%mfx, arrays%mfy, arrays%cx, arrays%cy, &
       scalars%remap_option)
  call cpu_time(finish)

  ! Compare output
  call get_output(output_file, output)
  ASSERT_(sum(arrays%pkz) == output%sum_pkz)
  ASSERT_(sum(arrays%pt) == output%sum_pt)
  
  ! Print results
  call MPI_Barrier(MPI_COMM_WORLD, mpierr)
  if (irank == 0) then
     print *, ''
     print *, '        rank      sum_pkz           sum_pt         time taken (s)'
     print *, '     --------------------------------------------------------------'
  end if
  call MPI_Barrier(MPI_COMM_WORLD, mpierr)
  print *, irank, '  ', output%sum_pkz, '  ', output%sum_pt, '  ', finish-start
  
  call MPI_Finalize(mpierr)

end program main
