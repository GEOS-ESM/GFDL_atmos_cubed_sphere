module input_mod

  use iso_fortran_env, only: real32
  use input_types_mod, only: ArrayExtents_T, Scalars_T, Arrays_T
  use fv_arrays_mod, only: fv_flags_type, fv_grid_type
  
  implicit none
  
  private

  public get_npx_npy, get_domain_layout
  public get_input_scalars, get_input_array_extents, get_input_arrays
  public get_input_flagstruct, get_input_gridstruct
  
contains

  subroutine get_domain_layout(filename, nx, ny)

    ! Arguments
    character(len=*), intent(in) :: filename
    integer, intent(out) :: nx
    integer, intent(out) :: ny

    ! Locals
    integer :: file_handle
    namelist /layout_nml/ nx, ny
    
    ! Start
    open(newunit = file_handle, file = filename, status = 'old')
    read(file_handle, nml = layout_nml)
    close(file_handle)
    
  end subroutine get_domain_layout
  
  subroutine get_npx_npy(filename, npx, npy)

    ! Arguments
    character(len=*) filename
    integer :: npx, npy

    ! Locals
    integer :: file_handle
    integer :: npz, n_sponge, fv_sg_adj, k_split, n_split
    integer :: hord_mt, hord_vt, hord_tm, hord_dp, hord_tr
    logical :: adiabatic, hydrostatic, make_nh, fv_debug, do_vort_damp
    real :: vtdm4, d_con
    namelist /fv_core_nml/ &
         npx, npy, npz, n_sponge, fv_sg_adj, adiabatic, hydrostatic, make_nh, &
         fv_debug, k_split, n_split, vtdm4, do_vort_damp, d_con, hord_mt, &
         hord_vt, hord_tm, hord_dp, hord_tr

    ! Start
    open(newunit = file_handle, file = filename, status = 'old')
    read(file_handle, nml = fv_core_nml)
    close(file_handle)
    
  end subroutine get_npx_npy
  
  subroutine get_input_scalars(filename, scalars)

    ! Arguments
    character(len=*), intent(in) :: filename
    type(Scalars_T), intent(out) :: scalars

    ! Locals
    integer :: file_handle
    
    ! Start
    print *, 'Reading scalar from file: ', trim(filename)

    open(newunit = file_handle, file = filename, form = 'unformatted', status = 'old')

    read(file_handle) scalars%last_step
    read(file_handle) scalars%consv
    read(file_handle) scalars%mdt
    read(file_handle) scalars%pdt
    read(file_handle) scalars%nwat
    read(file_handle) scalars%sphum
    read(file_handle) scalars%r_vir
    read(file_handle) scalars%cp
    read(file_handle) scalars%akap
    read(file_handle) scalars%kord_mt
    read(file_handle) scalars%kord_wz
    read(file_handle) scalars%kord_tm
    read(file_handle) scalars%ng
    read(file_handle) scalars%fill
    read(file_handle) scalars%reproduce_sum
    read(file_handle) scalars%out_dt
    read(file_handle) scalars%ptop
    read(file_handle) scalars%do_sat_adj
    read(file_handle) scalars%hydrostatic
    read(file_handle) scalars%hybrid_z
    read(file_handle) scalars%do_omega
    read(file_handle) scalars%adiabatic
    read(file_handle) scalars%do_adiabatic_init
    read(file_handle) scalars%remap_option

    close(file_handle)

  end subroutine get_input_scalars

  subroutine get_input_array_extents(filename, arrext)

    ! Arguments
    character(len=*), intent(in) :: filename
    type(ArrayExtents_T), intent(out) :: arrext

    ! Locals
    integer :: file_handle
    
    ! Start
    print *, 'Reading array extents from file: ', trim(filename)

    open(newunit = file_handle, file = filename, form = 'unformatted', status = 'old')

    read(file_handle) arrext%is
    read(file_handle) arrext%ie
    read(file_handle) arrext%js
    read(file_handle) arrext%je
    read(file_handle) arrext%isd
    read(file_handle) arrext%ied
    read(file_handle) arrext%jsd
    read(file_handle) arrext%jed
    read(file_handle) arrext%nq
    read(file_handle) arrext%km

  end subroutine get_input_array_extents

  subroutine allocate_memory_for_arrays_(arrext, arrays)

    ! Arguments
    type(ArrayExtents_T), intent(in) :: arrext
    type(Arrays_T), intent(inout) :: arrays
    integer :: is, ie, isd, ied, js, je, jsd, jed, km, nq

    ! Start

    ! Shorthands
    is = arrext%is; ie = arrext%ie
    isd = arrext%isd; ied = arrext%ied
    js = arrext%js; je = arrext%je
    jsd = arrext%jsd; jed = arrext%jed
    km = arrext%km
    nq = arrext%nq
    
    ! Now allocate memory
    allocate(arrays%ps(isd:ied, jsd:jed))
    allocate(arrays%pe(is-1:ie+1, km+1, js-1:je+1))
    allocate(arrays%delp(isd:ied, jsd:jed, km))
    allocate(arrays%pkz(is:ie, js:je, km))
    allocate(arrays%pk(is:ie, js:je, km+1))
    allocate(arrays%q_con(isd:isd, jsd:jsd, 1:1))
    allocate(arrays%u(isd:ied, jsd:jed+1, km))
    allocate(arrays%v(isd:ied+1, jsd:jed, km))
    allocate(arrays%w(isd:isd, jsd:jsd, 1:1))
    allocate(arrays%delz, mold = arrays%q_con)
    allocate(arrays%pt(isd:ied, jsd:jed, km))
    allocate(arrays%q(isd:ied, jsd:jed, km, nq))
    allocate(arrays%hs(isd:ied, jsd:jed))
    allocate(arrays%cappa, mold = arrays%q_con)
    allocate(arrays%kord_tr(nq))
    allocate(arrays%peln(is:ie, km+1, js:je))
    allocate(arrays%te0_2d(is:ie, js:je))
    allocate(arrays%ua(isd:ied, jsd:jed, km))
    allocate(arrays%va, mold = arrays%ua)
    allocate(arrays%omga, mold = arrays%ua)
    allocate(arrays%te(isd:ied, jsd:jed, km))
    allocate(arrays%ws(is:ie, js:je))
    allocate(arrays%dtdt(is:ie, js:je, km))
    allocate(arrays%ak(km+1))
    allocate(arrays%bk, mold = arrays%ak)
    allocate(arrays%pfull(km))
    allocate(arrays%mfx(is:ie+1, js:je, km))
    allocate(arrays%mfy(is:ie, js:je+1, km))
    allocate(arrays%cx(is:ie+1, jsd:jed, km))
    allocate(arrays%cy(isd:ied, js:je+1, km))
    
  end subroutine allocate_memory_for_arrays_
  
  subroutine get_input_arrays(filename, arrext, arrays)
  
    ! Arguments
    character(len=*), intent(in) :: filename
    type(ArrayExtents_T), intent(in) :: arrext
    type(Arrays_T), intent(out) :: arrays

    ! Locals
    integer :: is, ie, isd, ied, js, je, jsd, jed, km, nq
    integer :: file_handle
    
    ! Start
    print *, 'Reading array data from file: ', trim(filename)

    ! Allocate memory for all arrays
    call allocate_memory_for_arrays_(arrext, arrays)
    
    open(newunit = file_handle, file = filename, form = 'unformatted', status = 'old')

    read(file_handle) arrays%ps
    read(file_handle) arrays%pe
    read(file_handle) arrays%delp
    read(file_handle) arrays%pkz
    read(file_handle) arrays%pk
    read(file_handle) arrays%q_con
    read(file_handle) arrays%u
    read(file_handle) arrays%v
    read(file_handle) arrays%w
    read(file_handle) arrays%delz
    read(file_handle) arrays%pt
    read(file_handle) arrays%q(:,:,:,1)
    read(file_handle) arrays%hs
    read(file_handle) arrays%cappa
    read(file_handle) arrays%kord_tr
    read(file_handle) arrays%peln
    read(file_handle) arrays%te0_2d
    read(file_handle) arrays%ua
    read(file_handle) arrays%va
    read(file_handle) arrays%omga
    read(file_handle) arrays%te
    read(file_handle) arrays%ws
    read(file_handle) arrays%dtdt
    read(file_handle) arrays%ak
    read(file_handle) arrays%bk
    read(file_handle) arrays%pfull
    read(file_handle) arrays%mfx
    read(file_handle) arrays%mfy
    read(file_handle) arrays%cx
    read(file_handle) arrays%cy
    
    close(file_handle)

  end subroutine get_input_arrays

  subroutine get_input_flagstruct(filename, flagstruct)
  
    ! Arguments
    character(len=*), intent(in) :: filename
    type(fv_flags_type), intent(out) :: flagstruct

    ! Locals
    integer :: file_handle
    
    ! Start
    print *, 'Reading flagstruct data from file: ', trim(filename)

    open(newunit = file_handle, file = filename, form = 'unformatted', status = 'old')
    
    read(file_handle) flagstruct%fv_debug

    close(file_handle)
    
  end subroutine get_input_flagstruct

  subroutine get_input_gridstruct(filename, arrext, gridstruct)
  
    ! Arguments
    character(len=*), intent(in) :: filename
    type(ArrayExtents_T), intent(in) :: arrext
    type(fv_grid_type), intent(out) :: gridstruct

    ! Locals
    integer :: file_handle
    
    ! Start
    print *, 'Reading gridstruct data from file: ', trim(filename)

    allocate(gridstruct%rsin2(arrext%isd:arrext%ied, arrext%jsd:arrext%jed))
    allocate(gridstruct%cosa_s(arrext%isd:arrext%ied, arrext%jsd:arrext%jed))
    allocate(gridstruct%area_64(arrext%isd:arrext%ied, arrext%jsd:arrext%jed))
    
    open(newunit = file_handle, file = filename, form = 'unformatted', status = 'old')
    
    read(file_handle) gridstruct%rsin2
    read(file_handle) gridstruct%cosa_s
    read(file_handle) gridstruct%area_64

    close(file_handle)
    
  end subroutine get_input_gridstruct

end module input_mod
