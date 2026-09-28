!======================================================================
! cond_driver_mod.F90
!
! Thermal conduction driver (local-domain indexing):
!  - Build T = pt * pkz on the local owned interior (exclude halos via ng)
!  - Use halo-aware i,j for pt, but local ii,jj for pkz.
!  - Use gz interfaces to form altitude for MSIS sampling
!  - Restrict MSIS sampling and the solve to the caller's pressure mask
!  - Use an MSIS ghost layer to close the upper thermal-conduction boundary
!  - Compute dT/dt from thermal conduction and return in heat_tc (K/s)
!
!======================================================================

module cond_driver_mod
  implicit none
  private
  public :: cond_driver_from_msis
  public :: cond_driver_apply

  ! Sentinel threshold (huge_r=1e8 in R4 builds)
  real, parameter :: GZ_SENTINEL_THRESH = 1.0e7

  ! MSIS safety guard for GEOS-MLT thermal conduction.
  ! These guards prevent a bad geometric altitude from reaching NRLMSIS.
  real,    parameter :: GEOS_MLT_MSIS_ALT_MIN_KM = 0.0
  real,    parameter :: GEOS_MLT_MSIS_ALT_MAX_KM = 400.0

  integer, parameter :: GEOS_MLT_MSIS_WARN_LIMIT = 80

  ! The upper thermal-conduction boundary is controlled at runtime through
  ! FV3 fv_core_nml and passed into this module by dyn_core.
  ! Positive top flux means heat enters the GEOS-MLT column from above.

  integer, save :: msis_alt_warn_count = 0

  ! One-time prints per MPI rank (per process)
  logical, save :: printed_minmax = .false.

  ! Diagnostics captured in cond_driver_apply, printed in cond_driver_from_msis
  integer, save :: diag_T0_count   = -1
  integer, save :: diag_T_count    = -1
  integer, save :: diag_pkz_le0    = -1
  real,    save :: diag_pkz_min    =  huge(1.0)
  real,    save :: diag_pkz_max    = -huge(1.0)

  ! k-by-k counts of pkz<=0 on the local conduction domain (saved until printed)
  integer, save :: diag_ks = 0, diag_ke = -1
  integer, allocatable, save :: diag_pkz0_by_k(:)

contains

  !------------------------------------------------------------
  ! Compute the backward-Euler conduction temperature tendency (K/s).
  !------------------------------------------------------------
  subroutine cond_driver_from_msis(T, nO, nO2, nN2, dz, dz_if, dTdt, &
                                   T_ext_ref, K_ext_ref, dt, zero_top_flux, &
                                   ghost_dz_factor, top_flux_scale)
    use thermcond_mod,    only : tc_calc
    use cond_z_tend_mod, only : cond_z_tend
    implicit none

    real, intent(in)  :: T(:,:,:)
    real, intent(in)  :: nO(:,:,:), nO2(:,:,:), nN2(:,:,:)
    real, intent(in)  :: dz(:,:,:)
    real, intent(in)  :: dz_if(:,:,:)
    real, intent(out) :: dTdt(:,:,:)
    real, intent(in)  :: T_ext_ref(:,:)
    real, intent(in)  :: K_ext_ref(:,:)
    real, intent(in)  :: dt
    logical, intent(in), optional :: zero_top_flux
    real, intent(in), optional :: ghost_dz_factor
    real, intent(in), optional :: top_flux_scale

    integer :: is, ie, js, je, ks, ke
    integer :: i, j, kk

    real, allocatable :: K_tc(:,:,:)
    real, allocatable :: alpha(:,:,:)
    real, allocatable :: rho(:,:,:)
    real, allocatable :: cp(:,:,:)
    real, allocatable :: top_conductance_ij(:,:)

    real :: dz_top_m
    real :: dz_ghost_m
    real :: K_top
    real :: K_ghost
    real :: thermal_resistance
    logical :: zero_top_flux_use
    real :: ghost_dz_factor_use
    real :: top_flux_scale_use

    real, parameter :: amu = 1.66053906660e-27
    real, parameter :: mO  = 16.0 * amu
    real, parameter :: mO2 = 32.0 * amu
    real, parameter :: mN2 = 28.0 * amu

    is = lbound(T,1); ie = ubound(T,1)
    js = lbound(T,2); je = ubound(T,2)
    ks = lbound(T,3); ke = ubound(T,3)

    ! Backward-compatible defaults if this public helper is called directly.
    zero_top_flux_use = .true.
    ghost_dz_factor_use = 1.0
    top_flux_scale_use = 1.0
    if (present(zero_top_flux)) zero_top_flux_use = zero_top_flux
    if (present(ghost_dz_factor)) ghost_dz_factor_use = ghost_dz_factor
    if (present(top_flux_scale)) top_flux_scale_use = top_flux_scale

    allocate(K_tc(is:ie,js:je,ks:ke))
    allocate(alpha(is:ie,js:je,ks:ke))
    allocate(rho(is:ie,js:je,ks:ke))
    allocate(cp(is:ie,js:je,ks:ke))
    allocate(top_conductance_ij(is:ie,js:je))

    K_tc(:,:,:)     = 0.0
    alpha(:,:,:)    = 0.0
    rho(:,:,:)      = 0.0
    cp(:,:,:)       = 0.0
    top_conductance_ij(:,:) = 0.0
    dTdt(:,:,:)     = 0.0

    call tc_calc(T, nO, nO2, nN2, K_tc, alpha)

    ! Build the external top-boundary conductance.
    !
    ! The GEOS top cell and the MSIS ghost cell are treated as two half-cell
    ! thermal resistances in series:
    !
    !   R = 0.5*dz_top/K_top + 0.5*dz_ghost/K_ghost
    !
    ! so the boundary conductance is
    !
    !   G_top = top_flux_scale / R.
    !
    ! The MSIS ghost temperature and G_top are passed separately to
    ! cond_z_tend so the boundary flux is evaluated implicitly as
    ! G_top * (T_ghost - T_top_new).
    if (.not. zero_top_flux_use) then
      do j = js, je
        do i = is, ie
          dz_top_m = abs(dz(i,j,ks))
          dz_ghost_m = ghost_dz_factor_use * dz_top_m
          K_top = K_tc(i,j,ks)
          K_ghost = K_ext_ref(i,j)

          if (.not. is_finite_real(dz_top_m)) cycle
          if (.not. is_finite_real(dz_ghost_m)) cycle
          if (.not. is_finite_real(K_top)) cycle
          if (.not. is_finite_real(K_ghost)) cycle
          if (.not. is_finite_real(T(i,j,ks))) cycle
          if (.not. is_finite_real(T_ext_ref(i,j))) cycle

          if (dz_top_m <= 1.0 .or. dz_ghost_m <= 1.0) cycle
          if (K_top <= 0.0 .or. K_ghost <= 0.0) cycle

          thermal_resistance = 0.5*dz_top_m/K_top + &
                               0.5*dz_ghost_m/K_ghost

          if (.not. is_finite_real(thermal_resistance)) cycle
          if (thermal_resistance <= 0.0) cycle

          top_conductance_ij(i,j) = top_flux_scale_use / thermal_resistance
        end do
      end do
    end if

    do kk = ks, ke
      do j = js, je
        do i = is, ie
          rho(i,j,kk) = nO(i,j,kk)*mO + nO2(i,j,kk)*mO2 + nN2(i,j,kk)*mN2
          cp(i,j,kk)  = K_tc(i,j,kk) / &
                        max(1.0e-30, rho(i,j,kk)*alpha(i,j,kk))
        end do
      end do
    end do

    call cond_z_tend(T, K_tc, rho, cp, dz, dz_if, dTdt, dt, &
                     T_ext_ref, top_conductance_ij)

    deallocate(K_tc, alpha, rho, cp, top_conductance_ij)
  end subroutine cond_driver_from_msis


  !------------------------------------------------------------
  ! Full driver called from dyn_core.
  !------------------------------------------------------------
  subroutine cond_driver_apply(agrid, gz, pt, pkz, heat_tc, dt, ng, &
                               conduction_active, year_msis, doy_msis, &
                               ut_seconds_msis, zero_top_flux, &
                               ghost_dz_factor, top_flux_scale)
    use msis_wrapper, only : msis_point
    use thermcond_mod, only : tc_point
    implicit none

    real, intent(in)    :: agrid(:,:,:)          ! lon/lat radians (local storage)
    real, intent(in)    :: gz(:,:,:)             ! interface geopotential (m^2/s^2)
    real, intent(in)    :: pt(:,:,:)             ! temperature-like
    real, intent(in)    :: pkz(:,:,:)            ! Exner-like factor
    real, intent(out)   :: heat_tc(:,:,:)        ! dTdt (K/s)
    real, intent(in)    :: dt                     ! dynamics time step (s)
    integer, intent(in) :: ng                     ! halo width
    logical, intent(in) :: conduction_active(:,:,:) ! owned interior, no halos

    integer, intent(in), optional :: year_msis, doy_msis, ut_seconds_msis
    logical, intent(in), optional :: zero_top_flux
    real, intent(in), optional :: ghost_dz_factor
    real, intent(in), optional :: top_flux_scale

    integer :: ilb, iub, jlb, jub, klb, kub
    integer :: is, ie, js, je, ks, ke
    integer :: ni, nj, nk
    integer :: ii, jj, kkL
    integer :: i, j, kk
    integer :: y, doy, utsec

    real, parameter :: pi = 3.14159265358979323846
    real, parameter :: rad2deg = 180.0/pi
    real, parameter :: grav = 9.80665

    real, allocatable :: Tcol(:,:,:), dzcol(:,:,:), dzifcol(:,:,:)
    real, allocatable :: nO(:,:,:), nO2(:,:,:), nN2(:,:,:)
    real, allocatable :: T_ext_ref(:,:), K_ext_ref(:,:)

    real :: lon_deg, lat_deg, stl_hr, alt_km, raw_alt_km
    real :: O_cm3, N2_cm3, O2_cm3, Tmsis
    real :: T_here
    real :: top_interface_km
    real :: ghost_alt_km
    real :: ghost_dz_m
    real :: nO_ghost, nO2_ghost, nN2_ghost
    real :: alpha_ghost
    logical :: msis_ok, alt_ok
    logical :: zero_top_flux_use
    real :: ghost_dz_factor_use
    real :: top_flux_scale_use

    ! Always define output everywhere (including halos).
    heat_tc(:,:,:) = 0.0

    ilb = lbound(gz,1); iub = ubound(gz,1)
    jlb = lbound(gz,2); jub = ubound(gz,2)
    klb = lbound(gz,3); kub = ubound(gz,3)

    is = ilb + ng
    ie = iub - ng
    js = jlb + ng
    je = jub - ng

    ks = lbound(pt,3)
    ke = ubound(pt,3)

    ni = max(0, ie - is + 1)
    nj = max(0, je - js + 1)
    nk = max(0, ke - ks + 1)

    if (ni <= 0 .or. nj <= 0 .or. nk <= 0) return

    if (size(conduction_active,1) /= ni .or. &
        size(conduction_active,2) /= nj .or. &
        size(conduction_active,3) /= nk) then
      error stop 'GEOS-MLT conduction pressure-mask shape mismatch'
    end if

    y     = 2017
    doy   = 14
    utsec = 0

    if (present(year_msis))       y     = year_msis
    if (present(doy_msis))        doy   = doy_msis
    if (present(ut_seconds_msis)) utsec = ut_seconds_msis

    ! Runtime upper-boundary controls. Defaults preserve the historical
    ! zero-flux behavior if an older direct caller omits the new arguments.
    zero_top_flux_use = .true.
    ghost_dz_factor_use = 1.0
    top_flux_scale_use = 1.0
    if (present(zero_top_flux)) zero_top_flux_use = zero_top_flux
    if (present(ghost_dz_factor)) ghost_dz_factor_use = ghost_dz_factor
    if (present(top_flux_scale)) top_flux_scale_use = top_flux_scale

    if (ghost_dz_factor_use <= 0.0) then
      error stop 'GEOS-MLT ghost_dz_factor must be greater than zero'
    end if
    if (top_flux_scale_use < 0.0) then
      error stop 'GEOS-MLT top_flux_scale must be nonnegative'
    end if

    allocate(Tcol(1:ni,1:nj,1:nk))
    allocate(dzcol(1:ni,1:nj,1:nk))
    allocate(dzifcol(1:ni,1:nj,1:nk))
    allocate(nO(1:ni,1:nj,1:nk))
    allocate(nO2(1:ni,1:nj,1:nk))
    allocate(nN2(1:ni,1:nj,1:nk))

    Tcol(:,:,:)    = 0.0
    dzcol(:,:,:)   = 0.0
    dzifcol(:,:,:) = 0.0
    nO(:,:,:)      = 0.0
    nO2(:,:,:)     = 0.0
    nN2(:,:,:)     = 0.0

    allocate(T_ext_ref(is:ie,js:je))
    allocate(K_ext_ref(is:ie,js:je))
    T_ext_ref(:,:) = 0.0
    K_ext_ref(:,:) = 0.0

    ! Build geometric layer thicknesses.
    do kk = ks, ke
      kkL = kk - ks + 1
      do jj = 1, nj
        j = js + jj - 1
        do ii = 1, ni
          i = is + ii - 1

          if (kk+1 > ubound(gz,3)) cycle
          if (.not. is_finite_real(gz(i,j,kk))) cycle
          if (.not. is_finite_real(gz(i,j,kk+1))) cycle

          if (gz(i,j,kk) >= GZ_SENTINEL_THRESH) cycle
          if (gz(i,j,kk+1) >= GZ_SENTINEL_THRESH) cycle

          dzcol(ii,jj,kkL) = abs(gz(i,j,kk) - gz(i,j,kk+1)) / grav
        end do
      end do
    end do

    do kkL = 1, nk-1
      dzifcol(:,:,kkL) = 0.5*(dzcol(:,:,kkL) + dzcol(:,:,kkL+1))
    end do
    dzifcol(:,:,nk) = dzifcol(:,:,max(1,nk-1))

    ! Build physical temperature on the local owned interior.
    do kk = ks, ke
      kkL = kk - ks + 1
      do jj = 1, nj
        j = js + jj - 1
        do ii = 1, ni
          i = is + ii - 1

          if (is_finite_real(pt(i,j,kk)) .and. &
              is_finite_real(pkz(ii,jj,kk))) then
            T_here = pt(i,j,kk) * pkz(ii,jj,kk)
            if (is_finite_real(T_here)) then
              Tcol(ii,jj,kkL) = T_here
            end if
          end if
        end do
      end do
    end do

    ! Sample MSIS composition at each active GEOS layer center.
    do kk = ks, ke
      kkL = kk - ks + 1
      do jj = 1, nj
        j = js + jj - 1
        do ii = 1, ni
          i = is + ii - 1

          if (kk+1 > ubound(gz,3)) cycle
          if (.not. is_finite_real(gz(i,j,kk))) cycle
          if (.not. is_finite_real(gz(i,j,kk+1))) cycle
          if (gz(i,j,kk) >= GZ_SENTINEL_THRESH) cycle
          if (gz(i,j,kk+1) >= GZ_SENTINEL_THRESH) cycle

          ! The caller selects the physical pressure domain on any vertical grid.
          ! Keep densities zero outside it. cond_z_tend excludes zero-density
          ! layers and their interfaces, imposing zero flux at that boundary.
          if (.not. conduction_active(ii,jj,kkL)) cycle

          lon_deg = modulo(agrid(i,j,1) * rad2deg, 360.0)
          lat_deg = agrid(i,j,2) * rad2deg
          stl_hr  = modulo(real(utsec)/3600.0 + lon_deg/15.0, 24.0)

          raw_alt_km = 0.5*(gz(i,j,kk) + gz(i,j,kk+1)) / grav / 1000.0
          call sanitize_msis_alt(raw_alt_km, i, j, kk, gz(i,j,kk), &
                                 gz(i,j,kk+1), lat_deg, lon_deg, y, doy, &
                                 utsec, alt_km, alt_ok)
          if (.not. alt_ok) cycle

          if (.not. is_finite_real(lat_deg) .or. &
              .not. is_finite_real(lon_deg) .or. &
              .not. is_finite_real(stl_hr)) then
            call print_msis_alt_warning('bad lon/lat/stl for layer MSIS call', &
                                        i, j, kk, raw_alt_km, alt_km, &
                                        gz(i,j,kk), gz(i,j,kk+1), &
                                        lat_deg, lon_deg, y, doy, utsec)
            cycle
          end if

          call msis_point(y, doy, utsec, alt_km, lat_deg, lon_deg, stl_hr, &
                          O_cm3, N2_cm3, O2_cm3, Tmsis)

          msis_ok = (O_cm3 > 0.0 .or. O2_cm3 > 0.0 .or. N2_cm3 > 0.0)
          if (msis_ok) then
            nO(ii,jj,kkL)  = O_cm3  * 1.0e6
            nO2(ii,jj,kkL) = O2_cm3 * 1.0e6
            nN2(ii,jj,kkL) = N2_cm3 * 1.0e6
          end if
        end do
      end do
    end do

    ! Build the MSIS ghost state above the geometric model-top interface.
    do jj = 1, nj
      j = js + jj - 1
      do ii = 1, ni
        i = is + ii - 1

        ! An inactive top layer has no external conductive boundary flux.
        if (.not. conduction_active(ii,jj,1)) cycle

        ! Default fallback is zero top flux.
        T_ext_ref(i,j) = Tcol(ii,jj,1)
        K_ext_ref(i,j) = 0.0

        if (zero_top_flux_use) cycle

        if (ks+1 > ubound(gz,3)) cycle
        if (.not. is_finite_real(gz(i,j,ks))) cycle
        if (.not. is_finite_real(gz(i,j,ks+1))) cycle
        if (gz(i,j,ks) >= GZ_SENTINEL_THRESH) cycle
        if (gz(i,j,ks+1) >= GZ_SENTINEL_THRESH) cycle
        if (.not. is_finite_real(dzcol(ii,jj,1))) cycle
        if (dzcol(ii,jj,1) <= 1.0) cycle

        lon_deg = modulo(agrid(i,j,1) * rad2deg, 360.0)
        lat_deg = agrid(i,j,2) * rad2deg
        stl_hr  = modulo(real(utsec)/3600.0 + lon_deg/15.0, 24.0)

        if (.not. is_finite_real(lat_deg) .or. &
            .not. is_finite_real(lon_deg) .or. &
            .not. is_finite_real(stl_hr)) then
          call print_msis_alt_warning('bad lon/lat/stl for ghost MSIS call', &
                                      i, j, ks, -999.0, -999.0, &
                                      gz(i,j,ks), gz(i,j,ks+1), &
                                      lat_deg, lon_deg, y, doy, utsec)
          cycle
        end if

        ! Use the geometrically upper interface of the top GEOS layer.
        top_interface_km = max(gz(i,j,ks), gz(i,j,ks+1)) / grav / 1000.0
        ghost_dz_m = ghost_dz_factor_use * dzcol(ii,jj,1)
        raw_alt_km = top_interface_km + 0.5*ghost_dz_m/1000.0

        call sanitize_msis_alt(raw_alt_km, i, j, ks, gz(i,j,ks), &
                               gz(i,j,ks+1), lat_deg, lon_deg, y, doy, &
                               utsec, ghost_alt_km, alt_ok)
        if (.not. alt_ok) cycle

        call msis_point(y, doy, utsec, ghost_alt_km, lat_deg, lon_deg, &
                        stl_hr, O_cm3, N2_cm3, O2_cm3, Tmsis)

        msis_ok = is_finite_real(Tmsis) .and. Tmsis > 0.0 .and. &
                  is_finite_real(O_cm3) .and. O_cm3 >= 0.0 .and. &
                  is_finite_real(O2_cm3) .and. O2_cm3 >= 0.0 .and. &
                  is_finite_real(N2_cm3) .and. N2_cm3 >= 0.0 .and. &
                  (O_cm3 > 0.0 .or. O2_cm3 > 0.0 .or. N2_cm3 > 0.0)
        if (.not. msis_ok) cycle

        nO_ghost  = O_cm3  * 1.0e6
        nO2_ghost = O2_cm3 * 1.0e6
        nN2_ghost = N2_cm3 * 1.0e6

        call tc_point(Tmsis, nO_ghost, nO2_ghost, nN2_ghost, &
                      K_ext_ref(i,j), alpha_ghost)

        if (.not. is_finite_real(K_ext_ref(i,j))) then
          K_ext_ref(i,j) = 0.0
          cycle
        end if
        if (K_ext_ref(i,j) <= 0.0) cycle

        T_ext_ref(i,j) = Tmsis
      end do
    end do

    call cond_driver_from_msis(Tcol, nO, nO2, nN2, dzcol, dzifcol, &
                               heat_tc(is:ie,js:je,ks:ke), &
                               T_ext_ref, K_ext_ref, dt, zero_top_flux_use, &
                               ghost_dz_factor_use, top_flux_scale_use)

    deallocate(Tcol, dzcol, dzifcol, nO, nO2, nN2, T_ext_ref, K_ext_ref)

  end subroutine cond_driver_apply


  !------------------------------------------------------------
  ! Safety helpers for MSIS input.
  !------------------------------------------------------------
  logical function is_finite_real(x)
    implicit none
    real, intent(in) :: x

    ! This avoids an explicit IEEE module dependency and catches NaN/Inf.
    is_finite_real = (x == x) .and. (abs(x) < huge(x))
  end function is_finite_real


  subroutine sanitize_msis_alt(raw_alt_km, i, j, k, gz_top, gz_bot, lat_deg, &
                               lon_deg, year_msis, doy_msis, utsec_msis, &
                               alt_km, alt_ok)
    implicit none

    real,    intent(in)  :: raw_alt_km
    integer, intent(in)  :: i, j, k
    real,    intent(in)  :: gz_top, gz_bot
    real,    intent(in)  :: lat_deg, lon_deg
    integer, intent(in)  :: year_msis, doy_msis, utsec_msis
    real,    intent(out) :: alt_km
    logical, intent(out) :: alt_ok

    alt_km = raw_alt_km
    alt_ok = .true.

    if (.not. is_finite_real(raw_alt_km)) then
      ! Do not call NRLMSIS with NaN or Inf altitude.
      call print_msis_alt_warning('non-finite MSIS altitude: skipping MSIS call', &
                                  i, j, k, raw_alt_km, 0.0, gz_top, gz_bot, &
                                  lat_deg, lon_deg, year_msis, doy_msis, &
                                  utsec_msis)
      alt_km = 0.0
      alt_ok = .false.
      return
    end if

    if (raw_alt_km < GEOS_MLT_MSIS_ALT_MIN_KM) then
      ! Negative altitude is geometrically suspicious for an MSIS lookup.
      ! Skip the call instead of clamping to the surface.
      alt_km = GEOS_MLT_MSIS_ALT_MIN_KM
      alt_ok = .false.
      call print_msis_alt_warning('negative MSIS altitude: skipping MSIS call', &
                                  i, j, k, raw_alt_km, alt_km, gz_top, gz_bot, &
                                  lat_deg, lon_deg, year_msis, doy_msis, &
                                  utsec_msis)
      return
    end if

    if (raw_alt_km > GEOS_MLT_MSIS_ALT_MAX_KM) then
      alt_km = GEOS_MLT_MSIS_ALT_MAX_KM
      call print_msis_alt_warning('high MSIS altitude: clamped', &
                                  i, j, k, raw_alt_km, alt_km, gz_top, gz_bot, &
                                  lat_deg, lon_deg, year_msis, doy_msis, &
                                  utsec_msis)
      return
    end if
  end subroutine sanitize_msis_alt


  subroutine print_msis_alt_warning(reason, i, j, k, raw_alt_km, used_alt_km, &
                                    gz_top, gz_bot, lat_deg, lon_deg, year_msis, &
                                    doy_msis, utsec_msis)
    implicit none

    character(len=*), intent(in) :: reason
    integer, intent(in) :: i, j, k
    integer, intent(in) :: year_msis, doy_msis, utsec_msis
    real, intent(in) :: raw_alt_km, used_alt_km
    real, intent(in) :: gz_top, gz_bot
    real, intent(in) :: lat_deg, lon_deg

    msis_alt_warn_count = msis_alt_warn_count + 1

    if (msis_alt_warn_count <= GEOS_MLT_MSIS_WARN_LIMIT) then
      write(*,*) 'GEOS_MLT_MSIS_ALT_WARNING count=', msis_alt_warn_count, &
                 ' reason=', trim(reason)
      write(*,*) '  i=', i, ' j=', j, ' k=', k, &
                 ' raw_alt_km=', raw_alt_km, ' used_alt_km=', used_alt_km
      write(*,*) '  gz_top=', gz_top, ' gz_bot=', gz_bot, &
                 ' lat_deg=', lat_deg, ' lon_deg=', lon_deg
      write(*,*) '  year=', year_msis, ' doy=', doy_msis, ' utsec=', utsec_msis
    else if (msis_alt_warn_count == GEOS_MLT_MSIS_WARN_LIMIT + 1) then
      write(*,*) 'GEOS_MLT_MSIS_ALT_WARNING: further warnings suppressed on this rank.'
    end if
  end subroutine print_msis_alt_warning

end module cond_driver_mod
