!======================================================================
! cond_driver_mod.F90
!
! Thermal conduction driver (local-domain indexing):
!  - Build T = pt * pkz on the local owned interior (exclude halos via ng)
!  - Use halo-aware i,j for pt, but local ii,jj for pkz.
!  - Use gz interfaces to form altitude for MSIS sampling
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

  ! Do not call NRLMSIS below this altitude for thermal-conduction lookup.
  ! The GEOS-MLT conduction application is restricted to the upper atmosphere,
  ! so lower-atmosphere MSIS calls are unnecessary and can create noisy diagnostics.
  real,    parameter :: GEOS_MLT_MSIS_LOOKUP_ALT_MIN_KM = 70.0
  integer, parameter :: GEOS_MLT_MSIS_WARN_LIMIT = 80

  ! Only sample NRLMSIS in the upper column for thermal conduction.
  ! This avoids harmless near-surface or below-ground geometric-altitude
  ! warnings from lower layers that are not used by GEOS-MLT conduction.
  integer, parameter :: GEOS_MLT_MSIS_LOOKUP_KMAX = 40

  ! Diagnostic isolation switch for the top thermal-conduction boundary.
  ! .true. gives zero external heat flux at the model top.
  logical, parameter :: GEOS_MLT_TC_ZERO_TOP_FLUX = .true.

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
  ! Compute the backward-Euler conduction temperature tendency (K/s)
  !------------------------------------------------------------
  subroutine cond_driver_from_msis(T, nO, nO2, nN2, dz, dz_if, dTdt, T_ext_ref, z_top_km, dt)
    use thermcond_mod,   only : tc_calc
    use cond_z_tend_mod, only : cond_z_tend
    implicit none

    real, intent(in)  :: T(:,:,:)
    real, intent(in)  :: nO(:,:,:), nO2(:,:,:), nN2(:,:,:)
    real, intent(in)  :: dz(:,:,:)
    real, intent(in)  :: dz_if(:,:,:)
    real, intent(out) :: dTdt(:,:,:)
    real, intent(in) :: T_ext_ref(:,:), z_top_km(:,:)
    real, intent(in) :: dt

    integer :: is,ie,js,je,ks,ke
    integer :: i,j,kk
    integer :: k0

    real, allocatable :: K_tc(:,:,:)
    real, allocatable :: alpha(:,:,:)
    real, allocatable :: rho(:,:,:)
    real, allocatable :: cp(:,:,:)
    real, allocatable :: top_flux_ij(:,:)
    real :: dz_ref_m
    real, parameter :: z_ref_km = 220.0

    real, parameter :: amu = 1.66053906660e-27
    real, parameter :: mO  = 16.0 * amu
    real, parameter :: mO2 = 32.0 * amu
    real, parameter :: mN2 = 28.0 * amu

    is = lbound(T,1); ie = ubound(T,1)
    js = lbound(T,2); je = ubound(T,2)
    ks = lbound(T,3); ke = ubound(T,3)


    allocate(K_tc(is:ie,js:je,ks:ke))
    allocate(alpha(is:ie,js:je,ks:ke))
    allocate(rho(is:ie,js:je,ks:ke))
    allocate(cp(is:ie,js:je,ks:ke))
    allocate(top_flux_ij(is:ie, js:je))
    top_flux_ij(:,:) = 0.0

    K_tc(:,:,:)  = 0.0
    alpha(:,:,:) = 0.0
    rho(:,:,:)   = 0.0
    cp(:,:,:)    = 0.0
    dTdt(:,:,:)  = 0.0

    call tc_calc(T, nO, nO2, nN2, K_tc, alpha)

    if (GEOS_MLT_TC_ZERO_TOP_FLUX) then
      ! Diagnostic isolation: use a zero external heat flux at the model top.
      ! This prevents the MSIS 220-km temperature reservoir from forcing the
      ! top layer while we diagnose thermal/geometric stability.
      top_flux_ij(:,:) = 0.0
    else
      do j = js, je
        do i = is, ie
          ! dz between model top and reference altitude (meters)
          dz_ref_m = (z_ref_km - z_top_km(i,j)) * 1000.0
          if (is_finite_real(dz_ref_m) .and. is_finite_real(T_ext_ref(i,j)) .and. &
              dz_ref_m > 1.0) then
            ! F_up = -K_top * (T_top - T_ext)/dz
            top_flux_ij(i,j) = -K_tc(i,j,ks) * &
                 (T(i,j,ks) - T_ext_ref(i,j)) / dz_ref_m
          else
            top_flux_ij(i,j) = 0.0
          end if
        end do
      end do
    end if

    do kk = ks, ke
      do j = js, je
        do i = is, ie
          rho(i,j,kk) = nO(i,j,kk)*mO + nO2(i,j,kk)*mO2 + nN2(i,j,kk)*mN2
          cp(i,j,kk)  = K_tc(i,j,kk) / max(1.0e-30, (rho(i,j,kk)*alpha(i,j,kk)))
        end do
      end do
    end do

    call cond_z_tend(T, K_tc, rho, cp, dz, dz_if, dTdt, dt, top_flux_ij)

    deallocate(K_tc, alpha, rho, cp, top_flux_ij)
  end subroutine cond_driver_from_msis


  !------------------------------------------------------------
  ! Full driver called from dyn_core
  !------------------------------------------------------------
  subroutine cond_driver_apply(agrid, gz, pt, pkz, heat_tc, dt, ng, &
                               year_msis, doy_msis, ut_seconds_msis)
    use msis_wrapper, only : msis_point
    implicit none

    real, intent(in)    :: agrid(:,:,:)          ! lon/lat radians (local storage)
    real, intent(in)    :: gz(:,:,:)             ! interface geopotential (m^2/s^2)
    real, intent(in)    :: pt(:,:,:)             ! temperature-like
    real, intent(in)    :: pkz(:,:,:)            ! Exner-like factor
    real, intent(out)   :: heat_tc(:,:,:)        ! dTdt (K/s)
    real, intent(in)    :: dt                     ! dynamics time step (s)
    integer, intent(in) :: ng                    ! halo width

    integer, intent(in), optional :: year_msis, doy_msis, ut_seconds_msis

    integer :: ilb,iub,jlb,jub,klb,kub
    integer :: is,ie,js,je,ks,ke
    integer :: ni,nj,nk
    integer :: ii,jj,kkL
    integer :: i,j,kk
    integer :: y, doy, utsec

    real, parameter :: pi = 3.14159265358979323846
    real, parameter :: rad2deg = 180.0/pi
    real, parameter :: grav = 9.80665

    real, allocatable :: Tcol(:,:,:), dzcol(:,:,:), dzifcol(:,:,:)
    real, allocatable :: nO(:,:,:), nO2(:,:,:), nN2(:,:,:)
    real, allocatable :: T_ext_ref(:,:), z_top_km(:,:)

    real :: lon_deg, lat_deg, stl_hr, alt_km, raw_alt_km
    real :: O_cm3, N2_cm3, O2_cm3, Tmsis
    logical :: msis_ok, alt_ok
    real :: T_here

    ! Always define output everywhere (including halos)
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

    if (ni <= 0 .or. nj <= 0 .or. nk <= 0) then
      return
    end if

    y     = 2017
    doy   = 14
    utsec = 0
    
    if (present(year_msis))       y     = year_msis
    if (present(doy_msis))        doy   = doy_msis
    if (present(ut_seconds_msis)) utsec = ut_seconds_msis

    allocate(Tcol(1:ni, 1:nj, 1:nk))
    allocate(dzcol(1:ni, 1:nj, 1:nk))
    allocate(dzifcol(1:ni, 1:nj, 1:nk))
    allocate(nO(1:ni, 1:nj, 1:nk))
    allocate(nO2(1:ni, 1:nj, 1:nk))
    allocate(nN2(1:ni, 1:nj, 1:nk))

    Tcol(:,:,:)    = 0.0
    dzcol(:,:,:)   = 0.0
    dzifcol(:,:,:) = 0.0
    nO(:,:,:)      = 0.0
    nO2(:,:,:)     = 0.0
    nN2(:,:,:)     = 0.0

    allocate(T_ext_ref(is:ie, js:je))
    allocate(z_top_km(is:ie, js:je))
    T_ext_ref(:,:) = 0.0
    z_top_km(:,:)  = 0.0

    do kk = ks, ke
      kkL = kk - ks + 1
      do jj = 1, nj
        j = js + jj - 1
        do ii = 1, ni
          i = is + ii - 1
    
          ! Skip if this or next interface is invalid or sentinel.
          if (kk+1 > ubound(gz,3)) cycle
          if (.not. is_finite_real(gz(i,j,kk))) cycle
          if (.not. is_finite_real(gz(i,j,kk+1))) cycle
             if (gz(i,j,kk)   >= GZ_SENTINEL_THRESH) then
                dzcol(ii,jj,kkL) = 0.0  ! Mark as invalid
                cycle
             end if
          if (gz(i,j,kk+1) >= GZ_SENTINEL_THRESH) then
             dzcol(ii,jj,kkL) = 0.0  ! Mark as invalid
             cycle
          end if
    
          dzcol(ii,jj,kkL) = abs(gz(i,j,kk) - gz(i,j,kk+1)) / grav
        end do
      end do
    end do

    do kkL = 1, nk-1
      dzifcol(:,:,kkL) = 0.5*(dzcol(:,:,kkL) + dzcol(:,:,kkL+1))
    end do
    dzifcol(:,:,nk) = dzifcol(:,:,max(1,nk-1))
    
    do kk = ks, ke
      kkL = kk - ks + 1
      do jj = 1, nj
        j = js + jj - 1
        do ii = 1, ni
          i = is + ii - 1
          if (is_finite_real(pt(i,j,kk)) .and. is_finite_real(pkz(ii,jj,kk))) then
            T_here = pt(i,j,kk) * pkz(ii,jj,kk)
            if (is_finite_real(T_here)) then
              Tcol(ii,jj,kkL) = T_here
            else
              Tcol(ii,jj,kkL) = 0.0
            end if
          else
            Tcol(ii,jj,kkL) = 0.0
          end if
        end do
      end do
    end do


    do kk = ks, ke
      kkL = kk - ks + 1
      do jj = 1, nj
        j = js + jj - 1
        do ii = 1, ni
          i = is + ii - 1

          if (kk+1 > ubound(gz,3)) cycle
          if (.not. is_finite_real(gz(i,j,kk))) cycle
          if (.not. is_finite_real(gz(i,j,kk+1))) cycle
          if (gz(i,j,kk)   >= GZ_SENTINEL_THRESH) cycle
          if (gz(i,j,kk+1) >= GZ_SENTINEL_THRESH) cycle

          ! Thermal conduction/MSIS coupling is only active in the upper
          ! atmosphere.  Do not call NRLMSIS in lower model layers, because
          ! terrain-following or near-surface gz can be slightly negative and
          ! creates misleading low-altitude warnings.
          if (kkL > GEOS_MLT_MSIS_LOOKUP_KMAX) cycle

          lon_deg = modulo(agrid(i,j,1) * rad2deg, 360.0)
          lat_deg = agrid(i,j,2) * rad2deg
          stl_hr  = modulo(real(utsec)/3600.0 + lon_deg/15.0, 24.0)

          raw_alt_km = 0.5*(gz(i,j,kk) + gz(i,j,kk+1)) / grav / 1000.0
          call sanitize_msis_alt(raw_alt_km, i, j, kk, gz(i,j,kk), gz(i,j,kk+1), &
                                 lat_deg, lon_deg, y, doy, utsec, alt_km, alt_ok)
          if (.not. alt_ok) cycle

          if (.not. is_finite_real(lat_deg) .or. .not. is_finite_real(lon_deg) .or. &
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

    do jj = 1, nj
      j = js + jj - 1
      do ii = 1, ni
        i = is + ii - 1
    
        lon_deg = modulo(agrid(i,j,1) * rad2deg, 360.0)
        lat_deg = agrid(i,j,2) * rad2deg
        stl_hr  = modulo(real(utsec)/3600.0 + lon_deg/15.0, 24.0)
    
        ! Model top altitude from gz at k=ks.
        ! Use the same safety path so weird top altitude values are printed.
        if (ks+1 <= ubound(gz,3) .and. is_finite_real(gz(i,j,ks)) .and. &
            is_finite_real(gz(i,j,ks+1))) then
          raw_alt_km = 0.5*(gz(i,j,ks) + gz(i,j,ks+1)) / grav / 1000.0
          call sanitize_msis_alt(raw_alt_km, i, j, ks, gz(i,j,ks), gz(i,j,ks+1), &
                                 lat_deg, lon_deg, y, doy, utsec, alt_km, alt_ok)
          if (alt_ok) then
            z_top_km(i,j) = alt_km
          else
            z_top_km(i,j) = 220.0
          end if
        else
          raw_alt_km = -999.0
          z_top_km(i,j) = 220.0
          call print_msis_alt_warning('bad gz for model-top altitude', &
                                      i, j, ks, raw_alt_km, z_top_km(i,j), &
                                      gz(i,j,ks), gz(i,j,min(ks+1,ubound(gz,3))), &
                                      lat_deg, lon_deg, y, doy, utsec)
        end if

        if (GEOS_MLT_TC_ZERO_TOP_FLUX) then
          ! No external top flux is used, so do not call NRLMSIS for the
          ! 220-km reference temperature.  This keeps all MSIS calls restricted
          ! to layer densities needed by the conduction operator.
          T_ext_ref(i,j) = Tcol(ii,jj,1)
        else
          if (.not. is_finite_real(lat_deg) .or. .not. is_finite_real(lon_deg) .or. &
              .not. is_finite_real(stl_hr)) then
            call print_msis_alt_warning('bad lon/lat/stl for top MSIS call', &
                                        i, j, ks, raw_alt_km, z_top_km(i,j), &
                                        gz(i,j,ks), gz(i,j,ks+1), &
                                        lat_deg, lon_deg, y, doy, utsec)
            T_ext_ref(i,j) = Tcol(ii,jj,1)
          else
            call msis_point(y, doy, utsec, 220.0, lat_deg, lon_deg, stl_hr, &
                            O_cm3, N2_cm3, O2_cm3, Tmsis)
            T_ext_ref(i,j) = Tmsis
          end if
        end if
        !print *,'MSIS external temperature: ', T_ext_ref(i,j)
      end do
    end do


    call cond_driver_from_msis(Tcol, nO, nO2, nN2, dzcol, dzifcol, heat_tc(is:ie, js:je, ks:ke), &
                               T_ext_ref, z_top_km, dt)

    deallocate(Tcol, dzcol, dzifcol, nO, nO2, nN2, T_ext_ref, z_top_km)

  end subroutine cond_driver_apply


  !------------------------------------------------------------
  ! Safety helpers for MSIS input
  !------------------------------------------------------------
  logical function is_finite_real(x)
    implicit none
    real, intent(in) :: x

    ! This avoids an explicit IEEE module dependency and catches NaN/Inf.
    is_finite_real = (x == x) .and. (abs(x) < huge(x))
  end function is_finite_real


  subroutine sanitize_msis_alt(raw_alt_km, i, j, k, gz_top, gz_bot, lat_deg, lon_deg, &
                               year_msis, doy_msis, utsec_msis, alt_km, alt_ok)
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
                                  lat_deg, lon_deg, year_msis, doy_msis, utsec_msis)
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
                                  lat_deg, lon_deg, year_msis, doy_msis, utsec_msis)
      return
    end if

    if (raw_alt_km < GEOS_MLT_MSIS_LOOKUP_ALT_MIN_KM) then
      ! This layer is below the GEOS-MLT conduction/MSIS lookup region.
      ! Skip silently to avoid noisy lower-atmosphere diagnostics.
      alt_ok = .false.
      return
    end if

    if (raw_alt_km > GEOS_MLT_MSIS_ALT_MAX_KM) then
      alt_km = GEOS_MLT_MSIS_ALT_MAX_KM
      call print_msis_alt_warning('high MSIS altitude: clamped', &
                                  i, j, k, raw_alt_km, alt_km, gz_top, gz_bot, &
                                  lat_deg, lon_deg, year_msis, doy_msis, utsec_msis)
      return
    end if
  end subroutine sanitize_msis_alt


  subroutine print_msis_alt_warning(reason, i, j, k, raw_alt_km, used_alt_km, gz_top, &
                                    gz_bot, lat_deg, lon_deg, year_msis, doy_msis, &
                                    utsec_msis)
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
