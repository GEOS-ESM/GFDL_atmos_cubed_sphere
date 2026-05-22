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
  ! Compute conduction temperature tendency (K/s)
  !------------------------------------------------------------
  subroutine cond_driver_from_msis(T, nO, nO2, nN2, dz, dz_if, dTdt, T_ext_ref, z_top_km)
    use thermcond_mod,   only : tc_calc
    use cond_z_tend_mod, only : cond_z_tend
    implicit none

    real, intent(in)  :: T(:,:,:)
    real, intent(in)  :: nO(:,:,:), nO2(:,:,:), nN2(:,:,:)
    real, intent(in)  :: dz(:,:,:)
    real, intent(in)  :: dz_if(:,:,:)
    real, intent(out) :: dTdt(:,:,:)
    real, intent(in) :: T_ext_ref(:,:), z_top_km(:,:)

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

    do j = js, je
      do i = is, ie
        ! dz between model top and reference altitude (meters)
        dz_ref_m = (z_ref_km - z_top_km(i,j)) * 1000.0
!        print *, 'dz_ref_m', dz_ref_m
!        print *, 'z_top_km(i,j)', z_top_km(i,j), i, j
        if (dz_ref_m > 1.0) then
          ! F_up = -K_top * (T_top - T_ext)/dz
          top_flux_ij(i,j) = -K_tc(i,j,ks) * ( T(i,j,ks) - T_ext_ref(i,j) ) / dz_ref_m
        else
          top_flux_ij(i,j) = 0.0
        end if
      end do
    end do

    do kk = ks, ke
      do j = js, je
        do i = is, ie
          rho(i,j,kk) = nO(i,j,kk)*mO + nO2(i,j,kk)*mO2 + nN2(i,j,kk)*mN2
          cp(i,j,kk)  = K_tc(i,j,kk) / max(1.0e-30, (rho(i,j,kk)*alpha(i,j,kk)))
        end do
      end do
    end do

    call cond_z_tend(T, K_tc, rho, cp, dz, dz_if, dTdt, top_flux_ij)

    deallocate(K_tc, alpha, rho, cp, top_flux_ij)
  end subroutine cond_driver_from_msis


  !------------------------------------------------------------
  ! Full driver called from dyn_core
  !------------------------------------------------------------
  subroutine cond_driver_apply(agrid, gz, pt, pkz, heat_tc, ng, &
                               year_msis, mon_msis, day_msis, hour_msis)
    use msis_wrapper, only : msis_point
    implicit none

    real, intent(in)    :: agrid(:,:,:)          ! lon/lat radians (local storage)
    real, intent(in)    :: gz(:,:,:)             ! interface geopotential (m^2/s^2)
    real, intent(in)    :: pt(:,:,:)             ! temperature-like
    real, intent(in)    :: pkz(:,:,:)            ! Exner-like factor
    real, intent(out)   :: heat_tc(:,:,:)        ! dTdt (K/s)
    integer, intent(in) :: ng                    ! halo width

    integer, intent(in), optional :: year_msis, mon_msis, day_msis, hour_msis

    integer :: ilb,iub,jlb,jub,klb,kub
    integer :: is,ie,js,je,ks,ke
    integer :: ni,nj,nk
    integer :: ii,jj,kkL
    integer :: i,j,kk
    integer :: y,m,d,hh

    real, parameter :: pi = 3.14159265358979323846
    real, parameter :: rad2deg = 180.0/pi
    real, parameter :: grav = 9.80665

    real, allocatable :: Tcol(:,:,:), dzcol(:,:,:), dzifcol(:,:,:)
    real, allocatable :: nO(:,:,:), nO2(:,:,:), nN2(:,:,:)
    real, allocatable :: T_ext_ref(:,:), z_top_km(:,:)

    real :: lon_deg, lat_deg, stl_hr, alt_km
    real :: O_cm3, N2_cm3, O2_cm3, Tmsis
    logical :: msis_ok
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

    y  = 2017; m = 1; d = 14; hh = 0
    if (present(year_msis)) y  = year_msis
    if (present(mon_msis))  m  = mon_msis
    if (present(day_msis))  d  = day_msis
    if (present(hour_msis)) hh = hour_msis

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
    
          ! Skip if this or next interface is sentinel
          if (kk+1 > ubound(gz,3)) cycle
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
          T_here = pt(i,j,kk) * pkz(ii,jj,kk)
          Tcol(ii,jj,kkL) = T_here
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
          if (gz(i,j,kk)   >= GZ_SENTINEL_THRESH) cycle
          if (gz(i,j,kk+1) >= GZ_SENTINEL_THRESH) cycle

          lon_deg = modulo(agrid(i,j,1) * rad2deg, 360.0)
          lat_deg = agrid(i,j,2) * rad2deg
          stl_hr  = modulo(real(hh) + lon_deg/15.0, 24.0)

          alt_km = 0.5*(gz(i,j,kk) + gz(i,j,kk+1)) / grav / 1000.0
          if (alt_km < 0.0) alt_km = 0.0

          call msis_point(y, m, d, hh, alt_km, lat_deg, lon_deg, stl_hr, &
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
        stl_hr  = modulo(real(hh) + lon_deg/15.0, 24.0)
    
        ! model top altitude from gz at k=ks
        z_top_km(i,j) = 0.5*(gz(i,j,ks) + gz(i,j,ks+1)) / grav / 1000.0
    
        call msis_point(y, m, d, hh, 220.0, lat_deg, lon_deg, stl_hr, &
                        O_cm3, N2_cm3, O2_cm3, Tmsis)
        T_ext_ref(i,j) = Tmsis
        !print *,'MSIS external temperature: ', T_ext_ref(i,j)
      end do
    end do


    call cond_driver_from_msis(Tcol, nO, nO2, nN2, dzcol, dzifcol, heat_tc(is:ie, js:je, ks:ke), &
                               T_ext_ref, z_top_km)

    deallocate(Tcol, dzcol, dzifcol, nO, nO2, nN2, T_ext_ref, z_top_km)

  end subroutine cond_driver_apply

end module cond_driver_mod
