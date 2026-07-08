!======================================================================
! cond_driver_mod.F90
!
! Thermal conduction driver (local-domain indexing):
!  - Build T = pt * pkz on the local owned interior (exclude halos via ng)
!  - Use halo-aware i,j for pt/gz/MLT coefficients, but local ii,jj for pkz
!  - Use precomputed GEOS-MLT thermodynamic coefficients from geopk()
!  - Compute dT/dt from thermal conduction and return in heat_tc (K/s)
!
!======================================================================

module cond_driver_mod
  implicit none
  private
  public :: cond_driver_from_msis
  public :: cond_driver_from_coeffs
  public :: cond_driver_apply

  ! Sentinel threshold (huge_r=1e8 in R4 builds)
  real, parameter :: GZ_SENTINEL_THRESH = 1.0e7

  ! GEOS-MLT thermal conduction is active only in the upper atmosphere.
  ! The coefficient fields are precomputed by calc_mlt_thermo_state/geopk.
  real,    parameter :: GEOS_MLT_TC_ALT_MIN_KM = 70.0
  integer, parameter :: GEOS_MLT_TC_ACTIVE_KMAX = 40

  ! Diagnostic isolation switch for the top thermal-conduction boundary.
  ! .true. gives zero external heat flux at the model top.
  logical, parameter :: GEOS_MLT_TC_ZERO_TOP_FLUX = .true.

contains

  !------------------------------------------------------------
  ! Compute conduction temperature tendency from MSIS densities.
  !
  ! This legacy interface is kept for standalone tests. The main GEOS-MLT
  ! dycore path should call cond_driver_apply(), which reuses the centralized
  ! Lambda/Rho/Cp fields and avoids duplicate MSIS calls.
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
    real, intent(in)  :: T_ext_ref(:,:), z_top_km(:,:)

    integer :: is, ie, js, je, ks, ke
    integer :: i, j, kk

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
    allocate(top_flux_ij(is:ie,js:je))

    K_tc(:,:,:)     = 0.0
    alpha(:,:,:)    = 0.0
    rho(:,:,:)      = 0.0
    cp(:,:,:)       = 0.0
    top_flux_ij(:,:) = 0.0
    dTdt(:,:,:)     = 0.0

    call tc_calc(T, nO, nO2, nN2, K_tc, alpha)

    if (.not. GEOS_MLT_TC_ZERO_TOP_FLUX) then
      do j = js, je
        do i = is, ie
          dz_ref_m = (z_ref_km - z_top_km(i,j)) * 1000.0
          if (is_finite_real(dz_ref_m) .and. is_finite_real(T_ext_ref(i,j)) .and. &
              dz_ref_m > 1.0) then
            top_flux_ij(i,j) = -K_tc(i,j,ks) * &
                 (T(i,j,ks) - T_ext_ref(i,j)) / dz_ref_m
          endif
        enddo
      enddo
    endif

    do kk = ks, ke
      do j = js, je
        do i = is, ie
          rho(i,j,kk) = nO(i,j,kk)*mO + nO2(i,j,kk)*mO2 + nN2(i,j,kk)*mN2
          cp(i,j,kk)  = K_tc(i,j,kk) / max(1.0e-30, rho(i,j,kk)*alpha(i,j,kk))
        enddo
      enddo
    enddo

    call cond_z_tend(T, K_tc, rho, cp, dz, dz_if, dTdt, top_flux_ij)

    deallocate(K_tc, alpha, rho, cp, top_flux_ij)
  end subroutine cond_driver_from_msis


  !------------------------------------------------------------
  ! Compute conduction temperature tendency from precomputed
  ! Lambda/Rho/Cp coefficients.
  !------------------------------------------------------------
  subroutine cond_driver_from_coeffs(T, K_tc, rho, cp, dz, dz_if, dTdt)
    use cond_z_tend_mod, only : cond_z_tend
    implicit none

    real, intent(in)  :: T(:,:,:)
    real, intent(in)  :: K_tc(:,:,:)
    real, intent(in)  :: rho(:,:,:)
    real, intent(in)  :: cp(:,:,:)
    real, intent(in)  :: dz(:,:,:)
    real, intent(in)  :: dz_if(:,:,:)
    real, intent(out) :: dTdt(:,:,:)

    integer :: is, ie, js, je
    real, allocatable :: top_flux_ij(:,:)

    is = lbound(T,1); ie = ubound(T,1)
    js = lbound(T,2); je = ubound(T,2)

    allocate(top_flux_ij(is:ie,js:je))
    top_flux_ij(:,:) = 0.0
    dTdt(:,:,:) = 0.0

    ! The current GEOS-MLT configuration uses an isolated zero-flux top
    ! boundary. If a nonzero external top flux is reintroduced later, pass it
    ! as an explicit input here rather than calling MSIS again in this driver.
    call cond_z_tend(T, K_tc, rho, cp, dz, dz_if, dTdt, top_flux_ij)

    deallocate(top_flux_ij)
  end subroutine cond_driver_from_coeffs


  !------------------------------------------------------------
  ! Full driver called from dyn_core.
  !
  ! This path intentionally does not call NRLMSIS. It reuses the MLT
  ! thermodynamic fields already diagnosed by geopk()/calc_mlt_thermo_state.
  !------------------------------------------------------------
  subroutine cond_driver_apply(gz, pt, pkz, lambda_mlt, rho_mlt, cp_mlt, heat_tc, ng)
    implicit none

    real, intent(in)    :: gz(:,:,:)          ! interface geopotential (m^2/s^2)
    real, intent(in)    :: pt(:,:,:)          ! potential-temperature-like state
    real, intent(in)    :: pkz(:,:,:)         ! Exner-like factor on owned domain
    real, intent(in)    :: lambda_mlt(:,:,:)  ! thermal conductivity [W m-1 K-1]
    real, intent(in)    :: rho_mlt(:,:,:)     ! mass density [kg m-3]
    real, intent(in)    :: cp_mlt(:,:,:)      ! specific heat [J kg-1 K-1]
    real, intent(out)   :: heat_tc(:,:,:)     ! dTdt (K/s)
    integer, intent(in) :: ng                 ! halo width

    integer :: ilb, iub, jlb, jub
    integer :: is, ie, js, je, ks, ke
    integer :: ni, nj, nk
    integer :: ii, jj, kkL
    integer :: i, j, kk

    real, parameter :: grav = 9.80665

    real, allocatable :: Tcol(:,:,:)
    real, allocatable :: dzcol(:,:,:)
    real, allocatable :: dzifcol(:,:,:)
    real, allocatable :: Kcol(:,:,:)
    real, allocatable :: rhocol(:,:,:)
    real, allocatable :: cpcol(:,:,:)

    real :: T_here
    real :: raw_alt_km

    ! Always define output everywhere, including halos and inactive levels.
    heat_tc(:,:,:) = 0.0

    ilb = lbound(gz,1); iub = ubound(gz,1)
    jlb = lbound(gz,2); jub = ubound(gz,2)

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

    allocate(Tcol(1:ni,1:nj,1:nk))
    allocate(dzcol(1:ni,1:nj,1:nk))
    allocate(dzifcol(1:ni,1:nj,1:nk))
    allocate(Kcol(1:ni,1:nj,1:nk))
    allocate(rhocol(1:ni,1:nj,1:nk))
    allocate(cpcol(1:ni,1:nj,1:nk))

    Tcol(:,:,:)    = 0.0
    dzcol(:,:,:)   = 0.0
    dzifcol(:,:,:) = 0.0
    Kcol(:,:,:)    = 0.0
    rhocol(:,:,:)  = 0.0
    cpcol(:,:,:)   = 0.0

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
          if (gz(i,j,kk)   >= GZ_SENTINEL_THRESH) cycle
          if (gz(i,j,kk+1) >= GZ_SENTINEL_THRESH) cycle

          dzcol(ii,jj,kkL) = abs(gz(i,j,kk) - gz(i,j,kk+1)) / grav
        enddo
      enddo
    enddo

    do kkL = 1, nk-1
      dzifcol(:,:,kkL) = 0.5*(dzcol(:,:,kkL) + dzcol(:,:,kkL+1))
    enddo
    dzifcol(:,:,nk) = dzifcol(:,:,max(1,nk-1))

    do kk = ks, ke
      kkL = kk - ks + 1
      do jj = 1, nj
        j = js + jj - 1
        do ii = 1, ni
          i = is + ii - 1

          if (is_finite_real(pt(i,j,kk)) .and. is_finite_real(pkz(ii,jj,kk))) then
            T_here = pt(i,j,kk) * pkz(ii,jj,kk)
            if (is_finite_real(T_here)) Tcol(ii,jj,kkL) = T_here
          endif
        enddo
      enddo
    enddo

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

          ! Preserve the current thermal-conduction activation region while
          ! avoiding duplicate MSIS calls in this driver.
          if (kkL > GEOS_MLT_TC_ACTIVE_KMAX) cycle

          raw_alt_km = 0.5*(gz(i,j,kk) + gz(i,j,kk+1)) / grav / 1000.0
          if (.not. is_finite_real(raw_alt_km)) cycle
          if (raw_alt_km < GEOS_MLT_TC_ALT_MIN_KM) cycle

          if (.not. is_finite_real(lambda_mlt(i,j,kk))) cycle
          if (.not. is_finite_real(rho_mlt(i,j,kk))) cycle
          if (.not. is_finite_real(cp_mlt(i,j,kk))) cycle
          if (lambda_mlt(i,j,kk) < 0.0) cycle
          if (rho_mlt(i,j,kk) <= 0.0) cycle
          if (cp_mlt(i,j,kk) <= 0.0) cycle

          Kcol(ii,jj,kkL)   = lambda_mlt(i,j,kk)
          rhocol(ii,jj,kkL) = rho_mlt(i,j,kk)
          cpcol(ii,jj,kkL)  = cp_mlt(i,j,kk)
        enddo
      enddo
    enddo

    call cond_driver_from_coeffs(Tcol, Kcol, rhocol, cpcol, dzcol, dzifcol, &
                                 heat_tc(is:ie,js:je,ks:ke))

    deallocate(Tcol, dzcol, dzifcol, Kcol, rhocol, cpcol)
  end subroutine cond_driver_apply


  !------------------------------------------------------------
  ! Safety helper for local finite checks
  !------------------------------------------------------------
  logical function is_finite_real(x)
    implicit none
    real, intent(in) :: x

    ! This avoids an explicit IEEE module dependency and catches NaN/Inf.
    is_finite_real = (x == x) .and. (abs(x) < huge(x))
  end function is_finite_real

end module cond_driver_mod

