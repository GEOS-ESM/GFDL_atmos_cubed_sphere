module calc_gas_specific_heat_mlt_mod

    use fv_arrays_mod, only: fv_grid_type
    use msis_wrapper, only: msis_point
    use constants_mod, only: rdgas, cp_air
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

    implicit none
    private

    public :: calc_gas_specific_heat_mlt
    public :: calc_mlt_thermo_state
    public :: mlt_mixture_thermo_from_number_density

    real, parameter :: SAFE_MSIS_ALT_MIN_KM = 0.0
    real, parameter :: SAFE_MSIS_ALT_MAX_KM = 1000.0
    integer, parameter :: MAX_BAD_ALT_WARNINGS = 20
    integer, save :: bad_alt_warn_count = 0
    integer, save :: bad_msis_warn_count = 0

contains

    !====================================================================
    ! Compute mixture thermodynamics from MSIS number densities.
    !
    ! Notes:
    ! - MSIS densities are used diagnostically only.
    ! - The current GEOS-MLT branch uses O, N2, and O2 only.
    ! - Fractions returned here are mass fractions.
    !====================================================================
    subroutine mlt_mixture_thermo_from_number_density(Om, N2m, O2m, &
                                                      R_mix, Cp_mix, Cv_mix, Kappa_mix, &
                                                      phiO, phiN2, phiO2)
        implicit none

        real, intent(in)  :: Om, N2m, O2m
        real, intent(out) :: R_mix, Cp_mix, Cv_mix, Kappa_mix
        real, intent(out) :: phiO, phiN2, phiO2

        real :: rhoO, rhoN2, rhoO2, rhoTot

        ! SI constants
        real, parameter :: kB  = 1.380649e-23
        real, parameter :: amu = 1.66053906660e-27

        ! Molecular masses [kg]
        real, parameter :: mO  = 16.0 * amu
        real, parameter :: mO2 = 32.0 * amu
        real, parameter :: mN2 = 28.0 * amu

        ! Dry-air fallback. Fractions are not used when MSIS density is invalid.
        R_mix     = rdgas
        Cp_mix    = cp_air
        Cv_mix    = cp_air - rdgas
        Kappa_mix = rdgas / cp_air
        phiO      = 0.0
        phiN2     = 0.0
        phiO2     = 0.0

        rhoO   = Om  * mO
        rhoN2  = N2m * mN2
        rhoO2  = O2m * mO2
        rhoTot = rhoO + rhoN2 + rhoO2

        if (rhoTot > 0.0) then
           phiO  = rhoO  / rhoTot
           phiN2 = rhoN2 / rhoTot
           phiO2 = rhoO2 / rhoTot

           ! Mixture gas constant [J/(kg K)]
           R_mix = phiO  * (kB / mO)  + &
                   phiN2 * (kB / mN2) + &
                   phiO2 * (kB / mO2)

           ! Mixture specific heat at constant pressure [J/(kg K)].
           ! Monatomic O uses 5/2 R. Diatomic N2/O2 use 7/2 R.
           Cp_mix = phiO  * (2.5 * kB / mO)  + &
                    phiN2 * (3.5 * kB / mN2) + &
                    phiO2 * (3.5 * kB / mO2)

           Cv_mix = Cp_mix - R_mix

           if (Cp_mix > 0.0) then
              Kappa_mix = R_mix / Cp_mix
           else
              R_mix     = rdgas
              Cp_mix    = cp_air
              Cv_mix    = cp_air - rdgas
              Kappa_mix = rdgas / cp_air
              phiO      = 0.0
              phiN2     = 0.0
              phiO2     = 0.0
           endif
        endif

    end subroutine mlt_mixture_thermo_from_number_density


    !====================================================================
    ! Centralized diagnostic MSIS thermodynamic state.
    !
    ! This is the single source of MSIS-diagnosed thermodynamic properties
    ! for GEOS-MLT dycore calls. It does not transport composition. It only
    ! diagnoses the composition and mixture thermodynamics from the current
    ! grid, time, and optional height estimate.
    !====================================================================
    subroutine calc_mlt_thermo_state(is, ie, js, je, isd, ied, jsd, jed, km, pfull, &
                                     gridstruct, Cp_MLT, Kappa_MLT, year, doy, ut_seconds, &
                                     ifirst, ilast, jfirst, jlast, z_layer_km, &
                                     R_MLT, Cv_MLT, PhiO_MLT, PhiN2_MLT, PhiO2_MLT, &
                                     T_MSIS_MLT, Z_MSIS_MLT, Lambda_MLT, Rho_MLT, Alpha_MLT)

        implicit none

        ! --- Intent IN ---
        integer, intent(in) :: is, ie, js, je
        integer, intent(in) :: isd, ied, jsd, jed
        integer, intent(in) :: ifirst, ilast, jfirst, jlast
        integer, intent(in) :: km
        integer, intent(in) :: year, doy, ut_seconds
        real, intent(in) :: pfull(km)
        type(fv_grid_type), intent(in), target :: gridstruct
        real, intent(in), optional :: z_layer_km(isd:ied, jsd:jed, km)

        ! --- Required outputs ---
        real, intent(inout) :: Cp_MLT(isd:ied, jsd:jed, km)
        real, intent(inout) :: Kappa_MLT(isd:ied, jsd:jed, km)

        ! --- Optional outputs for shared diagnostic thermodynamic state ---
        real, intent(inout), optional :: R_MLT(isd:ied, jsd:jed, km)
        real, intent(inout), optional :: Cv_MLT(isd:ied, jsd:jed, km)
        real, intent(inout), optional :: PhiO_MLT(isd:ied, jsd:jed, km)
        real, intent(inout), optional :: PhiN2_MLT(isd:ied, jsd:jed, km)
        real, intent(inout), optional :: PhiO2_MLT(isd:ied, jsd:jed, km)
        real, intent(inout), optional :: T_MSIS_MLT(isd:ied, jsd:jed, km)
        real, intent(inout), optional :: Z_MSIS_MLT(isd:ied, jsd:jed, km)
        real, intent(inout), optional :: Lambda_MLT(isd:ied, jsd:jed, km)
        real, intent(inout), optional :: Rho_MLT(isd:ied, jsd:jed, km)
        real, intent(inout), optional :: Alpha_MLT(isd:ied, jsd:jed, km)

        ! --- Local scalars ---
        real :: lon_deg, lat_deg
        real :: R_mix, Cp_mix, Cv_mix, Kappa_mix
        real :: phiO, phiN2, phiO2
        real :: rhoO, rhoN2, rhoO2, rho_mix
        real :: lambda_mix, alpha_mix
        real :: ntot, xO, xN2, xO2
        real :: estz
        real :: stl
        real :: ut_hour
        logical :: valid_alt
        logical :: valid_msis
        logical :: use_geom_alt

        ! --- Local array ---
        real :: z_approx(km)

        ! --- Single-level MSIS outputs ---
        real(4) :: Om_k
        real(4) :: N2m_k
        real(4) :: O2m_k
        real(4) :: T_k

        ! --- Loop indices ---
        integer :: i, j, k

        ! --- Constants ---
        real, parameter :: rad2deg = 180.0 / 3.1415926535
        real, parameter :: kB  = 1.380649e-23
        real, parameter :: amu = 1.66053906660e-27
        real, parameter :: mO  = 16.0 * amu
        real, parameter :: mO2 = 32.0 * amu
        real, parameter :: mN2 = 28.0 * amu

        ! --- Fallback pressure-only height approximation ---
        real, parameter :: z_scale = 7.0
        real, parameter :: p_ref   = 1000.0

        Cp_MLT(:,:,:)    = cp_air
        Kappa_MLT(:,:,:) = rdgas / cp_air

        if (present(R_MLT))       R_MLT(:,:,:)       = rdgas
        if (present(Cv_MLT))      Cv_MLT(:,:,:)      = cp_air - rdgas
        if (present(PhiO_MLT))    PhiO_MLT(:,:,:)    = 0.0
        if (present(PhiN2_MLT))   PhiN2_MLT(:,:,:)   = 0.0
        if (present(PhiO2_MLT))   PhiO2_MLT(:,:,:)   = 0.0
        if (present(T_MSIS_MLT))  T_MSIS_MLT(:,:,:)  = 0.0
        if (present(Z_MSIS_MLT))  Z_MSIS_MLT(:,:,:)  = 0.0
        if (present(Lambda_MLT))  Lambda_MLT(:,:,:)  = 0.0
        if (present(Rho_MLT))     Rho_MLT(:,:,:)     = 0.0
        if (present(Alpha_MLT))   Alpha_MLT(:,:,:)   = 0.0

        do k = 1, km
           z_approx(k) = -z_scale * log(pfull(k)*0.01 / p_ref)
        enddo

        do j = jfirst, jlast
           do i = ifirst, ilast

              lon_deg = gridstruct%agrid(i, j, 1) * rad2deg
              lat_deg = gridstruct%agrid(i, j, 2) * rad2deg

              ut_hour = real(ut_seconds) / 3600.0
              stl = ut_hour + lon_deg/15.0

              if (stl < 0.0)  stl = stl + 24.0
              if (stl >= 24.0) stl = stl - 24.0

              do k = 1, km

                 ! Use geometric altitude everywhere that this routine is asked to
                 ! diagnose MSIS thermodynamics, including halo cells.  
                 use_geom_alt = present(z_layer_km)

                 if (use_geom_alt) then
                    estz = z_layer_km(i,j,k)
                 else
                    estz = z_approx(k)
                 endif

                 valid_alt = ieee_is_finite(estz) .and. &
                      estz >= SAFE_MSIS_ALT_MIN_KM .and. estz <= SAFE_MSIS_ALT_MAX_KM

                 if (.not. valid_alt) then
                    if (use_geom_alt .and. bad_alt_warn_count < MAX_BAD_ALT_WARNINGS) then
                       print *, 'GEOS_MLT_BAD_ALT_THERMO_FALLBACK: i,j,k,alt,z_approx,pfull=', &
                                i, j, k, estz, z_approx(k), pfull(k)
                    endif
                    if (use_geom_alt) bad_alt_warn_count = bad_alt_warn_count + 1

                    ! Fall back to pressure-based altitude only when the provided
                    ! geometric altitude is outside the safe MSIS range.
                    estz = z_approx(k)
                    valid_alt = ieee_is_finite(estz) .and. &
                         estz >= SAFE_MSIS_ALT_MIN_KM .and. estz <= SAFE_MSIS_ALT_MAX_KM
                 endif

                 if (.not. valid_alt) then
                    if (bad_alt_warn_count < MAX_BAD_ALT_WARNINGS) then
                       print *, 'GEOS_MLT_BAD_ALT_THERMO_AFTER_FALLBACK: i,j,k,alt,pfull=', &
                                i, j, k, estz, pfull(k)
                    endif
                    bad_alt_warn_count = bad_alt_warn_count + 1
                    cycle
                 endif

                 estz = min(max(estz, SAFE_MSIS_ALT_MIN_KM), SAFE_MSIS_ALT_MAX_KM)

                 call msis_point(year, doy, ut_seconds, estz, &
                      real(lat_deg, 4), real(lon_deg, 4), stl, &
                      Om_k, N2m_k, O2m_k, T_k)

                 valid_msis = ieee_is_finite(Om_k) .and. ieee_is_finite(N2m_k) .and. &
                      ieee_is_finite(O2m_k) .and. ieee_is_finite(T_k)
                 if (.not. valid_msis) then
                    if (bad_msis_warn_count < MAX_BAD_ALT_WARNINGS) then
                       print *, 'GEOS_MLT_BAD_MSIS_THERMO: i,j,k,alt,O,N2,O2,T=', &
                                i, j, k, estz, Om_k, N2m_k, O2m_k, T_k
                    endif
                    bad_msis_warn_count = bad_msis_warn_count + 1
                    Om_k = 0.0_4
                    N2m_k = 0.0_4
                    O2m_k = 0.0_4
                    T_k = 0.0_4
                 endif

                 call mlt_mixture_thermo_from_number_density(real(Om_k), real(N2m_k), real(O2m_k), &
                                                             R_mix, Cp_mix, Cv_mix, Kappa_mix, &
                                                             phiO, phiN2, phiO2)

                 ntot = real(Om_k) + real(N2m_k) + real(O2m_k)
                 xO  = real(Om_k)  / ntot
                 xN2 = real(N2m_k) / ntot
                 xO2 = real(O2m_k) / ntot

                 ! Convert MSIS number densities from cm-3 to m-3 for mass density.
                 rhoO    = real(Om_k)  * 1.0e6 * mO
                 rhoN2   = real(N2m_k) * 1.0e6 * mN2
                 rhoO2   = real(O2m_k) * 1.0e6 * mO2
                 rho_mix = rhoO + rhoN2 + rhoO2

                 ! Thermal-conductivity coefficient [W m-1 K-1 K^-0.69].
                 ! GEOS-MLT keeps the composition dependence here, but applies
                 ! the prognostic GEOS temperature dependence in the consumer:
                 !
                 !   lambda = lambda_mix * T_GEOS**0.69
                 !
                 ! Use mass fractions here to preserve the previous GEOS-MLT
                 ! conduction behavior during this diagnostic test.
                 if (rho_mix > 0.0 .and. Cp_mix > 0.0) then
                    lambda_mix = (56.0*(xO2+ xN2) + 75.9*xO) * 1.0e-5
                    !lambda_mix = (56.0*(phiO2 + phiN2) + 75.9*phiO) * 1.0e-5

                    ! Alpha_MLT is now a thermal-diffusivity coefficient:
                    !   alpha = alpha_mix * T_GEOS**0.69
                    alpha_mix  = lambda_mix / max(1.0e-30, rho_mix * Cp_mix)
                 else
                    lambda_mix = 0.0
                    alpha_mix  = 0.0
                 endif

                 Cp_MLT(i,j,k)    = Cp_mix
                 Kappa_MLT(i,j,k) = Kappa_mix

                 if (present(R_MLT))      R_MLT(i,j,k)      = R_mix
                 if (present(Cv_MLT))     Cv_MLT(i,j,k)     = Cv_mix
                 if (present(PhiO_MLT))   PhiO_MLT(i,j,k)   = phiO
                 if (present(PhiN2_MLT))  PhiN2_MLT(i,j,k)  = phiN2
                 if (present(PhiO2_MLT))  PhiO2_MLT(i,j,k)  = phiO2
                 if (present(T_MSIS_MLT)) T_MSIS_MLT(i,j,k) = real(T_k)
                 if (present(Z_MSIS_MLT)) Z_MSIS_MLT(i,j,k) = estz
                 if (present(Lambda_MLT)) Lambda_MLT(i,j,k) = lambda_mix
                 if (present(Rho_MLT))    Rho_MLT(i,j,k)    = rho_mix
                 if (present(Alpha_MLT))  Alpha_MLT(i,j,k)  = alpha_mix

              enddo
           enddo
        enddo

    end subroutine calc_mlt_thermo_state


    !====================================================================
    ! Backward-compatible wrapper used by existing GEOS-MLT call sites.
    ! New code should call calc_mlt_thermo_state directly when it also needs
    ! R, Cv, species fractions, MSIS temperature, or diagnostic height.
    !====================================================================
    subroutine calc_gas_specific_heat_MLT(is, ie, js, je, isd, ied, jsd, jed, km, pfull, &
                                      gridstruct, Cp_MLT, Kappa_MLT, year, doy, ut_seconds, &
                                      ifirst, ilast, jfirst, jlast, z_layer_km)

        implicit none

        integer, intent(in) :: is, ie, js, je
        integer, intent(in) :: isd, ied, jsd, jed
        integer, intent(in) :: ifirst, ilast, jfirst, jlast
        integer, intent(in) :: km
        integer, intent(in) :: year, doy, ut_seconds
        real, intent(in) :: pfull(km)
        type(fv_grid_type), intent(in), target :: gridstruct
        real, intent(in), optional :: z_layer_km(isd:ied, jsd:jed, km)
        real, intent(inout) :: Cp_MLT(isd:ied, jsd:jed, km)
        real, intent(inout) :: Kappa_MLT(isd:ied, jsd:jed, km)

        if (present(z_layer_km)) then
           call calc_mlt_thermo_state(is, ie, js, je, isd, ied, jsd, jed, km, pfull, &
                                      gridstruct, Cp_MLT, Kappa_MLT, year, doy, ut_seconds, &
                                      ifirst, ilast, jfirst, jlast, z_layer_km)
        else
           call calc_mlt_thermo_state(is, ie, js, je, isd, ied, jsd, jed, km, pfull, &
                                      gridstruct, Cp_MLT, Kappa_MLT, year, doy, ut_seconds, &
                                      ifirst, ilast, jfirst, jlast)
        endif

    end subroutine calc_gas_specific_heat_MLT

end module calc_gas_specific_heat_mlt_mod
