module calc_gas_specific_heat_mlt_mod

    use fv_arrays_mod, only: fv_grid_type
    use msis_wrapper, only: msis_point

    implicit none
    private

    public :: calc_gas_specific_heat_mlt

contains

    subroutine calc_gas_specific_heat_MLT(is, ie, js, je, isd, ied, jsd, jed, km, pfull, &
                                           gridstruct, Cp_MLT, Kappa_MLT, year, month, day, hour, minute, second, &
                                           ifirst, ilast, jfirst, jlast, z_layer_km)

        implicit none

        ! --- Intent IN ---
        integer, intent(in) :: is, ie, js, je
        integer, intent(in) :: isd, ied, jsd, jed
        integer, intent(in) :: ifirst, ilast, jfirst, jlast
        integer, intent(in) :: km
        real, intent(in) :: pfull(km)
        type(fv_grid_type), intent(in), target :: gridstruct
        real, intent(in), optional :: z_layer_km(isd:ied, jsd:jed, km)

        ! --- Intent OUT ---
        real, intent(inout) :: Cp_MLT(isd:ied, jsd:jed, km)      ! Variable specific heat [J/(kg·K)]
        real, intent(inout) :: Kappa_MLT(isd:ied, jsd:jed, km)   ! Variable kappa = Rg/Cp [-]

        ! --- Local scalars ---
        real :: lon_deg, lat_deg
        real :: Rg_MLT_k, Cp_MLT_k
        real :: estz
        real :: stl

        real :: rhoO_k, rhoN2_k, rhoO2_k, rhoTot_k
        real :: phiO_k, phiN2_k, phiO2_k

        ! --- Local array ---
        real :: z_approx(km)

        ! --- Single-level MSIS outputs ---
        real(4) :: Om_k
        real(4) :: N2m_k
        real(4) :: O2m_k
        real(4) :: T_k

        ! --- Loop indices ---
        integer :: i, j, k

        ! --- Time indices ---
        integer, intent(in) :: year, month, day, hour, minute, second

        ! --- Constants ---
        real, parameter :: rad2deg = 180.0 / 3.1415926535

        ! SI constants
        real, parameter :: kB  = 1.380649e-23
        real, parameter :: amu = 1.66053906660e-27

        ! Molecular masses [kg]
        real, parameter :: mO  = 16.0 * amu
        real, parameter :: mO2 = 32.0 * amu
        real, parameter :: mN2 = 28.0 * amu

        ! Fallback dry-air values
        real, parameter :: rdgas  = 287.0
        real, parameter :: cp_air = 1004.0

        ! --- Hacky way to calculate altitude grid. Should be fixed separately. ---
        real, parameter :: z_scale = 7.0
        real, parameter :: p_ref   = 1000.0

        Cp_MLT(:,:,:) = cp_air
        Kappa_MLT(:,:,:) = rdgas/cp_air       

        do k = 1, km
           z_approx(k) = -z_scale * log(pfull(k)*0.01 / p_ref)
        enddo

        do j = jfirst, jlast
           do i = ifirst, ilast

              lon_deg = gridstruct%agrid(i, j, 1) * rad2deg
              lat_deg = gridstruct%agrid(i, j, 2) * rad2deg

              stl = real(hour) + real(minute)/60.0 + real(second)/3600.0 + lon_deg/15.0

              if (stl < 0.0)  stl = stl + 24.0
              if (stl >= 24.0) stl = stl - 24.0

              do k = 1, km

                 if (present(z_layer_km)) then
                    estz = z_layer_km(i,j,k)
                 else
                    estz = z_approx(k)
                 endif
                 
                 call msis_point(year, month, day, hour, estz, &
                      real(lat_deg, 4), real(lon_deg, 4), stl, &
                      Om_k, N2m_k, O2m_k, T_k)

                 ! Convert MSIS number densities to mass densities
                 rhoO_k  = Om_k  * mO
                 rhoN2_k = N2m_k * mN2
                 rhoO2_k = O2m_k * mO2

                 rhoTot_k = rhoO_k + rhoN2_k + rhoO2_k

                 if (rhoTot_k > 0.0) then

                    phiO_k  = rhoO_k  / rhoTot_k
                    phiN2_k = rhoN2_k / rhoTot_k
                    phiO2_k = rhoO2_k / rhoTot_k

                    ! Mixture gas constant [J/(kg K)]
                    Rg_MLT_k = phiO_k  * (kB / mO)  + &
                               phiN2_k * (kB / mN2) + &
                               phiO2_k * (kB / mO2)

                    ! Mixture specific heat at constant pressure [J/(kg K)]
                    Cp_MLT_k = phiO_k  * (2.5 * kB / mO)  + &
                               phiN2_k * (3.5 * kB / mN2) + &
                               phiO2_k * (3.5 * kB / mO2)

                    Cp_MLT(i, j, k)    = Cp_MLT_k
                    Kappa_MLT(i, j, k) = Rg_MLT_k / Cp_MLT_k

                 else

                    Cp_MLT(i, j, k)    = cp_air
                    Kappa_MLT(i, j, k) = rdgas / cp_air

                 endif

              enddo

           enddo
        enddo

    end subroutine calc_gas_specific_heat_MLT

end module calc_gas_specific_heat_mlt_mod
