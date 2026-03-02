module calc_gas_specific_heat_mlt_mod

    use fv_arrays_mod, only: fv_grid_type
    use msis_wrapper, only: msis_point

    implicit none
    private

    public :: calc_gas_specific_heat_mlt

contains

    subroutine calc_gas_specific_heat_MLT(is, ie, js, je, isd, ied, jsd, jed, km, &
                                           gridstruct, Cp_MLT, Kappa_MLT)

        implicit none

        ! --- Intent IN ---
        integer, intent(in) :: is, ie, js, je
        integer, intent(in) :: isd, ied, jsd, jed
        integer, intent(in) :: km
        type(fv_grid_type), intent(in), target :: gridstruct

        ! --- Intent OUT ---
        real, intent(out) :: Cp_MLT(isd:ied, jsd:jed, km)      ! Variable specific heat [J/(kg·K)]
        real, intent(out) :: Kappa_MLT(isd:ied, jsd:jed, km)   ! Variable kappa = Rg/Cp [-]

        ! --- Local scalars ---
        real :: lon_deg, lat_deg
        real :: Rg_MLT_k, Cp_MLT_k
        real :: alt_km  ! altitude in km for MSIS call
        real :: stl     ! solar local time (you'll need to compute this)

        ! --- Single-level MSIS outputs (scalars, not arrays) ---
        real(4) :: Om_k    ! Atomic oxygen    at level k
        real(4) :: N2m_k   ! Diatomic N2      at level k
        real(4) :: O2m_k   ! Diatomic O2      at level k
        real(4) :: T_k     ! Temperature      at level k

        ! --- Loop indices ---
        integer :: i, j, k
        integer :: year, month, day, hour

        ! --- Constants ---
        real, parameter :: rad2deg     = 180.0 / 3.14159265358979
        real, parameter :: Rstar       = 8314.47  ! Universal gas constant [J/(kmol·K)]
        real, parameter :: Nmolar      = 14.0     ! Molar mass of N  [g/mol]
        real, parameter :: N2molar     = 28.0     ! Molar mass of N2 [g/mol]
        real, parameter :: Omolar      = 16.0     ! Molar mass of O  [g/mol]
        real, parameter :: O2molar     = 32.0     ! Molar mass of O2 [g/mol]
        real, parameter :: dof_diatomic = 7.0/2.0 ! N2, O2
        real, parameter :: dof_atomic   = 5.0/2.0 ! O, N

        ! Set date/time for MSIS (you'll need to get these from your model)
        year  = 2015
        month = 5    ! May (day 150 of year)
        day   = 30
        hour  = 12
        stl   = 12.0 ! Solar local time - you may want to compute this from lon

        do j = js, je
           do i = is, ie

              lon_deg = gridstruct%agrid(i, j, 1) * rad2deg
              lat_deg = gridstruct%agrid(i, j, 2) * rad2deg

              ! Loop over each vertical level
              do k = 1, km
                 ! Get altitude for this level (you need to compute this from your model)
                 ! This is a placeholder - replace with actual altitude calculation
                 alt_km = real(k * 2.0, 4)  ! Example: 2 km spacing

                 ! Call MSIS for THIS level only
                 call msis_point(year, month, day, hour, alt_km, &
                                 real(lat_deg, 4), real(lon_deg, 4), stl, &
                                 Om_k, N2m_k, O2m_k, T_k)

                 ! Gas constant for this composition
                 ! Note: MSIS doesn't return atomic N directly
                 ! You may need to derive it or assume negligible
                 Rg_MLT_k = Rstar / ((N2m_k * N2molar) + &
                                      (Om_k  * Omolar)  + &
                                      (O2m_k * O2molar))

                 ! Specific heat for this composition
                 Cp_MLT_k = Rg_MLT_k * ((dof_diatomic * (N2m_k + O2m_k)) + &
                                         (dof_atomic   * Om_k))

                 ! Store outputs
                 Cp_MLT(i, j, k)    = Cp_MLT_k
                 Kappa_MLT(i, j, k) = Rg_MLT_k / Cp_MLT_k

              enddo

           enddo
        enddo

    end subroutine calc_gas_specific_heat_MLT

end module calc_gas_specific_heat_mlt_mod
