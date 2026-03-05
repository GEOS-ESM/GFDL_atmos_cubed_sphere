module calc_gas_specific_heat_mlt_mod

    use fv_arrays_mod, only: fv_grid_type
    use msis_wrapper, only: msis_point
!    use time_manager_mod, only: time_type, get_date, get_time

    implicit none
    private

    public :: calc_gas_specific_heat_mlt

contains

    subroutine calc_gas_specific_heat_MLT(is, ie, js, je, isd, ied, jsd, jed, km, pfull, &
                                           gridstruct, Cp_MLT, Kappa_MLT, year, month, day, hour, minute, second)

        implicit none

        ! --- Intent IN ---
        integer, intent(in) :: is, ie, js, je
        integer, intent(in) :: isd, ied, jsd, jed
        integer, intent(in) :: km
        real, intent(in) :: pfull(km)  
        type(fv_grid_type), intent(in), target :: gridstruct
        
        ! --- Intent OUT ---
        real, intent(out) :: Cp_MLT(isd:ied, jsd:jed, km)      ! Variable specific heat [J/(kg·K)]
        real, intent(out) :: Kappa_MLT(isd:ied, jsd:jed, km)   ! Variable kappa = Rg/Cp [-]

        ! --- Local scalars ---
        real :: lon_deg, lat_deg
        real :: Rg_MLT_k, Cp_MLT_k
        real :: estz  ! altitude in km for MSIS call
        real :: stl

        ! --- Local array --- 
        real :: z_approx(km)

        ! --- Single-level MSIS outputs (scalars, not arrays) ---
        real(4) :: Om_k    ! Atomic oxygen    at level k
        real(4) :: N2m_k   ! Diatomic N2      at level k
        real(4) :: O2m_k   ! Diatomic O2      at level k
        real(4) :: T_k     ! Temperature      at level k

        ! --- Loop indices ---
        integer :: i, j, k
        
        ! --- Time indices ---
        integer, intent(in) :: year, month, day, hour, minute, second

        ! --- Constants ---
        real, parameter :: rad2deg     = 180.0 / 3.1415926535
        real, parameter :: Rstar       = 8314.47  ! Universal gas constant [J/(kmol·K)]
        real, parameter :: Nmolar      = 14.0     ! Molar mass of N  [g/mol]
        real, parameter :: N2molar     = 28.0     ! Molar mass of N2 [g/mol]
        real, parameter :: Omolar      = 16.0     ! Molar mass of O  [g/mol]
        real, parameter :: O2molar     = 32.0     ! Molar mass of O2 [g/mol]
        real, parameter :: dof_diatomic = 7.0/2.0 ! N2, O2
        real, parameter :: dof_atomic   = 5.0/2.0 ! O, N


        ! --- Hacky way to calculate altitude grid. Should be fixed. ---

        real, parameter :: z_scale = 7.0  ! scale height in km (approximate)
        real, parameter :: p_ref = 1000.0  ! reference pressure (hPa)

        do k=1,km
           z_approx(k) = -z_scale * log(pfull(k)*0.01 / p_ref)  ! Convert Pa to hPa
        enddo

        ! --- Get year, month, day, and hour from the model ---

        !call get_date(Time, year, month, day, hour, minute, second)
            

        do j = js, je
           do i = is, ie

              lon_deg = gridstruct%agrid(i, j, 1) * rad2deg
              lat_deg = gridstruct%agrid(i, j, 2) * rad2deg
            
              ! Calculate solar local time in hours (0-24)
              stl = real(hour) + real(minute)/60.0 + real(second)/3600.0 + lon_deg/15.0
              
              ! Normalize to 0-24 range
              if (stl < 0.0) stl = stl + 24.0
              if (stl >= 24.0) stl = stl - 24.0

              ! Loop over each vertical level
              do k = 1, km
                  
                 estz = z_approx(k)

!                 print *, 'About to call msis_point'
!                 print *, 'year, month, day, hour:', year, month, day, hour
!                 print *, 'estz:', estz
!                 print *, 'lat_deg, lon_deg:', lat_deg, lon_deg
!                 print *, 'stl:', stl
!                 print *, 'i, j, k:', i, j, k

                 call msis_point(year, month, day, hour, estz, &
                        real(lat_deg, 4), real(lon_deg, 4), stl, &
                        Om_k, N2m_k, O2m_k, T_k)
                
!                 print *, 'msis_point returned successfully'

                 ! Gas constant for this composition
                 ! Neglecting atomic nitrogen for now
                 Rg_MLT_k = Rstar / ((N2m_k * N2molar) + &
                                      (Om_k * Omolar)  + &
                                      (O2m_k * O2molar))

                 ! Specific heat for this composition
                 Cp_MLT_k = Rg_MLT_k * ((dof_diatomic * (N2m_k + O2m_k)) + &
                                         (dof_atomic * Om_k))

                 ! Outputs needed
                 Cp_MLT(i, j, k)    = Cp_MLT_k
                 Kappa_MLT(i, j, k) = Rg_MLT_k / Cp_MLT_k

              enddo

           enddo
        enddo

    end subroutine calc_gas_specific_heat_MLT

end module calc_gas_specific_heat_mlt_mod
