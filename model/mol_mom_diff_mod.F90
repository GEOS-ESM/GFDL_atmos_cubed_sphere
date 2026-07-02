! Molecular momentum diffusion for GEOS-MLT.
!
! This module computes a vertical molecular momentum-diffusion tendency
! for A-grid horizontal winds. The base diffusivity is diagnosed from
! thermal conductivity using
!
!   alpha = lambda / (rho * cp)
!   nu    = Pr * alpha
!

module mol_mom_diff_mod

  implicit none
  private

  public :: mol_mom_diff_compute_tend

contains

  subroutine mol_mom_diff_compute_tend(label, is, ie, js, je, isd, ied, jsd, jed, npz, &
                                      dt, pr_mol, pmax_pa, kmax_apply, rmax, nu_max, &
                                      nu_factor, pe, gz, ua, va, lambda_mlt, rho_mlt, &
                                      cp_mlt, u_tend, v_tend, ke_heat_tend, diag_enabled)

    implicit none

    character(len=*), intent(in) :: label
    integer, intent(in) :: is, ie, js, je
    integer, intent(in) :: isd, ied, jsd, jed
    integer, intent(in) :: npz
    integer, intent(in) :: kmax_apply
    real, intent(in) :: dt
    real, intent(in) :: pr_mol
    real, intent(in) :: pmax_pa
    real, intent(in) :: rmax
    real, intent(in) :: nu_max
    real, intent(in) :: nu_factor

    real, intent(in) :: pe(is-1:ie+1, npz+1, js-1:je+1)
    real, intent(in) :: gz(isd:ied, jsd:jed, npz+1)
    real, intent(in) :: ua(isd:ied, jsd:jed, npz)
    real, intent(in) :: va(isd:ied, jsd:jed, npz)
    real, intent(in) :: lambda_mlt(isd:ied, jsd:jed, npz)
    real, intent(in) :: rho_mlt(isd:ied, jsd:jed, npz)
    real, intent(in) :: cp_mlt(isd:ied, jsd:jed, npz)
    real, intent(out) :: u_tend(isd:ied, jsd:jed, npz)
    real, intent(out) :: v_tend(isd:ied, jsd:jed, npz)
    real, intent(out) :: ke_heat_tend(isd:ied, jsd:jed, npz)
    logical, intent(in) :: diag_enabled

    real, parameter :: GRAV0 = 9.80665
    real, parameter :: MIN_DZ = 10.0
    real, parameter :: MIN_RHO = 1.0e-30
    real, parameter :: MIN_CP = 1.0e-6
    real, parameter :: MIN_LAMBDA = 0.0
    real, parameter :: MIN_PRESSURE = 1.0e-30
    real, parameter :: SECONDS_PER_DAY = 86400.0
    ! Smooth pressure taper to avoid a binary on/off layer near pmax_pa.
    ! Width is in natural-log pressure units. log(10) gives about one decade.
    real, parameter :: TAPER_DLOGP = 2.302585093
    real, parameter :: MIN_TAPER = 1.0e-6

    integer :: i, j, k
    integer :: kend
    integer :: count
    integer :: limited_count
    logical :: limited_here
    real :: p_layer
    real :: taper_fac
    real :: log_ratio
    real :: z_here, z_above, z_below
    real :: dz_layer, dz_up, dz_down
    real :: nu_here, nu_above, nu_below
    real :: nu_base_here, nu_factored_here, nu_used_here
    real :: nu_base_tmp, nu_factored_tmp
    real :: rho_here, rho_above, rho_below
    real :: mu_if
    real :: du_dz, dv_dz
    real :: u_flux_up, u_flux_down
    real :: v_flux_up, v_flux_down
    real :: u_t, v_t
    real :: abs_t
    real :: r_here
    real :: ke_heat_rate
    real :: ke_heat_ks
    real :: work_rate
    real :: sum_abs_tend
    real :: sum_work
    real :: sum_ke_heat_ks
    real :: max_abs_tend
    real :: max_abs_update
    real :: max_r
    real :: min_nu_base
    real :: max_nu_base
    real :: max_nu_factored
    real :: max_nu_used
    real :: max_ke_heat_ks

    u_tend(:,:,:) = 0.0
    v_tend(:,:,:) = 0.0
    ke_heat_tend(:,:,:) = 0.0

    kend = min(npz, max(1, kmax_apply))

    count = 0
    limited_count = 0
    sum_abs_tend = 0.0
    sum_work = 0.0
    sum_ke_heat_ks = 0.0
    max_abs_tend = 0.0
    max_abs_update = 0.0
    max_r = 0.0
    min_nu_base = huge(1.0)
    max_nu_base = 0.0
    max_nu_factored = 0.0
    max_nu_used = 0.0
    max_ke_heat_ks = 0.0

    do k = 1, kend
      do j = js, je
        do i = is, ie

          if (pe(i,k,j) <= MIN_PRESSURE .or. pe(i,k+1,j) <= MIN_PRESSURE) cycle
          p_layer = sqrt(pe(i,k,j) * pe(i,k+1,j))
          if (pmax_pa > MIN_PRESSURE) then
             log_ratio = log(max(p_layer, MIN_PRESSURE) / max(pmax_pa, MIN_PRESSURE))
             taper_fac = 0.5 * (1.0 - tanh(log_ratio / max(TAPER_DLOGP, 1.0e-6)))
          else
             taper_fac = 1.0
          endif
          if (taper_fac < MIN_TAPER) cycle

          dz_layer = abs(gz(i,j,k) - gz(i,j,k+1)) / GRAV0
          if (dz_layer < MIN_DZ) cycle

          z_here = 0.5 * (gz(i,j,k) + gz(i,j,k+1)) / GRAV0
          call calc_nu(i, j, k, dz_layer, nu_base_here, nu_factored_here, &
                       nu_used_here, limited_here)
          if (nu_used_here <= 0.0) cycle

          nu_here = nu_used_here
          rho_here = max(rho_mlt(i,j,k), MIN_RHO)
          u_flux_up = 0.0
          u_flux_down = 0.0
          v_flux_up = 0.0
          v_flux_down = 0.0
          ke_heat_rate = 0.0

          if (k > 1) then
            z_above = 0.5 * (gz(i,j,k-1) + gz(i,j,k)) / GRAV0
            dz_up = abs(z_here - z_above)
            if (dz_up >= MIN_DZ) then
              call calc_nu(i, j, k-1, dz_up, nu_base_tmp, nu_factored_tmp, &
                           nu_above, limited_here)
              if (nu_above > 0.0 .and. rho_mlt(i,j,k-1) > MIN_RHO) then
                rho_above = max(rho_mlt(i,j,k-1), MIN_RHO)
                mu_if = 0.5 * (rho_here * nu_here + rho_above * nu_above)
                du_dz = (ua(i,j,k-1) - ua(i,j,k)) / dz_up
                dv_dz = (va(i,j,k-1) - va(i,j,k)) / dz_up
                u_flux_up = mu_if * du_dz
                v_flux_up = mu_if * dv_dz
                ke_heat_rate = ke_heat_rate + 0.5 * (mu_if / rho_here) * &
                               (du_dz*du_dz + dv_dz*dv_dz)
              endif
            endif
          endif

          if (k < npz) then
            z_below = 0.5 * (gz(i,j,k+1) + gz(i,j,k+2)) / GRAV0
            dz_down = abs(z_below - z_here)
            if (dz_down >= MIN_DZ) then
              call calc_nu(i, j, k+1, dz_down, nu_base_tmp, nu_factored_tmp, &
                           nu_below, limited_here)
              if (nu_below > 0.0 .and. rho_mlt(i,j,k+1) > MIN_RHO) then
                rho_below = max(rho_mlt(i,j,k+1), MIN_RHO)
                mu_if = 0.5 * (rho_here * nu_here + rho_below * nu_below)
                du_dz = (ua(i,j,k+1) - ua(i,j,k)) / dz_down
                dv_dz = (va(i,j,k+1) - va(i,j,k)) / dz_down
                u_flux_down = mu_if * du_dz
                v_flux_down = mu_if * dv_dz
                ke_heat_rate = ke_heat_rate + 0.5 * (mu_if / rho_here) * &
                               (du_dz*du_dz + dv_dz*dv_dz)
              endif
            endif
          endif

          u_t = taper_fac * (u_flux_up + u_flux_down) / (rho_here * dz_layer)
          v_t = taper_fac * (v_flux_up + v_flux_down) / (rho_here * dz_layer)
          ke_heat_ks = taper_fac * ke_heat_rate / max(MIN_CP, cp_mlt(i,j,k))

          u_tend(i,j,k) = u_t
          v_tend(i,j,k) = v_t
          ke_heat_tend(i,j,k) = ke_heat_ks

          count = count + 1
          if (nu_factored_here > nu_used_here * (1.0 + 1.0e-6)) limited_count = limited_count + 1

          abs_t = max(abs(u_t), abs(v_t))
          work_rate = ua(i,j,k) * u_t + va(i,j,k) * v_t

          sum_abs_tend = sum_abs_tend + 0.5 * (abs(u_t) + abs(v_t))
          sum_work = sum_work + work_rate
          sum_ke_heat_ks = sum_ke_heat_ks + ke_heat_ks
          max_abs_tend = max(max_abs_tend, abs_t)
          max_abs_update = max(max_abs_update, abs(dt) * abs_t)
          max_ke_heat_ks = max(max_ke_heat_ks, ke_heat_ks)

          r_here = nu_here * abs(dt) / max(MIN_DZ*MIN_DZ, dz_layer*dz_layer)
          max_r = max(max_r, r_here)
          min_nu_base = min(min_nu_base, nu_base_here)
          max_nu_base = max(max_nu_base, nu_base_here)
          max_nu_factored = max(max_nu_factored, nu_factored_here)
          max_nu_used = max(max_nu_used, nu_here)

        enddo
      enddo
    enddo

    if (diag_enabled) then
      if (count > 0) then
        write(*,*) 'GEOS_MLT_MOMDIFF ', trim(label), &
             ' COUNT=', count, &
             ' PR=', pr_mol, &
             ' NU_SCALE=', nu_factor, &
             ' PMAX_PA=', pmax_pa, &
             ' KMAX=', kend, &
             ' NU_BASE_MIN=', min_nu_base, &
             ' NU_BASE_MAX=', max_nu_base, &
             ' NU_FACTORED_MAX=', max_nu_factored, &
             ' NU_USED_MAX=', max_nu_used, &
             ' R_MAX=', max_r, &
             ' N_LIMITED=', limited_count, &
             ' MEAN_ABS_TEND=', sum_abs_tend / real(count), &
             ' MAX_ABS_TEND=', max_abs_tend, &
             ' MAX_ABS_UPDATE=', max_abs_update, &
             ' MEAN_WORK=', sum_work / real(count), &
             ' MEAN_KE2HEAT_KS=', sum_ke_heat_ks / real(count), &
             ' MAX_KE2HEAT_KS=', max_ke_heat_ks, &
             ' MAX_KE2HEAT_KDAY=', max_ke_heat_ks * SECONDS_PER_DAY
      else
        write(*,*) 'GEOS_MLT_MOMDIFF ', trim(label), ' COUNT=0'
      endif
    endif

  contains

    subroutine calc_nu(ii, jj, kk, dz_ref, nu_base, nu_factored, nu_used, limited)
      implicit none
      integer, intent(in) :: ii, jj, kk
      real, intent(in) :: dz_ref
      real, intent(out) :: nu_base
      real, intent(out) :: nu_factored
      real, intent(out) :: nu_used
      logical, intent(out) :: limited
      real :: nu_cap
      real :: safe_factor

      nu_base = 0.0
      nu_factored = 0.0
      nu_used = 0.0
      limited = .false.

      if (lambda_mlt(ii,jj,kk) <= MIN_LAMBDA) return
      if (rho_mlt(ii,jj,kk) <= MIN_RHO) return
      if (cp_mlt(ii,jj,kk) <= MIN_CP) return

      nu_base = pr_mol * lambda_mlt(ii,jj,kk) / (rho_mlt(ii,jj,kk) * cp_mlt(ii,jj,kk))
      if (nu_base <= 0.0) return

      safe_factor = max(0.0, nu_factor)
      nu_factored = safe_factor * nu_base
      if (nu_factored <= 0.0) return

      nu_cap = rmax * max(MIN_DZ*MIN_DZ, dz_ref*dz_ref) / max(abs(dt), 1.0e-6)
      nu_used = min(nu_factored, nu_cap, nu_max)
      limited = nu_used < nu_factored * (1.0 - 1.0e-6)
    end subroutine calc_nu

  end subroutine mol_mom_diff_compute_tend

end module mol_mom_diff_mod




