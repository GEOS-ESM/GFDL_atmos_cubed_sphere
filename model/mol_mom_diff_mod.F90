! Molecular momentum diffusion for GEOS-MLT.
!
! This module solves vertical molecular momentum diffusion for A-grid
! horizontal winds with a frozen-coefficient backward-Euler method.
! The kinematic viscosity is diagnosed from thermal conductivity using
!
!   lambda = lambda_coef * T_GEOS**0.69
!   alpha  = lambda / (rho * cp)
!   nu     = Pr * alpha
!
! The tridiagonal solve returns tendencies in the standard FV3 form:
!
!   tendency = (wind_new - wind_old) / dt
!
! Applying dt*tendency therefore reconstructs the implicitly solved wind.
!

module mol_mom_diff_mod

  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

  implicit none
  private

  public :: mol_mom_diff_compute_tend

contains

  subroutine mol_mom_diff_compute_tend(is, ie, js, je, isd, ied, jsd, jed, npz, dt, &
                                      pr_mol, pmax_pa, pe, gz, ua, va, lambda_mlt, &
                                      rho_mlt, cp_mlt, pt, pkz, u_tend, v_tend, &
                                      ke_heat_tend)

    implicit none

    integer, intent(in) :: is, ie, js, je
    integer, intent(in) :: isd, ied, jsd, jed
    integer, intent(in) :: npz
    real, intent(in) :: dt
    real, intent(in) :: pr_mol
    real, intent(in) :: pmax_pa

    real, intent(in) :: pe(is-1:ie+1, npz+1, js-1:je+1)
    real, intent(in) :: gz(isd:ied, jsd:jed, npz+1)
    real, intent(in) :: ua(isd:ied, jsd:jed, npz)
    real, intent(in) :: va(isd:ied, jsd:jed, npz)
    real, intent(in) :: lambda_mlt(isd:ied, jsd:jed, npz)
    real, intent(in) :: rho_mlt(isd:ied, jsd:jed, npz)
    real, intent(in) :: cp_mlt(isd:ied, jsd:jed, npz)
    real, intent(in) :: pt(isd:ied, jsd:jed, npz)
    real, intent(in) :: pkz(is:ie, js:je, npz)
    real, intent(out) :: u_tend(isd:ied, jsd:jed, npz)
    real, intent(out) :: v_tend(isd:ied, jsd:jed, npz)
    real, intent(out) :: ke_heat_tend(isd:ied, jsd:jed, npz)

    integer, parameter :: WORK_KIND = selected_real_kind(12, 100)
    real(kind=WORK_KIND), parameter :: GRAV0 = 9.80665_WORK_KIND
    real(kind=WORK_KIND), parameter :: MIN_DZ = 10.0_WORK_KIND
    real(kind=WORK_KIND), parameter :: MIN_RHO = 1.0e-30_WORK_KIND
    real(kind=WORK_KIND), parameter :: MIN_CP = 1.0e-6_WORK_KIND
    real(kind=WORK_KIND), parameter :: MIN_LAMBDA = 0.0_WORK_KIND
    real(kind=WORK_KIND), parameter :: MIN_PRESSURE = 1.0e-30_WORK_KIND
    real(kind=WORK_KIND), parameter :: DIAGONAL_MIN = 1.0e-20_WORK_KIND

    ! Smooth pressure taper to avoid a binary on/off layer near pmax_pa.
    ! The width is in natural-log pressure units. log(10) is one decade.
    real(kind=WORK_KIND), parameter :: TAPER_DLOGP = 2.302585093_WORK_KIND
    real(kind=WORK_KIND), parameter :: MIN_TAPER = 1.0e-6_WORK_KIND

    integer :: i, j, k
    real(kind=WORK_KIND) :: dt_work
    real(kind=WORK_KIND) :: pr_work
    real(kind=WORK_KIND) :: pmax_work
    real(kind=WORK_KIND) :: p_layer
    real(kind=WORK_KIND) :: log_ratio
    real(kind=WORK_KIND) :: t_geos
    real(kind=WORK_KIND) :: lambda_here
    real(kind=WORK_KIND) :: interface_taper
    real(kind=WORK_KIND) :: interface_viscosity
    real(kind=WORK_KIND) :: coupling
    real(kind=WORK_KIND) :: elimination_factor
    real(kind=WORK_KIND) :: grad_u
    real(kind=WORK_KIND) :: grad_v
    real(kind=WORK_KIND) :: dissipation_rate

    real(kind=WORK_KIND), allocatable :: lower(:)
    real(kind=WORK_KIND), allocatable :: diagonal(:)
    real(kind=WORK_KIND), allocatable :: upper(:)
    real(kind=WORK_KIND), allocatable :: rhs_u(:)
    real(kind=WORK_KIND), allocatable :: rhs_v(:)
    real(kind=WORK_KIND), allocatable :: u_new(:)
    real(kind=WORK_KIND), allocatable :: v_new(:)
    real(kind=WORK_KIND), allocatable :: z_center(:)
    real(kind=WORK_KIND), allocatable :: dz_layer(:)
    real(kind=WORK_KIND), allocatable :: mass_per_area(:)
    real(kind=WORK_KIND), allocatable :: rho_layer(:)
    real(kind=WORK_KIND), allocatable :: cp_layer(:)
    real(kind=WORK_KIND), allocatable :: dynamic_viscosity(:)
    real(kind=WORK_KIND), allocatable :: taper_factor(:)
    real(kind=WORK_KIND), allocatable :: interface_distance(:)
    real(kind=WORK_KIND), allocatable :: interface_mu(:)
    real(kind=WORK_KIND), allocatable :: ke_heat_work(:)
    logical, allocatable :: valid(:)

    u_tend(:,:,:) = 0.0
    v_tend(:,:,:) = 0.0
    ke_heat_tend(:,:,:) = 0.0

    if (.not. ieee_is_finite(dt) .or. dt <= 0.0) return
    if (.not. ieee_is_finite(pr_mol) .or. pr_mol <= 0.0) return
    if (npz <= 0) return

    dt_work = real(dt, kind=WORK_KIND)
    pr_work = real(pr_mol, kind=WORK_KIND)
    pmax_work = real(pmax_pa, kind=WORK_KIND)

    allocate(lower(1:npz))
    allocate(diagonal(1:npz))
    allocate(upper(1:npz))
    allocate(rhs_u(1:npz))
    allocate(rhs_v(1:npz))
    allocate(u_new(1:npz))
    allocate(v_new(1:npz))
    allocate(z_center(1:npz))
    allocate(dz_layer(1:npz))
    allocate(mass_per_area(1:npz))
    allocate(rho_layer(1:npz))
    allocate(cp_layer(1:npz))
    allocate(dynamic_viscosity(1:npz))
    allocate(taper_factor(1:npz))
    allocate(valid(1:npz))
    allocate(ke_heat_work(1:npz))

    if (npz > 1) then
      allocate(interface_distance(1:npz-1))
      allocate(interface_mu(1:npz-1))
    else
      allocate(interface_distance(1:1))
      allocate(interface_mu(1:1))
    end if

    do j = js, je
      do i = is, ie

        lower(:) = 0.0_WORK_KIND
        diagonal(:) = 1.0_WORK_KIND
        upper(:) = 0.0_WORK_KIND
        rhs_u(:) = 0.0_WORK_KIND
        rhs_v(:) = 0.0_WORK_KIND
        u_new(:) = 0.0_WORK_KIND
        v_new(:) = 0.0_WORK_KIND
        z_center(:) = 0.0_WORK_KIND
        dz_layer(:) = 0.0_WORK_KIND
        mass_per_area(:) = 0.0_WORK_KIND
        rho_layer(:) = 0.0_WORK_KIND
        cp_layer(:) = 0.0_WORK_KIND
        dynamic_viscosity(:) = 0.0_WORK_KIND
        taper_factor(:) = 0.0_WORK_KIND
        interface_distance(:) = 0.0_WORK_KIND
        interface_mu(:) = 0.0_WORK_KIND
        valid(:) = .false.
        ke_heat_work(:) = 0.0_WORK_KIND

        ! Diagnose layer properties and the frozen molecular viscosity.
        do k = 1, npz
          if (ieee_is_finite(ua(i,j,k))) then
            rhs_u(k) = real(ua(i,j,k), kind=WORK_KIND)
          end if
          if (ieee_is_finite(va(i,j,k))) then
            rhs_v(k) = real(va(i,j,k), kind=WORK_KIND)
          end if

          valid(k) = ieee_is_finite(ua(i,j,k)) .and. &
                     ieee_is_finite(va(i,j,k)) .and. &
                     ieee_is_finite(pe(i,k,j)) .and. pe(i,k,j) > MIN_PRESSURE .and. &
                     ieee_is_finite(pe(i,k+1,j)) .and. pe(i,k+1,j) > MIN_PRESSURE .and. &
                     ieee_is_finite(gz(i,j,k)) .and. &
                     ieee_is_finite(gz(i,j,k+1)) .and. &
                     ieee_is_finite(lambda_mlt(i,j,k)) .and. &
                     lambda_mlt(i,j,k) > MIN_LAMBDA .and. &
                     ieee_is_finite(rho_mlt(i,j,k)) .and. &
                     rho_mlt(i,j,k) > MIN_RHO .and. &
                     ieee_is_finite(cp_mlt(i,j,k)) .and. &
                     cp_mlt(i,j,k) > MIN_CP .and. &
                     ieee_is_finite(pt(i,j,k)) .and. &
                     ieee_is_finite(pkz(i,j,k))

          if (.not. valid(k)) cycle

          dz_layer(k) = abs(real(gz(i,j,k), kind=WORK_KIND) - &
                            real(gz(i,j,k+1), kind=WORK_KIND)) / GRAV0
          z_center(k) = 0.5_WORK_KIND * &
                        (real(gz(i,j,k), kind=WORK_KIND) + &
                         real(gz(i,j,k+1), kind=WORK_KIND)) / GRAV0

          if (.not. ieee_is_finite(dz_layer(k)) .or. dz_layer(k) < MIN_DZ) then
            valid(k) = .false.
            cycle
          end if

          t_geos = real(pt(i,j,k), kind=WORK_KIND) * &
                   real(pkz(i,j,k), kind=WORK_KIND)
          if (.not. ieee_is_finite(t_geos) .or. t_geos <= 0.0) then
            valid(k) = .false.
            cycle
          end if

          lambda_here = real(lambda_mlt(i,j,k), kind=WORK_KIND) * &
                        t_geos**0.69_WORK_KIND
          if (.not. ieee_is_finite(lambda_here) .or. lambda_here <= MIN_LAMBDA) then
            valid(k) = .false.
            cycle
          end if

          rho_layer(k) = real(rho_mlt(i,j,k), kind=WORK_KIND)
          cp_layer(k) = real(cp_mlt(i,j,k), kind=WORK_KIND)
          mass_per_area(k) = rho_layer(k) * dz_layer(k)

          if (.not. ieee_is_finite(mass_per_area(k)) .or. mass_per_area(k) <= 0.0) then
            valid(k) = .false.
            cycle
          end if

          ! mu = rho*nu = Pr*lambda/cp. Computing dynamic viscosity directly
          ! avoids overflow when nu becomes very large at low density.
          dynamic_viscosity(k) = pr_work * lambda_here / cp_layer(k)
          if (.not. ieee_is_finite(dynamic_viscosity(k)) .or. &
              dynamic_viscosity(k) <= 0.0) then
            valid(k) = .false.
            cycle
          end if

          p_layer = sqrt(real(pe(i,k,j), kind=WORK_KIND) * &
                         real(pe(i,k+1,j), kind=WORK_KIND))
          if (.not. ieee_is_finite(p_layer) .or. p_layer <= MIN_PRESSURE) then
            valid(k) = .false.
            cycle
          end if

          if (pmax_work > MIN_PRESSURE) then
            log_ratio = log(p_layer / max(pmax_work, MIN_PRESSURE))
            taper_factor(k) = 0.5_WORK_KIND * &
                              (1.0_WORK_KIND - tanh(log_ratio / TAPER_DLOGP))
          else
            taper_factor(k) = 1.0_WORK_KIND
          end if

          if (.not. ieee_is_finite(taper_factor(k)) .or. &
              taper_factor(k) < MIN_TAPER) then
            taper_factor(k) = 0.0_WORK_KIND
          end if
        end do

        ! Assemble the conservative frozen-coefficient backward-Euler system.
        ! Each interface coefficient is shared by its two adjacent layers.
        do k = 1, npz-1
          if (.not. valid(k) .or. .not. valid(k+1)) cycle

          interface_distance(k) = abs(z_center(k+1) - z_center(k))
          if (.not. ieee_is_finite(interface_distance(k)) .or. &
              interface_distance(k) < MIN_DZ) cycle

          interface_taper = min(taper_factor(k), taper_factor(k+1))
          if (interface_taper < MIN_TAPER) cycle

          interface_viscosity = interface_taper * &
                                0.5_WORK_KIND * (dynamic_viscosity(k) + &
                                       dynamic_viscosity(k+1))
          if (.not. ieee_is_finite(interface_viscosity) .or. &
              interface_viscosity <= 0.0) cycle

          interface_mu(k) = interface_viscosity

          coupling = dt_work * interface_viscosity / &
                     (interface_distance(k) * mass_per_area(k))
          if (.not. ieee_is_finite(coupling) .or. coupling < 0.0) cycle
          diagonal(k) = diagonal(k) + coupling
          upper(k) = upper(k) - coupling

          coupling = dt_work * interface_viscosity / &
                     (interface_distance(k) * mass_per_area(k+1))
          if (.not. ieee_is_finite(coupling) .or. coupling < 0.0) then
            diagonal(k) = diagonal(k) + upper(k)
            upper(k) = 0.0_WORK_KIND
            interface_mu(k) = 0.0_WORK_KIND
            cycle
          end if
          diagonal(k+1) = diagonal(k+1) + coupling
          lower(k+1) = lower(k+1) - coupling
        end do

        ! Factor the tridiagonal matrix once and apply it to both wind
        ! components. The top and bottom boundaries have zero diffusive flux.
        do k = 2, npz
          if (.not. ieee_is_finite(diagonal(k-1)) .or. &
              abs(diagonal(k-1)) <= DIAGONAL_MIN) then
            diagonal(k-1) = 1.0_WORK_KIND
            upper(k-1) = 0.0_WORK_KIND
            lower(k) = 0.0_WORK_KIND
          end if

          elimination_factor = lower(k) / diagonal(k-1)
          diagonal(k) = diagonal(k) - elimination_factor * upper(k-1)
          rhs_u(k) = rhs_u(k) - elimination_factor * rhs_u(k-1)
          rhs_v(k) = rhs_v(k) - elimination_factor * rhs_v(k-1)
        end do

        if (.not. ieee_is_finite(diagonal(npz)) .or. &
            abs(diagonal(npz)) <= DIAGONAL_MIN) then
          diagonal(npz) = 1.0_WORK_KIND
        end if

        u_new(npz) = rhs_u(npz) / diagonal(npz)
        v_new(npz) = rhs_v(npz) / diagonal(npz)

        do k = npz-1, 1, -1
          if (.not. ieee_is_finite(diagonal(k)) .or. &
              abs(diagonal(k)) <= DIAGONAL_MIN) then
            diagonal(k) = 1.0_WORK_KIND
            upper(k) = 0.0_WORK_KIND
          end if

          u_new(k) = (rhs_u(k) - upper(k) * u_new(k+1)) / diagonal(k)
          v_new(k) = (rhs_v(k) - upper(k) * v_new(k+1)) / diagonal(k)
        end do

        do k = 1, npz
          if (valid(k) .and. ieee_is_finite(u_new(k)) .and. &
              ieee_is_finite(v_new(k))) then
            u_tend(i,j,k) = real((u_new(k) - &
                                real(ua(i,j,k), kind=WORK_KIND)) / dt_work, &
                                kind=kind(u_tend(i,j,k)))
            v_tend(i,j,k) = real((v_new(k) - &
                                real(va(i,j,k), kind=WORK_KIND)) / dt_work, &
                                kind=kind(v_tend(i,j,k)))
          end if
        end do

        ! Diagnose viscous kinetic-energy loss from the implicitly solved
        ! wind gradients. Half of each interface loss is assigned to each
        ! adjacent layer. This remains optional in dyn_core.F90.
        do k = 1, npz-1
          if (interface_mu(k) <= 0.0) cycle
          if (interface_distance(k) < MIN_DZ) cycle
          if (.not. valid(k) .or. .not. valid(k+1)) cycle
          if (.not. ieee_is_finite(u_new(k)) .or. &
              .not. ieee_is_finite(u_new(k+1)) .or. &
              .not. ieee_is_finite(v_new(k)) .or. &
              .not. ieee_is_finite(v_new(k+1))) cycle

          grad_u = (u_new(k+1) - u_new(k)) / interface_distance(k)
          grad_v = (v_new(k+1) - v_new(k)) / interface_distance(k)
          dissipation_rate = interface_mu(k) * (grad_u*grad_u + grad_v*grad_v)

          if (.not. ieee_is_finite(dissipation_rate) .or. &
              dissipation_rate < 0.0) cycle

          ke_heat_work(k) = ke_heat_work(k) + &
                            0.5_WORK_KIND * dissipation_rate / &
                            (rho_layer(k) * cp_layer(k))
          ke_heat_work(k+1) = ke_heat_work(k+1) + &
                              0.5_WORK_KIND * dissipation_rate / &
                              (rho_layer(k+1) * cp_layer(k+1))
        end do

        do k = 1, npz
          if (ieee_is_finite(ke_heat_work(k)) .and. ke_heat_work(k) >= 0.0_WORK_KIND) then
            ke_heat_tend(i,j,k) = real(ke_heat_work(k), &
                                       kind=kind(ke_heat_tend(i,j,k)))
          end if
        end do

      end do
    end do

    deallocate(lower, diagonal, upper, rhs_u, rhs_v, u_new, v_new)
    deallocate(z_center, dz_layer, mass_per_area, rho_layer, cp_layer)
    deallocate(dynamic_viscosity, taper_factor, interface_distance)
    deallocate(interface_mu, valid, ke_heat_work)

  end subroutine mol_mom_diff_compute_tend

end module mol_mom_diff_mod
