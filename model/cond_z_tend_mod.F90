module cond_z_tend_mod
  implicit none
  private
  public :: cond_z_tend

contains

  subroutine cond_z_tend(T, K_tc, rho, cp, dz, dz_if, dTdt, dt, top_flux)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none

    real, intent(in)  :: T(:,:,:)
    real, intent(in)  :: K_tc(:,:,:)
    real, intent(in)  :: rho(:,:,:)
    real, intent(in)  :: cp(:,:,:)
    real, intent(in)  :: dz(:,:,:)
    real, intent(in)  :: dz_if(:,:,:)
    real, intent(out) :: dTdt(:,:,:)
    real, intent(in)  :: dt
    real, intent(in), optional :: top_flux(:,:)

    integer :: i, j, kk
    integer :: is, ie, js, je, ks, ke
    integer :: nk

    real, allocatable :: lower(:)
    real, allocatable :: diagonal(:)
    real, allocatable :: upper(:)
    real, allocatable :: rhs(:)
    real, allocatable :: T_new(:)
    real, allocatable :: heat_capacity(:)
    logical, allocatable :: valid(:)

    real :: interface_distance
    real :: interface_conductance
    real :: coupling
    real :: elimination_factor
    real :: applied_top_flux

    real, parameter :: dz_min = 100.0
    real, parameter :: diagonal_min = 1.0e-20

    is = lbound(T,1); ie = ubound(T,1)
    js = lbound(T,2); je = ubound(T,2)
    ks = lbound(T,3); ke = ubound(T,3)
    nk = ke - ks + 1

    dTdt(:,:,:) = 0.0

    if (.not. ieee_is_finite(dt) .or. dt <= 0.0 .or. nk <= 0) then
      return
    end if

    allocate(lower(ks:ke))
    allocate(diagonal(ks:ke))
    allocate(upper(ks:ke))
    allocate(rhs(ks:ke))
    allocate(T_new(ks:ke))
    allocate(heat_capacity(ks:ke))
    allocate(valid(ks:ke))

    do j = js, je
      do i = is, ie

        lower(:) = 0.0
        diagonal(:) = 1.0
        upper(:) = 0.0
        rhs(:) = 0.0
        T_new(:) = 0.0
        heat_capacity(:) = 0.0
        valid(:) = .false.

        do kk = ks, ke
          rhs(kk) = T(i,j,kk)

          valid(kk) = ieee_is_finite(T(i,j,kk)) .and. &
                      ieee_is_finite(K_tc(i,j,kk)) .and. K_tc(i,j,kk) >= 0.0 .and. &
                      ieee_is_finite(rho(i,j,kk)) .and. rho(i,j,kk) > 0.0 .and. &
                      ieee_is_finite(cp(i,j,kk)) .and. cp(i,j,kk) > 0.0 .and. &
                      ieee_is_finite(dz(i,j,kk)) .and. abs(dz(i,j,kk)) >= dz_min

          if (valid(kk)) then
            heat_capacity(kk) = rho(i,j,kk) * cp(i,j,kk) * abs(dz(i,j,kk))
            valid(kk) = ieee_is_finite(heat_capacity(kk)) .and. heat_capacity(kk) > 0.0
          end if
        end do

        ! Assemble the frozen-coefficient backward-Euler diffusion operator.
        do kk = ks, ke-1
          if (.not. valid(kk) .or. .not. valid(kk+1)) cycle

          interface_distance = abs(dz_if(i,j,kk))
          if (.not. ieee_is_finite(interface_distance) .or. &
              interface_distance < dz_min) cycle

          interface_conductance = 0.5 * (K_tc(i,j,kk) + K_tc(i,j,kk+1)) / &
                                  interface_distance
          if (.not. ieee_is_finite(interface_conductance) .or. &
              interface_conductance < 0.0) cycle

          coupling = dt * interface_conductance / heat_capacity(kk)
          diagonal(kk) = diagonal(kk) + coupling
          upper(kk) = upper(kk) - coupling

          coupling = dt * interface_conductance / heat_capacity(kk+1)
          diagonal(kk+1) = diagonal(kk+1) + coupling
          lower(kk+1) = lower(kk+1) - coupling
        end do

        if (valid(ks) .and. present(top_flux)) then
          applied_top_flux = top_flux(i,j)
          if (ieee_is_finite(applied_top_flux)) then
            rhs(ks) = rhs(ks) + dt * applied_top_flux / heat_capacity(ks)
          end if
        end if

        ! Thomas algorithm for the tridiagonal backward-Euler system.
        do kk = ks+1, ke
          if (abs(diagonal(kk-1)) <= diagonal_min) then
            diagonal(kk-1) = 1.0
            lower(kk) = 0.0
          end if

          elimination_factor = lower(kk) / diagonal(kk-1)
          diagonal(kk) = diagonal(kk) - elimination_factor * upper(kk-1)
          rhs(kk) = rhs(kk) - elimination_factor * rhs(kk-1)
        end do

        if (abs(diagonal(ke)) <= diagonal_min) diagonal(ke) = 1.0
        T_new(ke) = rhs(ke) / diagonal(ke)

        do kk = ke-1, ks, -1
          if (abs(diagonal(kk)) <= diagonal_min) diagonal(kk) = 1.0
          T_new(kk) = (rhs(kk) - upper(kk) * T_new(kk+1)) / diagonal(kk)
        end do

        do kk = ks, ke
          if (valid(kk) .and. ieee_is_finite(T_new(kk))) then
            dTdt(i,j,kk) = (T_new(kk) - T(i,j,kk)) / dt
          else
            dTdt(i,j,kk) = 0.0
          end if
        end do

      end do
    end do

    deallocate(lower, diagonal, upper, rhs, T_new, heat_capacity, valid)

  end subroutine cond_z_tend
end module cond_z_tend_mod
