module cond_z_tend_mod
  implicit none
  private
  public :: cond_z_tend

contains

  subroutine cond_z_tend(T, K_tc, rho, cp, dz, dz_if, dTdt, top_flux)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none

    real, intent(in)  :: T(:,:,:)
    real, intent(in)  :: K_tc(:,:,:)
    real, intent(in)  :: rho(:,:,:)
    real, intent(in)  :: cp(:,:,:)
    real, intent(in)  :: dz(:,:,:)
    real, intent(in)  :: dz_if(:,:,:)  ! Changed back to intent(in)
    real, intent(out) :: dTdt(:,:,:)
    real, intent(in), optional :: top_flux(:,:)

    integer :: i, j, kk
    integer :: is, ie, js, je, ks, ke

    real :: F_up, F_dn
    real :: K_if, dzif, dzk
    real :: rhok, cpk

    real, parameter :: eps = 1.0e-6
    real, parameter :: dz_min = 100.0  ! Minimum physical layer thickness
    logical :: ok_here, ok_above, ok_below

    is = lbound(T,1); ie = ubound(T,1)
    js = lbound(T,2); je = ubound(T,2)
    ks = lbound(T,3); ke = ubound(T,3)

    dTdt(:,:,:) = 0.0

    do j = js, je
      do i = is, ie

        ! Top layer (kk=ks)
        kk   = ks
        F_up = 0.0
        if (present(top_flux)) F_up = top_flux(i,j)
        F_dn = F_up

        rhok = rho(i,j,kk)
        cpk  = cp(i,j,kk)

        ok_here = ieee_is_finite(T(i,j,kk)) .and. &
                  ieee_is_finite(K_tc(i,j,kk)) .and. (K_tc(i,j,kk) >= 0.0) .and. &
                  ieee_is_finite(rhok) .and. (rhok > 0.0) .and. &
                  ieee_is_finite(cpk)  .and. (cpk  > 0.0)

        if (ok_here) then
          dzk = abs(dz(i,j,kk))
          
          ! Safety check for tiny layers
          if (dzk < dz_min) then
            dTdt(i,j,kk) = 0.0
            if (i==1 .and. j==1) then
              print *, 'WARNING: Top layer too thin at kk=', kk, ' dzk=', dzk
            end if
            cycle
          end if

          if (kk < ke) then
            dzif = 0.5 * (dz(i,j,kk) + dz(i,j,kk+1))  ! Distance between centers
            
            ! Safety check
            if (dzif < dz_min) then
              dTdt(i,j,kk) = 0.0
              if (i==1 .and. j==1) then
                print *, 'WARNING: Interface distance too small at kk=', kk, ' dzif=', dzif
              end if
              cycle
            end if
            
            ok_below = ieee_is_finite(T(i,j,kk+1)) .and. &
                       ieee_is_finite(K_tc(i,j,kk+1)) .and. (K_tc(i,j,kk+1) >= 0.0)

            if (ok_below) then
              K_if = 0.5*(K_tc(i,j,kk) + K_tc(i,j,kk+1))
              F_dn = -K_if * (T(i,j,kk+1) - T(i,j,kk)) / dzif
            else
              F_dn = F_up
            end if
          end if

          dTdt(i,j,kk) = -(F_dn - F_up) / (rhok*cpk*dzk)
          
          ! Diagnostic for first column, top layer
!          if (i==1 .and. j==1) then
!             print *, '=== Top layer diagnostics (kk=', kk, ') ==='
!             print *, '  T(kk), T(kk+1):', T(i,j,kk), T(i,j,kk+1)
!             print *, '  dz(kk), dz(kk+1):', dz(i,j,kk), dz(i,j,kk+1)
!             print *, '  dzk:', dzk
!             print *, '  dzif:', dzif
!             print *, '  K_tc(kk), K_tc(kk+1):', K_tc(i,j,kk), K_tc(i,j,kk+1)
!             print *, '  K_if:', K_if
!             print *, '  rho, cp:', rhok, cpk
!             print *, '  F_up, F_dn:', F_up, F_dn
!             print *, '  (F_dn - F_up):', (F_dn - F_up)
!             print *, '  (rhok*cpk*dzk):', (rhok*cpk*dzk)
!             print *, '  dTdt(i,j,kk):', dTdt(i,j,kk)
!             print *, '========================================'
!          end if
          
          if (.not. ieee_is_finite(dTdt(i,j,kk))) dTdt(i,j,kk) = 0.0
        else
          dTdt(i,j,kk) = 0.0
          if (i==1 .and. j==1) then
            print *, 'WARNING: Top layer failed ok_here check at kk=', kk
            print *, '  T, K_tc, rho, cp:', T(i,j,kk), K_tc(i,j,kk), rhok, cpk
          end if
        end if

        ! Interior layers (ks+1 .. ke-1)
        do kk = ks+1, ke-1

          rhok = rho(i,j,kk)
          cpk  = cp(i,j,kk)

          ok_here = ieee_is_finite(T(i,j,kk)) .and. &
                    ieee_is_finite(K_tc(i,j,kk)) .and. (K_tc(i,j,kk) >= 0.0) .and. &
                    ieee_is_finite(rhok) .and. (rhok > 0.0) .and. &
                    ieee_is_finite(cpk)  .and. (cpk  > 0.0)

          if (.not. ok_here) then
            dTdt(i,j,kk) = 0.0
            cycle
          end if

          dzk = abs(dz(i,j,kk))
          if (dzk < dz_min) then
            dTdt(i,j,kk) = 0.0
            cycle
          end if

          ! Flux from above interface (kk-1/kk)
          dzif = 0.5 * (dz(i,j,kk-1) + dz(i,j,kk))
          if (dzif < dz_min) then
            F_up = 0.0
          else
            ok_above = ieee_is_finite(T(i,j,kk-1)) .and. &
                       ieee_is_finite(K_tc(i,j,kk-1)) .and. (K_tc(i,j,kk-1) >= 0.0)

            if (ok_above) then
              K_if = 0.5*(K_tc(i,j,kk-1) + K_tc(i,j,kk))
              F_up = -K_if * (T(i,j,kk) - T(i,j,kk-1)) / dzif
            else
              F_up = 0.0
            end if
          end if

          ! Flux to below interface (kk/kk+1)
          dzif = 0.5 * (dz(i,j,kk) + dz(i,j,kk+1))
          if (dzif < dz_min) then
            F_dn = 0.0
          else
            ok_below = ieee_is_finite(T(i,j,kk+1)) .and. &
                       ieee_is_finite(K_tc(i,j,kk+1)) .and. (K_tc(i,j,kk+1) >= 0.0)

            if (ok_below) then
              K_if = 0.5*(K_tc(i,j,kk) + K_tc(i,j,kk+1))
              F_dn = -K_if * (T(i,j,kk+1) - T(i,j,kk)) / dzif
            else
              F_dn = 0.0
            end if
          end if

          dTdt(i,j,kk) = -(F_dn - F_up) / (rhok*cpk*dzk)
          if (.not. ieee_is_finite(dTdt(i,j,kk))) dTdt(i,j,kk) = 0.0

        end do

        ! Bottom layer (kk=ke)
        if (ke > ks) then
          kk = ke

          rhok = rho(i,j,kk)
          cpk  = cp(i,j,kk)

          ok_here = ieee_is_finite(T(i,j,kk)) .and. &
                    ieee_is_finite(K_tc(i,j,kk)) .and. (K_tc(i,j,kk) >= 0.0) .and. &
                    ieee_is_finite(rhok) .and. (rhok > 0.0) .and. &
                    ieee_is_finite(cpk)  .and. (cpk  > 0.0)

          if (ok_here) then
            dzk = abs(dz(i,j,kk))
            if (dzk < dz_min) then
              dTdt(i,j,kk) = 0.0
              cycle
            end if
            
            dzif = 0.5 * (dz(i,j,kk-1) + dz(i,j,kk))
            if (dzif < dz_min) then
              F_up = 0.0
            else
              ok_above = ieee_is_finite(T(i,j,kk-1)) .and. &
                         ieee_is_finite(K_tc(i,j,kk-1)) .and. (K_tc(i,j,kk-1) >= 0.0)

              if (ok_above) then
                K_if = 0.5*(K_tc(i,j,kk-1) + K_tc(i,j,kk))
                F_up = -K_if * (T(i,j,kk) - T(i,j,kk-1)) / dzif
              else
                F_up = 0.0
              end if
            end if

            F_dn = 0.0
            dTdt(i,j,kk) = -(F_dn - F_up) / (rhok*cpk*dzk)
            if (.not. ieee_is_finite(dTdt(i,j,kk))) dTdt(i,j,kk) = 0.0
          else
            dTdt(i,j,kk) = 0.0
          end if
        end if

      end do
    end do

  end subroutine cond_z_tend
end module cond_z_tend_mod
