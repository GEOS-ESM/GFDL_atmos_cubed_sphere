! Thermal conduction +++ geos_mlt

module thermcond_mod
  implicit none
  private
  public :: tc_calc

  real, parameter :: kB  = 1.380649e-23
  real, parameter :: amu = 1.66053906660e-27

  real, parameter :: mO  = 16.0 * amu
  real, parameter :: mO2 = 32.0 * amu
  real, parameter :: mN2 = 28.0 * amu

contains

  subroutine tc_calc(T, nO, nO2, nN2, K_tc, alpha)
    implicit none
    real, intent(in)  :: T(:,:,:)
    real, intent(in)  :: nO(:,:,:), nO2(:,:,:), nN2(:,:,:)
    real, intent(out) :: K_tc(:,:,:), alpha(:,:,:)

    integer :: i, j, kk
    integer :: is, ie, js, je, ks, ke

    real :: rho, rhoO, rhoO2, rhoN2
    real :: phiO, phiO2, phiN2
    real :: cp, Tloc, Kloc

    is = lbound(T,1); ie = ubound(T,1)
    js = lbound(T,2); je = ubound(T,2)
    ks = lbound(T,3); ke = ubound(T,3)

    do kk = ks, ke
      do j = js, je
        do i = is, ie

          Tloc = T(i,j,kk)

          ! If temperature is non-physical, turn off conduction here
          if (Tloc <= 0.0) then
            K_tc(i,j,kk) = 0.0
            alpha(i,j,kk) = 0.0
            cycle
          end if

          rhoO  = nO(i,j,kk)  * mO
          rhoO2 = nO2(i,j,kk) * mO2
          rhoN2 = nN2(i,j,kk) * mN2
          rho   = rhoO + rhoO2 + rhoN2

          ! If rho is zero, no conduction
          if (rho <= 0.0) then
            K_tc(i,j,kk)  = 0.0
            alpha(i,j,kk) = 0.0
            cycle
          end if

          phiO  = rhoO  / rho
          phiO2 = rhoO2 / rho
          phiN2 = rhoN2 / rho

          cp =  2.5 * (kB/mO)  * phiO  &
              + 3.5 * (kB/mO2) * phiO2 &
              + 3.5 * (kB/mN2) * phiN2

          Kloc = (56.0*(phiO2 + phiN2) + 75.9*phiO) * (Tloc**0.69) * 1.0e-5

          K_tc(i,j,kk)  = Kloc
          alpha(i,j,kk) = Kloc / max(1.0e-30, (rho * cp))

        end do
      end do
    end do

  end subroutine tc_calc

end module thermcond_mod


