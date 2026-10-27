! Thermal conduction +++ geos_mlt

module thermcond_mod
  implicit none
  private
  public :: tc_calc
  public :: tc_point

  real, parameter :: kB  = 1.380649e-23
  real, parameter :: amu = 1.66053906660e-27

  real, parameter :: mO  = 16.0 * amu
  real, parameter :: mO2 = 32.0 * amu
  real, parameter :: mN2 = 28.0 * amu

contains

  !------------------------------------------------------------
  ! Compute thermal-conduction properties at one point.
  !------------------------------------------------------------
  subroutine tc_point(Tloc, nOloc, nO2loc, nN2loc, Kloc, alpha_loc)
    implicit none

    real, intent(in)  :: Tloc
    real, intent(in)  :: nOloc, nO2loc, nN2loc
    real, intent(out) :: Kloc, alpha_loc

    real :: rho, rhoO, rhoO2, rhoN2
    real :: phiO, phiO2, phiN2
    real :: xO, xO2, xN2, ntot
    real :: cp

    Kloc     = 0.0
    alpha_loc = 0.0

    ! If temperature is non-physical, turn off conduction here.
    if (Tloc <= 0.0) return

    rhoO  = nOloc  * mO
    rhoO2 = nO2loc * mO2
    rhoN2 = nN2loc * mN2
    rho   = rhoO + rhoO2 + rhoN2

    ! If rho is zero, no conduction.
    if (rho <= 0.0) return

    ! Use mass fractions for Cp because Cp is mass-specific.
    phiO  = rhoO  / rho
    phiO2 = rhoO2 / rho
    phiN2 = rhoN2 / rho

    cp =  2.5 * (kB/mO)  * phiO  &
        + 3.5 * (kB/mO2) * phiO2 &
        + 3.5 * (kB/mN2) * phiN2

    ! Use number fractions for thermal conductivity because K is a
    ! molecular transport coefficient.
    ntot = nOloc + nO2loc + nN2loc
    if (ntot <= 0.0) return

    xO  = nOloc  / ntot
    xO2 = nO2loc / ntot
    xN2 = nN2loc / ntot

    Kloc = (75.9*xO + 56.0*(xO2 + xN2)) * (Tloc**0.69) * 1.0e-5
    alpha_loc = Kloc / max(1.0e-30, (rho * cp))

  end subroutine tc_point


  !------------------------------------------------------------
  ! Compute thermal-conduction properties on the 3-D grid.
  !------------------------------------------------------------
  subroutine tc_calc(T, nO, nO2, nN2, K_tc, alpha)
    implicit none

    real, intent(in)  :: T(:,:,:)
    real, intent(in)  :: nO(:,:,:), nO2(:,:,:), nN2(:,:,:)
    real, intent(out) :: K_tc(:,:,:), alpha(:,:,:)

    integer :: i, j, kk
    integer :: is, ie, js, je, ks, ke

    is = lbound(T,1); ie = ubound(T,1)
    js = lbound(T,2); je = ubound(T,2)
    ks = lbound(T,3); ke = ubound(T,3)

    do kk = ks, ke
      do j = js, je
        do i = is, ie
          call tc_point(T(i,j,kk), nO(i,j,kk), nO2(i,j,kk), nN2(i,j,kk), &
                        K_tc(i,j,kk), alpha(i,j,kk))
        end do
      end do
    end do

  end subroutine tc_calc

end module thermcond_mod
