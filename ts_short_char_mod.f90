!!!
! Elspeth KH Lee - May 2021
! Two-stream method following the short characteristics method (e.g. Helios-r2: Kitzmann et al. 2018)
! Uses the method of short characteristics (Olson & Kunasz 1987) with linear interpolants.
! Pros: Very fast, accurate at high optical depths, very stable
! Cons: No scattering
!!!

module ts_short_char_mod
  use, intrinsic :: iso_fortran_env
  implicit none

  !! Precision variables
  integer, parameter :: dp = REAL64

  !! Required constants
  real(dp), parameter :: pi = 4.0_dp * atan(1.0_dp)
  real(dp), parameter :: twopi = 2.0_dp * pi
  real(dp), parameter :: sb = 5.670374419e-8_dp

  !! Gauss quadrature variables, cosine angle values (uarr] and weights (w)
  !! here you can comment in/out groups of mu values for testing
  !! make sure to make clean and recompile if you change these

  !! single angle diffusion factor approximation - typically 1/1.66
  !integer, parameter :: nmu = 1
  !real(dp), dimension(nmu), parameter :: uarr = (/1.0_dp/1.66_dp/)
  !real(dp), dimension(nmu), parameter :: w = (/1.0_dp/)
  !real(dp), dimension(nmu), parameter :: wuarr = uarr * w


  !! Legendre quadrature for 2 nodes
  ! integer, parameter :: nmu = 2
  ! real(dp), dimension(nmu), parameter :: uarr = (/0.21132487_dp, 0.78867513_dp/)
  ! real(dp), dimension(nmu), parameter :: w = (/0.5_dp, 0.5_dp/)
  ! real(dp), dimension(nmu), parameter :: wuarr = uarr * w

  !! Lacis & Oinas (1991) 3 point numerical values - Does not work somehow, e-mail me if you know why :)
  ! integer, parameter :: nmu = 3
  ! real(dp), dimension(nmu), parameter :: uarr = (/0.1_dp, 0.5_dp, 1.0_dp/)
  ! real(dp), dimension(nmu), parameter :: wuarr = (/0.0433_dp, 0.5742_dp, 0.3825_dp/)

  !! Legendre quadrature for 4 nodes
  ! integer, parameter :: nmu = 4
  ! real(dp), dimension(nmu), parameter :: uarr = &
  !   & (/0.06943184_dp, 0.33000948_dp, 0.66999052_dp, 0.93056816_dp/)
  ! real(dp), dimension(nmu), parameter :: w = &
  !   & (/0.17392742_dp, 0.32607258_dp, 0.32607258_dp, 0.17392742_dp/)
  ! real(dp), dimension(nmu), parameter :: wuarr = uarr * w

  !! 5 point EGP quadrature values
  integer, parameter :: nmu = 5
  real(dp), dimension(nmu), parameter :: uarr = &
    &(/0.0985350858_dp, 0.3045357266_dp, 0.5620251898_dp, 0.8019865821_dp, 0.9601901429_dp/)
  real(dp), dimension(nmu), parameter :: wuarr = &
    & (/0.0157479145_dp, 0.0739088701_dp, 0.1463869871_dp, 0.1671746381_dp, 0.0967815902_dp/)

  public :: ts_short_char
  private :: lw_grey_updown_linear, sw_grey_down, linear_log_interp

contains

  subroutine ts_short_char(nlay, nlev, Tl, pl, pe, tau_Ve, tau_IRe, mu_z, F0, Tint, AB, Beta_V, Beta_IR, net_F)
    implicit none

    !! Input variables
    integer, intent(in) :: nlay, nlev
    real(dp), dimension(nlay), intent(in) :: Tl, pl
    real(dp), dimension(nlev), intent(in) :: pe
    real(dp), dimension(3,nlev), intent(in) :: tau_Ve
    real(dp), dimension(2,nlev), intent(in) :: tau_IRe
    real(dp), dimension(3), intent(in) :: Beta_V
    real(dp), dimension(2), intent(in) :: Beta_IR
    real(dp), intent(in) :: F0, mu_z, Tint, AB

    !! Output variables
    real(dp), dimension(nlev), intent(out) :: net_F

    !! Work variables
    integer :: i, b
    real(dp) :: Finc, be_int, Finc_b, be_int_b
    real(dp), dimension(nlev) :: Te, be, be_b
    real(dp), dimension(nlev) :: sw_down, sw_up, lw_down, lw_up
    real(dp), dimension(3,nlev) :: sw_down_b, sw_up_b
    real(dp), dimension(2,nlev) :: lw_down_b, lw_up_b
    real(dp), dimension(nlev) :: lw_net, sw_net

    !! Find temperature at layer edges through linear interpolation and extrapolation
    do i = 2, nlay
      call linear_log_interp(pe(i), pl(i-1), pl(i), Tl(i-1), Tl(i), Te(i))
      !print*, i, pl(i), pl(i-1), pe(i), Tl(i-1), Tl(i), Te(i)
    end do
    Te(1) = Tl(1) + (pe(1) - pe(2))/(pl(1) - pe(2)) * (Tl(1) - Te(2))
    Te(nlev) = Tl(nlay) + (pe(nlev) - pe(nlay))/(pl(nlay) - pe(nlay)) * (Tl(nlay) - Te(nlay))

    !! Shortwave flux calculation
    if (mu_z > 0.0_dp) then
      Finc = (1.0_dp - AB) * F0
      sw_down(:) = 0.0_dp
      do b = 1, 3
        Finc_b = Finc * Beta_V(b)
        call sw_grey_down(nlev, Finc_b, tau_Ve(b,:), mu_z, sw_down_b(b,:))
        sw_down(:) = sw_down(:) + sw_down_b(b,:)
      end do
    else
      sw_down(:) = 0.0_dp
    end if
    sw_up(:) = 0.0_dp ! sw_up is zero since we don't have shortwave scattering in this mode

    !! Longwave two-stream flux calculation
    be(:) = (sb * Te(:)**4)/pi  ! Integrated planck function intensity at levels
    be_int = (sb * Tint**4)/pi ! Integrated planck function intensity for internal temperature
    lw_up(:) = 0.0_dp
    lw_down(:) = 0.0_dp
    do b = 1, 2
      be_b(:) = be(:) * Beta_IR(b)
      be_int_b = be_int * Beta_IR(b)
      call lw_grey_updown_linear(nlay, nlev, be_b, be_int_b, tau_IRe(b,:), lw_up_b(b,:), lw_down_b(b,:))
      lw_up(:) = lw_up(:) + lw_up_b(b,:)
      lw_down(:) = lw_down(:) + lw_down_b(b,:)
    end do

    !! Net fluxes at each level
    lw_net(:) = lw_up(:) - lw_down(:)
    sw_net(:) = sw_up(:) - sw_down(:)
    net_F(:) = lw_net(:) + sw_net(:)

  end subroutine ts_short_char

  subroutine lw_grey_updown_linear(nlay, nlev, be, be_int, tau_IRe, lw_up, lw_down)
    implicit none

    !! Input variables
    integer, intent(in) :: nlay, nlev
    real(dp), dimension(nlev), intent(in) :: be, tau_IRe
    real(dp), intent(in) :: be_int

    !! Output variables
    real(dp), dimension(nlev), intent(out) :: lw_up, lw_down

    !! Work variables and arrays
    integer :: k, m
    real(dp), dimension(nlay) :: dtau, edel
    real(dp) :: del, e0i, e1i, e1i_del
    real(dp), dimension(nlay) :: Am, Bm, Gp, Bp
    real(dp), dimension(nlev) :: lw_up_g, lw_down_g

    !! Calculate dtau in each layer
    do k = 1, nlay
      dtau(k) = tau_IRe(k+1) - tau_IRe(k)
    end do

    ! Zero the total flux arrays
    lw_up(:) = 0.0_dp
    lw_down(:) = 0.0_dp

    !! Start loops to integrate in mu space
    do m = 1, nmu

      !! Prepare loop
      do k = 1, nlay
        ! Olson & Kunasz (1987) linear interpolant parameters
        del = dtau(k)/uarr(m)
        edel(k) = exp(-del)
        e0i = 1.0_dp - edel(k)
        e1i = del - e0i

        e1i_del = e1i/del ! The equivalent to the linear in tau term

        if (dtau(k) < 1.0e-6_dp) then
          ! If we are in very low optical depth regime, then use an isothermal approximation
          Am(k) = (0.5_dp*(be(k+1) + be(k)) * e0i)/be(k)
          Bm(k) = 0.0_dp
          Gp(k) = 0.0_dp
          Bp(k) = Am(k)
        else
          Am(k) = e0i - e1i_del ! Am(k) = Gp(k), just indexed differently
          Bm(k) = e1i_del ! Bm(k) = Bp(k), just indexed differently
          Gp(k) = Am(k)
          Bp(k) = Bm(k)
        end if
      end do

      !! Begin two-stream loops
      !! Perform downward loop first
      ! Top boundary condition - 0 flux downward from top boundary
      lw_down_g(1) = 0.0_dp
      do k = 1, nlay
        lw_down_g(k+1) = lw_down_g(k)*edel(k) + Am(k)*be(k) + Bm(k)*be(k+1) ! TS intensity
      end do

      !! Perform upward loop
      ! Lower boundary condition - internal heat definition Fint = F_down - F_up
      ! here we use the same condition but use intensity units to be consistent
      lw_up_g(nlev) = lw_down_g(nlev) + be_int
      do k = nlay, 1, -1
        lw_up_g(k) = lw_up_g(k+1)*edel(k) + Bp(k)*be(k) + Gp(k)*be(k+1) ! TS intensity
      end do

      !! Sum up flux arrays with Gaussian quadrature weights and points for this mu stream
      lw_down(:) = lw_down(:) + lw_down_g(:) * wuarr(m)
      lw_up(:) = lw_up(:) + lw_up_g(:) * wuarr(m)

    end do

    !! The flux is the intensity * 2pi
    lw_down(:) = twopi * lw_down(:)
    lw_up(:) = twopi * lw_up(:)

  end subroutine lw_grey_updown_linear

  subroutine sw_grey_down(nlev, Finc, tau_V, mu_z, sw_down)
    implicit none

    !! Input
    integer, intent(in) :: nlev
    real(dp), intent(in) :: Finc, mu_z
    real(dp), dimension(nlev), intent(in) :: tau_V

    !! Output
    real(dp), dimension(nlev), intent(out) :: sw_down

    sw_down(:) = Finc * mu_z * exp(-tau_V(:)/mu_z)

  end subroutine sw_grey_down

  ! Perform linear interpolation in log10 space
  subroutine linear_log_interp(xval, x1, x2, y1, y2, yval)
    implicit none

    real(dp), intent(in) :: xval, y1, y2, x1, x2
    real(dp) :: lxval, ly1, ly2, lx1, lx2
    real(dp), intent(out) :: yval
    real(dp) :: norm

    lxval = log10(xval)
    lx1 = log10(x1); lx2 = log10(x2)
    ly1 = log10(y1); ly2 = log10(y2)

    norm = 1.0_dp / (lx2 - lx1)

    yval = 10.0_dp**((ly1 * (lx2 - lxval) + ly2 * (lxval - lx1)) * norm)

  end subroutine linear_log_interp

end module ts_short_char_mod
