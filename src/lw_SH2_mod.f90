!!!
! Elspeth KH Lee - Oct 2026 : Initial version
!                - Nov 2026 : Reviewed against PICASO, provenance and differences documented;
!                             Hogan (2024) quadrature for the flux integration
!                - Nov 2026 : Optimisations - beam integrals from the stored transmissions (one
!                             exponential per layer and angle), compile-time band widths in the
!                             solver, reciprocals
!                - Nov 2026 : 2x2 block-tridiagonal Thomas solve (replaces the banded matrix solve)
!
! lw: Two-stream spherical harmonic method (SH2, Rooney et al. 2024), converted from PICASO,
!     with source function technique fluxes
!     Pros: Multiple scattering with source function technique accuracy
!     Cons: Two-stream scattering
!!!

module lw_SH2_mod
  use, intrinsic :: iso_fortran_env
  use WENO4_mod, only : interpolate_weno4
  implicit none

  !! Precision variables
  integer, parameter :: dp = REAL64

  !! Required constants
  real(dp), parameter :: pi = 4.0_dp * atan(1.0_dp)
  real(dp), parameter :: twopi = 2.0_dp * pi
  real(dp), parameter :: sb = 5.670374419e-8_dp

  !! Number of streams and the Legendre order of the delta-M scaling
  integer, parameter :: nstr = 2

  real(dp), parameter :: ubari = 0.5_dp

  !! Gauss quadrature variables, cosine angle values (uarr) and weights (wuarr)
  !! here you can comment in/out groups of mu values for testing
  !! make sure to make clean and recompile if you change these

  !! Gauss–Jacobi-5 quadrature for 4 nodes - 8 stream (Hogan 2024)
  !! Hogan's weights w are for the measure 2 mu dmu (sum 1, F = pi sum w I); here the 2 pi
  !! is carried in the intensities (F = 2 pi sum wuarr I), so wuarr = w/2 (sum 0.5)
  integer, parameter :: nmu = 4
  real(dp), dimension(nmu), parameter :: uarr = (/0.0454586727_dp, 0.2322334416_dp, 0.5740198775_dp, 0.9030775973_dp/)
  real(dp), dimension(nmu), parameter :: w = (/0.0092068785_dp, 0.1285704278_dp, 0.4323381850_dp, 0.4298845087_dp/)
  real(dp), dimension(nmu), parameter :: wuarr = 0.5_dp * w
  real(dp), dimension(nmu), parameter :: ruarr = 1.0_dp/uarr

  !! Band widths of the two-stream matrix

  private :: lw_SH_two_stream, sh2_solve, e_alp, e_bet
  public :: lw_SH2

contains

  subroutine lw_SH2(nlay, nlev, Tl, pl, pe, tau_e, ssa, gg, Tint, lw_up, lw_down, lw_net, olr)
    implicit none

    !! Input variables
    integer, intent(in) :: nlay, nlev
    real(dp), dimension(nlay), intent(in) :: Tl, pl
    real(dp), dimension(nlev), intent(in) :: pe
    real(dp), dimension(nlev), intent(in) :: tau_e
    real(dp), dimension(nlay), intent(in) :: ssa, gg
    real(dp), intent(in) :: Tint

    !! Output variables
    real(dp), dimension(nlev), intent(out) :: lw_up, lw_down, lw_net
    real(dp), intent(out) :: olr

    !! Work variables
    real(dp), dimension(nlev) :: Te
    real(dp), dimension(nlev) :: be
    real(dp) :: be_int

    !! Use WENO4 method to (smoothly) interpolate layers to levels, in ln p
    Te(:) = interpolate_weno4(log(pe), log(pl), Tl, .False.)

    !! Edges are linearly extrapolated in log(p)-log(T) space.
    Te(1) = 10.0_dp**(log10(Te(2)) + (log10(pe(1)/pe(2))/log10(pl(1)/pe(2))) * &
      & log10(Tl(1)/Te(2)))
    Te(nlev) = 10.0_dp**(log10(Te(nlay)) + (log10(pe(nlev)/pe(nlay))/ &
      & log10(pl(nlay)/pe(nlay))) * log10(Tl(nlay)/Te(nlay)))

    !! Find integrated planck function for each level
    be(:) = (sb * Te(:)**4)/pi
    be_int = (sb * Tint**4)/pi

    !! Longwave flux calculation
    call lw_SH_two_stream(nlay, nlev, be(:), be_int, tau_e(:), ssa(:), gg(:), &
      & lw_up(:), lw_down(:))

    !! Net lw fluxes at each level
    lw_net(:) = lw_up(:) - lw_down(:)

    !! Outgoing Longwave Radiation (olr)
    olr = lw_up(1)

  end subroutine lw_SH2

  subroutine lw_SH_two_stream(nlay, nlev, be, be_int, tau_in, w_in, g_in, flx_up, flx_down)
    implicit none

    !! Input variables
    integer, intent(in) :: nlay, nlev
    real(dp), dimension(nlev), intent(in) :: be, tau_in
    real(dp), dimension(nlay), intent(in) :: w_in, g_in
    real(dp), intent(in) :: be_int

    !! Output variables
    real(dp), dimension(nlev), intent(out) :: flx_up, flx_down

    !! Work variables
    integer :: k, m
    real(dp) :: tautop, Btop, u, ru
    real(dp), dimension(nlay) :: dtau, w0, ff, B0, B1
    real(dp), dimension(nlay) :: w_multi0, w_multi1, a0, a1
    real(dp), dimension(nlay) :: lam, q, Q1, Q2, exptrm
    real(dp), dimension(nlay) :: zmn_down, zmn_up, zpl_down, zpl_up
    real(dp), dimension(nlay) :: X1, X2
    real(dp), dimension(nlay) :: ex, c_ev, c_od, T_up, T_dn, Ip1, Nup, Ndn, Aup, Adn, ea, eb
    real(dp), dimension(nlay) :: src_up, src_dn, el, wa0
    real(dp), dimension(nlay) :: rhs_dn, rhs_up
    real(dp), dimension(nlev) :: lw_up_g, lw_down_g

    !! Delta-M scaling (PICASO convention, f = g**nstr)
    ff(:) = g_in(:)**nstr
    w0(:) = w_in(:)*(1.0_dp - ff(:))/(1.0_dp - w_in(:)*ff(:))
    w0(:) = min(w0(:), 1.0_dp - 1.0e-6_dp)
    dtau(:) = (tau_in(2:nlev) - tau_in(1:nlay))*(1.0_dp - w_in(:)*ff(:))

    !! Phase function Legendre coefficients and the PICASO 'a' coefficients
    w_multi0(:) = 1.0_dp
    w_multi1(:) = 3.0_dp*(g_in(:) - ff(:))/(1.0_dp - ff(:))
    a0(:) = 1.0_dp - w0(:)*w_multi0(:)
    a1(:) = 3.0_dp - w0(:)*w_multi1(:)

    do k = 1, nlay
      if (dtau(k) <= 1.0e-6_dp) then
        ! For low optical depths use the isothermal approimation
        B1(k) = 0.0_dp
        B0(k) = 0.5_dp*(be(k+1) + be(k))
      else
        B1(k) = (be(k+1) - be(k))/dtau(k) ! Linear in tau term
        B0(k) = be(k)
      endif
    end do

    !! Homogeneous solution: fluxes of mode 1 (decays with depth) and mode 2 (grows)
    lam(:) = sqrt(a0(:)*a1(:))
    q(:) = lam(:)/a1(:)
    Q1(:) = twopi*(0.5_dp + q(:))
    Q2(:) = twopi*(0.5_dp - q(:))
    exptrm(:) = exp(-min(lam(:)*dtau(:), 35.0_dp))
    ! uncapped exp(-lam dtau) for the beam integrals (= exptrm unless lam dtau > 35)
    el(:) = exptrm(:)
    where (lam(:)*dtau(:) > 35.0_dp) el(:) = exp(-lam(:)*dtau(:))

    !! Particular solution fluxes at layer top (down) and bottom (up), linear in tau Planck
    zmn_down(:) = twopi*(1.0_dp - w0(:))/a0(:) * (B0(:)/2.0_dp - B1(:)/a1(:))
    zmn_up(:) = twopi*(1.0_dp - w0(:))/a0(:) * (B0(:)/2.0_dp - B1(:)/a1(:) + B1(:)*dtau(:)/2.0_dp)
    zpl_down(:) = twopi*(1.0_dp - w0(:))/a0(:) * (B0(:)/2.0_dp + B1(:)/a1(:))
    zpl_up(:) = twopi*(1.0_dp - w0(:))/a0(:) * (B0(:)/2.0_dp + B1(:)/a1(:) + B1(:)*dtau(:)/2.0_dp)

    !! Top boundary, isothermal layer above the model top (as lw_Toon89)
    tautop = dtau(1)*exp(-1.0_dp)
    Btop = (1.0_dp - exp(-tautop/ubari))*be(1)

    !! Flux continuity right-hand sides at each interface (downward and upward rows)
    do k = 1, nlay-1
      rhs_dn(k) = zmn_down(k+1) - zmn_up(k)
      rhs_up(k) = zpl_down(k+1) - zpl_up(k)
    end do

    !! Block-tridiagonal solve: X1 = mode 1 at the layer top, X2 = mode 2 at the layer bottom
    !! Top: downward flux = incident isotropic flux. Bottom: gas giant, I_up = I_down + B_int
    call sh2_solve(nlay, Q1, Q2, exptrm, pi*Btop - zmn_down(1), rhs_dn, rhs_up, &
      & (Q2(nlay) - Q1(nlay))*exptrm(nlay), Q1(nlay) - Q2(nlay), pi*be_int - zpl_up(nlay) + zmn_up(nlay), X1, X2)

    !! Particular solution first moment (constant in tau)
    Ip1(:) = (1.0_dp - w0(:))/(a0(:)*a1(:))*B1(:)
    wa0(:) = w_multi0(:)*(1.0_dp - w0(:))/a0(:)

    ! Zero the total flux arrays
    flx_up(:) = 0.0_dp
    flx_down(:) = 0.0_dp

    do m = 1, nmu

      u = uarr(m)
      ru = ruarr(m)
      ex(:) = exp(-dtau(:)*ru)

      !! Scattering phase coefficients for each mode at +u; at -u the odd terms flip
      !! mode 1 Legendre moments (1, -q), mode 2 (1, q)
      c_ev(:) = w_multi0(:)
      c_od(:) = -w_multi1(:)*u*q(:)

      ea(:) = e_alp(nlay, ru + lam(:), ex(:), el(:))
      eb(:) = e_bet(nlay, lam(:), ru, dtau(:), ex(:), el(:))

      !! Thermal particular solution integrated along the beam
      T_up(:) = u*(B0(:)*(1.0_dp - ex(:)) + B1(:)*(u - (dtau(:) + u)*ex(:)))
      T_dn(:) = u*(B0(:)*(1.0_dp - ex(:)) + B1(:)*(dtau(:) - u + u*ex(:)))
      Nup(:) = wa0(:)*T_up(:) + w_multi1(:)*u*Ip1(:)*u*(1.0_dp - ex(:))
      Ndn(:) = wa0(:)*T_dn(:) - w_multi1(:)*u*Ip1(:)*u*(1.0_dp - ex(:))

      !! Homogeneous solution integrated along the beam
      Aup(:) = X1(:)*(c_ev(:) + c_od(:))*ea(:) + X2(:)*(c_ev(:) - c_od(:))*eb(:)
      Adn(:) = X1(:)*(c_ev(:) - c_od(:))*eb(:) + X2(:)*(c_ev(:) + c_od(:))*ea(:)

      src_up(:) = twopi*ru*(w0(:)*(Aup(:) + Nup(:)) + (1.0_dp - w0(:))*T_up(:))
      src_dn(:) = twopi*ru*(w0(:)*(Adn(:) + Ndn(:)) + (1.0_dp - w0(:))*T_dn(:))

      !! Perform downward loop first
      ! Top boundary condition - intensity downward from top boundary (tautop, assumed isothermal)
      lw_down_g(1) = twopi*(1.0_dp - exp(-tautop*ru))*be(1)
      do k = 1, nlay
        lw_down_g(k+1) = lw_down_g(k)*ex(k) + src_dn(k)
      end do

      lw_up_g(nlev) = lw_down_g(nlev) + twopi*be_int
      do k = nlay, 1, -1
        lw_up_g(k) = lw_up_g(k+1)*ex(k) + src_up(k)
      end do

      !! Sum up flux arrays with Gaussian quadrature weights and points for this mu stream
      flx_down(:) = flx_down(:) + lw_down_g(:) * wuarr(m)
      flx_up(:) = flx_up(:) + lw_up_g(:) * wuarr(m)

    end do

  end subroutine lw_SH_two_stream

  !! int_0^dtau exp(-alp t) dt, alp = 1/u + lam > 0, from ex = exp(-dtau/u) and el = exp(-lam dtau)
  function e_alp(nlay, alp, ex, el) result(res)
    implicit none
    integer, intent(in) :: nlay
    real(dp), dimension(nlay), intent(in) :: alp, ex, el
    real(dp), dimension(nlay) :: res

    res(:) = (1.0_dp - ex(:)*el(:))/alp(:)

  end function e_alp

  !! (exp(-lam dtau) - exp(-dtau/u))/(1/u - lam), with the series form as 1/u - lam -> 0,
  !! from ex = exp(-dtau/u) and el = exp(-lam dtau)
  function e_bet(nlay, lam, ru, dtau, ex, el) result(res)
    implicit none
    integer, intent(in) :: nlay
    real(dp), dimension(nlay), intent(in) :: lam, dtau, ex, el
    real(dp), intent(in) :: ru
    real(dp), dimension(nlay) :: res

    integer :: k
    real(dp) :: bet, x

    do k = 1, nlay
      bet = ru - lam(k)
      x = bet*dtau(k)
      if (abs(x) < 1.0e-6_dp) then
        res(k) = dtau(k)*ex(k)*(1.0_dp + x/2.0_dp + x*x/6.0_dp)
      else
        res(k) = (el(k) - ex(k))/bet
      end if
    end do

  end function e_bet

  !! Two-stream SH system as a 2x2 block-tridiagonal Thomas solve (no pivoting).
  !! Unknowns per layer: y1 = mode 1 amplitude at the layer top, y2 = mode 2 amplitude at the
  !! layer bottom, so every coefficient is bounded (|Q| at most, e = exp(-lam dtau) <= 1).
  !! Block row k pairs the downward-flux row of the interface above layer k (or the top
  !! boundary) with the upward-flux row of the interface below it (or the bottom boundary).
  !! The diagonal blocks then have |det| >= Q1^2 - Q2^2 = 8 pi^2 q > 0.
  subroutine sh2_solve(nlay, Q1, Q2, e, rtop, rhs_dn, rhs_up, cb1, cb2, rbot, y1, y2)
    implicit none

    integer, intent(in) :: nlay
    real(dp), dimension(nlay), intent(in) :: Q1, Q2, e
    real(dp), dimension(nlay), intent(in) :: rhs_dn, rhs_up   ! interface k (1..nlay-1): down and up rows
    real(dp), intent(in) :: rtop, cb1, cb2, rbot
    real(dp), dimension(nlay), intent(out) :: y1, y2

    integer :: k
    real(dp) :: d11, d12, d21, d22, rdet, v1, v2, r1, r2, t
    real(dp), dimension(nlay) :: i11, i12, i21, i22, s1, s2

    !! Forward elimination
    do k = 1, nlay
      ! row a: top boundary or downward flux continuity at the interface above
      if (k == 1) then
        d11 = Q1(1) ; d12 = Q2(1)*e(1)
        r1 = rtop
      else
        d11 = -Q1(k) ; d12 = -Q2(k)*e(k)
        r1 = rhs_dn(k-1)
        ! eliminate Y(k-1): coefficients l = (Q1 e, Q2) of layer k-1, v = l D'^-1
        v1 = Q1(k-1)*e(k-1)*i11(k-1) + Q2(k-1)*i21(k-1)
        v2 = Q1(k-1)*e(k-1)*i12(k-1) + Q2(k-1)*i22(k-1)
        ! U(k-1) = (0; -Q2(k), -Q1(k) e(k))
        d11 = d11 + v2*Q2(k)
        d12 = d12 + v2*Q1(k)*e(k)
        r1 = r1 - v1*s1(k-1) - v2*s2(k-1)
      end if
      ! row b: upward flux continuity at the interface below or bottom boundary
      if (k == nlay) then
        d21 = cb1 ; d22 = cb2
        r2 = rbot
      else
        d21 = Q2(k)*e(k) ; d22 = Q1(k)
        r2 = rhs_up(k)
      end if
      rdet = 1.0_dp/(d11*d22 - d12*d21)
      i11(k) = d22*rdet ; i12(k) = -d12*rdet
      i21(k) = -d21*rdet ; i22(k) = d11*rdet
      s1(k) = r1 ; s2(k) = r2
    end do

    !! Back substitution, Y(k) = D'^-1 (r' - (0; u.Y(k+1))), u = (-Q2(k+1), -Q1(k+1) e(k+1))
    y1(nlay) = i11(nlay)*s1(nlay) + i12(nlay)*s2(nlay)
    y2(nlay) = i21(nlay)*s1(nlay) + i22(nlay)*s2(nlay)
    do k = nlay-1, 1, -1
      t = s2(k) + Q2(k+1)*y1(k+1) + Q1(k+1)*e(k+1)*y2(k+1)
      y1(k) = i11(k)*s1(k) + i12(k)*t
      y2(k) = i21(k)*s1(k) + i22(k)*t
    end do

  end subroutine sh2_solve

end module lw_SH2_mod
