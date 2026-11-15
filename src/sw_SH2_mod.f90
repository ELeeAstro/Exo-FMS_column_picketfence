!!!
! Elspeth KH Lee - Jun 2024 : Initial version
!                - Nov 2026 : Rewritten following PICASO exactly, enabled
!                - Nov 2026 : Optimisations - compile-time band widths in the solver
!                - Nov 2026 : 2x2 block-tridiagonal Thomas solve (replaces the banded matrix solve)
!
! sw: Two-stream spherical harmonic method (SH2, Rooney et al. 2024), converted from PICASO
!     Pros: Fast multiple scattering method
!     Cons: Two-stream only, slower than sw_adding
!!!

module sw_SH2_mod
  use, intrinsic :: iso_fortran_env
  use sw_direct_mod, only : sw_direct
  implicit none

  !! Precision variables
  integer, parameter :: dp = REAL64

  !! Single scattering albedo cap - exactly conservative scattering (w = 1, eigenvalue k = 0) is singular.
  !! 1 - 1e-12 absorbs < 1e-8 of the incident flux even in a tau ~ 6000 column
  real(dp), parameter :: w_max = 1.0_dp - 1.0e-12_dp

  !! Required constants
  real(dp), parameter :: pi = 4.0_dp * atan(1.0_dp)
  real(dp), parameter :: twopi = 2.0_dp * pi
  real(dp), parameter :: fourpi = 4.0_dp * pi

  !! Number of streams
  integer, parameter :: nstr = 2

  !! Use Two-Term HG function (TTHG), parameters are the PICASO defaults
  !! forward fraction = frac_a + frac_b*g_back**frac_c, g_forward = c_forward*g, g_back = c_back*g
  logical, parameter :: TTHG = .False.
  real(dp), parameter :: tthg_forward = 1.0_dp, tthg_back = -0.5_dp
  real(dp), parameter :: tthg_frac_a = 1.0_dp, tthg_frac_b = -1.0_dp
  integer, parameter :: tthg_frac_c = 2

  !! Band widths of the two-stream matrix

  private :: sw_SH_two_stream, sh2_solve
  public :: sw_SH2

contains

  subroutine sw_SH2(nlay, nlev, tau_e, mu_z, Finc, ssa, gg, sw_up, sw_down, sw_net, asr)
    implicit none

    !! Input variables
    integer, intent(in) :: nlay, nlev
    real(dp), dimension(nlev), intent(in) :: tau_e
    real(dp), dimension(nlay), intent(in) :: ssa, gg
    real(dp), intent(in) :: mu_z, Finc

    !! Output variables
    real(dp), dimension(nlev), intent(out) :: sw_up, sw_down, sw_net
    real(dp), intent(out) :: asr

    if (all(ssa(:) <= 1.0e-12_dp)) then
      call sw_direct(nlev, tau_e, mu_z, Finc, sw_up, sw_down, sw_net, asr)
      return
    end if

    !! Shortwave flux calculation
    if ((mu_z > 0.0_dp) .and. (Finc > 0.0_dp)) then
      call sw_SH_two_stream(nlay, nlev, Finc, mu_z, tau_e, ssa, gg, sw_down, sw_up)
    else
      sw_down(:) = 0.0_dp
      sw_up(:) = 0.0_dp
    end if

    !! Net sw flux
    sw_net(:) = sw_up(:) - sw_down(:)

    !! Absorbed Stellar Radiation (ASR)
    asr = sw_down(1) - sw_up(1)

  end subroutine sw_SH2

  subroutine sw_SH_two_stream(nlay, nlev, F0PI, u0, tau_in, w_in, g_in, flx_down, flx_up)
    implicit none

    !! Input variables
    integer, intent(in) :: nlay, nlev
    real(dp), intent(in) :: F0PI, u0
    real(dp), dimension(nlev), intent(in) :: tau_in
    real(dp), dimension(nlay), intent(in) :: w_in, g_in

    !! Output variables
    real(dp), dimension(nlev), intent(out) :: flx_down, flx_up

    !! Work variables
    integer :: k, l
    real(dp) :: b_top, b_surface
    real(dp), dimension(nstr) :: Pu0
    real(dp), dimension(nlay) :: dtau, w0, f_deltaM, g_f, g_b, f_f
    real(dp), dimension(nstr,nlay) :: chi
    real(dp), dimension(nstr,nlay) :: w_multi, a, bb, eta
    real(dp), dimension(nlay) :: Del, lam, exptrm, q, Q1, Q2, rhs_dn, rhs_up, y1, y2
    real(dp), dimension(nlay) :: zmn, zpl, zmn_up, zpl_up, zmn_down, zpl_down
    real(dp), dimension(nlev) :: tau, expon

    !! If zero albedo across all atmospheric layers then return direct beam only
    if (all(w_in(:) <= 1.0e-12_dp)) then
      flx_down(:) = F0PI*u0 * exp(-(tau_in(:) - tau_in(1))/u0)
      flx_up(:) = 0.0_dp
      return
    end if

    !! Phase function Legendre moments chi_l (OTHG: chi_l = g**l, TTHG as PICASO)
    if (TTHG .eqv. .False.) then
      do l = 1, nstr
        chi(l,:) = g_in(:)**l
      end do
    else
      g_f(:) = tthg_forward*g_in(:)
      g_b(:) = tthg_back*g_in(:)
      f_f(:) = tthg_frac_a + tthg_frac_b*g_b(:)**tthg_frac_c
      do l = 1, nstr
        chi(l,:) = f_f(:)*g_f(:)**l + (1.0_dp - f_f(:))*g_b(:)**l
      end do
    end if

    !! Delta-M with f = chi_nstr (g**nstr for OTHG as PICASO, the TTHG moment for TTHG)
    f_deltaM(:) = chi(nstr,:)
    w0(:) = w_in(:)*(1.0_dp - f_deltaM(:))/(1.0_dp - w_in(:)*f_deltaM(:))
    w0(:) = min(w0(:), w_max)
    dtau(:) = (tau_in(2:nlev) - tau_in(1:nlay))*(1.0_dp - w_in(:)*f_deltaM(:))
    tau(1) = 0.0_dp
    do k = 1, nlay
      tau(k+1) = tau(k) + dtau(k)
    end do

    !! Legendre polynomials at -u0
    Pu0(1) = 1.0_dp
    Pu0(2) = -u0

    !! Phase function moments with delta-M
    w_multi(1,:) = 1.0_dp
    do l = 1, nstr-1
      w_multi(l+1,:) = (real(2*l+1,dp)*chi(l,:) - real(2*l+1,dp)*f_deltaM(:))/(1.0_dp - f_deltaM(:))
    end do

    do l = 0, nstr-1
      a(l+1,:) = real(2*l+1,dp) - w0(:)*w_multi(l+1,:)
      bb(l+1,:) = (F0PI*(w0(:)*w_multi(l+1,:)))*Pu0(l+1)/fourpi
    end do

    !! Boundary conditions - nothing enters at the top, gas giant lower boundary
    b_top = 0.0_dp
    b_surface = 0.0_dp

    !! Direct beam particular solution
    Del(:) = (1.0_dp/u0)**2 - a(1,:)*a(2,:)
    eta(1,:) = (bb(2,:)/u0 - a(2,:)*bb(1,:))/Del(:)
    eta(2,:) = (bb(1,:)/u0 - a(1,:)*bb(2,:))/Del(:)

    lam(:) = sqrt(a(1,:)*a(2,:))
    exptrm(:) = exp(-min(lam(:)*dtau(:), 35.0_dp))

    !! Parameters in the matrices
    q(:) = lam(:)/a(2,:)
    Q1(:) = (0.5_dp + q(:))*twopi
    Q2(:) = (0.5_dp - q(:))*twopi

    zmn(:) = (0.5_dp*eta(1,:) - eta(2,:))*twopi
    zpl(:) = (0.5_dp*eta(1,:) + eta(2,:))*twopi
    expon(:) = exp(-tau(:)/u0)
    zmn_up(:) = zmn(:)*expon(2:nlev)
    zpl_up(:) = zpl(:)*expon(2:nlev)
    zmn_down(:) = zmn(:)*expon(1:nlay)
    zpl_down(:) = zpl(:)*expon(1:nlay)

    !! Flux continuity right-hand sides at each interface (downward and upward rows)
    do k = 1, nlay-1
      rhs_dn(k) = zmn_down(k+1) - zmn_up(k)
      rhs_up(k) = zpl_down(k+1) - zpl_up(k)
    end do

    !! Block-tridiagonal solve: y1 = mode 1 at the layer top, y2 = mode 2 at the layer bottom
    !! (PICASO's X(2k) = y2 exptrm). Top: nothing enters. Bottom: nothing reflected (gas giant)
    call sh2_solve(nlay, Q1, Q2, exptrm, b_top - zmn_down(1), rhs_dn, rhs_up, &
      & Q2(nlay)*exptrm(nlay), Q1(nlay), b_surface - zpl_up(nlay), y1, y2)

    !! Level fluxes - top level from layer 1, level k+1 from the bottom of layer k
    flx_down(1) = Q1(1)*y1(1) + Q2(1)*exptrm(1)*y2(1) + zmn_down(1)
    flx_up(1) = Q2(1)*y1(1) + Q1(1)*exptrm(1)*y2(1) + zpl_down(1)
    do k = 1, nlay
      flx_down(k+1) = Q1(k)*exptrm(k)*y1(k) + Q2(k)*y2(k) + zmn_up(k)
      flx_up(k+1) = Q2(k)*exptrm(k)*y1(k) + Q1(k)*y2(k) + zpl_up(k)
    end do

    !! Add the direct beam to the downward flux
    flx_down(:) = flx_down(:) + F0PI*u0*expon(:)


  end subroutine sw_SH_two_stream

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

end module sw_SH2_mod
