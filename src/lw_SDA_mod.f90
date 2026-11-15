!!!
! Elspeth KH Lee - Nov 2026 : Initial version
!
! lw: Four-stream spherical harmonic adding method (delta-4SDA, Wu et al. 2016, Infrared Phys. Technol.
!     78, 254; solar counterpart Zhang & Li 2013), Planck function linear in tau as Rooney et al. (2023)
!     (linear_B = .False. for the exponential in tau form of Wu et al.), fluxes from the first moment
!     Pros: Four-stream multiple scattering accuracy, only 2x2 matrix inversions, no source function sweep
!     Cons: First moment fluxes less accurate than lw_SH4 (source function technique), mostly in clear sky
!!!

module lw_SDA_mod
  use, intrinsic :: iso_fortran_env
  use WENO4_mod, only : interpolate_weno4
  implicit none

  !! Precision variables
  integer, parameter :: dp = REAL64

  !! Required constants
  real(dp), parameter :: pi = 4.0_dp * atan(1.0_dp)
  real(dp), parameter :: twopi = 2.0_dp * pi
  real(dp), parameter :: sb = 5.670374419e-8_dp

  !! Single scattering albedo cap (conservative scattering, k2 = 0, is singular)
  real(dp), parameter :: w_max = 1.0_dp - 1.0e-6_dp

  !! Layers thinner than this use the isothermal (f0 = 0) Planck function
  real(dp), parameter :: dtau_min = 1.0e-6_dp

  !! Planck function in each layer: .True. = linear in tau (Rooney et al. 2023, Eqs. 9, 26-27),
  !! .False. = exponential in tau (Wu et al. 2016, Eq. 3)
  logical, parameter :: linear_B = .True.

  !! How close to the f0 = k resonance (exponential source = homogeneous mode) before f0 is nudged
  real(dp), parameter :: res_eps = 1.0e-6_dp

  private :: lw_SDA_method, inv2
  public :: lw_SDA

contains

  subroutine lw_SDA(nlay, nlev, Tl, pl, pe, tau_e, ssa, gg, Tint, lw_up, lw_down, lw_net, olr)
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
    real(dp), dimension(nlev) :: Te, be
    real(dp) :: be_int

    !! Use WENO4 method to (smoothly) interpolate layers to levels, in ln p
    Te(:) = interpolate_weno4(log(pe), log(pl), Tl, .False.)

    !! Edges are linearly extrapolated in log(p)-log(T) space.
    Te(1) = 10.0_dp**(log10(Te(2)) + (log10(pe(1)/pe(2))/log10(pl(1)/pe(2))) * log10(Tl(1)/Te(2)))
    Te(nlev) = 10.0_dp**(log10(Te(nlay)) + (log10(pe(nlev)/pe(nlay))/log10(pl(nlay)/pe(nlay))) * log10(Tl(nlay)/Te(nlay)))

    !! Find integrated planck function for each level
    be(:) = (sb * Te(:)**4)/pi
    be_int = (sb * Tint**4)/pi

    !! Longwave flux calculation
    call lw_SDA_method(nlay, nlev, be(:), be_int, tau_e(:), ssa(:), gg(:), lw_up(:), lw_down(:))

    !! Net lw fluxes at each level
    lw_net(:) = lw_up(:) - lw_down(:)

    !! Outgoing Longwave Radiation (olr)
    olr = lw_up(1)

  end subroutine lw_SDA

  subroutine lw_SDA_method(nlay, nlev, be, be_int, tau_in, w_in, g_in, flx_up, flx_down)
    implicit none

    !! Input variables
    integer, intent(in) :: nlay, nlev
    real(dp), dimension(nlev), intent(in) :: be, tau_in
    real(dp), dimension(nlay), intent(in) :: w_in, g_in
    real(dp), intent(in) :: be_int

    !! Output variables
    real(dp), dimension(nlev), intent(out) :: flx_up, flx_down

    !! Work variables
    integer :: k
    real(dp) :: dtau, w0, f, B0, rB, f0, b0s
    real(dp) :: om1, om2, om3, a0, a1, a2, a3, beta, gam, sq, k1, k2, e1, e2
    real(dp) :: R1, R2, P1, P2, Q1, Q2
    real(dp) :: c1, c2, c3, eta0, eta1, eta2, eta3
    real(dp), dimension(2) :: Zmt, Zpt, Zmb, Zpb, zs, zd, sv, dv
    real(dp) :: B1, fa
    real(dp), dimension(2,2) :: Am, Ap, S1, D1, S2, D2, Xm, Ym, Imat, Mi
    real(dp), dimension(2,2,nlay) :: R, T, Md
    real(dp), dimension(2,nlay) :: Eu, Ed
    real(dp), dimension(2,2,nlev) :: Rbot
    real(dp), dimension(2,nlev) :: Sup, Dn, Up

    Imat(:,:) = 0.0_dp ; Imat(1,1) = 1.0_dp ; Imat(2,2) = 1.0_dp

    !! Reflection, transmission and thermal sources of each layer (Wu et al. 2016, Sect. 2)
    do k = 1, nlay

      !! Delta-M scaling with f = g**4 (Henyey-Greenstein, chi_4)
      f = g_in(k)**4
      w0 = w_in(k)*(1.0_dp - f)/(1.0_dp - w_in(k)*f)
      w0 = min(max(w0, 0.0_dp), w_max)
      dtau = (tau_in(k+1) - tau_in(k))*(1.0_dp - w_in(k)*f)

      !! Phase function moments and a_l = (2l+1) - w omega_l
      om1 = 3.0_dp*(g_in(k) - f)/(1.0_dp - f)
      om2 = 5.0_dp*(g_in(k)**2 - f)/(1.0_dp - f)
      om3 = 7.0_dp*(g_in(k)**3 - f)/(1.0_dp - f)
      a0 = 1.0_dp - w0
      a1 = 3.0_dp - w0*om1
      a2 = 5.0_dp - w0*om2
      a3 = 7.0_dp - w0*om3

      !! Eigenvalues (Eq. 7), k2 in the cancellation free form k1^2 k2^2 = gamma
      beta = a0*a1 + 4.0_dp*a0*a3/9.0_dp + a2*a3/9.0_dp
      gam = a0*a1*a2*a3/9.0_dp
      sq = sqrt(max(beta**2 - 4.0_dp*gam, 0.0_dp))
      k1 = sqrt(0.5_dp*(beta + sq))
      k2 = sqrt(gam/(0.5_dp*(beta + sq)))
      e1 = exp(-k1*dtau)
      e2 = exp(-k2*dtau)

      !! Homogeneous mode coefficients (Eq. 7) and their hemispheric moments (Eq. 10)
      R1 = -a0/k1 ; R2 = -a0/k2
      P1 = 0.5_dp*(a0*a1/k1**2 - 1.0_dp) ; P2 = 0.5_dp*(a0*a1/k2**2 - 1.0_dp)
      Q1 = -1.5_dp*(a0*a1/k1 - k1)/a3 ; Q2 = -1.5_dp*(a0*a1/k2 - k2)/a3
      ! Am: downward (phi-, Phi-) of the mode decaying from the top; Ap: upward (phi+, Phi+)
      Am(1,1) = twopi*(0.5_dp - R1 + 0.625_dp*P1) ; Am(1,2) = twopi*(0.5_dp - R2 + 0.625_dp*P2)
      Am(2,1) = twopi*(-0.125_dp + 0.625_dp*P1 - Q1) ; Am(2,2) = twopi*(-0.125_dp + 0.625_dp*P2 - Q2)
      Ap(1,1) = twopi*(0.5_dp + R1 + 0.625_dp*P1) ; Ap(1,2) = twopi*(0.5_dp + R2 + 0.625_dp*P2)
      Ap(2,1) = twopi*(-0.125_dp + 0.625_dp*P1 + Q1) ; Ap(2,2) = twopi*(-0.125_dp + 0.625_dp*P2 + Q2)

      if (linear_B) then

        !! Planck function linear in tau, B = B0 + B1 tau (Rooney et al. 2023, Eqs. 9-10).
        !! Particular solution I0 = B, I1 = B1/a1, I2 = I3 = 0 (exact), moments (Eqs. 26-27):
        !! Z1+/- = pi (1 - w)/a0 (B +/- 2 B1/a1),  Z2+/- = -pi (1 - w)/(4 a0) B
        if (dtau <= dtau_min) then
          B0 = 0.5_dp*(be(k) + be(k+1))
          B1 = 0.0_dp
        else
          B0 = be(k)
          B1 = (be(k+1) - be(k))/dtau
        end if
        fa = pi*(1.0_dp - w0)/a0
        Zmt(1) = fa*(B0 - 2.0_dp*B1/a1) ; Zmt(2) = -0.25_dp*fa*B0
        Zpt(1) = fa*(B0 + 2.0_dp*B1/a1) ; Zpt(2) = Zmt(2)
        Zmb(1) = fa*(B0 + B1*dtau - 2.0_dp*B1/a1) ; Zmb(2) = -0.25_dp*fa*(B0 + B1*dtau)
        Zpb(1) = fa*(B0 + B1*dtau + 2.0_dp*B1/a1) ; Zpb(2) = Zmb(2)

      else

        !! Planck function exponential in tau, B = B0 exp(-f0 tau), f0 = -ln(B1/B0)/tau1 (Eq. 3)
        if (dtau <= dtau_min) then
          B0 = 0.5_dp*(be(k) + be(k+1))
          f0 = 0.0_dp
        else
          B0 = be(k)
          f0 = -log(be(k+1)/be(k))/dtau
          !! Nudge f0 off the f0^2 = k^2 resonance (singular particular solution)
          if (abs(1.0_dp - (f0/k1)**2) < res_eps) f0 = f0*(1.0_dp + 2.0_dp*res_eps)
          if (abs(1.0_dp - (f0/k2)**2) < res_eps) f0 = f0*(1.0_dp + 2.0_dp*res_eps)
        end if
        rB = exp(-f0*dtau)

        !! Particular solution I_l = eta_l exp(-f0 tau) of Eq. 5, b0 = (1 - w) B0:
        !!   a0 eta0 + f0 eta1 = b0,  f0 eta0 + a1 eta1 + 2 f0 eta2 = 0,
        !!   2 f0 eta1 + a2 eta2 + 3 f0 eta3 = 0,  3 f0 eta2 + a3 eta3 = 0
        !! (tridiagonal - eliminated from the bottom row; f0 = 0 gives eta0 = B0)
        b0s = a0*B0
        c3 = a2 - 9.0_dp*f0**2/a3
        c2 = a1 - 4.0_dp*f0**2/c3
        c1 = a0 - f0**2/c2
        eta0 = b0s/c1
        eta1 = -(f0/c2)*eta0
        eta2 = -(2.0_dp*f0/c3)*eta1
        eta3 = -(3.0_dp*f0/a3)*eta2

        !! Particular solution moments at the layer top (Eq. 10): Z- (downward), Z+ (upward);
        !! at the bottom they are scaled by exp(-f0 tau1)
        Zmt(1) = twopi*(0.5_dp*eta0 - eta1 + 0.625_dp*eta2)
        Zmt(2) = twopi*(-0.125_dp*eta0 + 0.625_dp*eta2 - eta3)
        Zpt(1) = twopi*(0.5_dp*eta0 + eta1 + 0.625_dp*eta2)
        Zpt(2) = twopi*(-0.125_dp*eta0 + 0.625_dp*eta2 + eta3)
        Zmb(:) = Zmt(:)*rB
        Zpb(:) = Zpt(:)*rB

      end if

      !! Layer diffuse reflection and transmission, A2 A1^-1 (Eqs. 11-12, Appendix A). With the
      !! columns (modes decaying from the top | from the bottom) and E = diag(e1, e2):
      !! A1 = [[Am, Ap E], [Ap E, Am]], A2 = [[Am E, Ap], [Ap, Am E]], which are block diagonal in the
      !! sum / difference basis: A2 A1^-1 = 1/2 [[X + Y, X - Y], [X - Y, X + Y]],
      !! X = S2 S1^-1, Y = D2 D1^-1, S1/D1 = Am +/- Ap E, S2/D2 = Am E +/- Ap
      S1(:,1) = Am(:,1) + Ap(:,1)*e1 ; S1(:,2) = Am(:,2) + Ap(:,2)*e2
      D1(:,1) = Am(:,1) - Ap(:,1)*e1 ; D1(:,2) = Am(:,2) - Ap(:,2)*e2
      S2(:,1) = Am(:,1)*e1 + Ap(:,1) ; S2(:,2) = Am(:,2)*e2 + Ap(:,2)
      D2(:,1) = Am(:,1)*e1 - Ap(:,1) ; D2(:,2) = Am(:,2)*e2 - Ap(:,2)
      Xm(:,:) = matmul(S2(:,:), inv2(S1(:,:)))
      Ym(:,:) = matmul(D2(:,:), inv2(D1(:,:)))
      T(:,:,k) = 0.5_dp*(Xm(:,:) + Ym(:,:))
      R(:,:,k) = 0.5_dp*(Xm(:,:) - Ym(:,:))

      !! Thermal emission of the layer (Eq. 12, zero incident diffuse flux): outgoing = H2 + A2 A1^-1 H1,
      !! H1 = -(Z- at top, Z+ at bottom), H2 = (Z- at bottom, Z+ at top)
      zs(:) = Zmt(:) + Zpb(:)
      zd(:) = Zmt(:) - Zpb(:)
      sv(:) = matmul(Xm(:,:), zs(:))
      dv(:) = matmul(Ym(:,:), zd(:))
      Ed(:,k) = Zmb(:) - 0.5_dp*(sv(:) + dv(:))   ! downward (F-, F-script-) out of the bottom
      Eu(:,k) = Zpt(:) - 0.5_dp*(sv(:) - dv(:))      ! upward (F+, F-script+) out of the top

    end do

    !! Adding (Wu et al. 2016, Sect. 3), one-sided form: composites from the bottom (Rbot, Sup) and
    !! Md_k = (E - R_k Rbot_{k+1})^-1, then one downward sweep (identical to Eqs. 22-24)
    !! Lower boundary - gas giant: I_up = I_down + B(Tint) in every direction, so the upward moments
    !! are the downward moments plus (pi B_int, -pi B_int/4)
    Rbot(:,:,nlev) = Imat(:,:)
    Sup(1,nlev) = pi*be_int ; Sup(2,nlev) = -0.25_dp*pi*be_int
    do k = nlay, 1, -1
      Mi(:,:) = inv2(Imat(:,:) - matmul(Rbot(:,:,k+1), R(:,:,k)))
      Rbot(:,:,k) = R(:,:,k) + matmul(T(:,:,k), matmul(Mi(:,:), matmul(Rbot(:,:,k+1), T(:,:,k))))
      Sup(:,k) = Eu(:,k) + matmul(T(:,:,k), matmul(Mi(:,:), Sup(:,k+1) + matmul(Rbot(:,:,k+1), Ed(:,k))))
      Md(:,:,k) = inv2(Imat(:,:) - matmul(R(:,:,k), Rbot(:,:,k+1)))
    end do

    !! Downward sweep - no downward flux at the top of the atmosphere (Eq. 26b)
    Dn(:,1) = 0.0_dp
    Up(:,1) = Sup(:,1)
    do k = 1, nlay
      Dn(:,k+1) = matmul(Md(:,:,k), matmul(T(:,:,k), Dn(:,k)) + Ed(:,k) + matmul(R(:,:,k), Sup(:,k+1)))
      Up(:,k+1) = Sup(:,k+1) + matmul(Rbot(:,:,k+1), Dn(:,k+1))
    end do

    !! Fluxes are the first (P1) moments (Eq. 25)
    flx_up(:) = Up(1,:)
    flx_down(:) = Dn(1,:)

  end subroutine lw_SDA_method

  !! Inverse of a 2x2 matrix
  pure function inv2(A) result(B)
    implicit none
    real(dp), dimension(2,2), intent(in) :: A
    real(dp), dimension(2,2) :: B
    real(dp) :: rdet
    rdet = 1.0_dp/(A(1,1)*A(2,2) - A(1,2)*A(2,1))
    B(1,1) = A(2,2)*rdet ; B(1,2) = -A(1,2)*rdet
    B(2,1) = -A(2,1)*rdet ; B(2,2) = A(1,1)*rdet
  end function inv2

end module lw_SDA_mod
