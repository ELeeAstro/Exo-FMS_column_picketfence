!!!
! Elspeth KH Lee - Oct 2026 : Initial version
!                - Nov 2026 : Reviewed against PICASO, provenance and differences documented;
!                             Hogan (2024) quadrature for the flux integration
!                - Nov 2026 : Optimisations - beam integrals from the stored transmissions (one
!                             exponential per layer and angle), compile-time band widths in the
!                             solver, reciprocals
!                - Nov 2026 : 4x4 block-tridiagonal Thomas solve (replaces the banded matrix solve)
!
! lw: Four-stream spherical harmonic method (SH4, Rooney et al. 2024), converted from PICASO,
!     with source function technique fluxes
!     Pros: Four-stream multiple scattering accuracy
!     Cons: Slower than the two-stream methods
!!!

module lw_SH4_mod
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
  integer, parameter :: nstr = 4

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

  private :: lw_SH_four_stream, e_alp, e_bet, sh4_block, solve4
  public :: lw_SH4

contains

  subroutine lw_SH4(nlay, nlev, Tl, pl, pe, tau_e, ssa, gg, Tint, lw_up, lw_down, lw_net, olr)
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
    call lw_SH_four_stream(nlay, nlev, be(:), be_int, tau_e(:), ssa(:), gg(:), &
      & lw_up(:), lw_down(:))

    !! Net lw fluxes at each level
    lw_net(:) = lw_up(:) - lw_down(:)

    !! Outgoing Longwave Radiation (olr)
    olr = lw_up(1)

  end subroutine lw_SH4

  subroutine lw_SH_four_stream(nlay, nlev, be, be_int, tau_in, w_in, g_in, flx_up, flx_down)
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
    real(dp) :: tautop, Btop, u, ru, P2u, P3u
    real(dp), dimension(nlay) :: dtau, w0, ff, B0, B1
    real(dp), dimension(nlay) :: w_multi0, w_multi1, w_multi2, w_multi3, a0, a1, a2, a3
    real(dp), dimension(nlay) :: beta, gama, disc, lam1, lam2, e1, e2
    real(dp), dimension(nlay) :: R1, R2, Q1, Q2, S1, S2
    real(dp), dimension(nlay) :: p1pl, p2pl, q1pl, q2pl, p1mn, p2mn, q1mn, q2mn
    real(dp), dimension(4,4,nlay) :: Ftop
    real(dp), dimension(4,nlay) :: zdown, zup
    real(dp), dimension(nlay) :: X1, X2, X3, X4
    real(dp), dimension(nlay) :: ex, ev1, od1, ev2, od2, T_up, T_dn, Ip1, Nup, Ndn, Aup, Adn
    real(dp), dimension(nlay) :: ea1, eb1, ea2, eb2
    real(dp), dimension(nlay) :: src_up, src_dn, el1, el2, wa0
    real(dp), dimension(nlev) :: lw_up_g, lw_down_g

    !! Delta-M scaling (PICASO convention, f = g**nstr)
    ff(:) = g_in(:)**nstr
    w0(:) = w_in(:)*(1.0_dp - ff(:))/(1.0_dp - w_in(:)*ff(:))
    w0(:) = min(w0(:), 1.0_dp - 1.0e-6_dp)
    dtau(:) = (tau_in(2:nlev) - tau_in(1:nlay))*(1.0_dp - w_in(:)*ff(:))

    !! Phase function Legendre coefficients and the PICASO 'a' coefficients
    w_multi0(:) = 1.0_dp
    w_multi1(:) = 3.0_dp*(g_in(:) - ff(:))/(1.0_dp - ff(:))
    w_multi2(:) = 5.0_dp*(g_in(:)**2 - ff(:))/(1.0_dp - ff(:))
    w_multi3(:) = 7.0_dp*(g_in(:)**3 - ff(:))/(1.0_dp - ff(:))
    a0(:) = 1.0_dp - w0(:)*w_multi0(:)
    a1(:) = 3.0_dp - w0(:)*w_multi1(:)
    a2(:) = 5.0_dp - w0(:)*w_multi2(:)
    a3(:) = 7.0_dp - w0(:)*w_multi3(:)

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

    !! Eigenvalues; lam2 in the cancellation free form lam1**2 * lam2**2 = gama
    beta(:) = a0(:)*a1(:) + 4.0_dp*a0(:)*a3(:)/9.0_dp + a2(:)*a3(:)/9.0_dp
    gama(:) = a0(:)*a1(:)*a2(:)*a3(:)/9.0_dp
    disc(:) = sqrt(max(beta(:)**2 - 4.0_dp*gama(:), 0.0_dp))
    lam1(:) = sqrt((beta(:) + disc(:))/2.0_dp)
    lam2(:) = sqrt(2.0_dp*gama(:)/(beta(:) + disc(:)))

    e1(:) = exp(-min(lam1(:)*dtau(:), 35.0_dp))
    e2(:) = exp(-min(lam2(:)*dtau(:), 35.0_dp))
    ! uncapped exp(-lam dtau) for the beam integrals (= e1, e2 unless lam dtau > 35)
    el1(:) = e1(:)
    where (lam1(:)*dtau(:) > 35.0_dp) el1(:) = exp(-lam1(:)*dtau(:))
    el2(:) = e2(:)
    where (lam2(:)*dtau(:) > 35.0_dp) el2(:) = exp(-lam2(:)*dtau(:))

    !! Legendre moments of the homogeneous modes
    R1(:) = -a0(:)/lam1(:) ; R2(:) = -a0(:)/lam2(:)
    Q1(:) = 0.5_dp*(a0(:)*a1(:)/lam1(:)**2 - 1.0_dp) ; Q2(:) = 0.5_dp*(a0(:)*a1(:)/lam2(:)**2 - 1.0_dp)
    S1(:) = -3.0_dp/(2.0_dp*a3(:))*(a0(:)*a1(:)/lam1(:) - lam1(:))
    S2(:) = -3.0_dp/(2.0_dp*a3(:))*(a0(:)*a1(:)/lam2(:) - lam2(:))

    !! Hemispheric flux (p) and second Marshak moment (q) of each mode
    p1pl(:) = twopi*(0.5_dp + R1(:) + 5.0_dp*Q1(:)/8.0_dp)
    p2pl(:) = twopi*(0.5_dp + R2(:) + 5.0_dp*Q2(:)/8.0_dp)
    q1pl(:) = twopi*(-1.0_dp/8.0_dp + 5.0_dp*Q1(:)/8.0_dp + S1(:))
    q2pl(:) = twopi*(-1.0_dp/8.0_dp + 5.0_dp*Q2(:)/8.0_dp + S2(:))
    p1mn(:) = twopi*(0.5_dp - R1(:) + 5.0_dp*Q1(:)/8.0_dp)
    p2mn(:) = twopi*(0.5_dp - R2(:) + 5.0_dp*Q2(:)/8.0_dp)
    q1mn(:) = twopi*(-1.0_dp/8.0_dp + 5.0_dp*Q1(:)/8.0_dp - S1(:))
    q2mn(:) = twopi*(-1.0_dp/8.0_dp + 5.0_dp*Q2(:)/8.0_dp - S2(:))

    !! Moments at the layer top: rows (down 1st, down 2nd, up 1st, up 2nd),
    !! columns modes (e^-lam1 t, e^+lam1 t, e^-lam2 t, e^+lam2 t)
    do k = 1, nlay
      Ftop(1,:,k) = (/p1mn(k), p1pl(k), p2mn(k), p2pl(k)/)
      Ftop(2,:,k) = (/q1mn(k), q1pl(k), q2mn(k), q2pl(k)/)
      Ftop(3,:,k) = (/p1pl(k), p1mn(k), p2pl(k), p2mn(k)/)
      Ftop(4,:,k) = (/q1pl(k), q1mn(k), q2pl(k), q2mn(k)/)
    end do

    !! Particular solution moments at layer top (down) and bottom (up), linear in tau Planck
    zdown(1,:) = twopi*(1.0_dp - w0(:))/a0(:) * (B0(:)/2.0_dp - B1(:)/a1(:))
    zdown(2,:) = -twopi*0.5_dp*(1.0_dp - w0(:))/(4.0_dp*a0(:)) * B0(:)
    zdown(3,:) = twopi*(1.0_dp - w0(:))/a0(:) * (B0(:)/2.0_dp + B1(:)/a1(:))
    zdown(4,:) = zdown(2,:)
    zup(1,:) = twopi*(1.0_dp - w0(:))/a0(:) * (B0(:)/2.0_dp - B1(:)/a1(:) + B1(:)*dtau(:)/2.0_dp)
    zup(2,:) = -twopi*0.5_dp*(1.0_dp - w0(:))/(4.0_dp*a0(:)) * (B0(:) + B1(:)*dtau(:))
    zup(3,:) = twopi*(1.0_dp - w0(:))/a0(:) * (B0(:)/2.0_dp + B1(:)/a1(:) + B1(:)*dtau(:)/2.0_dp)
    zup(4,:) = zup(2,:)

    !! Top boundary, isothermal layer above the model top (as lw_Toon89)
    tautop = dtau(1)*exp(-1.0_dp)
    Btop = (1.0_dp - exp(-tautop/ubari))*be(1)

    !! Block-tridiagonal solve: X1, X3 = decaying mode amplitudes at the layer top,
    !! X2, X4 = growing mode amplitudes at the layer bottom
    call sh4_block(nlay, Ftop, e1, e2, zdown, zup, Btop, be_int, X1, X2, X3, X4)

    !! Particular solution first moment (constant in tau)
    Ip1(:) = (1.0_dp - w0(:))/(a0(:)*a1(:))*B1(:)
    wa0(:) = w_multi0(:)*(1.0_dp - w0(:))/a0(:)

    ! Zero the total flux arrays
    flx_up(:) = 0.0_dp
    flx_down(:) = 0.0_dp

    do m = 1, nmu

      u = uarr(m)
      ru = ruarr(m)
      P2u = (3.0_dp*u**2 - 1.0_dp)/2.0_dp
      P3u = (5.0_dp*u**3 - 3.0_dp*u)/2.0_dp
      ex(:) = exp(-dtau(:)*ru)

      !! Scattering phase coefficients split into even and odd Legendre parts,
      !! at -u (and for the growing mode of each pair) the odd part flips sign
      ev1(:) = w_multi0(:) + w_multi2(:)*P2u*Q1(:)
      od1(:) = w_multi1(:)*u*R1(:) + w_multi3(:)*P3u*S1(:)
      ev2(:) = w_multi0(:) + w_multi2(:)*P2u*Q2(:)
      od2(:) = w_multi1(:)*u*R2(:) + w_multi3(:)*P3u*S2(:)

      ea1(:) = e_alp(nlay, ru + lam1(:), ex(:), el1(:))
      eb1(:) = e_bet(nlay, lam1(:), ru, dtau(:), ex(:), el1(:))
      ea2(:) = e_alp(nlay, ru + lam2(:), ex(:), el2(:))
      eb2(:) = e_bet(nlay, lam2(:), ru, dtau(:), ex(:), el2(:))

      !! Thermal particular solution integrated along the beam
      T_up(:) = u*(B0(:)*(1.0_dp - ex(:)) + B1(:)*(u - (dtau(:) + u)*ex(:)))
      T_dn(:) = u*(B0(:)*(1.0_dp - ex(:)) + B1(:)*(dtau(:) - u + u*ex(:)))
      Nup(:) = wa0(:)*T_up(:) + w_multi1(:)*u*Ip1(:)*u*(1.0_dp - ex(:))
      Ndn(:) = wa0(:)*T_dn(:) - w_multi1(:)*u*Ip1(:)*u*(1.0_dp - ex(:))

      !! Homogeneous solution integrated along the beam
      Aup(:) = X1(:)*(ev1(:) + od1(:))*ea1(:) + X2(:)*(ev1(:) - od1(:))*eb1(:) &
        & + X3(:)*(ev2(:) + od2(:))*ea2(:) + X4(:)*(ev2(:) - od2(:))*eb2(:)
      Adn(:) = X1(:)*(ev1(:) - od1(:))*eb1(:) + X2(:)*(ev1(:) + od1(:))*ea1(:) &
        & + X3(:)*(ev2(:) - od2(:))*eb2(:) + X4(:)*(ev2(:) + od2(:))*ea2(:)

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

  end subroutine lw_SH_four_stream

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

  !! Four-stream SH system as a 4x4 block-tridiagonal Thomas solve.
  !! Unknowns per layer y = (X1, X2, X3, X4): decaying modes at the layer top, growing modes at the
  !! layer bottom (Stamnes & Conklin 1984 scaling), so all coefficients are bounded.
  !! Block row k: rows 1:2 = downward moments at the interface above (or the top boundary),
  !! rows 3:4 = upward moments at the interface below (or the bottom boundary).
  subroutine sh4_block(nlay, Ftop, e1, e2, zdown, zup, Btop, be_int, X1, X2, X3, X4)
    implicit none
    integer, intent(in) :: nlay
    real(dp), dimension(4,4,nlay), intent(in) :: Ftop
    real(dp), dimension(nlay), intent(in) :: e1, e2
    real(dp), dimension(4,nlay), intent(in) :: zdown, zup
    real(dp), intent(in) :: Btop, be_int
    real(dp), dimension(nlay), intent(out) :: X1, X2, X3, X4

    integer :: k
    real(dp), dimension(4,4) :: Gt, Gb, D
    real(dp), dimension(4,3) :: Rh
    real(dp), dimension(2,4) :: Lr
    real(dp), dimension(2,4,nlay) :: Ur
    real(dp), dimension(2,2) :: M
    real(dp), dimension(4,nlay) :: wv, y
    real(dp), dimension(4,2,nlay) :: Cm
    real(dp), dimension(2) :: v

    do k = 1, nlay
      ! moments of layer k at its top (Gt) and bottom (Gb) per unit rescaled amplitude
      Gt(:,1) = Ftop(:,1,k) ; Gt(:,2) = Ftop(:,2,k)*e1(k) ; Gt(:,3) = Ftop(:,3,k) ; Gt(:,4) = Ftop(:,4,k)*e2(k)
      ! rows 1:2
      if (k == 1) then
        D(1:2,:) = Gt(1:2,:)
        Rh(1,1) = pi*Btop - zdown(1,1)
        Rh(2,1) = -pi*Btop/4.0_dp - zdown(2,1)
      else
        D(1:2,:) = -Gt(1:2,:)
        Rh(1:2,1) = zdown(1:2,k) - zup(1:2,k-1)
        ! eliminate y(k-1): Lr = Gb(k-1)(1:2,:)
        Lr(:,1) = Ftop(1:2,1,k-1)*e1(k-1) ; Lr(:,2) = Ftop(1:2,2,k-1)
        Lr(:,3) = Ftop(1:2,3,k-1)*e2(k-1) ; Lr(:,4) = Ftop(1:2,4,k-1)
        M(:,:) = matmul(Lr(:,:), Cm(:,:,k-1))
        D(1:2,:) = D(1:2,:) - matmul(M(:,:), Ur(:,:,k-1))
        Rh(1:2,1) = Rh(1:2,1) - matmul(Lr(:,:), wv(:,k-1))
      end if
      ! rows 3:4
      Gb(:,1) = Ftop(:,1,k)*e1(k) ; Gb(:,2) = Ftop(:,2,k) ; Gb(:,3) = Ftop(:,3,k)*e2(k) ; Gb(:,4) = Ftop(:,4,k)
      if (k == nlay) then
        D(3:4,:) = Gb(3:4,:) - Gb(1:2,:)
        Rh(3,1) = pi*be_int - zup(3,k) + zup(1,k)
        Rh(4,1) = -pi*be_int/4.0_dp - zup(4,k) + zup(2,k)
      else
        D(3:4,:) = Gb(3:4,:)
        Rh(3:4,1) = zdown(3:4,k+1) - zup(3:4,k)
        ! Ur = -Gt(k+1)(3:4,:)
        Ur(:,1,k) = -Ftop(3:4,1,k+1) ; Ur(:,2,k) = -Ftop(3:4,2,k+1)*e1(k+1)
        Ur(:,3,k) = -Ftop(3:4,3,k+1) ; Ur(:,4,k) = -Ftop(3:4,4,k+1)*e2(k+1)
      end if
      ! D [w C] = [r e3 e4]
      Rh(:,2:3) = 0.0_dp ; Rh(3,2) = 1.0_dp ; Rh(4,3) = 1.0_dp
      call solve4(D, Rh)
      wv(:,k) = Rh(:,1)
      Cm(:,:,k) = Rh(:,2:3)
    end do

    y(:,nlay) = wv(:,nlay)
    do k = nlay-1, 1, -1
      v(:) = matmul(Ur(:,:,k), y(:,k+1))
      y(:,k) = wv(:,k) - matmul(Cm(:,:,k), v(:))
    end do
    X1(:) = y(1,:) ; X2(:) = y(2,:) ; X3(:) = y(3,:) ; X4(:) = y(4,:)

  end subroutine sh4_block

  !! Solve A X = B in place (4x4, 3 right-hand sides), Gaussian elimination with partial pivoting
  subroutine solve4(A, B)
    implicit none
    real(dp), dimension(4,4), intent(inout) :: A
    real(dp), dimension(4,3), intent(inout) :: B
    integer :: i, j, p
    real(dp) :: fac
    real(dp), dimension(4) :: tr
    real(dp), dimension(3) :: tb
    do j = 1, 3
      p = j - 1 + maxloc(abs(A(j:4,j)), 1)
      if (p /= j) then
        tr(:) = A(j,:) ; A(j,:) = A(p,:) ; A(p,:) = tr(:)
        tb(:) = B(j,:) ; B(j,:) = B(p,:) ; B(p,:) = tb(:)
      end if
      do i = j+1, 4
        fac = A(i,j)/A(j,j)
        A(i,j+1:4) = A(i,j+1:4) - fac*A(j,j+1:4)
        B(i,:) = B(i,:) - fac*B(j,:)
      end do
    end do
    do i = 4, 1, -1
      do j = i+1, 4
        B(i,:) = B(i,:) - A(i,j)*B(j,:)
      end do
      B(i,:) = B(i,:)/A(i,i)
    end do
  end subroutine solve4

end module lw_SH4_mod
