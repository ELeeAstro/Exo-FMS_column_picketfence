!!!
! Elspeth KH Lee - Jun 2024 : Initial version
!                - Nov 2026 : Rewritten following PICASO exactly, enabled
!                - Nov 2026 : Optimisations - compile-time band widths in the solver
!                - Nov 2026 : 4x4 block-tridiagonal Thomas solve (replaces the banded matrix solve), direct
!                             beam exponentials computed once
!
! sw: Four-stream spherical harmonic method (SH4, Rooney et al. 2024), converted from PICASO
!     Pros: Four-stream accuracy, works well for clouds and hazes
!     Cons: Slower than the two-stream methods and sw_SDA (same answer)
!!!

module sw_SH4_mod
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
  integer, parameter :: nstr = 4

  !! Use Two-Term HG function (TTHG), parameters are the PICASO defaults
  !! forward fraction = frac_a + frac_b*g_back**frac_c, g_forward = c_forward*g, g_back = c_back*g
  logical, parameter :: TTHG = .False.
  real(dp), parameter :: tthg_forward = 1.0_dp, tthg_back = -0.5_dp
  real(dp), parameter :: tthg_frac_a = 1.0_dp, tthg_frac_b = -1.0_dp
  integer, parameter :: tthg_frac_c = 2

  !! Band widths of the four-stream matrix

  private :: sw_SH_four_stream, sh4_block, solve4
  public :: sw_SH4

contains

  subroutine sw_SH4(nlay, nlev, tau_e, mu_z, Finc, ssa, gg, sw_up, sw_down, sw_net, asr)
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
      call sw_SH_four_stream(nlay, nlev, Finc, mu_z, tau_e, ssa, gg, sw_down, sw_up)
    else
      sw_down(:) = 0.0_dp
      sw_up(:) = 0.0_dp
    end if

    !! Net sw flux
    sw_net(:) = sw_up(:) - sw_down(:)

    !! Absorbed Stellar Radiation (ASR)
    asr = sw_down(1) - sw_up(1)

  end subroutine sw_SH4

  subroutine sw_SH_four_stream(nlay, nlev, F0PI, u0, tau_in, w_in, g_in, flx_down, flx_up)
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
    real(dp) :: Del, b_top, b_surface, b_surface_SH4
    real(dp), dimension(nstr) :: Pu0
    real(dp), dimension(nlay) :: dtau, w0, f_deltaM, g_f, g_b, f_f
    real(dp), dimension(nstr,nlay) :: chi
    real(dp), dimension(nstr,nlay) :: w_multi, a, bb, Dels, eta
    real(dp), dimension(nlay) :: beta, gama, disc, lam1, lam2, exptrm1, exptrm2
    real(dp), dimension(nlay) :: R1, R2, Q1, Q2, S1, S2
    real(dp), dimension(nlay) :: p1pl, p2pl, q1pl, q2pl, p1mn, p2mn, q1mn, q2mn
    real(dp), dimension(nlay) :: z1pl, z1mn, z2pl, z2mn
    real(dp), dimension(4,4,nlay) :: Ftop
    real(dp), dimension(4,nlay) :: zdown, zup
    real(dp), dimension(nlay) :: X1, X2, X3, X4
    real(dp), dimension(nlay) :: z1mn_up, z2mn_up, z1pl_up, z2pl_up
    real(dp), dimension(nlay) :: z1mn_down, z2mn_down, z1pl_down, z2pl_down
    real(dp), dimension(nlev) :: tau, expon, cum_trans, trans

    !! If zero albedo across all atmospheric layers then return direct beam only
    if (all(w_in(:) <= 1.0e-12_dp)) then
      cum_trans(:) = (tau_in(:) - tau_in(1))/u0
      flx_down(:) = F0PI * u0 * exp(-cum_trans(:))
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
    Pu0(3) = (3.0_dp*u0**2 - 1.0_dp)/2.0_dp
    Pu0(4) = (5.0_dp*(-u0)**3 - 3.0_dp*(-u0))/2.0_dp

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
    b_surface_SH4 = 0.0_dp

    !! Eigenvalues (lam2 from lam1**2 lam2**2 = gama to avoid cancellation)
    beta(:) = a(1,:)*a(2,:) + 4.0_dp*a(1,:)*a(4,:)/9.0_dp + a(3,:)*a(4,:)/9.0_dp
    gama(:) = a(1,:)*a(2,:)*a(3,:)*a(4,:)/9.0_dp
    disc(:) = sqrt(max(beta(:)**2 - 4.0_dp*gama(:), 0.0_dp))
    lam1(:) = sqrt((beta(:) + disc(:))/2.0_dp)
    lam2(:) = sqrt(2.0_dp*gama(:)/(beta(:) + disc(:)))

    !! Direct beam particular solution
    do k = 1, nlay
      Del = 9.0_dp*((1.0_dp/u0)**4 - beta(k)*(1.0_dp/u0)**2 + gama(k))
      Dels(1,k) = (a(2,k)*bb(1,k) - bb(2,k)/u0)*(a(3,k)*a(4,k) - 9.0_dp/u0**2) &
        & + 2.0_dp*(a(4,k)*bb(3,k) - 2.0_dp*a(4,k)*bb(1,k) - 3.0_dp*bb(4,k)/u0)/u0**2
      Dels(2,k) = (a(1,k)*bb(2,k) - bb(1,k)/u0)*(a(3,k)*a(4,k) - 9.0_dp/u0**2) &
        & - 2.0_dp*a(1,k)*(a(4,k)*bb(3,k) - 3.0_dp*bb(4,k)/u0)/u0
      Dels(3,k) = (a(4,k)*bb(3,k) - 3.0_dp*bb(4,k)/u0)*(a(1,k)*a(2,k) - 1.0_dp/u0**2) &
        & - 2.0_dp*a(4,k)*(a(1,k)*bb(2,k) - bb(1,k)/u0)/u0
      Dels(4,k) = (a(3,k)*bb(4,k) - 3.0_dp*bb(3,k)/u0)*(a(1,k)*a(2,k) - 1.0_dp/u0**2) &
        & + 2.0_dp*(3.0_dp*a(1,k)*bb(2,k) - 2.0_dp*a(1,k)*bb(4,k) - 3.0_dp*bb(1,k)/u0)/u0**2
      eta(:,k) = Dels(:,k)/Del
    end do

    z1pl(:) = (eta(1,:)/2.0_dp + eta(2,:) + 5.0_dp*eta(3,:)/8.0_dp)*twopi
    z1mn(:) = (eta(1,:)/2.0_dp - eta(2,:) + 5.0_dp*eta(3,:)/8.0_dp)*twopi
    z2pl(:) = (-eta(1,:)/8.0_dp + 5.0_dp*eta(3,:)/8.0_dp + eta(4,:))*twopi
    z2mn(:) = (-eta(1,:)/8.0_dp + 5.0_dp*eta(3,:)/8.0_dp - eta(4,:))*twopi

    exptrm1(:) = exp(-min(lam1(:)*dtau(:), 35.0_dp))
    exptrm2(:) = exp(-min(lam2(:)*dtau(:), 35.0_dp))

    !! Legendre moments of the homogeneous modes
    R1(:) = -a(1,:)/lam1(:) ; R2(:) = -a(1,:)/lam2(:)
    Q1(:) = 0.5_dp*(a(1,:)*a(2,:)/(lam1(:)**2) - 1.0_dp) ; Q2(:) = 0.5_dp*(a(1,:)*a(2,:)/(lam2(:)**2) - 1.0_dp)
    S1(:) = -3.0_dp/(2.0_dp*a(4,:))*(a(1,:)*a(2,:)/lam1(:) - lam1(:))
    S2(:) = -3.0_dp/(2.0_dp*a(4,:))*(a(1,:)*a(2,:)/lam2(:) - lam2(:))

    !! Hemispheric flux (p) and second Marshak moment (q) of each mode
    p1pl(:) = (0.5_dp + R1(:) + 5.0_dp*Q1(:)/8.0_dp)*twopi
    p2pl(:) = (0.5_dp + R2(:) + 5.0_dp*Q2(:)/8.0_dp)*twopi
    q1pl(:) = (-1.0_dp/8.0_dp + 5.0_dp*Q1(:)/8.0_dp + S1(:))*twopi
    q2pl(:) = (-1.0_dp/8.0_dp + 5.0_dp*Q2(:)/8.0_dp + S2(:))*twopi
    p1mn(:) = (0.5_dp - R1(:) + 5.0_dp*Q1(:)/8.0_dp)*twopi
    p2mn(:) = (0.5_dp - R2(:) + 5.0_dp*Q2(:)/8.0_dp)*twopi
    q1mn(:) = (-1.0_dp/8.0_dp + 5.0_dp*Q1(:)/8.0_dp - S1(:))*twopi
    q2mn(:) = (-1.0_dp/8.0_dp + 5.0_dp*Q2(:)/8.0_dp - S2(:))*twopi

    !! Moments at the layer top: rows (down 1st, down 2nd, up 1st, up 2nd),
    !! columns modes (e^-lam1 t, e^+lam1 t, e^-lam2 t, e^+lam2 t)
    do k = 1, nlay
      Ftop(1,:,k) = (/p1mn(k), p1pl(k), p2mn(k), p2pl(k)/)
      Ftop(2,:,k) = (/q1mn(k), q1pl(k), q2mn(k), q2pl(k)/)
      Ftop(3,:,k) = (/p1pl(k), p1mn(k), p2pl(k), p2mn(k)/)
      Ftop(4,:,k) = (/q1pl(k), q1mn(k), q2pl(k), q2mn(k)/)
    end do

    !! Direct beam transmission, and the capped form for the particular solution:
    !! exp(-min(x, 35)) = max(exp(-x), exp(-35)) (one exponential per level)
    trans(:) = exp(-tau(:)/u0)
    expon(:) = max(trans(:), exp(-35.0_dp))
    z1mn_up(:) = z1mn(:)*expon(2:nlev)
    z2mn_up(:) = z2mn(:)*expon(2:nlev)
    z1pl_up(:) = z1pl(:)*expon(2:nlev)
    z2pl_up(:) = z2pl(:)*expon(2:nlev)
    z1mn_down(:) = z1mn(:)*expon(1:nlay)
    z2mn_down(:) = z2mn(:)*expon(1:nlay)
    z1pl_down(:) = z1pl(:)*expon(1:nlay)
    z2pl_down(:) = z2pl(:)*expon(1:nlay)

    !! Particular solution moments at the layer top (down) and bottom (up)
    zdown(1,:) = z1mn_down(:) ; zdown(2,:) = z2mn_down(:) ; zdown(3,:) = z1pl_down(:) ; zdown(4,:) = z2pl_down(:)
    zup(1,:) = z1mn_up(:) ; zup(2,:) = z2mn_up(:) ; zup(3,:) = z1pl_up(:) ; zup(4,:) = z2pl_up(:)

    !! Block-tridiagonal solve: X1, X3 = decaying mode amplitudes at the layer top,
    !! X2, X4 = growing mode amplitudes at the layer bottom (PICASO's X(4k-2) = X2 exptrm1, X(4k) = X4 exptrm2)
    !! Top: nothing enters (b_top). Bottom: nothing reflected (gas giant, b_surface)
    call sh4_block(nlay, Ftop, exptrm1, exptrm2, zdown, zup, b_top - z1mn_down(1), -b_top/4.0_dp - z2mn_down(1), &
      & b_surface - z1pl_up(nlay), b_surface_SH4 - z2pl_up(nlay), X1, X2, X3, X4)

    !! Level fluxes - top level from layer 1, level k+1 from the bottom of layer k
    flx_down(1) = p1mn(1)*X1(1) + p1pl(1)*exptrm1(1)*X2(1) + p2mn(1)*X3(1) + p2pl(1)*exptrm2(1)*X4(1) + z1mn_down(1)
    flx_up(1) = p1pl(1)*X1(1) + p1mn(1)*exptrm1(1)*X2(1) + p2pl(1)*X3(1) + p2mn(1)*exptrm2(1)*X4(1) + z1pl_down(1)
    do k = 1, nlay
      flx_down(k+1) = p1mn(k)*exptrm1(k)*X1(k) + p1pl(k)*X2(k) + p2mn(k)*exptrm2(k)*X3(k) + p2pl(k)*X4(k) + z1mn_up(k)
      flx_up(k+1) = p1pl(k)*exptrm1(k)*X1(k) + p1mn(k)*X2(k) + p2pl(k)*exptrm2(k)*X3(k) + p2mn(k)*X4(k) + z1pl_up(k)
    end do

    !! Add the direct beam to the downward flux
    flx_down(:) = flx_down(:) + F0PI*u0*trans(:)


  end subroutine sw_SH_four_stream

  !! Four-stream SH system as a 4x4 block-tridiagonal Thomas solve.
  !! Unknowns per layer y = (X1, X2, X3, X4): decaying modes at the layer top, growing modes at the
  !! layer bottom (Stamnes & Conklin 1984 scaling), so all coefficients are bounded.
  !! Block row k: rows 1:2 = downward moments at the interface above (or the top boundary),
  !! rows 3:4 = upward moments at the interface below (or the bottom boundary, no reflection).
  subroutine sh4_block(nlay, Ftop, e1, e2, zdown, zup, rt1, rt2, rb1, rb2, X1, X2, X3, X4)
    implicit none
    integer, intent(in) :: nlay
    real(dp), dimension(4,4,nlay), intent(in) :: Ftop
    real(dp), dimension(nlay), intent(in) :: e1, e2
    real(dp), dimension(4,nlay), intent(in) :: zdown, zup
    real(dp), intent(in) :: rt1, rt2, rb1, rb2
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
        Rh(1,1) = rt1
        Rh(2,1) = rt2
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
        D(3:4,:) = Gb(3:4,:)
        Rh(3,1) = rb1
        Rh(4,1) = rb2
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

end module sw_SH4_mod
