!!!
! Elspeth KH Lee - Nov 2026 : Initial version
!                - Nov 2026 : Optimisations - quadrature tables once per call, even/odd Legendre
!                             sums, sum/difference layer mode outputs, angle independent source
!                             function terms once per layer, thin-layer kernels from the
!                             transmissions, reciprocal angles
!
! lw: Four-stream discrete ordinates adding method (4DDA, Zhang et al. 2013, 2016)
!     Pros: Four-stream multiple scattering at near two-stream cost, only small matrix inversions
!     Cons: Four streams only - lw_AD is better for very strong forward scattering
!!!

module lw_DDA_mod
  use, intrinsic :: iso_fortran_env
  use WENO4_mod, only : interpolate_weno4
  implicit none

  !! Precision variables
  integer, parameter :: dp = REAL64

  !! Required constants
  real(dp), parameter :: pi = 4.0_dp * atan(1.0_dp)
  real(dp), parameter :: sb = 5.670374419e-8_dp

  !! Number of streams per hemisphere (4-stream method) and retained phase function moments
  integer, parameter :: nmu = 2
  integer, parameter :: nmom = 2*nmu

  !! Double-Gauss nodes and weights on [0,1] (sum(w) = 1, sum(w*mu) = 1/2)
  real(dp), dimension(nmu), parameter :: uarr = (/0.2113248654051871_dp, 0.7886751345948129_dp/)
  real(dp), dimension(nmu), parameter :: w = (/0.5_dp, 0.5_dp/)
  real(dp), dimension(nmu), parameter :: ruarr = 1.0_dp/uarr

  !! Single scattering albedo cap and minimum layer optical depth
  real(dp), parameter :: w_max = 1.0_dp - 1.0e-6_dp
  real(dp), parameter :: dtau_min = 1.0e-10_dp

  private :: lw_DDA_method, layer_RT, inv2
  public :: lw_DDA

contains

  subroutine lw_DDA(nlay, nlev, Tl, pl, pe, tau_e, ssa, gg, Tint, lw_up, lw_down, lw_net, olr)
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

    !! Quadrature constants 
    integer :: l, i, j
    real(dp), dimension(nmom, nmu) :: Plg
    real(dp), dimension(nmom) :: c2l
    real(dp), dimension(nmom,nmu,nmu) :: PP

    !! Legendre polynomials P_l(mu_i), l = 0..3 (stored at l+1)
    Plg(1,:) = 1.0_dp
    Plg(2,:) = uarr(:)
    Plg(3,:) = 0.5_dp*(3.0_dp*uarr(:)**2 - 1.0_dp)
    Plg(4,:) = 0.5_dp*(5.0_dp*uarr(:)**3 - 3.0_dp*uarr(:))

    !! (2l+1) factors of the phase function expansion
    do l = 0, nmom-1
      c2l(l+1) = real(2*l+1,dp)
    end do

    !! Quadrature-weighted Legendre products, PP(l,i,j) = w_j (2l+1) P_l(mu_i) P_l(mu_j)/2
    do j = 1, nmu
      do i = 1, nmu
        PP(:,i,j) = 0.5_dp*w(j)*c2l(:)*Plg(:,i)*Plg(:,j)
      end do
    end do

    !! Use WENO4 method to (smoothly) interpolate layers to levels, in ln p
    Te(:) = interpolate_weno4(log(pe), log(pl), Tl, .False.)

    !! Edges are linearly extrapolated in log(p)-log(T) space.
    Te(1) = 10.0_dp**(log10(Te(2)) + (log10(pe(1)/pe(2))/log10(pl(1)/pe(2))) * log10(Tl(1)/Te(2)))
    Te(nlev) = 10.0_dp**(log10(Te(nlay)) + (log10(pe(nlev)/pe(nlay))/log10(pl(nlay)/pe(nlay))) * log10(Tl(nlay)/Te(nlay)))

    !! Find integrated planck function for each level
    be(:) = (sb * Te(:)**4)/pi
    be_int = (sb * Tint**4)/pi

    !! Longwave flux calculation
    call lw_DDA_method(nlay, nlev, be(:), be_int, tau_e(:), ssa(:), gg(:), PP, lw_up(:), lw_down(:))

    !! Net lw fluxes at each level
    lw_net(:) = lw_up(:) - lw_down(:)

    !! Outgoing Longwave Radiation (olr)
    olr = lw_up(1)

  end subroutine lw_DDA

  subroutine lw_DDA_method(nlay, nlev, be, be_int, tau_in, w_in, g_in, PP, flx_up, flx_down)
    implicit none

    !! Input variables
    integer, intent(in) :: nlay, nlev
    real(dp), dimension(nlev), intent(in) :: be, tau_in
    real(dp), dimension(nlay), intent(in) :: w_in, g_in
    real(dp), intent(in) :: be_int
    real(dp), dimension(nmom,nmu,nmu), intent(in) :: PP

    !! Output variables
    real(dp), dimension(nlev), intent(out) :: flx_up, flx_down

    !! Work variables and arrays
    integer :: k, l, i, j
    real(dp), dimension(nmom) :: chi
    real(dp), dimension(nlay) :: dtau, w0
    real(dp), dimension(nmu,nmu) :: Ws, Wd, Imat, Mi
    real(dp) :: pe_, po_
    real(dp), dimension(nmu,nmu,nlay) :: R, T, Md
    real(dp), dimension(nmu,nlay) :: Eu, Ed
    real(dp), dimension(nmu,nmu,nlev) :: Rbot
    real(dp), dimension(nmu,nlev) :: Sup
    real(dp), dimension(nmu) :: Iup, Idn

    real(dp) :: f

    Imat(:,:) = 0.0_dp
    do i = 1, nmu
      Imat(i,i) = 1.0_dp
    end do

    !! Layer optical depths and single scattering albedo
    dtau(:) = tau_in(2:nlev) - tau_in(1:nlay)
    w0(:) = min(max(w_in(:), 0.0_dp), w_max)

    !! Reflection, transmission and thermal sources of each layer
    do k = 1, nlay

      !! Henyey-Greenstein moments, delta-M scaled with the first truncated moment chi_4 = g^4
      do l = 0, nmom-1
        chi(l+1) = g_in(k)**l
      end do
      f = g_in(k)**nmom
      if (w0(k) > 0.0_dp .and. f > 1.0e-12_dp .and. f < 1.0_dp - 1.0e-12_dp) then
        chi(:) = (chi(:) - f)/(1.0_dp - f)
        dtau(k) = (1.0_dp - w0(k)*f)*dtau(k)
        w0(k) = w0(k)*(1.0_dp - f)/(1.0_dp - w0(k)*f)
      end if

      !! Quadrature-weighted phase matrices, normalised so (1/2) int P dmu' = 1:
      !! Ws(i,j) = w_j P(mu_i, mu_j)/2 (same hemisphere), Wd(i,j) = w_j P(mu_i,-mu_j)/2 (opposite)
      !! P_l(-mu) = (-1)^l P_l(mu): Ws = even + odd l terms, Wd = even - odd l terms
      do j = 1, nmu
        do i = 1, nmu
          pe_ = sum(chi(1:nmom:2)*PP(1:nmom:2,i,j))
          po_ = sum(chi(2:nmom:2)*PP(2:nmom:2,i,j))
          Ws(i,j) = pe_ + po_
          Wd(i,j) = pe_ - po_
        end do
      end do

      call layer_RT(dtau(k), w0(k), be(k), be(k+1), Ws, Wd, Imat, R(:,:,k), T(:,:,k), Eu(:,k), Ed(:,k))

    end do

    !! Upward pass - composites from the bottom: Rbot = reflection of the slab below a level (for
    !! downward light), Sup = upward intensity emitted by the slab below, and Md_k = (I - R_k Rbot_{k+1})^-1
    !! for the downward pass. Lower boundary - gas giant: I_up = I_down + B(Tint) in every stream
    Rbot(:,:,nlev) = Imat(:,:)
    Sup(:,nlev) = be_int
    do k = nlay, 1, -1
      Mi(:,:) = inv2(Imat(:,:) - matmul(Rbot(:,:,k+1), R(:,:,k)))
      Rbot(:,:,k) = R(:,:,k) + matmul(T(:,:,k), matmul(Mi(:,:), matmul(Rbot(:,:,k+1), T(:,:,k))))
      Sup(:,k) = Eu(:,k) + matmul(T(:,:,k), matmul(Mi(:,:), Sup(:,k+1) + matmul(Rbot(:,:,k+1), Ed(:,k))))
      Md(:,:,k) = inv2(Imat(:,:) - matmul(R(:,:,k), Rbot(:,:,k+1)))
    end do

    !! Downward pass - no incident intensity at the top
    Idn(:) = 0.0_dp
    do k = 1, nlev
      if (k > 1) then
        Idn(:) = matmul(Md(:,:,k-1), matmul(T(:,:,k-1), Idn(:)) + Ed(:,k-1) + matmul(R(:,:,k-1), Sup(:,k)))
      end if
      Iup(:) = Sup(:,k) + matmul(Rbot(:,:,k), Idn(:))
      flx_down(k) = 2.0_dp*pi*sum(w(:)*uarr(:)*Idn(:))
      flx_up(k) = 2.0_dp*pi*sum(w(:)*uarr(:)*Iup(:))
    end do

  end subroutine lw_DDA_method

  !! Analytic 4-stream discrete ordinates solution of one homogeneous layer:
  !! R, T (2x2, incident from either side) and thermal sources Eu (up, top), Ed (down, bottom)
  subroutine layer_RT(dtau, w0, B_top, B_bot, Ws, Wd, Imat, R, T, Eu, Ed)
    implicit none

    real(dp), intent(in) :: dtau, w0, B_top, B_bot
    real(dp), dimension(nmu,nmu), intent(in) :: Ws, Wd, Imat
    real(dp), dimension(nmu,nmu), intent(out) :: R, T
    real(dp), dimension(nmu), intent(out) :: Eu, Ed

    integer :: i, m
    real(dp), dimension(nmu,nmu) :: al, bt, amb, apb, Q, Gp, Gm, GmE, GpE, P1, P2, M1o, M2o
    real(dp), dimension(nmu) :: kk, ex, X, Y, v1, v2, d, one, sum_o, dif_o
    real(dp) :: tr, det, disc, lam1, lam, B1

    !! Optically negligible layer - transparent
    if (dtau <= dtau_min) then
      R(:,:) = 0.0_dp
      T(:,:) = Imat(:,:)
      Eu(:) = 0.0_dp
      Ed(:) = 0.0_dp
      return
    end if

    !! mu dI+/dtau = alpha I+ - beta I-,  dI-/dtau = beta I+ - alpha I-  (per unit tau, divided by mu)
    do i = 1, nmu
      al(i,:) = (Imat(i,:) - w0*Ws(i,:))*ruarr(i)
      bt(i,:) = w0*Wd(i,:)*ruarr(i)
    end do
    amb(:,:) = al(:,:) - bt(:,:)
    apb(:,:) = al(:,:) + bt(:,:)
    Q(:,:) = matmul(apb(:,:), amb(:,:))

    !! Eigenvalues k^2 of Q (quadratic in k^2); smaller root via det/lam1 to avoid cancellation
    tr = Q(1,1) + Q(2,2)
    det = Q(1,1)*Q(2,2) - Q(1,2)*Q(2,1)
    disc = sqrt(max(tr**2 - 4.0_dp*det, 0.0_dp))
    lam1 = 0.5_dp*(tr + disc)
    kk(1) = sqrt(lam1)
    kk(2) = sqrt(max(det/lam1, tiny(1.0_dp)))

    !! Eigenvectors X of Q, then the decaying (e^{-k tau}) modes:
    !! I+ = Gp e^{-k tau}, I- = Gm e^{-k tau}, with Gp = (X - Y)/2, Gm = (X + Y)/2, Y = (alpha-beta) X/k
    !! The growing mode is the mirror image (I+ = Gm, I- = Gp), written as e^{-k (dtau - tau)}
    do m = 1, nmu
      lam = kk(m)**2
      v1(:) = (/Q(1,2), lam - Q(1,1)/)
      v2(:) = (/lam - Q(2,2), Q(2,1)/)
      if (sum(v1**2) >= sum(v2**2)) then
        X(:) = v1(:)*(1.0_dp/sqrt(sum(v1**2)))
      else
        X(:) = v2(:)*(1.0_dp/sqrt(sum(v2**2)))
      end if
      Y(:) = matmul(amb(:,:), X(:))*(1.0_dp/kk(m))
      Gp(:,m) = 0.5_dp*(X(:) - Y(:))
      Gm(:,m) = 0.5_dp*(X(:) + Y(:))
      ex(m) = exp(-kk(m)*dtau)
    end do

    !! Thermal particular solution, B(tau) = B_top + B1 tau: I+/- = B(tau) +/- d B1
    B1 = (B_bot - B_top)/dtau
    d(:) = matmul(inv2(Imat(:,:) - w0*(Ws(:,:) - Wd(:,:))), uarr(:))

    !! Boundary value problem for the coefficients (A_m, C_m):
    !! I+ = sum_m [A_m Gp_m e^{-k_m tau} + C_m Gm_m e^{-k_m (dtau - tau)}] + Ip+
    !! I- = sum_m [A_m Gm_m e^{-k_m tau} + C_m Gp_m e^{-k_m (dtau - tau)}] + Ip-
    !! with I-(0) = a and I+(dtau) = b given. With u = A + C, s = A - C:
    !!   a + b - (Ip-(0) + Ip+(D)) = (Gm + Gp E) u,   a - b - (Ip-(0) - Ip+(D)) = (Gm - Gp E) s
    !!   I+(0) + I-(D) = (Gp + Gm E) u + Ip+(0) + Ip-(D),   I+(0) - I-(D) = (Gp - Gm E) s + Ip+(0) - Ip-(D)
    do m = 1, nmu
      GmE(:,m) = Gm(:,m)*ex(m)
      GpE(:,m) = Gp(:,m)*ex(m)
    end do
    M1o(:,:) = inv2(Gm(:,:) + GpE(:,:))
    M2o(:,:) = inv2(Gm(:,:) - GpE(:,:))
    P1(:,:) = matmul(Gp(:,:) + GmE(:,:), M1o(:,:))   ! R + T
    P2(:,:) = matmul(Gp(:,:) - GmE(:,:), M2o(:,:))   ! R - T

    R(:,:) = 0.5_dp*(P1(:,:) + P2(:,:))
    T(:,:) = 0.5_dp*(P1(:,:) - P2(:,:))

    !! Thermal sources (a = b = 0): Ip-(0) + Ip+(D) = Ip+(0) + Ip-(D) = (B_top + B_bot),
    !! Ip-(0) - Ip+(D) = (B_top - B_bot) - 2 d B1,  Ip+(0) - Ip-(D) = (B_top - B_bot) + 2 d B1
    one(:) = 1.0_dp
    sum_o(:) = (B_top + B_bot)*(one(:) - matmul(P1(:,:), one(:)))
    dif_o(:) = (B_top - B_bot)*(one(:) - matmul(P2(:,:), one(:))) + 2.0_dp*B1*(d(:) + matmul(P2(:,:), d(:)))
    Eu(:) = 0.5_dp*(sum_o(:) + dif_o(:))
    Ed(:) = 0.5_dp*(sum_o(:) - dif_o(:))

  end subroutine layer_RT

  !! Inverse of a 2x2 matrix
  pure function inv2(A) result(B)
    implicit none

    real(dp), dimension(2,2), intent(in) :: A
    real(dp), dimension(2,2) :: B
    real(dp) :: rdet

    rdet = 1.0_dp/(A(1,1)*A(2,2) - A(1,2)*A(2,1))
    B(1,1) = A(2,2)*rdet
    B(1,2) = -A(1,2)*rdet
    B(2,1) = -A(2,1)*rdet
    B(2,2) = A(1,1)*rdet

  end function inv2

end module lw_DDA_mod