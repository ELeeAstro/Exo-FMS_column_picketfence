!!!
! Elspeth KH Lee - Oct 2026 : Riccati (invariant-imbedding) solution of the linear sc scattering scheme
!
! lw: Short characteristics discrete ordinates with anisotropic multiple scattering, solved by
!     Riccati invariant imbedding
!     Pros: Accurate heating rates, any number of streams
!     Cons: Not well tested yet - use with caution, inaccurate for optically thick near-conservative
!           scattering unless the vertical resolution is increased
!!!

module lw_sc_scat_ric_mod
  use, intrinsic :: iso_fortran_env
  use WENO4_mod, only : interpolate_weno4
  implicit none

  !! Precision variables
  integer, parameter :: dp = REAL64

  !! Required constants
  real(dp), parameter :: pi = 4.0_dp * atan(1.0_dp)
  real(dp), parameter :: twopi = 2.0_dp * pi
  real(dp), parameter :: sb = 5.670374419e-8_dp

  !! Double-Gauss quadrature: Gauss-Legendre on each hemisphere [0,1], sum(w) = 1, sum(w*uarr) = 1/2
  !! here you can comment in/out groups of mu values for testing
  !! make sure to make clean and recompile if you change these

  !! 1 node per hemisphere - 2 stream
  ! integer, parameter :: nmu = 1
  ! real(dp), dimension(nmu), parameter :: uarr = (/0.5_dp/)
  ! real(dp), dimension(nmu), parameter :: w = (/1.0_dp/)

  !! 2 nodes per hemisphere - 4 stream
  ! integer, parameter :: nmu = 2
  ! real(dp), dimension(nmu), parameter :: uarr = (/0.2113248654051871_dp, 0.7886751345948129_dp/)
  ! real(dp), dimension(nmu), parameter :: w = (/0.5_dp, 0.5_dp/)

  !! 3 nodes per hemisphere - 6 stream
  ! integer, parameter :: nmu = 3
  ! real(dp), dimension(nmu), parameter :: uarr = &
  !   & (/0.1127016653792583_dp, 0.5_dp, 0.8872983346207417_dp/)
  ! real(dp), dimension(nmu), parameter :: w = &
  !   & (/0.2777777777777778_dp, 0.4444444444444444_dp, 0.2777777777777778_dp/)

  !! 4 nodes per hemisphere - 8 stream
  integer, parameter :: nmu = 4
  real(dp), dimension(nmu), parameter :: uarr = &
    & (/0.0694318442029737_dp, 0.3300094782075719_dp, 0.6699905217924281_dp, 0.9305681557970263_dp/)
  real(dp), dimension(nmu), parameter :: w = &
    & (/0.1739274225687269_dp, 0.3260725774312731_dp, 0.3260725774312731_dp, 0.1739274225687269_dp/)

  !! Total number of streams
  integer, parameter :: nstr = 2*nmu

  private :: solve_nopiv, lu_factor_nopiv, lu_solve_nopiv, lu_solve_transpose_nopiv
  private :: solve4_nopiv_5, lu_factor4_nopiv, lu_solve4_5_nopiv, lu_solve4_transpose_4_nopiv
  public :: lw_sc_scat_ric, lw_sc_scat_ric_solve

contains

  subroutine lw_sc_scat_ric(nlay, nlev, Tl, pl, pe, tau_e, ssa, gg, Tint, lw_up, lw_down, lw_net, olr)
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
    call lw_sc_scat_ric_solve(nlay, nlev, be(:), be_int, tau_e(:), ssa(:), gg(:), lw_up(:), lw_down(:))

    !! Net lw fluxes at each level
    lw_net(:) = lw_up(:) - lw_down(:)

    !! Outgoing Longwave Radiation (olr)
    olr = lw_up(1)

  end subroutine lw_sc_scat_ric

  subroutine lw_sc_scat_ric_solve(nlay, nlev, be, be_int, tau_in, w_in, g_in, flx_up, flx_down)
    implicit none

    !! Input variables
    integer, intent(in) :: nlay, nlev
    real(dp), dimension(nlev), intent(in) :: be, tau_in
    real(dp), dimension(nlay), intent(in) :: w_in, g_in
    real(dp), intent(in) :: be_int

    !! Output variables
    real(dp), dimension(nlev), intent(out) :: flx_up, flx_down

    !! Work variables and arrays
    integer :: k, l, i, j, m, p
    real(dp) :: del, sw, fc, se, so

    real(dp), dimension(nlay) :: dtau, w0
    real(dp), dimension(0:nstr-1, nlay) :: chi
    real(dp), dimension(nlev) :: ext, sca, ssa_lev, q
    real(dp), dimension(0:nstr-1, nlev) :: chi_lev
    logical, dimension(nlev) :: noscat

    real(dp), dimension(0:nstr-1, nmu) :: Pl            ! P_l(u_m), upward directions
    real(dp), dimension(0:nstr-1, nmu, nmu) :: PPw      ! 1/2 (2l+1) P_l(u_i) P_l(u_j) w_j
    real(dp), dimension(nmu, nmu, nlev) :: Fm, Gm       ! forward / backward scattering blocks

    real(dp), dimension(nmu, nlay) :: edel, psi_u, psi_0

    real(dp), dimension(nmu, nmu, nlev) :: Rr, Tt
    real(dp), dimension(nmu, nlev) :: ss, tv
    real(dp), dimension(nmu, nmu) :: Wm, Wi, Vm, Zm, Aa, Bg, Xm, GR, FR
    real(dp), dimension(nmu, nmu+1) :: Ym
    real(dp), dimension(nmu) :: c1, c2, Gs, Fs, e, pu, p0
    real(dp), dimension(nmu) :: dn, up

    !! Delta-M scaled optical depths, albedo and Henyey-Greenstein Legendre moments, f = g**nstr
    do k = 1, nlay
      if (g_in(k) >= 1.0e-6_dp) then
        fc = g_in(k)**nstr
      else
        fc = 0.0_dp
      end if
      w0(k) = w_in(k)*(1.0_dp - fc)/(1.0_dp - fc*w_in(k))
      dtau(k) = (1.0_dp - w_in(k)*fc)*(tau_in(k+1) - tau_in(k))
      chi(0,k) = 1.0_dp
      do l = 1, nstr-1
        chi(l,k) = (g_in(k)**l - fc)/(1.0_dp - fc)
      end do
    end do

    !! Layer -> level single scattering properties (identical to lw_sc_scat)
    ext(:) = 0.0_dp
    sca(:) = 0.0_dp
    chi_lev(:,:) = 0.0_dp
    do k = 1, nlay
      sw = w0(k)*dtau(k)
      ext(k) = ext(k) + dtau(k)
      ext(k+1) = ext(k+1) + dtau(k)
      sca(k) = sca(k) + sw
      sca(k+1) = sca(k+1) + sw
      chi_lev(:,k) = chi_lev(:,k) + sw*chi(:,k)
      chi_lev(:,k+1) = chi_lev(:,k+1) + sw*chi(:,k)
    end do
    do k = 1, nlev
      if (ext(k) > 0.0_dp) then
        ssa_lev(k) = sca(k)/ext(k)
      else
        ssa_lev(k) = 0.0_dp
      end if
      if (sca(k) > 0.0_dp) then
        chi_lev(:,k) = chi_lev(:,k)/sca(k)
      else
        chi_lev(:,k) = 0.0_dp
      end if
      chi_lev(0,k) = 1.0_dp
      q(k) = (1.0_dp - ssa_lev(k))*be(k)
      noscat(k) = (ssa_lev(k) <= 0.0_dp)
    end do

    !! Legendre polynomials at the upward directions (downward: P_l(-u) = (-1)^l P_l(u))
    Pl(0,:) = 1.0_dp
    if (nstr > 1) Pl(1,:) = uarr(:)
    do l = 1, nstr-2
      Pl(l+1,:) = (real(2*l+1,dp)*uarr(:)*Pl(l,:) - real(l,dp)*Pl(l-1,:))/real(l+1,dp)
    end do
    do j = 1, nmu
      do i = 1, nmu
        do l = 0, nstr-1
          PPw(l,i,j) = 0.5_dp*real(2*l+1,dp)*Pl(l,i)*Pl(l,j)*w(j)
        end do
      end do
    end do

    !! Forward (F) and backward (G) scattering blocks: even/odd Legendre split
    !!   F = sum_even + sum_odd,  G = sum_even - sum_odd
    do k = 1, nlev
      if (noscat(k)) then
        Fm(:,:,k) = 0.0_dp
        Gm(:,:,k) = 0.0_dp
        cycle
      end if
      do j = 1, nmu
        do i = 1, nmu
          se = 0.0_dp
          so = 0.0_dp
          do l = 0, nstr-1, 2
            se = se + chi_lev(l,k)*PPw(l,i,j)
          end do
          do l = 1, nstr-1, 2
            so = so + chi_lev(l,k)*PPw(l,i,j)
          end do
          Fm(i,j,k) = ssa_lev(k)*(se + so)
          Gm(i,j,k) = ssa_lev(k)*(se - so)
        end do
      end do
    end do

    !! Linear short characteristic coefficients per layer and stream (as lw_sc_scat)
    do k = 1, nlay
      do m = 1, nmu
        del = dtau(k)/uarr(m)
        edel(m,k) = exp(-del)
        if (del < 1.0e-6_dp) then
          psi_u(m,k) = del*(0.5_dp - del*(1.0_dp/3.0_dp - del*0.125_dp))
          psi_0(m,k) = del*(0.5_dp - del*(1.0_dp/6.0_dp - del/24.0_dp))
        else
          psi_u(m,k) = (1.0_dp - (1.0_dp + del)*edel(m,k))/del
          psi_0(m,k) = (1.0_dp - edel(m,k)) - psi_u(m,k)
        end if
      end do
    end do

    !! ---- Upward Riccati sweep: u_k = R_k d_k + s_k ----
    !! Lower boundary condition - internal heat definition: I_up(mu) = I_down(mu) + be_int
    Rr(:,:,nlev) = 0.0_dp
    do m = 1, nmu
      Rr(m,m,nlev) = 1.0_dp
    end do
    ss(:,nlev) = be_int

    do k = nlev-1, 1, -1
      e(:) = edel(:,k)
      pu(:) = psi_u(:,k)
      p0(:) = psi_0(:,k)

      if (noscat(k) .and. noscat(k+1)) then
        !! Clear layer: F = G = 0 -> W = I, Z = E R, X = I, R_k = E R_{k+1} E (diagonal scaling)
        do i = 1, nmu
          c1(i) = pu(i)*q(k) + p0(i)*q(k+1)
          c2(i) = e(i)*ss(i,k+1) + pu(i)*q(k+1) + p0(i)*q(k)
        end do
        do j = 1, nmu
          do i = 1, nmu
            Rr(i,j,k) = e(i)*Rr(i,j,k+1)*e(j)
            Tt(i,j,k) = 0.0_dp
          end do
          Tt(j,j,k) = e(j)
        end do
        do i = 1, nmu
          sw = 0.0_dp
          do p = 1, nmu
            sw = sw + Rr(i,p,k+1)*c1(p)
          end do
          ss(i,k) = e(i)*sw + c2(i)
          tv(i,k) = c1(i)
        end do
        cycle
      end if

      !! GR = G_{k+1} R_{k+1},  FR = F_{k+1} R_{k+1},  Gs = G_{k+1} s_{k+1},  Fs = F_{k+1} s_{k+1}
      do j = 1, nmu
        do i = 1, nmu
          sw = 0.0_dp
          se = 0.0_dp
          do p = 1, nmu
            sw = sw + Gm(i,p,k+1)*Rr(p,j,k+1)
            se = se + Fm(i,p,k+1)*Rr(p,j,k+1)
          end do
          GR(i,j) = sw
          FR(i,j) = se
        end do
      end do
      do i = 1, nmu
        sw = 0.0_dp
        se = 0.0_dp
        do p = 1, nmu
          sw = sw + Gm(i,p,k+1)*ss(p,k+1)
          se = se + Fm(i,p,k+1)*ss(p,k+1)
        end do
        Gs(i) = sw
        Fs(i) = se
      end do

      !! W, V, c1, c2
      do j = 1, nmu
        do i = 1, nmu
          Wm(i,j) = -p0(i)*(Fm(i,j,k+1) + GR(i,j))
          Vm(i,j) = e(i)*Rr(i,j,k+1) + pu(i)*(FR(i,j) + Gm(i,j,k+1))
        end do
        Wm(j,j) = Wm(j,j) + 1.0_dp
      end do
      do i = 1, nmu
        c1(i) = pu(i)*q(k) + p0(i)*q(k+1) + p0(i)*Gs(i)
        c2(i) = e(i)*ss(i,k+1) + pu(i)*Fs(i) + pu(i)*q(k+1) + p0(i)*q(k)
      end do

      !! Factor W once. Obtain Z = V W^-1 from W^T Z^T = V^T,
      !! retaining the factors for the later T and t solves.
      if (nmu == 4) then
        call lu_factor4_nopiv(Wm)
      else
        call lu_factor_nopiv(Wm)
      end if
      do j = 1, nmu
        do i = 1, nmu
          Wi(i,j) = Vm(j,i)
        end do
      end do
      if (nmu == 4) then
        call lu_solve4_transpose_4_nopiv(Wm, Wi)
      else
        call lu_solve_transpose_nopiv(nmu, Wm, Wi)
      end if

      !! Transpose Z^T into Z; Aa = E + P F_k, Bg = P G_k
      do j = 1, nmu
        do i = 1, nmu
          Zm(i,j) = Wi(j,i)
          Aa(i,j) = pu(i)*Fm(i,j,k)
          Bg(i,j) = pu(i)*Gm(i,j,k)
        end do
        Aa(j,j) = Aa(j,j) + e(j)
      end do

      !! X = I - Q F_k - Z Bg,  [Y | rhs] = [Q G_k + Z Aa | Z c1 + c2]
      do j = 1, nmu
        do i = 1, nmu
          sw = 0.0_dp
          se = 0.0_dp
          do p = 1, nmu
            sw = sw + Zm(i,p)*Bg(p,j)
            se = se + Zm(i,p)*Aa(p,j)
          end do
          Xm(i,j) = -p0(i)*Fm(i,j,k) - sw
          Ym(i,j) = p0(i)*Gm(i,j,k) + se
        end do
        Xm(j,j) = Xm(j,j) + 1.0_dp
      end do
      do i = 1, nmu
        sw = 0.0_dp
        do p = 1, nmu
          sw = sw + Zm(i,p)*c1(p)
        end do
        Ym(i,nmu+1) = sw + c2(i)
      end do

      !! R_k = X^-1 Y,  s_k = X^-1 (Z c1 + c2)
      if (nmu == 4) then
        call solve4_nopiv_5(Xm, Ym)
      else
        call solve_nopiv(nmu+1, Xm, Ym)
      end if
      Rr(:,:,k) = Ym(:,1:nmu)
      ss(:,k) = Ym(:,nmu+1)

      !! T_k = W^-1 (Aa + Bg R_k),  t_k = W^-1 (Bg s_k + c1)
      do j = 1, nmu
        do i = 1, nmu
          sw = 0.0_dp
          do p = 1, nmu
            sw = sw + Bg(i,p)*Rr(p,j,k)
          end do
          GR(i,j) = Aa(i,j) + sw
        end do
      end do
      do i = 1, nmu
        sw = 0.0_dp
        do p = 1, nmu
          sw = sw + Bg(i,p)*ss(p,k)
        end do
        Gs(i) = sw + c1(i)
      end do
      Ym(:,1:nmu) = GR(:,:)
      Ym(:,nmu+1) = Gs(:)
      if (nmu == 4) then
        call lu_solve4_5_nopiv(Wm, Ym)
      else
        call lu_solve_nopiv(nmu+1, Wm, Ym)
      end if
      Tt(:,:,k) = Ym(:,1:nmu)
      tv(:,k) = Ym(:,nmu+1)

    end do

    !! ---- Downward sweep ----
    !! Top boundary condition - 0 intensity downward from top boundary
    dn(:) = 0.0_dp
    do k = 1, nlev
      if (k > 1) then
        ! d_k = T_{k-1} d_{k-1} + t_{k-1}
        do i = 1, nmu
          sw = tv(i,k-1)
          do p = 1, nmu
            sw = sw + Tt(i,p,k-1)*up(p)     ! up() holds d_{k-1} here
          end do
          c1(i) = sw
        end do
        dn(:) = c1(:)
      end if
      ! u_k = R_k d_k + s_k
      do i = 1, nmu
        sw = ss(i,k)
        do p = 1, nmu
          sw = sw + Rr(i,p,k)*dn(p)
        end do
        c2(i) = sw
      end do
      flx_down(k) = 0.0_dp
      flx_up(k) = 0.0_dp
      do m = 1, nmu
        flx_down(k) = flx_down(k) + dn(m)*w(m)*uarr(m)
        flx_up(k) = flx_up(k) + c2(m)*w(m)*uarr(m)
      end do
      flx_down(k) = twopi*flx_down(k)
      flx_up(k) = twopi*flx_up(k)
      up(:) = dn(:)
    end do

  end subroutine lw_sc_scat_ric_solve

  !! Fully unrolled kernels for the default four ordinates per hemisphere.
  !! The generic routines below remain as a fallback if nmu is changed.
  subroutine solve4_nopiv_5(A, B)
    implicit none
    real(dp), dimension(:,:), intent(inout) :: A, B
    real(dp) :: fac

    fac = A(2,1)/A(1,1)
    A(2,2) = A(2,2) - fac*A(1,2)
    A(2,3) = A(2,3) - fac*A(1,3)
    A(2,4) = A(2,4) - fac*A(1,4)
    B(2,1:5) = B(2,1:5) - fac*B(1,1:5)
    fac = A(3,1)/A(1,1)
    A(3,2) = A(3,2) - fac*A(1,2)
    A(3,3) = A(3,3) - fac*A(1,3)
    A(3,4) = A(3,4) - fac*A(1,4)
    B(3,1:5) = B(3,1:5) - fac*B(1,1:5)
    fac = A(4,1)/A(1,1)
    A(4,2) = A(4,2) - fac*A(1,2)
    A(4,3) = A(4,3) - fac*A(1,3)
    A(4,4) = A(4,4) - fac*A(1,4)
    B(4,1:5) = B(4,1:5) - fac*B(1,1:5)

    fac = A(3,2)/A(2,2)
    A(3,3) = A(3,3) - fac*A(2,3)
    A(3,4) = A(3,4) - fac*A(2,4)
    B(3,1:5) = B(3,1:5) - fac*B(2,1:5)
    fac = A(4,2)/A(2,2)
    A(4,3) = A(4,3) - fac*A(2,3)
    A(4,4) = A(4,4) - fac*A(2,4)
    B(4,1:5) = B(4,1:5) - fac*B(2,1:5)

    fac = A(4,3)/A(3,3)
    A(4,4) = A(4,4) - fac*A(3,4)
    B(4,1:5) = B(4,1:5) - fac*B(3,1:5)

    B(4,1:5) = B(4,1:5)/A(4,4)
    B(3,1:5) = B(3,1:5) - A(3,4)*B(4,1:5)
    B(3,1:5) = B(3,1:5)/A(3,3)
    B(2,1:5) = B(2,1:5) - A(2,3)*B(3,1:5)
    B(2,1:5) = B(2,1:5) - A(2,4)*B(4,1:5)
    B(2,1:5) = B(2,1:5)/A(2,2)
    B(1,1:5) = B(1,1:5) - A(1,2)*B(2,1:5)
    B(1,1:5) = B(1,1:5) - A(1,3)*B(3,1:5)
    B(1,1:5) = B(1,1:5) - A(1,4)*B(4,1:5)
    B(1,1:5) = B(1,1:5)/A(1,1)

  end subroutine solve4_nopiv_5

  subroutine lu_factor4_nopiv(A)
    implicit none
    real(dp), dimension(:,:), intent(inout) :: A

    A(2,1) = A(2,1)/A(1,1)
    A(2,2:4) = A(2,2:4) - A(2,1)*A(1,2:4)
    A(3,1) = A(3,1)/A(1,1)
    A(3,2:4) = A(3,2:4) - A(3,1)*A(1,2:4)
    A(4,1) = A(4,1)/A(1,1)
    A(4,2:4) = A(4,2:4) - A(4,1)*A(1,2:4)

    A(3,2) = A(3,2)/A(2,2)
    A(3,3:4) = A(3,3:4) - A(3,2)*A(2,3:4)
    A(4,2) = A(4,2)/A(2,2)
    A(4,3:4) = A(4,3:4) - A(4,2)*A(2,3:4)

    A(4,3) = A(4,3)/A(3,3)
    A(4,4) = A(4,4) - A(4,3)*A(3,4)

  end subroutine lu_factor4_nopiv

  subroutine lu_solve4_5_nopiv(A, B)
    implicit none
    real(dp), dimension(:,:), intent(in) :: A
    real(dp), dimension(:,:), intent(inout) :: B

    B(2,1:5) = B(2,1:5) - A(2,1)*B(1,1:5)
    B(3,1:5) = B(3,1:5) - A(3,1)*B(1,1:5)
    B(3,1:5) = B(3,1:5) - A(3,2)*B(2,1:5)
    B(4,1:5) = B(4,1:5) - A(4,1)*B(1,1:5)
    B(4,1:5) = B(4,1:5) - A(4,2)*B(2,1:5)
    B(4,1:5) = B(4,1:5) - A(4,3)*B(3,1:5)

    B(4,1:5) = B(4,1:5)/A(4,4)
    B(3,1:5) = B(3,1:5) - A(3,4)*B(4,1:5)
    B(3,1:5) = B(3,1:5)/A(3,3)
    B(2,1:5) = B(2,1:5) - A(2,3)*B(3,1:5)
    B(2,1:5) = B(2,1:5) - A(2,4)*B(4,1:5)
    B(2,1:5) = B(2,1:5)/A(2,2)
    B(1,1:5) = B(1,1:5) - A(1,2)*B(2,1:5)
    B(1,1:5) = B(1,1:5) - A(1,3)*B(3,1:5)
    B(1,1:5) = B(1,1:5) - A(1,4)*B(4,1:5)
    B(1,1:5) = B(1,1:5)/A(1,1)

  end subroutine lu_solve4_5_nopiv

  subroutine lu_solve4_transpose_4_nopiv(A, B)
    implicit none
    real(dp), dimension(:,:), intent(in) :: A
    real(dp), dimension(:,:), intent(inout) :: B

    B(1,1:4) = B(1,1:4)/A(1,1)
    B(2,1:4) = B(2,1:4) - A(1,2)*B(1,1:4)
    B(2,1:4) = B(2,1:4)/A(2,2)
    B(3,1:4) = B(3,1:4) - A(1,3)*B(1,1:4)
    B(3,1:4) = B(3,1:4) - A(2,3)*B(2,1:4)
    B(3,1:4) = B(3,1:4)/A(3,3)
    B(4,1:4) = B(4,1:4) - A(1,4)*B(1,1:4)
    B(4,1:4) = B(4,1:4) - A(2,4)*B(2,1:4)
    B(4,1:4) = B(4,1:4) - A(3,4)*B(3,1:4)
    B(4,1:4) = B(4,1:4)/A(4,4)

    B(3,1:4) = B(3,1:4) - A(4,3)*B(4,1:4)
    B(2,1:4) = B(2,1:4) - A(3,2)*B(3,1:4)
    B(2,1:4) = B(2,1:4) - A(4,2)*B(4,1:4)
    B(1,1:4) = B(1,1:4) - A(2,1)*B(2,1:4)
    B(1,1:4) = B(1,1:4) - A(3,1)*B(3,1:4)
    B(1,1:4) = B(1,1:4) - A(4,1)*B(4,1:4)

  end subroutine lu_solve4_transpose_4_nopiv

  !! Solve A X = B in place for nr right-hand sides (B overwritten by X), Gaussian elimination,
  !! no pivoting (W and X are diagonally dominant M-matrix Schur complements for p >= 0)
  subroutine solve_nopiv(nr, A, B)
    implicit none
    integer, intent(in) :: nr
    real(dp), dimension(nmu, nmu), intent(inout) :: A
    real(dp), dimension(nmu, nr), intent(inout) :: B

    integer :: i, j, c
    real(dp) :: fac, piv

    do j = 1, nmu
      piv = 1.0_dp/A(j,j)
      do i = j+1, nmu
        fac = A(i,j)*piv
        if (fac /= 0.0_dp) then
          do c = j+1, nmu
            A(i,c) = A(i,c) - fac*A(j,c)
          end do
          do c = 1, nr
            B(i,c) = B(i,c) - fac*B(j,c)
          end do
        end if
      end do
    end do
    do j = nmu, 1, -1
      piv = 1.0_dp/A(j,j)
      do c = 1, nr
        fac = B(j,c)
        do i = j+1, nmu
          fac = fac - A(j,i)*B(i,c)
        end do
        B(j,c) = fac*piv
      end do
    end do

  end subroutine solve_nopiv

  !! In-place no-pivot LU factorization, A = L U. The strict lower triangle
  !! stores the unit-diagonal L multipliers and the upper triangle stores U.
  subroutine lu_factor_nopiv(A)
    implicit none
    real(dp), dimension(nmu, nmu), intent(inout) :: A

    integer :: i, j, c
    real(dp) :: fac, piv

    do j = 1, nmu-1
      piv = 1.0_dp/A(j,j)
      do i = j+1, nmu
        fac = A(i,j)*piv
        A(i,j) = fac
        if (fac /= 0.0_dp) then
          do c = j+1, nmu
            A(i,c) = A(i,c) - fac*A(j,c)
          end do
        end if
      end do
    end do

  end subroutine lu_factor_nopiv

  !! Solve L U X = B using a previously computed no-pivot LU factorization.
  subroutine lu_solve_nopiv(nr, A, B)
    implicit none
    integer, intent(in) :: nr
    real(dp), dimension(nmu, nmu), intent(in) :: A
    real(dp), dimension(nmu, nr), intent(inout) :: B

    integer :: i, j, c
    real(dp) :: fac, piv

    do i = 2, nmu
      do c = 1, nr
        fac = B(i,c)
        do j = 1, i-1
          fac = fac - A(i,j)*B(j,c)
        end do
        B(i,c) = fac
      end do
    end do
    do i = nmu, 1, -1
      piv = 1.0_dp/A(i,i)
      do c = 1, nr
        fac = B(i,c)
        do j = i+1, nmu
          fac = fac - A(i,j)*B(j,c)
        end do
        B(i,c) = fac*piv
      end do
    end do

  end subroutine lu_solve_nopiv

  !! Solve U^T L^T X = B using a previously computed no-pivot LU factorization.
  subroutine lu_solve_transpose_nopiv(nr, A, B)
    implicit none
    integer, intent(in) :: nr
    real(dp), dimension(nmu, nmu), intent(in) :: A
    real(dp), dimension(nmu, nr), intent(inout) :: B

    integer :: i, j, c
    real(dp) :: fac, piv

    do i = 1, nmu
      piv = 1.0_dp/A(i,i)
      do c = 1, nr
        fac = B(i,c)
        do j = 1, i-1
          fac = fac - A(j,i)*B(j,c)
        end do
        B(i,c) = fac*piv
      end do
    end do
    do i = nmu-1, 1, -1
      do c = 1, nr
        fac = B(i,c)
        do j = i+1, nmu
          fac = fac - A(j,i)*B(j,c)
        end do
        B(i,c) = fac
      end do
    end do

  end subroutine lu_solve_transpose_nopiv

end module lw_sc_scat_ric_mod
