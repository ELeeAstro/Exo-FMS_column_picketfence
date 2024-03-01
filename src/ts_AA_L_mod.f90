!!!
! Elspeth KH Lee - Aug 2023 : Initial version
! sw: alpha-4SDA - 4 stream analytic spherical harmonic method with doubling adding method (Zhang & Li 2013)
! lw: Absorption Approximation following Li (2002) - This is the linear in tau method
!     Pros: Very fast method with LW scattering approximation, no matrix inversions
!     Cons: Not technically multiple scattering (can be quite innacurate at even moderate albedo)
!!!

module ts_AA_L_mod
  use, intrinsic :: iso_fortran_env
  implicit none

  !! Precision variables
  integer, parameter :: dp = REAL64

  !! Required constants
  real(dp), parameter :: pi = 4.0_dp * atan(1.0_dp)
  real(dp), parameter :: twopi = 2.0_dp * pi
  real(dp), parameter :: sb = 5.670374419e-8_dp

  !! Legendre quadrature for 1 nodes
  ! integer, parameter :: nmu = 1
  ! real(dp), dimension(nmu), parameter :: uarr = (/1.0_dp/1.6487213_dp/)
  ! real(dp), dimension(nmu), parameter :: w = (/1.0_dp/)
  ! real(dp), dimension(nmu), parameter :: wuarr = uarr * w

  !! Legendre quadrature for 2 nodes
  integer, parameter :: nmu = 2
  real(dp), dimension(nmu), parameter :: uarr = (/0.21132487_dp, 0.78867513_dp/)
  real(dp), dimension(nmu), parameter :: w = (/0.5_dp, 0.5_dp/)
  real(dp), dimension(nmu), parameter :: wuarr = uarr * w

  !! Use Two-Term HG function for sw
  logical, parameter :: TTHG = .False.

  public :: ts_AA_L
  private :: lw_AA_L, sw_SDA, linear_log_interp, bezier_interp, &
      & matinv4, ludcmp, lubksb, inv_LU

contains

  subroutine ts_AA_L(Bezier, nlay, nlev, Tl, pl, pe, tau_Ve, tau_IRe, mu_z, F0, Tint, AB, Beta_V, Beta_IR, &
    & sw_a, sw_g, lw_a, lw_g, sw_a_surf, lw_a_surf, net_F, olr, asr)
    implicit none

    !! Input variables
    logical, intent(in) :: Bezier
    integer, intent(in) :: nlay, nlev
    real(dp), dimension(nlay), intent(in) :: Tl, pl
    real(dp), dimension(nlev), intent(in) :: pe, mu_z
    real(dp), dimension(nlev,3), intent(in) :: tau_Ve
    real(dp), dimension(nlev,2), intent(in) :: tau_IRe
    real(dp), dimension(nlay,3), intent(in) :: sw_a, sw_g, lw_a, lw_g
    real(dp), dimension(3), intent(in) :: Beta_V
    real(dp), dimension(2), intent(in) :: Beta_IR
    real(dp), intent(in) :: F0, Tint, AB, sw_a_surf, lw_a_surf

    !! Output variables
    real(dp), intent(out) :: olr, asr
    real(dp), dimension(nlev), intent(out) :: net_F

    !! Work variables
    integer :: i, b
    real(dp) :: Finc, be_int, Finc_b, be_int_b
    real(dp), dimension(nlev) :: lpe
    real(dp), dimension(nlay) :: lTl, lpl
    real(dp), dimension(nlev) :: Te, be, be_b
    real(dp), dimension(nlev) :: sw_down, sw_up, lw_down, lw_up
    real(dp), dimension(nlev,3) :: sw_down_b, sw_up_b
    real(dp), dimension(nlev,2) :: lw_down_b, lw_up_b
    real(dp), dimension(nlev) :: lw_net, sw_net

    !! Find temperature at layer edges through interpolation and extrapolation
    if (Bezier .eqv. .True.) then

      ! Log the layer values and pressure edges for more accurate interpolation
      lTl(:) = log10(Tl(:))
      lpl(:) = log10(pl(:))
      lpe(:) = log10(pe(:))

      ! Perform interpolation using Bezier peicewise polynomial interpolation
      do i = 2, nlay-1
        call bezier_interp(lpl(i-1:i+1), lTl(i-1:i+1), 3, lpe(i), Te(i))
        Te(i) = 10.0_dp**(Te(i))
        !print*, i, pl(i), pl(i-1), pe(i), Tl(i-1), Tl(i), Te(i)
      end do
      call bezier_interp(lpl(nlay-2:nlay), lTl(nlay-2:nlay), 3, lpe(nlay), Te(nlay))
      Te(nlay) = 10.0_dp**(Te(nlay))
    else
      ! Perform interpolation using linear interpolation
      do i = 2, nlay
        call linear_log_interp(pe(i), pl(i-1), pl(i), Tl(i-1), Tl(i), Te(i))
        !print*, i, pl(i), pl(i-1), pe(i), Tl(i-1), Tl(i), Te(i)
      end do
    end if

    ! Edges are linearly interpolated
    Te(1) = 10.0_dp**(log10(Tl(1)) + (log10(pe(1)/pe(2))/log10(pl(1)/pe(2))) * log10(Tl(1)/Te(2)))
    Te(nlev) = 10.0_dp**(log10(Tl(nlay)) + (log10(pe(nlev)/pe(nlay))/log10(pl(nlay)/pe(nlay))) * log10(Tl(nlay)/Te(nlay)))

    !! Shortwave flux calculation
    if (mu_z(nlev) > 0.0_dp) then
      Finc = (1.0_dp - AB) * F0
      sw_down(:) = 0.0_dp
      sw_up(:) = 0.0_dp
      do b = 1, 3
        Finc_b = Finc * Beta_V(b)
        call sw_SDA(nlay, nlev, Finc_b, tau_Ve(:,b), mu_z(nlev), sw_a(:,b), sw_g(:,b), sw_a_surf, &
          & sw_down_b(:,b), sw_up_b(:,b))
        sw_down(:) = sw_down(:) + sw_down_b(:,b)
        sw_up(:) = sw_up(:) + sw_up_b(:,b)
      end do
    else
      sw_down(:) = 0.0_dp
      sw_up(:) = 0.0_dp
    end if

    !! Longwave two-stream flux calculation
    be(:) = (sb * Te(:)**4)/pi  ! Integrated planck function intensity at levels
    be_int = (sb * Tint**4)/pi ! Integrated planck function intensity for internal temperature
    lw_up(:) = 0.0_dp
    lw_down(:) = 0.0_dp
    do b = 1, 2
      be_b(:) = be(:) * Beta_IR(b)
      be_int_b = be_int * Beta_IR(b)
      call lw_AA_L(nlay, nlev, be_b, be_int_b, tau_IRe(:,b), lw_a(:,b), lw_g(:,b), lw_a_surf, &
        & lw_up_b(:,b), lw_down_b(:,b))
      lw_up(:) = lw_up(:) + lw_up_b(:,b)
      lw_down(:) = lw_down(:) + lw_down_b(:,b)
    end do

    !! Net fluxes at each level
    lw_net(:) = lw_up(:) - lw_down(:)
    sw_net(:) = sw_up(:) - sw_down(:)
    net_F(:) = lw_net(:) + sw_net(:)

    !! Output olr
    olr = lw_up(1)

    !! Output asr
    asr = sw_down(1) - sw_up(1)

  end subroutine ts_AA_L

  subroutine lw_AA_L(nlay, nlev, be, be_int, tau_IRe, w0, g0, lw_a_surf, lw_up, lw_down)
    implicit none

    !! Input variables
    integer, intent(in) :: nlay, nlev
    real(dp), dimension(nlev), intent(in) :: be, tau_IRe
    real(dp), dimension(nlay), intent(in) :: w0, g0
    real(dp), intent(in) :: be_int, lw_a_surf

    !! Output variables
    real(dp), dimension(nlev), intent(out) :: lw_up, lw_down

    !! Work variables and arrays
    integer :: k, m
    real(dp), dimension(nlay) :: w0_d, g0_d, fc
    real(dp), dimension(nlay) :: dtau, Blin, eps, dtau_a
    real(dp), dimension(nlay) :: T, zeta
    real(dp), dimension(nlev) :: lw_up_g, lw_down_g

    real(dp), dimension(nlay) :: sigma_sq, pmom2, c
    integer, parameter :: nstr = nmu*2

    !! Calculate dtau in each layer
    dtau(:) = tau_IRe(2:nlev) - tau_IRe(1:nlay)

    !! Delta-M+ scaling (Following DISORT: Lin et al. 2018)
    !! Assume HG phase function for scaling
    fc(:) = g0(:)**(nstr)
    pmom2(:) = g0(:)**(nstr+1)

    where (fc(:) /=  pmom2(:))
      sigma_sq(:) = real((nstr+1)**2 - nstr**2,dp) / &
      & ( log(fc(:)**2/pmom2(:)**2) )
      c(:) = exp(real(nstr**2,dp)/(2.0_dp*sigma_sq(:)))
      fc(:) = c(:)*fc(:)

      w0_d(:) = w0(:)*((1.0_dp - fc(:))/(1.0_dp - fc(:)*w0(:)))
      dtau(:) = (1.0_dp - w0(:)*fc(:))*dtau(:)

    elsewhere
      w0_d(:) = w0(:)
    end where

    g0_d(:) = g0(:)

    !! Linear B with tau function
    Blin(:) = be(2:nlev)-be(1:nlay)

    !! modified co-albedo epsilon
    eps(:) = sqrt((1.0_dp - w0_d(:))*(1.0_dp - g0_d(:)*w0_d(:)))
    !eps(:) = (1.0_dp - w0(:))

    !! Absorption/modified optical depth for transmission function
    dtau_a(:) = eps(:)*dtau(:)

    ! Zero the total flux arrays
    lw_up(:) = 0.0_dp
    lw_down(:) = 0.0_dp

    !! Start loops to integrate in mu space
    do m = 1, nmu

      !! Efficency variables
      T(:) = exp(-dtau_a(:)/uarr(m))
      zeta(:) = (Blin(:)*uarr(m)/dtau_a(:))

      !! Begin two-stream loops
      !! Perform downward loop first
      ! Top boundary condition - 0 flux downward from top boundary
      lw_down_g(1) = 0.0_dp
      do k = 1, nlay
        if (dtau_a(k) > 1e-4_dp) then
          lw_down_g(k+1) = lw_down_g(k)*T(k) + &
            & be(k+1) - zeta(k) - (be(k) - zeta(k))*T(k)
        else
          lw_down_g(k+1) = lw_down_g(k)*T(k) + be(k)*(1.0_dp - T(k))
        end if
      end do

      !! Perform upward loop
      ! Lower boundary condition - internal heat definition Fint = F_down - F_up
      ! here the lw_a_surf is assumed to be = 1 as per the definition
      ! here we use the same condition but use intensity units to be consistent
      lw_up_g(nlev) = lw_down_g(nlev) + be_int

      do k = nlay, 1, -1
        if (dtau_a(k) > 1e-4_dp) then
          lw_up_g(k) = lw_up_g(k+1)*T(k) + &
            & be(k) + zeta(k) - (be(k+1) + zeta(k))*T(k)
        else
          lw_up_g(k) = lw_up_g(k+1)*T(k) + be(k+1)*(1.0_dp - T(k))
        end if
      end do

      !! Sum up flux arrays with Gaussian quadrature weights and points for this mu stream
      lw_down(:) = lw_down(:) + lw_down_g(:) * wuarr(m)
      lw_up(:) = lw_up(:) + lw_up_g(:) * wuarr(m)

    end do

    !! The flux is the integrated intensity * 2pi
    lw_down(:) = twopi * lw_down(:)
    lw_up(:) = twopi * lw_up(:)

  end subroutine lw_AA_L

  subroutine sw_SDA(nlay, nlev, Finc, tau_Ve, mu_z, ww, gg, w_surf, sw_down, sw_up)
    implicit none

    !! Input variables
    integer, intent(in) :: nlay, nlev
    real(dp), intent(in) :: Finc, mu_z, w_surf
    real(dp), dimension(nlev), intent(in) :: tau_Ve
    real(dp), dimension(nlay), intent(in) :: ww, gg

    !! Output variables
    real(dp), dimension(nlev), intent(out) :: sw_down, sw_up

    !! Work variables
    integer :: k
    real(dp), dimension(nlay) :: w0, dtau, hg
    real(dp), dimension(nlev) :: tau, T
    real(dp) :: mu1, mu2, mu3, f0
    real(dp) :: om0, om1, om2, om3
    real(dp) :: a0, a1, a2, a3, b0, b1, b2, b3
    real(dp) :: e1, e2
    real(dp) :: beta, gam, k1, k2, R1, R2, P1, P2, Q1, Q2
    real(dp) :: eta0, eta1, eta2, eta3, del0, del1, del2, del3, delta
    real(dp) :: z1p, z1m, z2p, z2m
    real(dp) :: phi1p, phi1m, phi2p, phi2m
    real(dp) :: Cphi1p, Cphi1m, Cphi2p, Cphi2m

    real(dp), dimension(4) :: H1, H2, H3, H4
    real(dp), dimension(4) :: AA12H1
    real(dp), dimension(4,4) :: AA1, AA2, AA1_i, AA12
    real(dp), dimension(4) :: Fdir, Fdiffa, Fdiffb

    real(dp), dimension(nlay,2) :: Rdir, Tdir
    real(dp), dimension(nlay,2,2) :: Rdiff, Tdiff 

    real(dp), dimension(nlev,2) :: T1k, R1k
    real(dp), dimension(nlev,2,2) :: Rst1k, Rd1k

    real(dp), dimension(nlay) :: fc, sigma_sq, pmom2, c
    integer, parameter :: nstr = 4
    real(dp), parameter :: eps_20 = 1.0e-20_dp

    integer :: l, km1, lp1
    real(dp)  :: c11, c12, c21, c22, pmod, t11, t12, t21, t22, b11, b12, &
           &   b21, b22, uu1, uu2, du1, du2, a11, a12, a21, a22, r11, r12, r21, r22, d1, d2

    real(dp), dimension(nlev) :: scumdtr, reflt, trant
    real(dp), dimension(nlay) :: dtr

    real(dp) :: hg2, alp

    !! If zero albedo across all atmospheric layers then return direct beam only
    if (all(ww(:) <= 1.0e-6_dp)) then
      sw_down(:) = Finc * mu_z * exp(-tau_Ve(:)/mu_z)
      sw_down(nlev) = sw_down(nlev) * (1.0_dp - w_surf) ! The surface flux for surface heating is the amount of flux absorbed by surface
      sw_up(:) = 0.0_dp ! We assume no upward flux here even if surface albedo
      return
    end if

    !! Calculate dtau in each layer
    dtau(:) = tau_Ve(2:) - tau_Ve(1:nlay)

    !! Delta-M+ scaling (Following DISORT: Lin et al. 2018)
    !! Assume HG phase function for scaling
    fc(:) = gg(:)**(nstr)
    pmom2(:) = gg(:)**(nstr+1)

    where (fc(:) /=  pmom2(:))
      sigma_sq(:) = real((nstr+1)**2 - nstr**2,dp) / &
      & ( log(fc(:)**2/pmom2(:)**2) )
      c(:) = exp(real(nstr**2,dp)/(2.0_dp*sigma_sq(:)))
      fc(:) = c(:)*fc(:)

      w0(:) = ww(:)*((1.0_dp - fc(:))/(1.0_dp - fc(:)*ww(:)))
      dtau(:) = (1.0_dp - ww(:)*fc(:))*dtau(:)

    elsewhere
      w0(:) = ww(:)
      fc(:) = 0.0_dp
    end where

    hg(:) = gg(:)

    !! Reform edge optical depths
    tau(1) = tau_Ve(1)
    do k = 1, nlay
      tau(k+1) = tau(k) + dtau(k)
    end do

    !! Start SDA calculation

    !! First find the Reflection and Transmission coefficents (direct and diffuse) for each layer
    do k = 1, nlay

      !! Layer transmission
      dtr(k) = exp(-dtau(k)/mu_z)

      !! Mu moments
      mu1 = mu_z
      mu2 = mu_z**2
      mu3 = mu_z**3

      !! Inverse zenith angle
      f0 = 1.0_dp/mu_z

      !! Omega Legendre polynomial coefficents - scale with delta-M+
      if (hg(k) /= 0.0_dp) then
        if (TTHG .eqv. .False.) then
          ! Use HG phase function
          om0 = 1.0_dp
          om1 = 3.0_dp * (hg(k) - fc(k))/(1.0_dp - fc(k))
          om2 = 5.0_dp * (hg(k)**2 - fc(k))/(1.0_dp - fc(k))
          om3 = 7.0_dp * (hg(k)**3 - fc(k))/(1.0_dp - fc(k))
        else
          ! Use TTHG phase function with default parameters
          hg2 = hg(k)/2.0_dp
          alp = 1.0_dp - hg2**2
          om0 = 1.0_dp
          om1 = 3.0_dp * ((alp*hg(k) + (1.0_dp - alp)*hg2) - fc(k))/(1.0_dp - fc(k))
          om2 = 5.0_dp * ((alp*hg(k)**2 + (1.0_dp - alp)*hg2**2) - fc(k))/(1.0_dp - fc(k))
          om3 = 7.0_dp * ((alp*hg(k)**3 + (1.0_dp - alp)*hg2**3) - fc(k))/(1.0_dp - fc(k))
        end if
      else
        ! Use Rayleigh scattering phase function for isotropic scattering
        om0 = 1.0_dp
        om1 = 0.0_dp
        om2 = 0.5_dp
        om3 = 0.0_dp
      end if

      !! Find the a coefficents
      a0 =  1.0_dp - w0(k)*om0 + eps_20
      a1 =  3.0_dp - w0(k)*om1 + eps_20
      a2 =  5.0_dp - w0(k)*om2 + eps_20
      a3 =  7.0_dp - w0(k)*om3 + eps_20

      !! Find the b coefficents
      b0 = 0.25_dp * w0(k)*om0 * Finc/pi
      b1 = 0.25_dp * w0(k)*om1 * -(mu1) * Finc/pi
      b2 = 0.125_dp * w0(k)*om2 * (3.0_dp * mu2 - 1.0_dp) * Finc/pi
      b3 = 0.125_dp * w0(k)*om3 * (5.0_dp * -mu3 - 3.0_dp*-(mu_z)) * Finc/pi

      !! Find beta and gamma
      beta = a0*a1 + (4.0_dp/9.0_dp)*a0*a3 + (1.0_dp/9.0_dp)*a2*a3
      gam = (1.0_dp/9.0_dp)*a0*a1*a2*a3

      !! Find k values - lambda in Rooney
      k1 = (beta + sqrt((beta**2 - 4.0_dp*gam)))/2.0_dp
      k2 = (beta - sqrt((beta**2 - 4.0_dp*gam)))/2.0_dp

      k1 = sqrt(k1)
      k2 = sqrt(k2)

      !! Find e values
      e1 = exp(-k1*dtau(k))
      e2 = exp(-k2*dtau(k))      
      
      !! Find R, P and Q coefficents 
      !! NOTE: Zhang et al. (2013) has the wrong coefficent definitions
      !! Rooney et al. (2023) has the correct definitions and order in the matrix
      !! So we we use the Rooney definitions, but keep the Zhang notation
      Q1 = -3.0_dp/2.0_dp * (a0*a1/k1 - k1)/a3
      Q2 = -3.0_dp/2.0_dp * (a0*a1/k2 - k2)/a3
      R1 = -a0/k1
      R2 = -a0/k2
      P1 = 0.3125_dp * (a0*a1/k1**2 - 1.0_dp)
      P2 = 0.3125_dp * (a0*a1/k2**2 - 1.0_dp)

      !! Find the delta values
      delta = 9.0_dp*(f0**4 - beta*f0**2 + gam)
      del0 = (a1*b0 - b1*f0)*(a2*a3 - 9.0_dp*f0**2) + 2.0_dp*f0**2*(a3*b2 - 2.0_dp*a3*b0 - 3.0_dp*b3*f0)
      del1 = (a0*b1 - b0*f0)*(a2*a3 - 9.0_dp*f0**2) - 2.0_dp*a0*f0*(a3*b2 - 3.0_dp*b3*f0)
      del2 = (a3*b2 - 3.0_dp*b3*f0)*(a0*a1 - f0**2) - 2.0_dp*a3*f0*(a0*b1 - b0*f0)
      del3 = (a2*b3 - 3.0_dp*b2*f0)*(a0*a1 - f0**2) + f0**2*(6.0_dp*a0*b1 - 4.0_dp*a0*b3 - 6.0_dp*b0*f0)

      !! Find the eta values
      eta0 = del0/delta
      eta1 = del1/delta
      eta2 = 0.625_dp * del2/delta
      eta3 = del3/delta

      !! Find the phi values
      phi1p = 2.0_dp*pi*(0.5_dp + R1 + 5.0_dp/8.0_dp*P1)
      phi1m = 2.0_dp*pi*(0.5_dp - R1 + 5.0_dp/8.0_dp*P1)
      phi2p = 2.0_dp*pi*(0.5_dp + R2 + 5.0_dp/8.0_dp*P2)
      phi2m = 2.0_dp*pi*(0.5_dp - R2 + 5.0_dp/8.0_dp*P2)

      !! Find the Phi values
      Cphi1p = 2.0_dp*pi*(-1.0_dp/8.0_dp + 5.0_dp/8.0_dp*P1 + Q1)
      Cphi1m = 2.0_dp*pi*(-1.0_dp/8.0_dp + 5.0_dp/8.0_dp*P1 - Q1)
      Cphi2p = 2.0_dp*pi*(-1.0_dp/8.0_dp + 5.0_dp/8.0_dp*P2 + Q2)
      Cphi2m = 2.0_dp*pi*(-1.0_dp/8.0_dp + 5.0_dp/8.0_dp*P2 - Q2)

      !! Find the Z values
      z1p = 2.0_dp*pi*(0.5_dp*eta0 + eta1 + eta2)
      z1m = 2.0_dp*pi*(0.5_dp*eta0 - eta1 + eta2)
      z2p = 2.0_dp*pi*(-1.0_dp/8.0_dp*eta0 + eta2 + eta3)
      z2m = 2.0_dp*pi*(-1.0_dp/8.0_dp*eta0 + eta2 - eta3)

      !! Find A1 matrix
      AA1(1,1) = phi1m ; AA1(1,2) = phi1p*e1 ; AA1(1,3) = phi2m ; AA1(1,4) = phi2p*e2
      AA1(2,1) = Cphi1m ; AA1(2,2) = Cphi1p*e1 ; AA1(2,3) = Cphi2m; AA1(2,4) = Cphi2p*e2
      AA1(3,1) = phi1p*e1 ; AA1(3,2) = phi1m ; AA1(3,3) = phi2p*e2 ; AA1(3,4) = phi2m    
      AA1(4,1) = Cphi1p*e1 ; AA1(4,2) = Cphi1m ; AA1(4,3) = Cphi2p*e2 ; AA1(4,4) = Cphi2m

      !! Find H1 vector
      H1(1) = -z1m ;  H1(2) = -z2m ;  H1(3) = -z1p * dtr(k) ; H1(4) = -z2p * dtr(k)

      !! Find A2 matrix
      AA2(1,1) = phi1m*e1 ; AA2(1,2) = phi1p ; AA2(1,3) = phi2m*e2 ; AA2(1,4) = phi2p
      AA2(2,1) = Cphi1m*e1 ; AA2(2,2) = Cphi1p ; AA2(2,3) = Cphi2m*e2; AA2(2,4) = Cphi2p
      AA2(3,1) = phi1p ; AA2(3,2) = phi1m*e1 ; AA2(3,3) = phi2p ; AA2(3,4) = phi2m*e2    
      AA2(4,1) = Cphi1p ; AA2(4,2) = Cphi1m*e1 ; AA2(4,3) = Cphi2p ; AA2(4,4) = Cphi2m*e2

      !! Find H2 vector
      H2(1) = z1m * dtr(k) ;  H2(2) = z2m * dtr(k) ;  H2(3) = z1p ; H2(4) = z2p

      !! Now we need to invert the A1 matrix - Zhang and Li (2013) use a adjugate matrix method 
      !! with some reduction in the matrix order (not 100% sure what they did)
      !! We primarily use the same method but keep the 4x4 layout, with code from the Fortran wiki (matinv4)
      !! We have LU decomposition here as an alternative in case of numerical instability (and for testing)

      AA1_i(:,:) = matinv4(AA1(:,:)) ! Use matrix determinant and adjugate method (faster but can be numerically unstable)
      !call inv_LU(AA1,4,4,AA1_i)    ! Use LU decomposition (slower but probably more stable)

      !! Multiply the AA1_i and AA2 matrices
      AA12(:,:) = matmul(AA2(:,:),AA1_i(:,:))

      !! Multiply the AA12 and H2 matrix
      AA12H1(:) = matmul(AA12(:,:),H1(:))

      !! Firect flux component - now we have the flux array (F(1)(2) = neg flux at lower,F(3)(4) = pos flux at upper)
      Fdir(:) = AA12H1(:) + H2(:)

      !! Store the direct beam reflection and transmission for this layer - normalised by the beam flux
      Rdir(k,1) = Fdir(3)/(Finc*mu_z)
      Rdir(k,2) = Fdir(4)/(Finc*mu_z)

      Tdir(k,1) = Fdir(1)/(Finc*mu_z)
      Tdir(k,2) = Fdir(2)/(Finc*mu_z)

      !! Now find the diffusive flux component

      !! Vector H3 is
      H3(1) = 1.0_dp; H3(2) = 0.0_dp; H3(3) = 0.0_dp; H3(4) = 0.0_dp 

      !! Multiply the AA12 and H3 matrix to get `a' boundary fluxes at layer edges (levels)
      Fdiffa(:) = matmul(AA12(:,:),H3(:))

      !! Vector H4 is
      H4(1) = 0.0_dp; H4(2) = 1.0_dp; H4(3) = 0.0_dp; H4(4) = 0.0_dp 

      !! Multiply the AA12 and H4 matrix to get `b' boundary fluxes at layer edges (levels)
      Fdiffb(:) = matmul(AA12(:,:),H4(:))

      !! Store the diffuse reflection and transmission for this layer - no normalisation
      Rdiff(k,1,1) = Fdiffa(3)
      Rdiff(k,1,2) = Fdiffb(3)
      Rdiff(k,2,2) = Fdiffa(4)
      Rdiff(k,2,1) = Fdiffb(4)
  
      Tdiff(k,1,1) = Fdiffa(1)
      Tdiff(k,1,2) = Fdiffb(1)
      Tdiff(k,2,1) = Fdiffa(2)
      Tdiff(k,2,2) = Fdiffb(2)

    end do

    !! We now have the transmission and reflection coefficents for both the direct and diffuse components
    !! Now we perform the doubling-adding method accros multiple layers

    !! Here we directly copy the code from J. Li 
    !! we can probably make an improvement on this at some point through vectorisation

    !! Do boundary conditons first
    ! Upper
    T1k(1,:) = 0.0_dp
    Rd1k(1,:,:) = 0.0_dp

    ! Lower
    R1k(nlev,1) = w_surf ; R1k(nlev,2) = -w_surf/4.0_dp 
    Rst1k(nlev,1,1) = w_surf ; Rst1k(nlev,1,2) = 0.0_dp
    Rst1k(nlev,2,1) = -w_surf/4.0_dp ; Rst1k(nlev,2,2) = 0.0_dp

    !! Direct beam transmission to level
    T(:) = exp(-tau(:)/mu_z)

    !! Cumulative transmission
    scumdtr(1) = T(1)

    do k = 2, nlev
      km1 = k - 1        
        c11                =  1.0d0 - Rdiff(km1,1,1) * Rd1k(km1,1,1) - Rdiff(km1,1,2) * Rd1k(km1,2,1)
        c12                =      - Rdiff(km1,1,1) * Rd1k(km1,1,2) - Rdiff(km1,1,2) * Rd1k(km1,2,2)
        c21                =      - Rdiff(km1,2,1) * Rd1k(km1,1,1) - Rdiff(km1,2,2) * Rd1k(km1,2,1)
        c22                =  1.0d0 - Rdiff(km1,2,1) * Rd1k(km1,1,2) - Rdiff(km1,2,2) * Rd1k(km1,2,2)
        pmod             =  c11 * c22 - c12 * c21
        t11              =  c22 / pmod
        t12              = -c12 / pmod
        t21              = -c21 / pmod
        t22              =  c11 / pmod

        b11              =  t11 * Tdiff(km1,1,1) + t12 * Tdiff(km1,2,1)
        b12              =  t11 * Tdiff(km1,1,2) + t12 * Tdiff(km1,2,2)
        b21              =  t21 * Tdiff(km1,1,1) + t22 * Tdiff(km1,2,1)
        b22              =  t21 * Tdiff(km1,1,2) + t22 * Tdiff(km1,2,2)
!
        scumdtr(k)   =  scumdtr(km1)*dtr(km1)
        d1               =  Rdir(km1,1) * scumdtr(km1) + Rdiff(km1,1,1) * T1k(km1,1) + &
                            Rdiff(km1,1,2) * T1k(km1,2)
        d2               =  Rdir(km1,2) * scumdtr(km1) + Rdiff(km1,2,1) * T1k(km1,1) + &
                            Rdiff(km1,2,2) * T1k(km1,2)
        uu1              =  d1 * t11 + d2 * t12
        uu2              =  d1 * t21 + d2 * t22
        du1              =  T1k(km1,1) + uu1 * Rd1k(km1,1,1) + uu2 * Rd1k(km1,1,2)
        du2              =  T1k(km1,2) + uu1 * Rd1k(km1,2,1) + uu2 * Rd1k(km1,2,2)
        T1k(k,1)   =  Tdir(km1,1) * scumdtr(km1) + du1 * Tdiff(km1,1,1) + du2 * Tdiff(km1,1,2)
        T1k(k,2)   =  Tdir(km1,2) * scumdtr(km1) + du1 * Tdiff(km1,2,1) + du2 * Tdiff(km1,2,2)
!
        a11              =  Tdiff(km1,1,1) * Rd1k(km1,1,1) + Tdiff(km1,1,2) * Rd1k(km1,2,1)
        a12              =  Tdiff(km1,1,1) * Rd1k(km1,1,2) + Tdiff(km1,1,2) * Rd1k(km1,2,2)
        a21              =  Tdiff(km1,2,1) * Rd1k(km1,1,1) + Tdiff(km1,2,2) * Rd1k(km1,2,1)
        a22              =  Tdiff(km1,2,1) * Rd1k(km1,1,2) + Tdiff(km1,2,2) * Rd1k(km1,2,2)
        Rd1k(k,1,1) =  Rdiff(km1,1,1) + a11 * b11 + a12 * b21
        Rd1k(k,1,2) =  Rdiff(km1,1,2) + a11 * b12 + a12 * b22
        Rd1k(k,2,1) =  Rdiff(km1,2,1) + a21 * b11 + a22 * b21
        Rd1k(k,2,2) =  Rdiff(km1,2,2) + a21 * b12 + a22 * b22
    end do
  

    do l = nlay, 1, -1
      lp1 = l + 1
        c11              =  1.0d0 - Rst1k(lp1,1,1) * Rdiff(l,1,1) - Rst1k(lp1,1,2) * Rdiff(l,2,1)
        c12              =      - Rst1k(lp1,1,1) * Rdiff(l,1,2) - Rst1k(lp1,1,2) * Rdiff(l,2,2)
        c21              =      - Rst1k(lp1,2,1) * Rdiff(l,1,1) - Rst1k(lp1,2,2) * Rdiff(l,2,1)
        c22              =  1.0d0 - Rst1k(lp1,2,1) * Rdiff(l,1,2) - Rst1k(lp1,2,2) * Rdiff(l,2,2)
        pmod             =  c11 * c22 - c12 * c21
        t11              =  c22 / pmod
        t12              = -c12 / pmod
        t21              = -c21 / pmod
        t22              =  c11 / pmod
        d1               =  R1k(lp1,1) * dtr(l) + Rst1k(lp1,1,1) * Tdir(l,1) + &
                            Rst1k(lp1,1,2) * Tdir(l,2)
        d2               =  R1k(lp1,2) * dtr(l) + Rst1k(lp1,2,1) * Tdir(l,1) + &
                            Rst1k(lp1,2,2) * Tdir(l,2)
        uu1              =  d1 * t11 + d2 * t12
        uu2              =  d1 * t21 + d2 * t22
        R1k(l,1)   =  Rdir(l,1) + uu1 * Tdiff(l,1,1) + uu2 * Tdiff(l,1,2)
        R1k(l,2)   =  Rdir(l,2) + uu1 * Tdiff(l,2,1) + uu2 * Tdiff(l,2,2)
!
        c11              =  1.0d0 - Rdiff(l,1,1) * Rst1k(lp1,1,1) - Rdiff(l,1,2) * Rst1k(lp1,2,1)
        c12              =      - Rdiff(l,1,1) * Rst1k(lp1,1,2) - Rdiff(l,1,2) * Rst1k(lp1,2,2)
        c21              =      - Rdiff(l,2,1) * Rst1k(lp1,1,1) - Rdiff(l,2,2) * Rst1k(lp1,2,1)
        c22              =  1.0d0 - Rdiff(l,2,1) * Rst1k(lp1,1,2) - Rdiff(l,2,2) * Rst1k(lp1,2,2)
        pmod             =  c11 * c22 - c12 * c21
        t11              =  c22 / pmod
        t12              = -c12 / pmod
        t21              = -c21 / pmod
        t22              =  c11 / pmod
!
        r11              =  Tdiff(l,1,1) * (t11 * Rst1k(lp1,1,1) +  t21 * Rst1k(lp1,1,2)) + &
                           & Tdiff(l,1,2) * (t11 * Rst1k(lp1,2,1) +  t21 * Rst1k(lp1,2,2))
        r12              =  Tdiff(l,1,1) * (t12 * Rst1k(lp1,1,1) +  t22 * Rst1k(lp1,1,2)) + &
                          & Tdiff(l,1,2) * (t12 * Rst1k(lp1,2,1) +  t22 * Rst1k(lp1,2,2))
        r21              =  Tdiff(l,2,1) * (t11 * Rst1k(lp1,1,1) +  t21 * Rst1k(lp1,1,2)) + &
                          &  Tdiff(l,2,2) * (t11 * Rst1k(lp1,2,1) +  t21 * Rst1k(lp1,2,2))
        r22              =  Tdiff(l,2,1) * (t12 * Rst1k(lp1,1,1) +  t22 * Rst1k(lp1,1,2)) + &
                          &  Tdiff(l,2,2) * (t12 * Rst1k(lp1,2,1) +  t22 * Rst1k(lp1,2,2))
!
        Rst1k(l,1,1) =  Rdiff(l,1,1) + r11 * Tdiff(l,1,1) + r12 * Tdiff(l,2,1)
        Rst1k(l,1,2) =  Rdiff(l,1,2) + r11 * Tdiff(l,1,2) + r12 * Tdiff(l,2,2)
        Rst1k(l,2,1) =  Rdiff(l,2,1) + r21 * Tdiff(l,1,1) + r22 * Tdiff(l,2,1)
        Rst1k(l,2,2) =  Rdiff(l,2,2) + r21 * Tdiff(l,1,2) + r22 * Tdiff(l,2,2)
    end do

    do k = 1, nlev
        c11              =  1.0d0 - Rst1k(k,1,1) * Rd1k(k,1,1) - Rst1k(k,1,2) * Rd1k(k,2,1)
        c12              =      - Rst1k(k,1,1) * Rd1k(k,1,2) - Rst1k(k,1,2) * Rd1k(k,2,2)
        c21              =      - Rst1k(k,2,1) * Rd1k(k,1,1) - Rst1k(k,2,2) * Rd1k(k,2,1)
        c22              =  1.0d0 - Rst1k(k,2,1) * Rd1k(k,1,2) - Rst1k(k,2,2) * Rd1k(k,2,2)
        pmod             =  c11 * c22 - c12 * c21
        t11              =  c22 / pmod
        t12              = -c12 / pmod
        t21              = -c21 / pmod
        t22              =  c11 / pmod
!
        d1               =  R1k(k,1) * scumdtr(k) + Rst1k(k,1,1) * T1k(k,1) + &
                           & Rst1k(k,1,2) * T1k(k,2)
        d2               =  R1k(k,2) * scumdtr(k) + Rst1k(k,2,1) * T1k(k,1) + &
                           & Rst1k(k,2,2) * T1k(k,2)
        uu1              =  d1 * t11 + d2 * t12
        du1              =  T1k(k,1) + Rd1k(k,1,1) * uu1 + Rd1k(k,1,2) * (d1 * t21 + d2 * t22)
        reflt(k)   =  uu1
        trant(k)   =  du1 + scumdtr(k)

    end do

    !! Down and up fluxes are multiplied by the incident flux
    !! up is defined as negative in the adding method, so we make it positive here
    sw_down(:) = trant(:)*mu_z*Finc
    sw_up(:) = -reflt(:)*mu_z*Finc

  end subroutine sw_SDA

  ! Perform linear interpolation in log10 space
  subroutine linear_log_interp(xval, x1, x2, y1, y2, yval)
    implicit none

    real(dp), intent(in) :: xval, y1, y2, x1, x2
    real(dp) :: ly1, ly2
    real(dp), intent(out) :: yval
    real(dp) :: norm

    ly1 = log10(y1); ly2 = log10(y2)

    norm = 1.0_dp / log10(x2/x1)

    yval = 10.0_dp**((ly1 * log10(x2/xval) + ly2 * log10(xval/x1)) * norm)

  end subroutine linear_log_interp

  subroutine bezier_interp(xi, yi, ni, x, y)
    implicit none

    integer, intent(in) :: ni
    real(dp), dimension(ni), intent(in) :: xi, yi
    real(dp), intent(in) :: x
    real(dp), intent(out) :: y

    real(dp) :: dx, dx1, dy, dy1, wh, yc, t, wlim, wlim1

    !xc = (xi(1) + xi(2))/2.0_dp ! Control point (no needed here, implicitly included)
    dx = xi(2) - xi(1)
    dx1 = xi(3) - xi(2)
    dy = yi(2) - yi(1)
    dy1 = yi(3) - yi(2)

    if (x > xi(1) .and. x < xi(2)) then
      ! left hand side interpolation
      !print*,'left'
      wh = dx1/(dx + dx1)
      wlim = 1.0_dp + 1.0_dp/(1.0_dp - (dy1/dy) * (dx/dx1))
      wlim1 = 1.0_dp/(1.0_dp - (dy/dy1) * (dx1/dx))
      if (wh <= min(wlim,wlim1) .or. wh >= max(wlim,wlim1)) then
        wh = 1.0_dp
      end if
      yc = yi(2) - dx/2.0_dp * (wh*dy/dx + (1.0_dp - wh)*dy1/dx1)
      t = (x - xi(1))/dx
      y = (1.0_dp - t)**2 * yi(1) + 2.0_dp*t*(1.0_dp - t)*yc + t**2*yi(2)
    else ! (x > xi(2) and x < xi(3)) then
      ! right hand side interpolation
      !print*,'right'
      wh = dx/(dx + dx1)
      wlim = 1.0_dp/(1.0_dp - (dy1/dy) * (dx/dx1))
      wlim1 = 1.0_dp + 1.0_dp/(1.0_dp - (dy/dy1) * (dx1/dx))
      if (wh <= min(wlim,wlim1) .or. wh >= max(wlim,wlim1)) then
        wh = 1.0_dp
      end if
      yc = yi(2) + dx1/2.0_dp * (wh*dy1/dx1 + (1.0_dp - wh)*dy/dx)
      t = (x - xi(2))/(dx1)
      y = (1.0_dp - t)**2 * yi(2) + 2.0_dp*t*(1.0_dp - t)*yc + t**2*yi(3)
    end if

  end subroutine bezier_interp

  pure function matinv4(A) result(B)
    !! Performs a direct calculation of the inverse of a 4×4 matrix.
    real(dp), intent(in) :: A(4,4)   !! Matrix
    real(dp)             :: B(4,4)   !! Inverse matrix
    real(dp)             :: detinv, s0, s1, s2, s3, s4, s5, c5, c4, c3, c2, c1, c0

    s0 = A(1,1) * A(2,2) - A(2,1) * A(1,2)
    s1 = A(1,1) * A(2,3) - A(2,1) * A(1,3)
    s2 = A(1,1) * A(2,4) - A(2,1) * A(1,4)
    s3 = A(1,2) * A(2,3) - A(2,2) * A(1,3)
    s4 = A(1,2) * A(2,4) - A(2,2) * A(1,4)
    s5 = A(1,3) * A(2,4) - A(2,3) * A(1,4)

    c5 = A(3,3) * A(4,4) - A(4,3) * A(3,4)
    c4 = A(3,2) * A(4,4) - A(4,2) * A(3,4)
    c3 = A(3,2) * A(4,3) - A(4,2) * A(3,3)
    c2 = A(3,1) * A(4,4) - A(4,1) * A(3,4)
    c1 = A(3,1) * A(4,3) - A(4,1) * A(3,3)
    c0 = A(3,1) * A(4,2) - A(4,1) * A(3,2)

    detinv = 1.0_dp / (s0 * c5 - s1 * c4 + s2 * c3 + s3 * c2 - s4 * c1 + s5 * c0)

    B(1,1) = ( A(2,2) * c5 - A(2,3) * c4 + A(2,4) * c3) * detinv
    B(1,2) = (-A(1,2) * c5 + A(1,3) * c4 - A(1,4) * c3) * detinv
    B(1,3) = ( A(4,2) * s5 - A(4,3) * s4 + A(4,4) * s3) * detinv
    B(1,4) = (-A(3,2) * s5 + A(3,3) * s4 - A(3,4) * s3) * detinv

    B(2,1) = (-A(2,1) * c5 + A(2,3) * c2 - A(2,4) * c1) * detinv
    B(2,2) = ( A(1,1) * c5 - A(1,3) * c2 + A(1,4) * c1) * detinv
    B(2,3) = (-A(4,1) * s5 + A(4,3) * s2 - A(4,4) * s1) * detinv
    B(2,4) = ( A(3,1) * s5 - A(3,3) * s2 + A(3,4) * s1) * detinv

    B(3,1) = ( A(2,1) * c4 - A(2,2) * c2 + A(2,4) * c0) * detinv
    B(3,2) = (-A(1,1) * c4 + A(1,2) * c2 - A(1,4) * c0) * detinv
    B(3,3) = ( A(4,1) * s4 - A(4,2) * s2 + A(4,4) * s0) * detinv
    B(3,4) = (-A(3,1) * s4 + A(3,2) * s2 - A(3,4) * s0) * detinv

    B(4,1) = (-A(2,1) * c3 + A(2,2) * c1 - A(2,3) * c0) * detinv
    B(4,2) = ( A(1,1) * c3 - A(1,2) * c1 + A(1,3) * c0) * detinv
    B(4,3) = (-A(4,1) * s3 + A(4,2) * s1 - A(4,3) * s0) * detinv
    B(4,4) = ( A(3,1) * s3 - A(3,2) * s1 + A(3,3) * s0) * detinv

  end function matinv4

  subroutine ludcmp(A,n,np,indx,D)
    implicit none

    integer, intent(in) :: n, np
    real(dp), dimension(np,np), intent(inout) :: A

    integer, dimension(n), intent(out) :: indx
    real(dp), intent(out) :: D

    integer, parameter :: nmax = 100
    real(dp), parameter :: tiny = 1.0e-20_dp
    real(dp), dimension(nmax) :: vv

    integer :: i, j, k, imax
    real(dp) :: aamax, dum, sum

    D = 1.0_dp

    do i = 1, n
      aamax = 0.0_dp
      do j = 1, n
        if (abs(A(i,j)) > aamax) then
          aamax = abs(A(i,j))
        end if
      end do
      if (aamax == 0.0_dp) then
        print*, 'singualr matrix in LU decomp!'
        stop
      end if
      vv(i) = 1.0_dp/aamax
    end do

  
    do j = 1, n
      do i = 1, j-1
        sum = A(i,j)
        do k = 1, i-1
          sum = sum  - A(i,k)*A(k,j)
        end do
        A(i,j) = sum
      end do
      aamax = 0.0_dp
      do i = j, n
        sum = A(i,j)
        do k = 1, j-1
          sum = sum  - A(i,k)*A(k,j)
        end do
        A(i,j) = sum
        dum = vv(i)*abs(sum)
        if (dum >= aamax) then
          imax = i
          aamax = dum
        end if
      end do
      if (j /= imax) then
        do k = 1, n
          dum = A(imax,k)
          A(imax,k) = A(j,k)
          A(j,k) = dum
        end do 
        D = -D
        vv(imax) = vv(j)
      end if
      indx(j) = imax
      if (A(j,j) == 0.0_dp) then
        A(j,j) = tiny
      end if
      if (j /= n) then
        dum = 1.0_dp/A(j,j)
        do i = j+1, n
          A(i,j) = A(i,j)*dum
        end do
      end if
    end do

  end subroutine ludcmp

  subroutine lubksb(A, n, np, indx, B)
    implicit none

    integer, intent(in) :: n, np
    integer, dimension(n), intent(in) :: indx
    real(dp), dimension(np,np), intent(in) :: A

    real(dp), dimension(n), intent(out) :: B

    integer :: i, j, ii, ll
    real(dp) :: sum

    ii = 0

    do i = 1, n
      ll = indx(i)
      sum = B(ll)
      b(ll) = b(i)
      if (ii /= 0) then
        do j = ii,i-1
          sum = sum - A(i,j)*B(j)
        end do
      else if (sum /= 0.0_dp) then
        ii = i
      end if
      B(i) = sum
    end do

    do i = n, 1, -1
      sum = B(i)
      if (i < n) then
        do j = i+1, n
          sum = sum - A(i,j)*B(j)
        end do
      end if
      B(i) = sum/A(i,i)
    end do
    
  end subroutine lubksb

  subroutine inv_LU(A,n,np,Y)
    implicit none

    integer, intent(in) :: n, np
    real(dp), dimension(np,np), intent(inout) :: A

    integer, dimension(n) :: indx
    real(dp), dimension(np,np), intent(out) :: Y

    real(dp) :: D
    integer :: i, j

    do i = 1, n
      do j = 1, n
        Y(i,j) = 0.0_dp
      end do
      Y(i,i) = 1.0_dp
    end do

    call ludcmp(A,n,np,indx,D)

    do j = 1, n
      call lubksb(A,n,np,indx,Y(1,j))
    end do

  end subroutine inv_LU

end module ts_AA_L_mod
