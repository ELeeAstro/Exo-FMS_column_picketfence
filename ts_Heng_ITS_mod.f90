
!!! Work in progress - do not use!

module ts_Heng_ITS_mod
  use, intrinsic :: iso_fortran_env
  implicit none

  !! Precision variables
  integer, parameter :: dp = REAL64

  !! Required constants
  real(dp), parameter :: pi = 4.0_dp * atan(1.0_dp)
  real(dp), parameter :: twopi = 2.0_dp * pi
  real(dp), parameter :: sb = 5.670374419e-8_dp

  real(dp), parameter :: e1 = 1.0_dp/2.0_dp ! First Eddington coefficent
  real(dp), parameter :: e2 = 2.0_dp/3.0_dp ! Second Eddington coefficent

  integer, parameter :: nit = 3

  public :: ts_Heng_ITS
  private :: lw_Heng_ITS, sw_Heng_ITS, linear_log_interp, bezier_interp

contains
 subroutine ts_Heng_ITS(surf, Bezier, nlay, nlev, Ts, Tl, pl, pe, tau_Ve, tau_IRe, mu_z, F0, Tint, AB, &
    & sw_a, sw_g, lw_a, lw_g, sw_a_surf, lw_a_surf, net_F, olr, asr, net_Fs)
    implicit none

    !! Input variables
    logical, intent(in) :: surf, Bezier
    integer, intent(in) :: nlay, nlev
    real(dp), intent(in) :: F0, Tint, AB, sw_a_surf, lw_a_surf, Ts
    real(dp), dimension(nlay), intent(in) :: Tl, pl
    real(dp), dimension(nlev), intent(in) :: pe, mu_z
    real(dp), dimension(nlev), intent(in) :: tau_Ve, tau_IRe
    real(dp), dimension(nlay), intent(in) :: sw_a, sw_g, lw_a, lw_g

    !! Output variables
    real(dp), intent(out) :: olr, asr, net_Fs
    real(dp), dimension(nlev), intent(out) :: net_F

    !! Work variables
    integer :: i
    real(dp) :: Finc, be_int
    real(dp), dimension(nlev) :: Te, be
    real(dp), dimension(nlev) :: sw_down, sw_up, lw_down, lw_up
    real(dp), dimension(nlev) :: lw_net, sw_net

    !! Find temperature at layer edges through interpolation and extrapolation
    if (Bezier .eqv. .True.) then
      ! Perform interpolation using Bezier peicewise polynomial interpolation
      do i = 2, nlay-1
        call bezier_interp(pl(i-1:i+1), Tl(i-1:i+1), 3, pe(i), Te(i))
        !print*, i, pl(i), pl(i-1), pe(i), Tl(i-1), Tl(i), Te(i)
      end do
      call bezier_interp(pl(nlay-2:nlay), Tl(nlay-2:nlay), 3, pe(nlay), Te(nlay))
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
      call sw_Heng_ITS(nlay, nlev, Finc, tau_Ve(:), mu_z(nlev), sw_a, sw_g, sw_a_surf, sw_down(:), sw_up(:))
    else
      sw_down(:) = 0.0_dp
      sw_up(:) = 0.0_dp
    end if

    !! Longwave two-stream flux calculation
    be(:) = (sb * Te(:)**4)/pi  ! Integrated planck function intensity at levels
    if (surf .eqv. .True.) then
      be_int = (sb * Ts**4)/pi
    else
      be_int = (sb * Tint**4)/pi ! Integrated planck function intensity for internal temperature
    end if

    call lw_Heng_ITS(surf, nlay, nlev, be, be_int, tau_IRe(:), lw_a, lw_g, lw_a_surf, lw_up(:), lw_down(:))

    !! Net fluxes at each level
    lw_net(:) = lw_up(:) - lw_down(:)
    sw_net(:) = sw_up(:) - sw_down(:)
    net_F(:) = lw_net(:) + sw_net(:)

    !! Net surface flux (for surface temperature evolution)
    !! We have to define positive as downward (heating) and cooling (upward) in this case
    net_Fs = sw_down(nlev) + lw_down(nlev) - lw_up(nlev)

    !! Output the olr
    olr = lw_up(1)

    !! Output asr
    asr = sw_down(1) - sw_up(1)

  end subroutine ts_Heng_ITS

  subroutine lw_Heng_ITS(surf, nlay, nlev, be, be_int, tau_IRe, w0, g0, lw_a_surf, lw_up, lw_down)
    implicit none

    !! Input variables
    logical, intent(in) :: surf 
    integer, intent(in) :: nlay, nlev
    real(dp), dimension(nlev), intent(in) :: be, tau_IRe
    real(dp), dimension(nlay), intent(in) :: w0, g0
    real(dp), intent(in) :: be_int, lw_a_surf

    !! Output variables
    real(dp), dimension(nlev), intent(out) :: lw_up, lw_down

    integer :: k, n
    real(dp) :: D, BB_term
    real(dp), dimension(nlay) :: chi, eta, phi, zetap, zetam, pie
    real(dp), dimension(nlay) :: dtau, T, Blin, E
    real(dp), external :: eone
    real(dp), dimension(nlev) :: Fup, Fdown

    !! Calculate dtau in each layer
    dtau(:) = tau_IRe(2:) - tau_IRe(1:nlay)

    do k = 1, nlay

      D = -log((1.0_dp - dtau(k))*exp(-dtau(k)) + dtau(k)**2*eone(dtau(k)))/dtau(k)
      T(k) = exp(-D*sqrt((1.0_dp - w0(k))*(1.0_dp - w0(k)*g0(k)))*dtau(k))

      !! Linear B with tau function
      if (dtau(k) < 1.0e-6_dp) then
        ! Too low optical depth for numerical stability, Blin = 0
        Blin(k) = 0.0_dp
      else
        ! Linear with tau Planck function 
        Blin(k) = (be(k+1)-be(k))/dtau(k)
      end if

    end do

    !! First do a non-scattering sweep using the absorption only part of the ITS

    !! Do downward sweep
    lw_down(1) = 0.0_dp
    do k = 1, nlay
      lw_down(k+1) = lw_down(k)*T(k) + (1.0_dp - w0(k))*pi*be(k)*(1.0_dp - T(k)) - &
        & (1.0_dp - w0(k))*pi*Blin(k)* &
        & (2.0_dp/3.0_dp * (1.0_dp - exp(-dtau(k))) - dtau(k)*(1.0_dp - T(k)/3.0_dp))
    end do

    !! Perform upward loop
    if (surf .eqv. .True.) then
      ! Surface boundary condition given by surface temperature + reflected longwave radiaiton
      lw_up(nlev) = lw_down(nlev)*lw_a_surf + (1.0_dp - lw_a_surf)*pi*be_int
    else
      ! Lower boundary condition - internal heat definition Fint = F_down - F_up
      ! here the lw_a_surf is assumed to be = 1 as per the definition
      ! here we use the same condition but use intensity units to be consistent
      lw_up(nlev) = lw_down(nlev) + pi*be_int
    end if

    do k = nlay, 1, -1
      lw_up(k) = lw_up(k+1)*T(k) + (1.0_dp - w0(k))*pi*be(k+1)*(1.0_dp - T(k)) + &
        & (1.0_dp - w0(k))*pi*Blin(k)* &
        & (2.0_dp/3.0_dp * (1.0_dp - exp(-dtau(k))) - dtau(k)*(1.0_dp - T(k)/3.0_dp))
    end do

    !! If no scattering component in profile, then just return two-stream flux
    !! no need to perform any scattering calculations
    if (all(w0(:) <= 1.0e-6_dp)) then
      return
    end if

    !! E values
    E(:) = 1.225_dp - 0.1582_dp*g0(:) - 0.1777_dp*w0(:) - &
      & 0.07465_dp*g0(:)**2 + 0.2351_dp*w0(:)*g0(:) - 0.05582_dp*w0(:)**2

    !! Recalculate transmission function with 1st Eddington coefficent
    T(:) = exp(-1.0_dp/e1*sqrt(E(:)*(E(:) - w0(:))*(1.0_dp - w0(:)*g0(:)))*dtau(:))

    !! Coupling coefficents
    zetap(:) = 0.5_dp * (1.0_dp + sqrt((E(:) - w0(:))/(E(:)*(1.0_dp - w0(:)*g0(:)))))
    zetam(:) = 0.5_dp * (1.0_dp - sqrt((E(:) - w0(:))/(E(:)*(1.0_dp - w0(:)*g0(:)))))

    chi(:) = zetam(:)**2*T(:)**2 - zetap(:)**2
    eta(:) = zetap(:)*zetam(:)*(1.0_dp - T(:)**2)
    phi(:) = (zetam(:)**2 - zetap(:)**2)*T(:)
    pie(:) = pi * ((1.0_dp - w0(:))/(E(:) - w0(:)))

    do n = 1, nit

      !! Now we have first order boundary conditions, we can calculate the two-stream scattering coefficents
      do k = 1, nlay

        BB_term = pie(k) * (be(k)*(chi(k) + eta(k)) - phi(k)*be(k+1) + &
          & Blin(k)/(2.0_dp*E(k)*(1.0_dp - w0(k)*g0(k))) * (chi(k) - phi(k) - eta(k)))

        Fup(k) = 1.0_dp/chi(k) * (phi(k)*lw_up(k+1) - eta(k)*lw_down(k) + BB_term)

        BB_term = pie(k) * (be(k+1)*(chi(k) + eta(k)) - phi(k)*be(k) + &
          & Blin(k)/(2.0_dp*E(k)*(1.0_dp - w0(k)*g0(k))) * (eta(k) + phi(k) - chi(k)))

        Fdown(k) = 1.0_dp/chi(k) * (phi(k)*lw_down(k) - eta(k)*lw_up(k+1) + BB_term)   
        
      end do

      Fup(nlev) = 0.0_dp
      Fdown(nlev) = 0.0_dp

      !! Then we can do the two-stream fluxes with the scattering coefficents

      !! Do downward sweep
      lw_down(1) = 0.0_dp
      do k = 1, nlay
        lw_down(k+1) = lw_down(k)*T(k) + (1.0_dp - w0(k))*pi*be(k)*(1.0_dp - T(k)) - &
          & (1.0_dp - w0(k))*pi*Blin(k)* &
          & (2.0_dp/3.0_dp * (1.0_dp - exp(-dtau(k))) - dtau(k)*(1.0_dp - T(k)/3.0_dp)) + &
          & w0(k)/2.0_dp * ((1.0_dp - g0(k))*Fup(k) + (1.0_dp + g0(k))*Fdown(k+1)) * (1.0_dp - T(k)) 
      end do

      !! Perform upward loop
      if (surf .eqv. .True.) then
        ! Surface boundary condition given by surface temperature + reflected longwave radiaiton
        lw_up(nlev) = lw_down(nlev)*lw_a_surf + (1.0_dp - lw_a_surf)*pi*be_int
      else
        ! Lower boundary condition - internal heat definition Fint = F_down - F_up
        ! here the lw_a_surf is assumed to be = 1 as per the definition
        ! here we use the same condition but use intensity units to be consistent
        lw_up(nlev) = lw_down(nlev) + pi*be_int
      end if

      do k = nlay, 1, -1
        lw_up(k) = lw_up(k+1)*T(k) + (1.0_dp - w0(k))*pi*be(k+1)*(1.0_dp - T(k)) + &
          & (1.0_dp - w0(k))*pi*Blin(k)* &
          & (2.0_dp/3.0_dp * (1.0_dp - exp(-dtau(k))) - dtau(k)*(1.0_dp - T(k)/3.0_dp)) + &
          & w0(k)/2.0_dp * ((1.0_dp + g0(k))*Fup(k) + (1.0_dp - g0(k))*Fdown(k+1)) * (1.0_dp - T(k)) 
      end do


      !! This process of iteration can be repeated nit times until convergence, or just left as first order
      !! Best way would be to perform a matrix inversion across all layers once which would allow multiple scattering

    end do

  end subroutine lw_Heng_ITS

  subroutine sw_Heng_ITS(nlay, nlev, Finc, tau_Ve, mu_z, w0, g0, w_surf, sw_down, sw_up)
    implicit none

    !! Input variables
    integer, intent(in) :: nlay, nlev
    real(dp), dimension(nlev), intent(in) :: tau_Ve
    real(dp), dimension(nlay), intent(in) :: w0, g0
    real(dp), intent(in) :: w_surf, mu_z, Finc

    !! Output variables
    real(dp), dimension(nlev), intent(out) :: sw_up, sw_down

    integer :: k, ki
    real(dp), dimension(nlay) :: chi, eta, phi, zetap, zetam
    real(dp), dimension(nlay) :: dtau, T,  E, Gp, Gm
    real(dp), dimension(nlev) :: Fbeam
    real(dp), dimension(nlay,2,2) :: A, B, C
    real(dp), dimension(nlay,2) :: X, D

    !! If zero albedo across all atmospheric layers then return direct beam only
    if (all(w0(:) <= 1.0e-6_dp)) then
      sw_down(:) = Finc * mu_z * exp(-tau_Ve(:)/mu_z)
      sw_down(nlev) = sw_down(nlev) * (1.0_dp - w_surf) ! The surface flux for surface heating is the amount of flux absorbed by surface
      sw_up(:) = 0.0_dp ! We assume no upward flux here even if surface albedo
      return
    end if

    !! Direct beam flux
    Fbeam(:) = Finc * mu_z * exp(-tau_Ve(:)/mu_z) 

    !! Calculate dtau in each layer
    dtau(:) = tau_Ve(2:) - tau_Ve(1:nlay)

    !! E values
    E(:) = 1.225_dp - 0.1582_dp*g0(:) - 0.1777_dp*w0(:) - &
      & 0.07465_dp*g0(:)**2 + 0.2351_dp*w0(:)*g0(:) - 0.05582_dp*w0(:)**2

    !! Transmission function with 1st Eddington coefficent
    T(:) = exp(-1.0_dp/e1*sqrt(E(:)*(E(:) - w0(:))*(1.0_dp - w0(:)*g0(:)))*dtau(:))

    !! Coupling coefficents
    zetap(:) = 0.5_dp * (1.0_dp + sqrt((E(:) - w0(:))/(E(:)*(1.0_dp - w0(:)*g0(:)))))
    zetam(:) = 0.5_dp * (1.0_dp - sqrt((E(:) - w0(:))/(E(:)*(1.0_dp - w0(:)*g0(:)))))

    chi(:) = zetam(:)**2*T(:)**2 - zetap(:)**2
    eta(:) = zetap(:)*zetam(:)*(1.0_dp - T(:)**2)
    phi(:) = (zetam(:)**2 - zetap(:)**2)*T(:)

    !! Find downward and upward scattering coefficents
    Gp(:) = 0.5_dp * (((w0(:) * (2.0_dp*E(:)*(1.0_dp - w0(:)*g0(:)) + g0(:)/e2))/ &
      & (4.0_dp*E(:)*mu_z**2*(E(:) - w0(:))*(1.0_dp - w0(:)*g0(:)) - 1.0_dp)) * &
      & (-mu_z + 1.0_dp/(2.0_dp*E(:)*(1.0_dp - w0(:)*g0(:)))) + &
      & ((w0(:)*g0(:))/(2.0_dp*e2*E(:)*(1.0_dp - w0(:)*g0(:)))))
    Gm(:) = 0.5_dp * (((w0(:) * (2.0_dp*E(:)*(1.0_dp - w0(:)*g0(:)) + g0(:)/e2))/ &
      & (4.0_dp*E(:)*mu_z**2*(E(:) - w0(:))*(1.0_dp - w0(:)*g0(:)) - 1.0_dp)) * &
      & (-mu_z - 1.0_dp/(2.0_dp*E(:)*(1.0_dp - w0(:)*g0(:)))) - &
      & ((w0(:)*g0(:))/(2.0_dp*e2*E(:)*(1.0_dp - w0(:)*g0(:)))))

    !! Find upward and downward scattered fluxes
    do k = 1, nlay
      sw_up(k) = 1.0_dp/chi(k) * (phi(k)*Gp(k)*Fbeam(k+1) - (eta(k)*Gm(k) + chi(k)*Gp(k))*Fbeam(k))
      sw_down(k) = Fbeam(k) + 1.0_dp/chi(k) * (phi(k)*Gm(k)*Fbeam(k) - (eta(k)*Gp(k) + chi(k)*Gm(k))*Fbeam(k+1))
    end do

    sw_down(nlev) = Fbeam(nlev) 
    sw_up(nlev) = w_surf*sw_down(nlev)


  end subroutine sw_Heng_ITS

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

    if ((x > xi(1)) .and. (x < xi(2))) then
      ! left hand side interpolation
      !print*,'left'
      wh = dx1/(dx + dx1)
      wlim = 1.0_dp + 1.0_dp/(1.0_dp - (dy1/dy) * (dx/dx1))
      wlim1 = 1.0_dp/(1.0_dp - (dy/dy1) * (dx1/dx))
      if ((wh <= min(wlim,wlim1)) .or. (wh >= max(wlim,wlim1))) then
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
      if ((wh <= min(wlim,wlim1)) .or. (wh >= max(wlim,wlim1))) then
        wh = 1.0_dp
      end if
      yc = yi(2) + dx1/2.0_dp * (wh*dy1/dx1 + (1.0_dp - wh)*dy/dx)
      t = (x - xi(2))/(dx1)
      y = (1.0_dp - t)**2 * yi(2) + 2.0_dp*t*(1.0_dp - t)*yc + t**2*yi(3)
    end if

  end subroutine bezier_interp

end module ts_Heng_ITS_mod
