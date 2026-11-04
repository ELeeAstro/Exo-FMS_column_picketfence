!!!
! Elspeth KH Lee - May 2021
! A simple program that emulates a column inside the Exo-FMS GCM
! This is useful for testing different RT solutions
! This version uses three visible and two infrared picket-fence bands
!
! Input parameters are set via the namelist (FMS_RC.nml) - hopefully paramaters are
! somewhat self explanatory. See the github readme for more information.
! Not guarenteed to be bug free
! NOTE: Indexing starts at 1, where 1 is the Top Of Atmosphere level or layer
!!!

program Exo_FMS_RC
  use, intrinsic :: iso_fortran_env
  use picket_radiation_mod, only: solve_sw_band, solve_lw_band
  use Manabe_dry_adjust_mod, only: Manabe_dry_adjust
  use Pierrehumbert_dry_adjust_mod, only: Pierrehumbert_dry_adjust
  use MLT_mod, only: MLT
  use k_Rosseland_mod, only : k_Ross_Freedman, k_Ross_Valencia, gam_Parmentier, Bond_Parmentier
  use IC_mod, only : IC_profile
  use ieee_arithmetic
  implicit none

  ! Precision variable
  integer, parameter :: dp = REAL64

  ! Constants
  real(dp), parameter :: sb = 5.670374419e-8_dp

  integer :: n, i, k, u, j, ib
  integer :: nstep, nlay, nlev
  integer :: table_num
  real(dp) :: t_step, t_tot
  real(dp) :: mu_z, Tirr, Tint, F0, Fint, pref, pu, met
  real(dp) :: Teff, gam_1, gam_2, tau_lim, gam_P
  real(dp), allocatable, dimension(:) :: Tl, pl, pe, dpe
  real(dp), allocatable, dimension(:) :: Beta_V, Beta_IR, gam_V
  real(dp), allocatable, dimension(:,:) :: k_Vl, k_IRl
  real(dp), allocatable, dimension(:,:) :: tau_Ve, tau_IRe
  real(dp), allocatable, dimension(:) :: dT_rad, dT_conv, net_F

  real(dp), allocatable, dimension(:,:) :: sw_a, sw_g, lw_a, lw_g

  real(dp) :: cp_air, grav, k_IR, k_V, kappa_air, Rd_air
  real(dp) :: olr, asr, band_olr, band_asr, Finc
  real(dp), allocatable :: sw_up(:), sw_down(:), sw_net(:), lw_up(:), lw_down(:), lw_net(:)
  real(dp), allocatable :: band_up(:), band_down(:), band_net(:)
  real(dp), allocatable :: sw_up_band(:,:), sw_down_band(:,:), sw_net_band(:,:)
  real(dp), allocatable :: lw_up_band(:,:), lw_down_band(:,:), lw_net_band(:,:)
  real(dp), allocatable :: Tl_rad(:), Kzz(:), Kzz_e(:), Rd_bar(:), cp_bar(:)
  logical :: use_bond_fit = .true.

  real(dp), dimension(3) :: sw_ac, sw_gc
  real(dp), dimension(2) :: lw_ac, lw_gc

  integer :: iIC
  logical :: corr
  real(dp) :: prc
  real(dp) :: Tdeep_init = 2000.0_dp, Tstrat_init = 500.0_dp

  real(dp) :: fl, AB

  integer :: ua, ub
  character(len=1024) :: a_sh, b_sh
  real(dp), allocatable, dimension(:) :: a_hs, b_hs

  real(dp) :: start, finish

  character(len=50) :: sw_scheme = 'sw_adding', lw_scheme = 'lw_Toon89', opac_scheme, adj_scheme

  integer :: u_nml


  namelist /FMS_RC_nml/ sw_scheme, lw_scheme, opac_scheme, adj_scheme, nlay, a_sh, b_sh, pref, &
          & t_step, nstep, Rd_air, cp_air, grav, mu_z, Tirr, Tint, k_V, k_IR, AB, fl, met, &
          & iIC, corr, table_num, sw_ac, sw_gc, lw_ac, lw_gc, use_bond_fit, Tdeep_init, Tstrat_init

  !! Read input variables from namelist
  open(newunit=u_nml, file='FMS_RC.nml', status='old', action='read')
  read(u_nml, nml=FMS_RC_nml)
  close(u_nml)

  if (nlay < 3) error stop 'WENO interpolation requires at least three layers'
  if (nstep < 0 .or. t_step <= 0.0_dp) error stop 'Require nstep >= 0 and t_step > 0'

  !! Number of layer edges (levels)
  nlev = nlay + 1

  !! Read in hybrid sigma grid values
  open(newunit=ua,file=trim(a_sh), action='read', status='old')
  open(newunit=ub,file=trim(b_sh), action='read', status='old')
  allocate(a_hs(nlev),b_hs(nlev))
  do k = 1, nlev
    read(ua,*) a_hs(k)
    read(ub,*) b_hs(k)
  end do
  close(ua); close(ub)

  !! Contruct pressure array [pa] at the levels using the hybrid sigma formula
  ! Reference surface pressure [pa] is pref
  allocate(pe(nlev))
  do k = 1, nlev
    pe(k) = a_hs(k) + b_hs(k)*pref
  end do
  pu = pe(1)

  !! Pressure at the layers
  allocate(pl(nlay),dpe(nlay))
  do k = 1, nlay
    dpe(k) = pe(k+1) - pe(k)
    pl(k) = dpe(k) / log(pe(k+1)/pe(k))
  end do

  !! Allocate other arrays we need
  allocate(Tl(nlay), dT_rad(nlay), dT_conv(nlay), net_F(nlev))
  allocate(tau_Ve(nlev,3),tau_IRe(nlev,2), k_Vl(nlay,3), k_IRl(nlay,2))
  allocate(Beta_V(3), Beta_IR(2), gam_V(3))

  allocate(sw_up(nlev), sw_down(nlev), sw_net(nlev), lw_up(nlev), lw_down(nlev), lw_net(nlev))
  allocate(band_up(nlev), band_down(nlev), band_net(nlev))
  allocate(sw_up_band(nlev,3), sw_down_band(nlev,3), sw_net_band(nlev,3))
  allocate(lw_up_band(nlev,2), lw_down_band(nlev,2), lw_net_band(nlev,2))
  allocate(Tl_rad(nlay), Kzz(nlay), Kzz_e(nlev), Rd_bar(nlay), cp_bar(nlay))
  dT_rad = 0.0_dp
  dT_conv = 0.0_dp
  net_F = 0.0_dp
  tau_Ve = 0.0_dp
  tau_IRe = 0.0_dp
  sw_up = 0.0_dp; sw_down = 0.0_dp; sw_net = 0.0_dp
  lw_up = 0.0_dp; lw_down = 0.0_dp; lw_net = 0.0_dp
  sw_up_band = 0.0_dp; sw_down_band = 0.0_dp; sw_net_band = 0.0_dp
  lw_up_band = 0.0_dp; lw_down_band = 0.0_dp; lw_net_band = 0.0_dp
  Kzz = 0.0_dp
  asr = 0.0_dp; olr = 0.0_dp

  ! Allocate cloud properties (constant all layers for now)
  allocate(sw_a(nlay,3), sw_g(nlay,3), lw_a(nlay,2), lw_g(nlay,2))
  sw_a(:,1) = sw_ac(1) ; sw_a(:,2) = sw_ac(2) ; sw_a(:,3) = sw_ac(3)
  sw_g(:,1) = sw_gc(1) ; sw_g(:,2) = sw_gc(2) ; sw_g(:,3) = sw_gc(3)
  lw_a(:,1) = lw_ac(1) ; lw_a(:,2) = lw_ac(2)
  lw_g(:,1) = lw_gc(1) ; lw_g(:,2) = lw_gc(2)

  !! Calculate the adiabatic coefficent
  kappa_air = Rd_air/cp_air   ! kappa = Rd/cp
  Rd_bar = Rd_air
  cp_bar = cp_air

  F0 = sb * Tirr**4 ! Substellar point irradiation flux
  Fint = sb * Tint**4 ! Internal flux

  print*, 'Tint ', 'Tirr ', 'pref ', 'pu ', 'mu_z ', 'grav '
  print*, Tint, Tirr, pref/1e5_dp, pu/1e5_dp, mu_z, grav
  print*, '--------------'

  ! Semi-grey atmosphere values (here they are not used, but just need to be passed to IC routine)
  k_Vl(:,:) = k_V
  k_IRl(:,:) = k_IR

  !! Initial condition T-p profile - see the routine for options
  call IC_profile(iIC,corr,nlay,pref,pl,k_Vl(:,1),k_IRl(:,1),Tint,mu_z,Tirr,grav,fl,Tl,prc,table_num,met, &
    & Tdeep_init,Tstrat_init,kappa_air)

  !! Parmentier opacity profile parameters - first get Bond albedo
  Teff = (Tint**4 + (1.0_dp/sqrt(3.0_dp)) * Tirr**4)**(0.25_dp)
  if (use_bond_fit) call Bond_Parmentier(Teff, grav, AB)
  if (AB < 0.0_dp .or. AB >= 1.0_dp) error stop 'AB must be in [0,1)'
  Finc = (1.0_dp - AB) * F0
  !! Recalculate Teff and then find parameters
  Teff = (Tint**4 + (1.0_dp - AB) * max(0.0_dp,mu_z) * Tirr**4)**(0.25_dp)
  call gam_Parmentier(Teff, table_num, gam_V, Beta_V, Beta_IR, gam_1, gam_2, gam_P, tau_lim)

  !! Print variables from Parmentier non-grey scheme
  print*, 'Teff ', 'AB ', 'gam_V ', 'Beta_V ', 'Beta_IR ', 'gam_1 ', 'gam_2 ', 'gam_P ','tau_lim ','prc'
  print*, Teff, AB, gam_V, Beta_V, Beta_IR, gam_1, gam_2, gam_P, tau_lim, prc/1e5_dp
  print*, '--------------'

  !! Print initial T-p profile
  do i = 1, nlay
    print*, i, pl(i)/1e5_dp, Tl(i)
  end do
  print*, '--------------'


  !! Write out initial conditions
  open(newunit=u,file='FMS_RC_ic.out',status='replace',action='write')
  do i = 1, nlay
    write(u,*) i, pl(i), Tl(i)
  end do
  close(u)

  !! Time stepping loop
  print*, 'Start timestepping, using: ', trim(sw_scheme), ' + ', trim(lw_scheme)

  t_tot = 0.0_dp
  n = 0

  !! cpu timer start
  call cpu_time(start)

  do n = 1, nstep

    net_F(:) = 0.0_dp
    dT_conv(:) = 0.0_dp

    select case(opac_scheme)

    case('Freedman')
      ! Calculate optical depth structure for Freedman et al. (2014) Rosseland mean fitting function scheme

      ! Include optical depth component from 0 pressure, assuming constant T and p at boundary
      call k_Ross_Freedman(Tl(1), pe(1), met, k_IRl(1,1))
      k_Vl(1,:) = k_IRl(1,1) * gam_V(:)
      k_IRl(1,2) = k_IRl(1,1) * gam_2
      k_IRl(1,1) = k_IRl(1,1) * gam_1

      tau_Ve(1,:) = (k_Vl(1,:) * pe(1)) / grav
      tau_IRe(1,:) = (k_IRl(1,:) * pe(1)) / grav

      do k = 1, nlay
        call k_Ross_Freedman(Tl(k), pl(k), met, k_IRl(k,1))
        k_Vl(k,:) = k_IRl(k,1) * gam_V(:)
        k_IRl(k,2) = k_IRl(k,1) * gam_2
        k_IRl(k,1) = k_IRl(k,1) * gam_1

        tau_Ve(k+1,:) = tau_Ve(k,:) + (k_Vl(k,:) * dpe(k)) / grav
        tau_IRe(k+1,:) = tau_IRe(k,:) + (k_IRl(k,:) * dpe(k)) / grav
      end do

    case('Valencia')
      ! Calculate optical depth structure for Valencia et al. (2013) Rosseland mean fitting function scheme

      ! Include optical depth component from 0 pressure, assuming constant T and p at boundary
      call k_Ross_Valencia(Tl(1), pe(1), met, k_IRl(1,1))
      k_Vl(1,:) = k_IRl(1,1) * gam_V(:)
      k_IRl(1,2) = k_IRl(1,1) * gam_2
      k_IRl(1,1) = k_IRl(1,1) * gam_1

      tau_Ve(1,:) = (k_Vl(1,:) * pe(1)) / grav
      tau_IRe(1,:) = (k_IRl(1,:) * pe(1)) / grav

      do k = 1, nlay
        call k_Ross_Valencia(Tl(k), pl(k), met, k_IRl(k,1))
        k_Vl(k,:) = k_IRl(k,1) * gam_V(:)
        k_IRl(k,2) = k_IRl(k,1) * gam_2
        k_IRl(k,1) = k_IRl(k,1) * gam_1

        tau_Ve(k+1,:) = tau_Ve(k,:) + (k_Vl(k,:) * dpe(k)) / grav
        tau_IRe(k+1,:) = tau_IRe(k,:) + (k_IRl(k,:) * dpe(k)) / grav
      end do

    case default
      print*, 'Invalid opac_scheme: ', trim(opac_scheme)
      stop
    end select


    !! Independent shortwave solves: the incident stellar flux is split once.
    do ib = 1, 3
      call solve_sw_band(sw_scheme, nlay, nlev, tau_Ve(:,ib), mu_z, Finc*Beta_V(ib), &
        & sw_a(:,ib), sw_g(:,ib), sw_up_band(:,ib), sw_down_band(:,ib), sw_net_band(:,ib), band_asr)
    end do
    sw_up = sum(sw_up_band, dim=2)
    sw_down = sum(sw_down_band, dim=2)
    sw_net = sum(sw_net_band, dim=2)
    asr = -sw_net(1)

    !! A grey solve is linear in all thermal sources. Weight the entire
    !! solution, including its internal flux, by the infrared band fraction.
    do ib = 1, 2
      call solve_lw_band(lw_scheme, nlay, nlev, Tl, pl, pe, tau_IRe(:,ib), lw_a(:,ib), lw_g(:,ib), Tint, &
        & band_up, band_down, band_net, band_olr)
      lw_up_band(:,ib) = Beta_IR(ib)*band_up
      lw_down_band(:,ib) = Beta_IR(ib)*band_down
      lw_net_band(:,ib) = Beta_IR(ib)*band_net
    end do
    lw_up = sum(lw_up_band, dim=2)
    lw_down = sum(lw_down_band, dim=2)
    lw_net = sum(lw_net_band, dim=2)
    olr = lw_up(1)
    net_F = sw_net + lw_net

    if (abs(sw_up(nlev)) > 1.0e-10_dp*max(1.0_dp, Finc)) &
      & error stop 'Invalid shortwave lower-boundary upward flux'
    if (lw_scheme /= 'none' .and. lw_scheme /= 'None') then
      if (abs(lw_net(nlev) - Fint) > 1.0e-7_dp*max(1.0_dp, Fint)) &
        & error stop 'Invalid longwave lower-boundary internal flux'
    end if

    !! Calculate the temperature tendency due to radiation
    do i = 1, nlay
      dT_rad(i) = (grav/cp_air) * (net_F(i+1)-net_F(i))/(dpe(i))
    end do

    !! Convective adjustment scheme
    !! The schemes act on the column after this step's radiative heating,
    !! so the column is convectively stable at the end of the step
    Tl_rad(:) = Tl(:) + t_step * dT_rad(:)
    select case(adj_scheme)
    case('Manabe')
      ! Dry convective adjustment after Manabe, Smagorinsky & Strickler (1965), no iteration
      call Manabe_dry_adjust(nlay, Tl_rad, pl, pe, kappa_air, t_step, dT_conv)
      Kzz(:) = 1e1_dp
    case('Pierrehumbert','Ray_dry')
      ! Pairwise dry convective adjustment following Ray Pierrehumbert's python script
      call Pierrehumbert_dry_adjust(nlay, Tl_rad, pl, pe, kappa_air, t_step, dT_conv)
      Kzz(:) = 1e1_dp
    case('MLT')
      ! Use mixing length theory (MLT) to time dependently adjust the adiabat and estimate Kzz
      call MLT(nlay, Tl_rad, pl, pe, Rd_bar, cp_bar, grav, t_step, dT_conv, Kzz=Kzz_e)
      ! MLT gives Kzz at the levels, output the layer average
      Kzz(:) = 0.5_dp * (Kzz_e(1:nlay) + Kzz_e(2:nlev))
    case('none','None')
      dT_conv(:) = 0.0_dp
      Kzz(:) = 1e1_dp
    case default
      print*, 'Invalid adj_scheme: ', trim(adj_scheme)
      error stop
    end select

    !! Forward march the temperature change in each layer from convection and radiation
    Tl(:) = Tl(:) + t_step * (dT_conv(:) + dT_rad(:))

    !! Fail visibly on non-finite or nonphysical temperature profiles.
    if (any(.not. ieee_is_finite(Tl)) .or. any(Tl <= 0.0_dp)) then
      print*, 'Invalid temperature at timestep: ', n
      do j = 1, nlay
        print*, j, Tl(j), net_F(j), dT_rad(j), dT_conv(j)
      end do
      error stop 'Invalid temperature profile'
    end if

    !! Increase the total time simulated
    t_tot = t_tot + t_step

  end do

  !! cpu timer end
  call cpu_time(finish)

  !! Output the results
  print*, 'sec: ', 'hours: ', 'days: '
  print*, t_tot, t_tot/60.0_dp/60.0_dp,t_tot/60.0_dp/60.0_dp/24.0_dp

  print*, 'For profile properties: '
  print*, Tint, Tirr, pref, mu_z

  print*, 'OLR [W m-2]: '
  print*, olr

  print*, 'ASR [W m-2]: '
  print*, asr

  print*, 'Outputting results: '
  open(newunit=u,file='FMS_RC_pp.out', status='replace', action='write')
  do i = 1, nlay
    write(u,*) i, pl(i), Tl(i), dT_rad(i), dT_conv(i), 0.5_dp*(tau_Ve(i+1,:)+tau_Ve(i,:)), 0.5_dp*(tau_IRe(i+1,:)+tau_IRe(i,:)), &
      & k_Vl(i,:), k_IRl(i,:), Kzz(i)
  end do
  close(u)

  ! Broadband fluxes followed by all three visible/two IR optical depths.
  open(newunit=u,file='FMS_RC_flx.out',status='replace',action='write')
  do i = 1, nlev
    write(u,*) i, pe(i), sw_up(i), sw_down(i), sw_net(i), lw_up(i), lw_down(i), lw_net(i), &
      & tau_Ve(i,:), tau_IRe(i,:)
  end do
  close(u)

  ! Band-resolved fluxes: SW up/down/net (3 each), then LW up/down/net (2 each).
  open(newunit=u,file='FMS_RC_band_flx.out',status='replace',action='write')
  do i = 1, nlev
    write(u,*) i, pe(i), sw_up_band(i,:), sw_down_band(i,:), sw_net_band(i,:), &
      & lw_up_band(i,:), lw_down_band(i,:), lw_net_band(i,:)
  end do
  close(u)

  print*, min(n,nstep), 'steps took: '
  print '("Time = ",f8.3," seconds.")', finish-start

end program Exo_FMS_RC
