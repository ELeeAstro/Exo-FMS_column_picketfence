!!!
! Elspeth KH Lee - May 2021
! A simple program that emulates a column inside the Exo-FMS GCM
! This is useful for testing different RT solutions
! This version is for semi-grey radiation only
!
! Input parameters are set via the namelist (FMS_RC.nml) - hopefully paramaters are
! somewhat self explanatory. See the github readme for more information.
! Not guarenteed to be bug free
! NOTE: Indexing starts at 1, where 1 is the Top Of Atmosphere level or layer
!!!

program Exo_FMS_RC
  use, intrinsic :: iso_fortran_env
  use ts_isothermal_mod, only : ts_isothermal
  use ts_isothermal_2_mod, only : ts_isothermal_2
  use ts_Toon_mod, only : ts_Toon
  use ts_Toon_scatter_mod, only : ts_Toon_scatter
  use ts_Heng_mod, only : ts_Heng
  use ts_short_char_mod, only : ts_short_char
  use ts_Lewis_scatter_mod, only : ts_Lewis_scatter
  use ts_disort_scatter_mod, only : ts_disort_scatter
  use ts_Mendonca_mod, only : ts_Mendonca
  use k_Rosseland_mod, only : k_Ross_Freedman, k_Ross_Valencia, gam_Parmentier, Bond_Parmentier
  use IC_mod, only : IC_profile
  use dry_conv_adj_mod, only : Ray_dry_adj
  use ieee_arithmetic
  implicit none

  ! Precision variable
  integer, parameter :: dp = REAL64

  ! Constants
  real(dp), parameter :: sb = 5.670374419e-8_dp

  integer :: n, i, k, u, j, b, inan
  integer :: nstep, nlay, nlev
  integer :: table_num
  real(dp) :: t_step, t_tot
  real(dp) :: mu_z, Tirr, Tint, F0, Fint, pref, pu, met
  real(dp) :: Teff, gam_1, gam_2, tau_lim, gam_P
  real(dp), allocatable, dimension(:) :: Tl, pl, pe, dpe
  real(dp), allocatable, dimension(:) :: Beta_V, Beta_IR, gam_V
  real(dp), allocatable, dimension(:,:) :: k_Vl, k_IRl
  real(dp), allocatable, dimension(:,:) :: tau_Ve, tau_IRe, tau_IRl
  real(dp), allocatable, dimension(:) :: dT_rad, dT_conv, net_F

  real(dp), allocatable, dimension(:,:) :: sw_a, sw_g, lw_a, lw_g

  real(dp) :: cp_air, grav, k_IR, k_V, kappa_air, Rd_air
  real(dp) :: olr, asr

  real(dp), dimension(3) :: sw_ac, sw_gc
  real(dp), dimension(2) :: lw_ac, lw_gc

  integer :: iIC
  logical :: corr
  real(dp) :: prc

  real(dp) :: fl, AB

  integer :: ua, ub
  character(len=50) :: a_sh, b_sh
  real(dp), allocatable, dimension(:) :: a_hs, b_hs

  real(dp) :: start, finish

  character(len=50) :: ts_scheme, opac_scheme, adj_scheme

  integer :: u_nml

  logical :: Bezier

  namelist /FMS_RC_nml/ ts_scheme, opac_scheme, adj_scheme, nlay, a_sh, b_sh, pref, &
          & t_step, nstep, Rd_air, cp_air, grav, mu_z, Tirr, Tint, k_V, k_IR, AB, fl, met, &
          & iIC, corr, table_num, sw_ac, sw_gc, lw_ac, lw_gc, Bezier

  !! Read input variables from namelist
  open(newunit=u_nml, file='FMS_RC.nml', status='old', action='read')
  read(u_nml, nml=FMS_RC_nml)
  close(u_nml)

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
  allocate(tau_Ve(3,nlev),tau_IRe(2,nlev), k_Vl(3,nlay), k_IRl(2,nlay))
  allocate(Beta_V(3), Beta_IR(2), gam_V(3))

  if (ts_scheme == 'Heng') then
    allocate(tau_IRl(2,nlay))
  end if

  ! Allocate cloud properties (constant all layers for now)
  allocate(sw_a(3,nlay), sw_g(3,nlay), lw_a(2,nlay), lw_g(2,nlay))
  sw_a(1,:) = sw_ac(1) ; sw_a(2,:) = sw_ac(2) ; sw_a(3,:) = sw_ac(3)
  sw_g(1,:) = sw_gc(1) ; sw_g(2,:) = sw_gc(2) ; sw_g(3,:) = sw_gc(3)
  lw_a(1,:) = lw_ac(1) ; lw_a(2,:) = lw_ac(2)
  lw_g(1,:) = lw_gc(1) ; lw_g(2,:) = lw_gc(2)

  !! Calculate the adiabatic coefficent
  kappa_air = Rd_air/cp_air   ! kappa = Rd/cp

  F0 = sb * Tirr**4 ! Substellar point irradiation flux
  Fint = sb * Tint**4 ! Internal flux

  print*, 'Tint ', 'Tirr ', 'pref ', 'pu ', 'mu_z ', 'grav '
  print*, Tint, Tirr, pref/1e5_dp, pu/1e5_dp, mu_z, grav
  print*, '--------------'

  ! Semi-grey atmosphere values (here they are not used, but just need to be passed to IC routine)
  k_Vl(1,:) = k_V
  k_IRl(1,:) = k_IR

  !! Initial condition T-p profile - see the routine for options
  call IC_profile(iIC,corr,nlay,pref,pl,k_Vl(1,:),k_IRl(1,:),Tint,mu_z,Tirr,grav,fl,Tl,prc,table_num,met)

  !! Parmentier opacity profile parameters - first get Bond albedo
  Teff = (Tint**4 + (1.0_dp/sqrt(3.0_dp)) * Tirr**4)**(0.25_dp)
  call Bond_Parmentier(Teff, grav,  AB)
  !! Recalculate Teff and then find parameters
  Teff = (Tint**4 + (1.0_dp - AB) * mu_z * Tirr**4)**(0.25_dp)
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
  open(newunit=u,file='FMS_RC_ic.out',action='readwrite')
  do i = 1, nlay
    write(u,*) i, pl(i), Tl(i)
  end do
  close(u)

  !! Time stepping loop
  print*, 'Start timestepping'

  t_tot = 0.0_dp
  inan = 0

  !! cpu timer start
  call cpu_time(start)

  do n = 1, nstep

    net_F(:) = 0.0_dp
    dT_conv(:) = 0.0_dp

    select case(opac_scheme)
    case('Freedman')
      ! Calculate optical depth structure for Freedman et al. (2014) Rosseland mean fitting function scheme
      tau_Ve(:,1) = 0.0_dp
      tau_IRe(:,1) = 0.0_dp
      do k = 1, nlay
        call k_Ross_Freedman(Tl(k), pl(k), met, k_IRl(1,k))
        k_Vl(:,k) = k_IRl(1,k) * gam_V(:)
        k_IRl(2,k) = k_IRl(1,k) * gam_2
        k_IRl(1,k) = k_IRl(1,k) * gam_1

        tau_Ve(:,k+1) = tau_Ve(:,k) + (k_Vl(:,k) * dpe(k)) / grav
        tau_IRe(:,k+1) = tau_IRe(:,k) + (k_IRl(:,k) * dpe(k)) / grav
      end do
    case('Valencia')
      ! Calculate optical depth structure for Valencia et al. (2013) Rosseland mean fitting function scheme
      tau_Ve(:,1) = 0.0_dp
      tau_IRe(:,1) = 0.0_dp
      do k = 1, nlay
        call k_Ross_Valencia(Tl(k), pl(k), met, k_IRl(1,k))
        k_Vl(:,k) = k_IRl(1,k) * gam_V(:)
        k_IRl(2,k) = k_IRl(1,k) * gam_2
        k_IRl(1,k) = k_IRl(1,k) * gam_1

        tau_Ve(:,k+1) = tau_Ve(:,k) + (k_Vl(:,k) * dpe(k)) / grav
        tau_IRe(:,k+1) = tau_IRe(:,k) + (k_IRl(:,k) * dpe(k)) / grav
      end do
   case default
      print*, 'Invalid opac_scheme: ', trim(opac_scheme)
      stop
    end select


    !! Two stream radiative transfer step
    select case(ts_scheme)
    case('Isothermal')
      ! Isothermal layers approximation
      call ts_isothermal(nlay, nlev, Tl, tau_Ve, tau_IRe, mu_z, F0, Tint, AB, Beta_V, Beta_IR, &
      & sw_a, sw_g, 0.0_dp, net_F, olr, asr)
    case('Isothermal_2')
      ! Isothermal layers approximation - first order fix for high optical depths
      call ts_isothermal_2(nlay, nlev, Tl, tau_Ve, tau_IRe, mu_z, F0, Tint, AB, Beta_V, Beta_IR, &
      & sw_a, sw_g, 0.0_dp, net_F, olr, asr)
    case('Toon')
      ! Toon method without IR scattering
      call ts_Toon(Bezier, nlay, nlev, Tl, pl, pe, tau_Ve, tau_IRe, mu_z, F0, Tint, AB, Beta_V, Beta_IR, &
      & sw_a, sw_g, 0.0_dp, net_F, olr, asr)
    case("Toon_scatter")
      !! In development !!
      ! Toon method with scattering
      call ts_Toon_scatter(Bezier, nlay, nlev, Tl, pl, pe, tau_Ve, tau_IRe, mu_z, F0, Tint, AB,  Beta_V, Beta_IR, &
      & sw_a, sw_g, lw_a, lw_g, net_F)
    case('Shortchar')
      ! Short characteristics method without IR scattering
      call ts_short_char(Bezier, nlay, nlev, Tl, pl, pe, tau_Ve, tau_IRe, mu_z, F0, Tint, AB, Beta_V, Beta_IR, &
      & sw_a, sw_g, 0.0_dp, net_F, olr, asr)
    case('Heng')
      ! Heng flux method without IR scattering
      do b = 1, 2
        tau_IRl(b,:) = (tau_IRe(b,1:nlay) + tau_IRe(b,2:nlev)) / 2.0_dp
      end do
      call ts_Heng(Bezier, nlay, nlev, Tl, pl, pe, tau_Ve, tau_IRe, tau_IRl, mu_z, F0, Tint, AB, Beta_V, Beta_IR, net_F, olr)
    case('Lewis_scatter')
      ! Neil Lewis's code with scattering
      call ts_Lewis_scatter(nlay, nlev, Tl, tau_Ve, tau_IRe, mu_z, F0, Tint, AB,  Beta_V, Beta_IR, &
      & sw_a, sw_g, lw_a, lw_g, net_F, olr, asr)
    case('Disort_scatter')
      call ts_disort_scatter(Bezier, nlay, nlev, Tl, pl, pe, tau_Ve, tau_IRe, mu_z, F0, Tint, AB,  Beta_V, Beta_IR, &
      & sw_a, sw_g, lw_a, lw_g, net_F, olr, asr)
    case('Mendonca')
      ! Mendonca method without IR scattering
      call ts_Mendonca(Bezier, nlay, nlev, Tl, pl, pe, tau_Ve, tau_IRe, mu_z, F0, Tint, AB, Beta_V, Beta_IR, &
      & sw_a, sw_g, 0.0_dp, net_F, olr, asr)
    case('None')
    case default
      print*, 'Invalid ts_scheme: ', trim(ts_scheme)
      stop
    end select

    !! Calculate the temperature tendency due to radiation
    do i = 1, nlay
      dT_rad(i) = (grav/cp_air) * (net_F(i+1)-net_F(i))/(dpe(i))
    end do

    !! Convective adjustment scheme
    select case(adj_scheme)
    case('Ray_dry')
      ! Dry convective adjustment following Ray Pierrehumbert's python script
      call Ray_dry_adj(nlay, nlev, t_step, kappa_air, Tl, pl, pe, dT_conv)
    case('None')
    case default
      print*, 'Invalid adj_scheme: ', trim(adj_scheme)
      stop
    end select

    !! Forward march the temperature change in each layer from convection and radiation
    Tl(:) = Tl(:) + t_step * (dT_conv(:) + dT_rad(:))

    !! Check for NaN's in the temperature and exit the simulation if detected
    do i = 1, nlay
      if (ieee_is_nan(Tl(i)) .eqv. .True.) then
        do j = 1, nlay
          print*, j, Tl(j), net_F(j), dT_rad(j), dT_conv(j)
        end do
        print*, nlev, net_F(nlev)
        inan = 1
        exit
      end if
    end do
    if (inan == 1) then
      exit
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
  open(newunit=u,file='FMS_RC_pp.out', action='readwrite')
  do i = 1, nlay
    write(u,*) i, pl(i), Tl(i), dT_rad(i), dT_conv(i), 0.5_dp*(tau_Ve(:,i+1)+tau_Ve(:,i)), 0.5_dp*(tau_IRe(:,i+1)+tau_IRe(:,i)), &
      & k_Vl(:,i), k_IRl(:,i)
  end do
  close(u)

  print*, nstep, 'steps took: '
  print '("Time = ",f8.3," seconds.")', finish-start

end program Exo_FMS_RC
